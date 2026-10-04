""" Initial velocity of a decelerating meteor from a fit of a single-body drag and ablation model to the lengths
    along the first part of its trajectory.

The trajectory solver estimates the initial velocity as the slope of a straight line fitted to the first part
of the meteor (at least the first 25% of the points). A meteor that already decelerates there is slower on
average over that part than at its first point, so the slope underestimates the velocity at the reference
point. For a fireball first seen at 45-60 km this can be hundreds of m/s, comparable to or larger than the
deceleration in the atmosphere above the first point, and it carries into the orbit.

Here the lengths of the points within a time from the first one (1 s by default), and optionally above a height,
are fitted with a MetSim single body, in 3D with gravity and the Coriolis acceleration:
dv/dt = -B rho(h) v^2 (m0/m)^(1/3) plus gravity, where rho is the atmosphere density and the mass ablates at the
rate set by sigma, as m/m0 = exp(sigma (v^2 - v0^2)/2) without gravity. The fitted parameters are the velocity at
t = 0, sigma and B = Gamma A rho_m^(-2/3) m0^(-1/3), so neither the mass nor the bulk density are needed: MetSim's
drag only takes that combination, and the mass it is given follows from B. The measured deceleration constrains B
and sigma, and with them the velocity at the first point. The time offsets of the stations other than the
reference one are fitted again, since the solver estimates them with a lag model that absorbs part of the
deceleration; they are only used for this fit.

MetSim erodes mass with the same law as it ablates it, so a constant erosion coefficient eta is absorbed by the
fitted sigma, which is then sigma + eta; the eroded grains are not followed. Erosion that starts within the fitted
part, or a fragmentation, is not described, and the fit does not detect it: its RMS stays below the straight
line's. In synthetic tests (20 and 30 km/s, first seen at 55 and 70 km), the main body losing 50-80% of its mass
10-20 km below the first point biased the velocity by 40-480 m/s, 5-29 times its uncertainty, and erosion with
eta = 0.1-0.3 s^2/km^2 starting 10 km below by up to 110 m/s, 5 times; fitting only the points above the event,
all were within their uncertainty. Hence the fitted part ends at a time or a height, which should be before the
first fragmentation, e.g. from the light curve. Without fragmentation, on the true lengths of 35 synthetic
meteoroids, the uncertainty of the velocity had a median of 97, 33 and 9 m/s fitting 0.5 s, 1 s and all points.

Fitting all the points of 50 synthetic meteoroids without fragmentation (15-40 km/s, B from 1.1e-3 to 2.5e-2
m^2/kg, sigma from 0.005 to 0.05 s^2/km^2, first seen at 50 and 70 km, observed while they keep 1e-3 of their mass
and 40% of their speed, with NRLMSISE-00 as the atmosphere) the velocity error over its formal uncertainty had an
RMS of 1.02 and stayed within 2.7. The fitted velocity depends on the shape of the density profile with height,
not on its scale, which B absorbs: in those cases a density 6% higher at 40 km than NRLMSISE-00, growing linearly
from 60 km, moved the velocity by up to 42 m/s, beyond 3 times its uncertainty in 6 of them. The polynomial MetSim
takes, fitted over the observed heights, was within 2.3% of NRLMSISE-00. The fit also inherits the errors of the
solver's lengths: through the solver, the synthetic fireballs of the tests come out 19-56 m/s high fitting all
their points, within 2.3 times their uncertainty.
"""

import math

import numpy as np
import scipy.optimize

from wmpl.MetSim.BackwardAtmIntegration import backwardConstants
from wmpl.MetSim.MetSimErosion import runSimulation
from wmpl.Utils.AtmosphereDensity import fitAtmPoly


# Default time limit of the fitted points from the reference time (s)
DEFAULT_TIME_LIMIT = 1.0

class DragVelocityFit(object):
    def __init__(self, v_init, v_init_stddev, intercept, drag_coeff, sigma, sigma_stddev, time_offsets, rms,
            v_init_linear, rms_linear, n_points, ht_range, t_range):
        """ Result of fitDragInitialVelocity().

        Arguments:
            v_init: [float] Velocity at the reference time, t = 0 (m/s).
            v_init_stddev: [float] Its formal 1-sigma uncertainty (m/s).
            intercept: [float] Length at t = 0 (m).
            drag_coeff: [float] B = Gamma A rho_m^(-2/3) m0^(-1/3) at t = 0 (m^2/kg).
            sigma: [float] Ablation coefficient (s^2/km^2).
            sigma_stddev: [float] Its formal 1-sigma uncertainty (s^2/km^2). Large when the ablation does not
                change the deceleration enough to be measured, which is when it does not matter for v_init.
            time_offsets: [dict] Time offsets added to each station's times for this fit (s), by station ID.
            rms: [float] RMS of the length residuals (m).
            v_init_linear: [float] The solver's straight-line initial velocity (m/s).
            rms_linear: [float] RMS of the length residuals of that straight line (m).
            n_points: [int] Number of fitted points.
            ht_range: [tuple] Lowest and highest height of the fitted points (m).
            t_range: [tuple] Earliest and latest time of the fitted points from the reference time (s).
        """

        self.v_init = v_init
        self.v_init_stddev = v_init_stddev
        self.intercept = intercept
        self.drag_coeff = drag_coeff
        self.sigma = sigma
        self.sigma_stddev = sigma_stddev
        self.time_offsets = time_offsets
        self.rms = rms
        self.v_init_linear = v_init_linear
        self.rms_linear = rms_linear
        self.n_points = n_points
        self.ht_range = ht_range
        self.t_range = t_range


def _metsimLengths(const, v0, drag_coeff, sigma, t_lo, t_hi, v_rotation):
    """ Run MetSim from t = 0 forwards to t_hi, and backwards to t_lo if it is negative, and return the times and
        the lengths along the path from t = 0 in the solver's inertial frame.

    MetSim follows the motion relative to the ground, whose speed along the path differs from the inertial one by
    the Earth's rotation, v_rotation, practically constant over the few seconds of a meteor: MetSim starts at
    v0 - v_rotation, and v_rotation*t is added back to its lengths.
    """

    # MetSim's drag is gamma*shape_factor*rho^(-2/3)*m^(-1/3), so the drag coefficient sets the mass
    const.v_init, const.sigma = v0 - v_rotation, sigma
    const.m_init = (const.gamma*const.shape_factor*const.rho**(-2/3.0)/drag_coeff)**3

    times, lengths = [np.zeros(1)], [np.zeros(1)]
    for sign, t_kill in [(1, t_hi), (-1, -t_lo)]:
        if t_kill <= 0:
            continue
        const.dt, const.t_kill = sign*abs(const.dt), t_kill
        const.h_kill = 1.0 if (sign > 0) else const.h_init + 1e5
        results = np.array(runSimulation(const)[1])
        times.append(results[:, 0])
        lengths.append(results[:, 18])

    times, lengths = np.concatenate(times), np.concatenate(lengths)
    order = np.argsort(times)

    return times[order], lengths[order] + v_rotation*times[order]


def fitDragInitialVelocity(traj, ht_min=None, t_max=None, fine_dt=0.001):
    """ Fit the single-body drag and ablation model to the lengths of the non-ignored points of a solved
        trajectory above a height and before a time (see the module docstring).

    Arguments:
        traj: [Trajectory] Solved trajectory, with time_data, state_vect_dist and model_ht of its observations.

    Keyword arguments:
        ht_min: [float] Only points above this height are fitted (m). None by default, for no limit.
        t_max: [float] Only points before this time from the reference time, the first point in a solved
            trajectory, are fitted (s). None by default, for no limit.
        fine_dt: [float] The larger of the two MetSim time steps whose fits are extrapolated to a zero step (s).
            The starting fit uses MetSim's default step.

    Return:
        [DragVelocityFit] or None if there are not more points than parameters, or the fit did not converge or
            does not fit better than the straight line.
    """

    ref_id = traj.observations[traj.t_ref_station].station_id

    observations, times, lengths, heights, station_index = [], [], [], [], []
    for obs in traj.observations:
        good = (obs.ignore_list == 0) & (obs.model_ht > (-np.inf if ht_min is None else ht_min)) \
            & (obs.time_data < (np.inf if t_max is None else t_max))
        if obs.ignore_station or (not np.any(good)):
            continue
        station_index.append(np.full(np.count_nonzero(good), len(observations)))
        observations.append(obs)
        times.append(obs.time_data[good])
        lengths.append(obs.state_vect_dist[good])
        heights.append(obs.model_ht[good])
    if not observations:
        return None
    times, lengths, heights = np.concatenate(times), np.concatenate(lengths), np.concatenate(heights)
    station_index = np.concatenate(station_index)

    # The stations whose time offset is fitted, all but the reference one, or but the first one if the reference
    #   station has no points in the fitted part, since a common offset is the intercept
    offset_stations = [i for i, obs in enumerate(observations) if obs.station_id != ref_id]
    if len(offset_stations) == len(observations):
        offset_stations = offset_stations[1:]
    if len(times) <= 4 + len(offset_stations):
        return None

    # MetSim starts from the solver's state vector, relative to the ground, with the air density over the heights
    #   the meteor reaches. The bulk density only scales the mass that B gives, so its value does not matter
    ht_lo, ht_hi = np.min(heights) - 5000.0, np.max(heights) + 5000.0
    const = backwardConstants(traj.jdt_ref, np.concatenate([traj.state_vect_mini, traj.v_init*traj.radiant_eci_mini]),
        1.0, h_kill=ht_hi)
    const.dens_co = fitAtmPoly(traj.rbeg_lat, traj.rbeg_lon, ht_lo, ht_hi, traj.jdt_ref)
    const.v_kill = 100.0
    v_rotation = traj.v_init - const.v_init

    v_lin, intercept_lin = traj.velocity_fit
    res_lin = lengths - (v_lin*times + intercept_lin)
    rms_linear = np.sqrt(np.mean(res_lin**2))

    def residuals(params):
        v0, log_b, intercept, sigma = params[:4]
        offsets = np.zeros(len(observations))
        offsets[offset_stations] = params[-len(offset_stations):] if offset_stations else []
        t_model = times + offsets[station_index]
        try:
            t_grid, l_grid = _metsimLengths(const, v0, math.exp(log_b), sigma*1e-6, np.min(t_model), np.max(t_model),
                v_rotation)
        except (ValueError, OverflowError, ZeroDivisionError):
            return np.full_like(lengths, 1e9)
        res = np.interp(t_model, t_grid, l_grid) + intercept - lengths
        return res if np.all(np.isfinite(res)) else np.full_like(res, 1e9)

    n_off = len(offset_stations)
    # The solver's own offsets can be off by as much as it allows them, so they can be corrected as much
    lb = [0.5*v_lin, -25.0, intercept_lin - 1e4, 0.0] + [-traj.max_toffset]*n_off
    ub = [1.5*v_lin, 3.0, intercept_lin + 1e4, 0.5] + [traj.max_toffset]*n_off
    x_scale = [100.0, 1.0, 10.0, 0.005] + [1e-3]*n_off
    # The straight line's residuals are dominated by the deceleration, so this robust scale is loose and only
    #   downweights gross outliers
    f_scale = max(1.4826*np.median(np.abs(res_lin - np.median(res_lin))), 1.0)

    def fit(p0):
        return scipy.optimize.least_squares(residuals, p0, bounds=(lb, ub), x_scale=x_scale, loss="soft_l1",
            f_scale=f_scale)

    # MetSim advances the speed and the mass one after the other, so its lengths are first order in the time
    #   step, which biases the velocity of a strongly ablating meteoroid high. The fit starts with MetSim's default
    #   step, then fits with fine_dt and half of it are extrapolated to a zero step
    start = fit([v_lin, math.log(1e-3), intercept_lin, 0.01] + [0.0]*n_off)
    const.dt = fine_dt
    full = fit(start.x)
    const.dt = fine_dt/2
    best = fit(full.x)
    params = np.clip(2*best.x - full.x, lb, ub)

    # The straight line is the model without drag, so a worse fit than it means the fit went wrong
    rms = np.sqrt(np.mean(best.fun**2))
    if (best.status <= 0) or (not np.all(np.isfinite(params))) or (rms > rms_linear):
        return None

    # Formal uncertainties from the Jacobian, scaled by the residual variance
    dof = max(len(lengths) - len(params), 1)
    cov = np.linalg.pinv(best.jac.T.dot(best.jac))*np.sum(best.fun**2)/dof
    stddev = np.sqrt(np.abs(np.diag(cov)))

    time_offsets = {obs.station_id: 0.0 for obs in observations}
    for k, i in enumerate(offset_stations):
        time_offsets[observations[i].station_id] = params[4 + k]

    return DragVelocityFit(params[0], stddev[0], params[2], math.exp(params[1]), params[3], stddev[3],
        time_offsets, rms, v_lin, rms_linear, len(times), (np.min(heights), np.max(heights)),
        (np.min(times), np.max(times)))
