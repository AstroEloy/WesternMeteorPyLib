""" Initial velocity of a decelerating meteor from a fit of a single-body drag and ablation model to the lengths
    along its whole trajectory.

The trajectory solver estimates the initial velocity as the slope of a straight line fitted to the first part
of the meteor (at least the first 25% of the points). A meteor that already decelerates there is slower on
average over that part than at its first point, so the slope underestimates the velocity at the reference
point. For a fireball first seen at 45-60 km this can be hundreds of m/s, comparable to or larger than the
deceleration in the atmosphere above the first point, and it carries into the orbit.

Here the lengths of all points are fitted with a MetSim single body, in 3D with gravity and the Coriolis
acceleration: dv/dt = -B rho(h) v^2 (m0/m)^(1/3) plus gravity, where rho is the atmosphere density and the mass
ablates as m/m0 = exp(sigma (v^2 - v0^2)/2). The fitted parameters are the velocity at t = 0, sigma and
B = Gamma A rho_m^(-2/3) m0^(-1/3), so neither the mass nor the bulk density are needed: MetSim's drag only takes
that combination, and the mass it is given follows from B. The deceleration
measured further down, where it is large, constrains B and sigma, and with them the velocity at the first point.
The time offsets of the stations other than the reference one are fitted again, since the solver estimates them
with a lag model that absorbs part of the deceleration; they are only used for this fit.

The model is a single body, so it does not describe fragmentation. The fitted velocity depends on the shape of
the density profile with height, not on its scale, which B absorbs: in synthetic tests a 6% error in the shape
between 40 and 60 km moved sigma by up to 65% and the velocity by up to 100 m/s. It was tested for ablation
coefficients up to 0.05 s^2/km^2, and it inherits the errors of the solver's lengths.
"""

import math

import numpy as np
import scipy.optimize

from wmpl.MetSim.BackwardAtmIntegration import backwardConstants
from wmpl.MetSim.MetSimErosion import runSimulation
from wmpl.Utils.AtmosphereDensity import fitAtmPoly


class DragVelocityFit(object):
    def __init__(self, v_init, v_init_stddev, intercept, drag_coeff, sigma, sigma_stddev, time_offsets, rms,
            v_init_linear, rms_linear):
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


def fitDragInitialVelocity(traj, fine_dt=0.001):
    """ Fit the single-body drag and ablation model to the lengths of all non-ignored points of a solved
        trajectory (see the module docstring).

    Arguments:
        traj: [Trajectory] Solved trajectory, with time_data, state_vect_dist and model_ht of its observations.

    Keyword arguments:
        fine_dt: [float] MetSim time step of the final fit (s). The starting fits use MetSim's default step.

    Return:
        [DragVelocityFit] or None if the fit did not converge or does not fit better than the straight line.
    """

    observations = [obs for obs in traj.observations if not obs.ignore_station]
    ref_id = traj.observations[traj.t_ref_station].station_id

    times, lengths, heights, station_index = [], [], [], []
    for i, obs in enumerate(observations):
        good = obs.ignore_list == 0
        times.append(obs.time_data[good])
        lengths.append(obs.state_vect_dist[good])
        heights.append(obs.model_ht[good])
        station_index.append(np.full(np.count_nonzero(good), i))
    times, lengths, heights = np.concatenate(times), np.concatenate(lengths), np.concatenate(heights)
    station_index = np.concatenate(station_index)

    # The stations whose time offset is fitted, all but the reference one
    offset_stations = [i for i, obs in enumerate(observations) if obs.station_id != ref_id]

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

    def residuals(params, sigma_fixed=None):
        v0, log_b, intercept = params[:3]
        sigma = sigma_fixed if (sigma_fixed is not None) else params[3]
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
    lb = [0.5*v_lin, -25.0, intercept_lin - 1e4] + [-traj.max_toffset]*n_off
    ub = [1.5*v_lin, 3.0, intercept_lin + 1e4] + [traj.max_toffset]*n_off
    x_scale = [100.0, 1.0, 10.0] + [1e-3]*n_off
    # The straight line's residuals are dominated by the deceleration, so this robust scale is loose and only
    #   downweights gross outliers
    f_scale = max(1.4826*np.median(np.abs(res_lin - np.median(res_lin))), 1.0)

    def fit(p0, lower, upper, scale, sigma_fixed=None):
        return scipy.optimize.least_squares(residuals, p0, bounds=(lower, upper), x_scale=scale, loss="soft_l1",
            f_scale=f_scale, kwargs={"sigma_fixed": sigma_fixed})

    # Without ablation first, then with the ablation coefficient free from three starting values, keeping the
    #   lowest cost: started at 0.05 s^2/km^2 alone, fits of meteoroids with that coefficient ended 7 km/s off
    first = fit([v_lin, math.log(1e-3), intercept_lin] + [0.0]*n_off, lb, ub, x_scale, sigma_fixed=0.0)
    v0, log_b, intercept = first.x[:3]
    offsets = list(first.x[3:])
    best = min((fit([v0, log_b, intercept, s] + offsets, lb[:3] + [0.0] + lb[3:], ub[:3] + [0.5] + ub[3:],
        x_scale[:3] + [0.005] + x_scale[3:]) for s in (0.001, 0.01, 0.05)), key=lambda r: r.cost)

    # MetSim advances the speed and the mass one after the other, to first order in the time step, which biases
    #   the velocity of a strongly ablating meteoroid high; the final fit refines the best start with a finer step
    const.dt = fine_dt
    best = fit(best.x, lb[:3] + [0.0] + lb[3:], ub[:3] + [0.5] + ub[3:], x_scale[:3] + [0.005] + x_scale[3:])

    # The straight line is the model without drag, so a worse fit than it means the fit went wrong
    rms = np.sqrt(np.mean(best.fun**2))
    if (best.status <= 0) or (not np.all(np.isfinite(best.x))) or (rms > rms_linear):
        return None

    # Formal uncertainties from the Jacobian, scaled by the residual variance
    dof = max(len(lengths) - len(best.x), 1)
    cov = np.linalg.pinv(best.jac.T.dot(best.jac))*np.sum(best.fun**2)/dof
    stddev = np.sqrt(np.abs(np.diag(cov)))

    time_offsets = {obs.station_id: 0.0 for obs in observations}
    for k, i in enumerate(offset_stations):
        time_offsets[observations[i].station_id] = best.x[4 + k]

    return DragVelocityFit(best.x[0], stddev[0], best.x[2], math.exp(best.x[1]), best.x[3], stddev[3],
        time_offsets, rms, v_lin, rms_linear)
