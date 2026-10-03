""" Initial velocity of a decelerating meteor from a fit of a single-body drag and ablation model to the lengths
    along its whole trajectory.

The trajectory solver estimates the initial velocity as the slope of a straight line fitted to the first part
of the meteor (at least the first 25% of the points). A meteor that already decelerates there is slower on
average over that part than at its first point, so the slope underestimates the velocity at the reference
point. For a fireball first seen at 45-60 km this can be hundreds of m/s, comparable to or larger than the
deceleration in the atmosphere above the first point, and it carries into the orbit.

Here the lengths of all points are fitted with dv/dt = -B rho(h) v^2 (m0/m)^(1/3) + g cos(z), where rho is the
atmosphere density, B = Gamma A rho_m^(-2/3) m0^(-1/3) is fitted, so neither the mass nor the bulk density are
needed, and the mass ablates as m/m0 = exp(sigma (v^2 - v0^2)/2) with sigma fitted too. The deceleration
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

from wmpl.Utils.AtmosphereDensity import atmDensPoly, fitAtmPoly


# Integration step of the model (s). In synthetic tests the fitted velocity did not change between 2 and 0.5 ms
DRAG_FIT_DT = 0.002

# Earth's surface gravity (m/s^2) and mean radius (m), for the gravity along the path
G0 = 9.81
R_EARTH = 6371008.7714


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


def _integrateLengths(v0, drag_coeff, sigma, t_lo, t_hi, ht_poly, log_dens, ht_lo, ht_step):
    """ Integrate the model from t = 0 forwards to t_hi and backwards to t_lo (RK2), and return the times and the
        lengths from t = 0. Pure Python floats, as this runs for every residual evaluation. """

    h0, h1, h2 = ht_poly
    n_dens = len(log_dens) - 1

    def accel(length, v):
        ht = h0 + h1*length + h2*length**2
        x = min(max((ht - ht_lo)/ht_step, 0.0), n_dens - 1e-9)
        i = int(x)
        dens = 10**(log_dens[i] + (x - i)*(log_dens[i + 1] - log_dens[i]))
        cos_z = -(h1 + 2*h2*length)
        g = G0*(R_EARTH/(R_EARTH + ht))**2
        return -drag_coeff*math.exp(-sigma*(v*v - v0*v0)/6)*dens*v*v + g*cos_z

    def run(t_end, dt):
        times, lengths = [0.0], [0.0]
        length, v = 0.0, v0
        for _ in range(int(math.ceil(abs(t_end)/DRAG_FIT_DT)) + 1):
            a1 = accel(length, v)
            vm, lm = v + a1*dt/2, length + v*dt/2
            v, length = v + accel(lm, vm)*dt, length + vm*dt
            times.append(times[-1] + dt)
            lengths.append(length)
        return times, lengths

    t_fwd, l_fwd = run(t_hi, DRAG_FIT_DT)
    if t_lo < 0:
        t_bwd, l_bwd = run(t_lo, -DRAG_FIT_DT)
        return np.array(t_bwd[:0:-1] + t_fwd), np.array(l_bwd[:0:-1] + l_fwd)

    return np.array(t_fwd), np.array(l_fwd)


def fitDragInitialVelocity(traj):
    """ Fit the single-body drag and ablation model to the lengths of all non-ignored points of a solved
        trajectory (see the module docstring).

    Arguments:
        traj: [Trajectory] Solved trajectory, with time_data, state_vect_dist and model_ht of its observations.

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

    # Height along the straight path, with its curvature over the Earth, and the air density over the heights the
    #   model can reach, tabulated every 100 m
    ht_poly = np.polyfit(lengths, heights, 2)[::-1]
    ht_lo, ht_hi, ht_step = np.min(heights) - 5000.0, np.max(heights) + 5000.0, 100.0
    dens_co = fitAtmPoly(traj.rbeg_lat, traj.rbeg_lon, ht_lo, ht_hi, traj.jdt_ref)
    log_dens = list(np.log10([atmDensPoly(ht, dens_co) for ht in np.arange(ht_lo, ht_hi + ht_step, ht_step)]))

    v_lin, intercept_lin = traj.velocity_fit
    res_lin = lengths - (v_lin*times + intercept_lin)
    rms_linear = np.sqrt(np.mean(res_lin**2))

    def residuals(params, sigma_fixed=None):
        v0, log_b, intercept = params[:3]
        sigma = sigma_fixed if (sigma_fixed is not None) else params[3]
        offsets = np.zeros(len(observations))
        offsets[offset_stations] = params[-len(offset_stations):] if offset_stations else []
        t_model = times + offsets[station_index]
        t_grid, l_grid = _integrateLengths(v0, math.exp(log_b), sigma*1e-6, min(np.min(t_model), 0.0),
            np.max(t_model), ht_poly, log_dens, ht_lo, ht_step)
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
