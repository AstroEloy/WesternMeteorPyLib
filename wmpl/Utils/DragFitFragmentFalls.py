""" Falls of the fragments of a fireball with OpenDarkflight: where the fragments which broke off at given
heights fall (direct), and where the meteorites which were found started their dark flight and broke off
(inverse).

A fragment breaks off the leading body with its speed, and slows down by its own ballistic coefficient K,
faster than the leading body as it is smaller. The fragment can be given by its mass, with its density, shape
and the drag model of OpenDarkflight, and then K comes from the same drag model at every step, or by a fixed
K.

Without ablation (--sigma 0, the default), OpenDarkflight flies the fragment from the breakup, with the speed
and the direction of the leading body there: the whole flight in three dimensions, the path bent by gravity
and pushed by the whole wind, exact within the drag model. A run from the breakup takes some tens of seconds,
as OpenDarkflight takes steps of 0.1 m in height above 5.8 km/s.

OpenDarkflight does not ablate, so with ablation (--sigma above 0) the flight down to where the fragment has
slowed down to a few km/s (--v_dark, 3.5 km/s by default) is computed here instead
(wmpl.Utils.LuminousFlight), in the air density and the winds of the same atmosphere profile, and
OpenDarkflight takes it over from there. That flight keeps to the ground track of the leading body: it takes
the gravity and the wind along the motion, but not the bending of the path by gravity nor the wind across it.
On a synthetic fireball with a 20 m/s wind, against OpenDarkflight from the breakup, it moved the falls of 4 g
fragments by 0.4 to 16 m and of 60 g ones by 40 to 85 m along the motion and 15 to 48 m across, less for a
higher --v_dark, as a large fragment flies far below the breakup before it is slow. Two fragments of similar K
are moved alike. The dark flight never ablates.

The sideways kick of the breakup can be added (--kicks) after Passey and Melosh (1980), Eq. 24, the speed
across U = V*sqrt(3/2*C*(R1/R2)*rho_air/rho_met), with C = 1 (they fit 0.02 to 1.52), V and rho_air at the
breakup and R1/R2 the ratio of the radii of the leading body and the fragment, taken from their K as if they
had the same shape and density. It turns the direction of motion by atan(U/V) in a random direction around it,
which drag keeps.

Direct: for every breakup height and every fragment, OpenDarkflight flies it from the breakup (or from where
it has slowed down, with ablation), with the kicks drawn and the Monte Carlo of OpenDarkflight for the errors
of the wind, and the falls are collected by breakup in fragment_falls.csv and fragment_falls.png.

Inverse: for every meteorite found, with its mass, density and shape, OpenDarkflight flies it from breakups on
the line at heights every --scan_step, the breakup whose fall is as far along the trajectory as the meteorite
is found, and the distance across the trajectory left over tells how well the wind and the line agree with it.
The errors of the wind (and the kicks and the lines, if drawn) give the range of breakup heights. The breakup
assumes that the meteorite flew alone from the leading body, with its speed. With ablation, the dark flight is
started at --v_dark on the line instead, the start whose fall is as far along is found, and flying back up
from it gives the breakup, where the fragment had the speed of the leading body, and the mass it had there. As
the wind is the same for all the meteorites, their differences are better known than where each broke off. The
results go to inverse.json, inverse.txt and inverse.png.

WMPL does not import OpenDarkflight: it writes its inputs and runs it from its own environment, and also runs
the flight down to the dark flight there when K comes from its drag models.

On the command line, with the heights of the flares (-f), the atmosphere profile and OpenDarkflight:

    python -m wmpl.Utils.DragFitFragmentFalls traj.pickle out_dir -f 35.0 30.0 --atm_profile profile.csv \\
        --odf_python <environment of OpenDarkflight>/bin/python --odf_dir <OpenDarkflight directory> \\
        --direct 35.0 30.0 --masses 0.001 0.01 0.1 1 --density 3500 --shape sphere --drag_model dfn \\
        --kicks 10 --mc_runs 5 --end_ht 0.25 --run

    python -m wmpl.Utils.DragFitFragmentFalls traj.pickle out_dir -f 35.0 30.0 --atm_profile profile.csv \\
        --odf_python ... --odf_dir ... --inverse meteorites.csv --mc_runs 30

The meteorites file is a CSV with the columns name, lat, lon (deg), mass_kg and density (kg/m^3), and
optionally ground_ht_m, shape, drag_model, sphericity, flatness, elongation and volume_factor as in
OpenDarkflight, or k_m2_kg for a fixed K.

The flight down to the dark flight can ablate (--sigma, s^2/km^2, 0 by default), with the mass following the
speed through the air as m = m0*exp(sigma/2*(v^2 - v0^2)): sigma is the fraction of the mass lost for every
unit of kinetic energy per kg dissipated. The fireball types of Ceplecha et al. (1998, Table XVII) give
apparent values, 0.014 s^2/km^2 for the type I, and the alpha-beta fit of the leading body
(wmpl.Utils.AlphaBeta) gives sigma = 2*beta/((1 - mu)*V_e^2) (Gritsevich 2009), which holds its fragmentation
too and is an upper value for a fragment. In the inverse with ablation, the start of the dark flight does not
depend on it, and several values give the breakup for each.

    References:
        Passey, Q. R., & Melosh, H. J. (1980). Effects of atmospheric breakup on crater field formation.
            Icarus, 42(2), 211-233.
        Ceplecha, Z., Borovicka, J., Elford, W. G., ReVelle, D. O., Hawkes, R. L., Porubcan, V., & Simek, M.
            (1998). Meteor phenomena and bodies. Space Science Reviews, 84(3-4), 327-471.
        Gritsevich, M. I. (2009). Determination of parameters of meteor bodies based on flight observational
            data. Advances in Space Research, 44(3), 323-334.
"""

from __future__ import print_function, division, absolute_import

import argparse
import csv
import json
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np

from wmpl.Utils.AtmosphereProfile import addProfileArguments, profileFromArguments
from wmpl.Utils.DragFit import DragFit
from wmpl.Utils.DragFitDarkFlight import _cleanRunDir, alongCross, checkDarkFlightInputs, drawLine, \
    findMcUncertainties, handoverState, loadMcLines, readOdfRun, runOdf, writeOdfInput
from wmpl.Utils.Pickling import loadPickle
from wmpl.Utils.TrajConversions import jd2Date


# Speed relative to the ground where the dark flight of a fragment starts (m/s)
V_DARK = 3500.0

# Constant C of Passey and Melosh (1980), who fit 0.02 to 1.52 and adopt 1
PM_C = 1.0

# Step of the heights scanned by the inverse (m)
SCAN_STEP = 250.0

# How far above the highest fragmentation the inverse scans the breakups by default (m)
BREAKUP_MARGIN = 5000.0

# Shape and drag model of a meteorite which gives none, for both the flight down and OpenDarkflight
SHAPE = 'sphere'
DRAG_MODEL = 'dfn'

# Options of the meteorite which go to OpenDarkflight as they are
SHAPE_OPTIONS = ['sphericity', 'flatness', 'elongation', 'volume_factor']

# Longest time the flights of the fragments may take (s)
LUMINOUS_TIMEOUT = 3600.0

# Where the lines drawn around the fitted one come from
LINES_FROM = {'monte_carlo': "the Monte Carlo of the solution", \
    'covariance': "the covariance of the line-of-sight fit", 'none': "nowhere"}

# Flight down to the dark flight, run with the Python of OpenDarkflight
LUMINOUS_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'LuminousFlight.py')



def fragmentLine(drag_fit):
    """ Take the ground track of the leading body from a drag fit, along its descending modelled part, for the
        flights of the fragments (see wmpl.Utils.LuminousFlight.lineAt()). Drag keeps the direction of motion
        relative to the ground, so a fragment flies along the same track, and the lengths along it are those
        the leading body flies relative to the ground, its inertial lengths less the rotation of the Earth
        along the motion times the time.

    Arguments:
        drag_fit: [DragFit object] Drag fit of the leading body.

    Return:
        [dict] Lengths along the ground track (m), the inertial lengths of the drag fit there (m) and the
            times of the leading body (s), heights (m), air density (kg/m^3), the wind along the motion (m/s,
            the speed relative to the ground less that relative to the air), the speed of the leading body
            relative to the ground (m/s) and the gravity along the motion (m/s^2), as lists.
    """

    len_arr = drag_fit.len_grid[drag_fit._valid]
    ht_arr = drag_fit.ht_grid[drag_fit._valid]
    vel_ground = drag_fit.vel_grid[drag_fit._valid] - drag_fit.v_rotation

    # Keep the descending part, down to the lowest point, while the leading body moves on over the ground
    n = np.argmin(ht_arr) + 1
    if np.any(vel_ground[:n] <= 0):
        n = np.argmax(vel_ground <= 0)
    time_arr = drag_fit.time_grid[drag_fit._valid][:n]

    line = dict(len=len_arr[:n] - drag_fit.v_rotation*(time_arr - time_arr[0]), len_inertial=len_arr[:n], \
        time=time_arr, ht=ht_arr[:n], rho=drag_fit.airDensityAtLength(len_arr[:n]), \
        air_offset=drag_fit.v_air_offset[drag_fit._valid][:n] - drag_fit.v_rotation, \
        vel_parent=vel_ground[:n], g_along=drag_fit.g_grid[drag_fit._valid][:n])


    return {key: np.asarray(value, dtype=float).tolist() for key, value in line.items()}



def groundLength(line, len_inertial):
    """ Length along the ground track (m) at an inertial length of the drag fit (m), NaN outside it. """

    return float(np.interp(len_inertial, line['len_inertial'], line['len'], left=np.nan, right=np.nan))



def inertialLength(line, len_ground):
    """ Inertial length of the drag fit (m) at a length along the ground track (m), NaN outside it. """

    return float(np.interp(len_ground, line['len'], line['len_inertial'], left=np.nan, right=np.nan))



def trackState(drag_fit, line, len_ground, k_val, vel_ground, time=None, line_drawn=None, kick=None):
    """ Compute the state of a dark flight on the ground track of the leading body (see handoverState()).

    Arguments:
        drag_fit: [DragFit object] Drag fit of the leading body.
        line: [dict] Its ground track (see fragmentLine()).
        len_ground: [float] Length along the ground track (m).
        k_val: [float] Ballistic coefficient (m^2/kg).
        vel_ground: [float] Speed relative to the ground (m/s).

    Keyword arguments:
        time: [float] Time of the fragment there (s). None by default, for that of the leading body.
        line_drawn: [tuple] Line drawn around the fitted one (see wmpl.Utils.DragFitDarkFlight.drawLine()).
            None by default.
        kick: [tuple] Sideways kick of the breakup (see handoverState()). None by default.

    Return:
        [dict] The state, at the place on the ground where the leading body was at that length.
    """

    len_val = inertialLength(line, len_ground)
    time_parent = float(np.interp(len_val, line['len_inertial'], line['time']))

    state = handoverState(drag_fit, len_val, k_val, vel_ground + drag_fit.v_rotation, time_parent, \
        line=line_drawn, kick=kick)
    if time is not None:
        state['time'] = float(time)


    return state



def leadingBallisticCoeff(drag_fit, len_val):
    """ Ballistic coefficient of the leading body at an inertial length (m^2/kg), taken within the observed
        part, as the drag fit has no drag above it.
    """

    len_obs = np.clip(len_val, drag_fit.len_first, drag_fit.len_last)

    return float(drag_fit.ballisticCoeffAtLength(np.array([len_obs]))[0][0])



def pmKickAngle(vel, rho_air, rho_met, radius_ratio, pm_c=PM_C):
    """ Compute the angle by which the sideways kick of a breakup turns a fragment, after Passey and Melosh
        (1980), Eq. 24.

    Arguments:
        vel: [float] Speed at the breakup (m/s), which the angle does not depend on.
        rho_air: [float] Air density at the breakup (kg/m^3).
        rho_met: [float] Bulk density of the fragment (kg/m^3).
        radius_ratio: [float] Ratio of the radii of the larger and the smaller body, at least 1.

    Keyword arguments:
        pm_c: [float] Constant C. PM_C by default.

    Return:
        [float] Angle (radians), atan(U/V).
    """

    return float(np.arctan(np.sqrt(1.5*pm_c*max(radius_ratio, 1.0)*rho_air/rho_met)))



def _runLuminous(work_dir, name, line, flights, atm_profile, odf_python=None, odf_dir=None):
    """ Run the flights of the fragments along the line (see wmpl.Utils.LuminousFlight.runJobs()), with the
        Python of OpenDarkflight if any K comes from its drag models, and else with this one.
    """

    jobs = dict(line=line, profile=dict(path=os.path.abspath(atm_profile.path), \
        type=atm_profile.profile_type), flights=flights)

    jobs_path = os.path.join(work_dir, name + '_jobs.json')
    results_path = os.path.join(work_dir, name + '_results.json')
    with open(jobs_path, 'w') as f:
        json.dump(jobs, f)

    # Fixed values of K need nothing from OpenDarkflight
    if all(flight.get('k') is not None for flight in flights):
        python, cwd = sys.executable, None

    else:

        if (odf_python is None) or (odf_dir is None):
            raise ValueError("K comes from the drag models of OpenDarkflight, so its Python and directory "
                "are needed (odf_python, odf_dir)")

        python, cwd = odf_python, odf_dir

    try:
        result = subprocess.run([python, LUMINOUS_SCRIPT, jobs_path, results_path], cwd=cwd, \
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, \
            timeout=LUMINOUS_TIMEOUT)
    except subprocess.TimeoutExpired:
        raise RuntimeError("The flights of the fragments took over {:.0f} s".format(LUMINOUS_TIMEOUT))

    if result.returncode != 0:
        raise RuntimeError("The flights of the fragments failed:\n{:s}".format(result.stdout[-2000:]))

    with open(results_path) as f:
        return json.load(f)



def _fragmentFlight(meteorite, mass, sigma):
    """ The part of a flight which describes the fragment: its fixed K, or its meteorite for the drag models
        of OpenDarkflight.
    """

    if meteorite.get('k') is not None:
        return dict(k=float(meteorite['k']), mass=float(mass), sigma=sigma)

    return dict(meteorite=meteorite, mass=float(mass), sigma=sigma)



def resolveMeteorite(meteorite):
    """ Give a meteorite the shape and the drag model it does not give, SHAPE and DRAG_MODEL, so that the
        flight down to the dark flight and OpenDarkflight, whose own default shape is a brick, take the same.

    Arguments:
        meteorite: [dict] Density (kg/m^3), and shape, drag_model and the SHAPE_OPTIONS of OpenDarkflight.

    Return:
        [dict] A copy with the shape and the drag model.
    """

    meteorite = dict(meteorite)
    if meteorite.get('shape') is None:
        meteorite['shape'] = SHAPE
    if meteorite.get('drag_model') is None:
        meteorite['drag_model'] = DRAG_MODEL


    return meteorite



def _odfMeteorite(meteorite, state, mass):
    """ Keyword arguments of writeOdfInput() for the meteorite: its mass, or its K for a fixed K, and its
        shape.
    """

    return dict(shape=meteorite['shape'], drag_model=meteorite['drag_model'], \
        mass=None if (meteorite.get('k') is not None) else mass, \
        shape_options={key: meteorite[key] for key in SHAPE_OPTIONS if meteorite.get(key) is not None})



def _setup(traj, knot_hts, atm_profile):
    """ Check the inputs, fit the drag of the leading body with its fragmentations as jumps, and take its
        line, its state at the last observed point and the reference time.
    """

    checkDarkFlightInputs(traj, atm_profile)

    drag_fit = DragFit(traj, knot_hts=knot_hts, atm_profile=atm_profile)
    line = fragmentLine(drag_fit)

    len_last = np.array([drag_fit.len_last])
    last = handoverState(drag_fit, drag_fit.len_last, drag_fit.ballisticCoeffAtLength(len_last)[0][0], \
        drag_fit.velocityAtLength(len_last)[0], float(drag_fit.timeAtLength(len_last)[0]))
    reference = dict(lat=last['lat'], lon=last['lon'], ht=last['ht'], azim_motion=(last['azim'] + 180.0)%360)

    ref_time = jd2Date(traj.jdt_ref, dt_obj=True).strftime('%Y-%m-%d %H:%M:%S.%f')


    return drag_fit, line, reference, ref_time



def _odfOptions(odf_options, mc_runs):
    """ Give the Monte Carlo of OpenDarkflight its statistical wind errors, unless the options do. """

    odf_options = list(odf_options or [])
    if mc_runs > 0:
        odf_options = ['Monte Carlo.wind_error_mode=statistical'] + odf_options

    return odf_options



def prepareFragmentFalls(traj, out_dir, knot_hts, atm_profile, sep_hts, masses, meteorite, v_dark=V_DARK, \
    n_kicks=0, pm_c=PM_C, sigma=0.0, n_lines=0, mc_lines=None, end_ht=None, mc_runs=0, seed=0, \
    odf_options=None, odf_python=None, odf_dir=None):
    """ Prepare the flights of the fragments which broke off at the given heights: from the breakup
        without ablation, or, with ablation, from where the fragment has slowed down to v_dark, flying it down
        to there; draw its kicks, and write an OpenDarkflight input for each, with fragment_falls.json. Run
        them with runOdf() and collect the falls with summarizeFragmentFalls().

    Arguments:
        traj: [Trajectory object] Solved trajectory, with its orbit.
        out_dir: [str] Output directory.
        knot_hts: [list] Heights of the fragmentations of the drag fit (m), e.g. those of the flares.
        atm_profile: [AtmosphereProfile object] Atmosphere profile, with its winds, which has to cover the
            breakup heights.
        sep_hts: [list] Heights where the fragments broke off (m).
        masses: [list] Masses of the fragments (kg), or the reference masses of fixed K values.
        meteorite: [dict] Density (kg/m^3), and shape, drag_model and the SHAPE_OPTIONS of OpenDarkflight, or
            'k', a list of fixed K values (m^2/kg) at the masses.

    Keyword arguments:
        v_dark: [float] Speed relative to the ground where the dark flight starts with ablation (m/s). V_DARK
            by default.
        n_kicks: [int] Kicks of Passey and Melosh drawn for every fragment, on top of its run without a kick.
            0 by default.
        pm_c: [float] Constant C of the kicks. PM_C by default.
        sigma: [float] Ablation coefficient in the flight down to the dark flight (s^2/m^2). 0 by default, and
            then OpenDarkflight flies the fragment from the breakup.
        n_lines: [int] Runs of every fragment on a line drawn around the fitted one, for the error of the
            line. 0 by default.
        mc_lines: [list] Lines of the Monte Carlo of the solution to draw them from (see
            wmpl.Utils.DragFitDarkFlight.loadMcLines()). None by default, for the covariance of the
            line-of-sight fit.
        end_ht: [float or str] Height where the dark flight ends (km), or 'auto'. None by default.
        mc_runs: [int] Monte Carlo runs of OpenDarkflight for every run, with its statistical wind errors.
            0 by default.
        seed: [int] Seed of the kicks and of the Monte Carlo of OpenDarkflight. 0 by default.
        odf_options: [list] Further OpenDarkflight options as 'Section.option=value'. None by default.
        odf_python: [str] Python of the environment of OpenDarkflight, needed for K from its drag models.
            None by default.
        odf_dir: [str] Directory of OpenDarkflight. None by default.

    Return:
        (setup, ini_paths):
            - setup: [dict] The fragments, their runs and the settings, as saved to fragment_falls.json.
            - ini_paths: [list] Paths of the OpenDarkflight inputs.
    """

    meteorite = resolveMeteorite(meteorite)
    ks = meteorite.get('k')
    if (ks is not None) and (len(ks) != len(masses)):
        raise ValueError("Give one fixed K for every mass")
    if seed < 0:
        raise ValueError("The seed has to be positive, not {:d}".format(seed))

    drag_fit, line, reference, ref_time = _setup(traj, knot_hts, atm_profile)

    # Without ablation OpenDarkflight flies the whole way, from the breakup
    from_breakup = (sigma == 0)

    # The breakups have to be in the air of the profile, where OpenDarkflight gives the flow
    try:
        atm_profile.checkCoverage(min(sep_hts), max(sep_hts))
    except ValueError as e:
        raise ValueError("The breakups are not covered by the atmosphere profile: {:s}".format(str(e)))

    os.makedirs(out_dir, exist_ok=True)

    # Fly every fragment down from its breakup, with the speed of the leading body there
    fragments, flights = [], []
    for i, sep_ht in enumerate(sep_hts):

        len_sep = float(drag_fit.lengthAtHeight(np.array([sep_ht]))[0])
        if np.isnan(len_sep):
            raise ValueError("The breakup at {:.2f} km is outside the drag fit".format(sep_ht/1000))

        sep_state = handoverState(drag_fit, len_sep, 1.0, drag_fit.velocityAtLength(np.array([len_sep]))[0], \
            float(drag_fit.timeAtLength(np.array([len_sep]))[0]))
        ground_sep = groundLength(line, len_sep)

        for j, mass in enumerate(masses):

            frag_meteorite = dict(meteorite, k=None if ks is None else ks[j])
            fragments.append(dict(block=i, sep_ht=float(sep_ht), number=j, mass=float(mass), \
                k_fixed=frag_meteorite['k'], len_sep=len_sep, ground_sep=ground_sep, \
                vel_sep=sep_state['vel'], \
                time_sep=sep_state['time'], zenith=90.0 - sep_state['elev']))
            flights.append(dict(_fragmentFlight(frag_meteorite, mass, sigma), direction='down', \
                len_start=ground_sep, vel_start=sep_state['vel'], vel_end=v_dark, \
                zenith=90.0 - sep_state['elev']))

    results = _runLuminous(out_dir, 'luminous', line, flights, atm_profile, odf_python, odf_dir)

    odf_options = _odfOptions(odf_options, mc_runs)
    rng = np.random.default_rng(seed)
    ini_paths, runs, notes = [], [], []
    for fragment, result in zip(fragments, results):

        # With ablation, a fragment which has not slowed down by the end of the line has no dark flight from
        #   here
        if (not from_breakup) and (not result['ok']):
            notes.append("The fragment {:d} of the breakup at {:.2f} km has not slowed down to {:.1f} km/s "
                "by the end of the drag fit".format(fragment['number'], fragment['sep_ht']/1000, v_dark/1000))
            continue

        # The kick of the breakup, with the leading body there
        ground_sep = fragment['ground_sep']
        k_main = leadingBallisticCoeff(drag_fit, fragment['len_sep'])
        rho_air = float(np.interp(ground_sep, line['len'], line['rho']))
        angle = pmKickAngle(fragment['vel_sep'], rho_air, meteorite['density'], result['k_start']/k_main, \
            pm_c=pm_c)

        fragment.update(ground_dark=result['len_end'], len_dark=inertialLength(line, result['len_end']), \
            vel_dark=result['vel_end'], k_sep=result['k_start'], k_dark=result['k_end'], \
            mass_dark=result['mass_end'], time_dark=fragment['time_sep'] + result['time'], kick_angle=angle)

        # Without ablation the dark flight starts at the breakup, where the kick turns the motion
        if from_breakup:
            ground_start, k_start, vel_start, time_start, mass_start = ground_sep, result['k_start'], \
                fragment['vel_sep'], fragment['time_sep'], fragment['mass']
        else:
            ground_start, k_start, vel_start, time_start, mass_start = result['len_end'], result['k_end'], \
                result['vel_end'], fragment['time_dark'], result['mass_end']

        # The run without a kick nor a line drawn, and the runs which draw a kick, a line, or both
        n_draws = max(n_kicks, n_lines)
        for draw in range(n_draws + 1):

            kick = draw if (0 < draw <= n_kicks) else 0
            line_num = draw if (0 < draw <= n_lines) else 0
            phase = float(rng.uniform(0.0, 2*np.pi)) if kick else 0.0
            line_drawn = drawLine(traj, rng, mc_lines=mc_lines) if line_num else None
            if line_num and (line_drawn is None):
                line_num = 0
                if not kick:
                    continue

            state = trackState(drag_fit, line, ground_start, k_start, vel_start, time=time_start, \
                line_drawn=line_drawn, kick=(angle if kick else 0.0, phase, ground_start - ground_sep))

            name = 'b{:02d}_f{:03d}_r{:03d}'.format(fragment['block'], fragment['number'], draw)
            run_dir = os.path.join(out_dir, 'fragments', name)
            _cleanRunDir(run_dir)
            ini_paths.append(os.path.join(run_dir, name + '_df_params.ini'))
            writeOdfInput(ini_paths[-1], state, ref_time, atm_profile, meteorite['density'], end_ht=end_ht, \
                mc_runs=mc_runs, seed=seed + len(runs), odf_options=odf_options, \
                **_odfMeteorite(dict(meteorite, k=fragment['k_fixed']), state, mass_start))

            runs.append(dict(name=name, block=fragment['block'], number=fragment['number'], draw=draw, \
                kick=kick, line=line_num, phase=phase, state=state))

    # Remove the falls of an earlier run
    for file_name in ('fragment_falls.csv', 'fragment_falls.png'):
        if os.path.isfile(os.path.join(out_dir, file_name)):
            os.remove(os.path.join(out_dir, file_name))

    setup = dict(reference=reference, fragments=fragments, runs=runs, notes=notes, \
        atm_profile=os.path.abspath(atm_profile.path), atm_profile_sha256=atm_profile.sha256, \
        settings=dict(knot_hts=list(map(float, knot_hts)), sep_hts=list(map(float, sep_hts)), \
        masses=list(map(float, masses)), meteorite=meteorite, v_dark=v_dark, n_kicks=n_kicks, pm_c=pm_c, \
        sigma=sigma, start='breakup' if from_breakup else 'v_dark', n_lines=n_lines, \
        lines_from=('monte_carlo' if mc_lines else 'covariance') if n_lines else 'none', end_ht=end_ht, \
        mc_runs=mc_runs, seed=seed, odf_options=odf_options))

    with open(os.path.join(out_dir, 'fragment_falls.json'), 'w') as f:
        json.dump(setup, f, indent=1)


    return setup, ini_paths



def _runFalls(out_dir, run):
    """ Read the nominal and the Monte Carlo falls of a run, None if it has no output. """

    return readOdfRun(os.path.join(out_dir, 'fragments', run['name']))



def _isNominal(row):
    """ Tell whether a fall is the nominal one of its fragment, without a kick, a line drawn or wind
        errors.
    """

    return (row['draw'] == 0) and (row['mc_run'] == 0)



def spreadFalls(rows):
    """ Take the falls of a fragment which hold all its errors: those of the runs with the kicks and the
        lines drawn, and of their Monte Carlo runs if any, so that every fall has every error.

    Arguments:
        rows: [list] Falls of one fragment (see summarizeFragmentFalls()).

    Return:
        [list] The falls of its spread, none if it has only the nominal one.
    """

    drawn = [row for row in rows if row['draw'] > 0] or rows
    with_mc = [row for row in drawn if row['mc_run'] > 0]

    return with_mc or [row for row in drawn if row['draw'] > 0]



def summarizeFragmentFalls(out_dir):
    """ Collect the falls of the fragments prepared by prepareFragmentFalls() after OpenDarkflight has run,
        save them to fragment_falls.csv and plot them to fragment_falls.png, by breakup, and return a summary.

    Arguments:
        out_dir: [str] Output directory.

    Return:
        [str] Summary of the falls, along and across the motion from the last observed point of the leading
            body.
    """

    with open(os.path.join(out_dir, 'fragment_falls.json')) as f:
        setup = json.load(f)

    reference = setup['reference']
    fragments = {(frag['block'], frag['number']): frag for frag in setup['fragments']}

    rows, n_missing = [], 0
    for run in setup['runs']:

        falls = _runFalls(out_dir, run)
        if falls is None:
            n_missing += 1
            continue

        frag = fragments[(run['block'], run['number'])]
        lat_arr = [falls['lat']] + falls['mc_lat']
        lon_arr = [falls['lon']] + falls['mc_lon']
        along, cross = alongCross(lat_arr, lon_arr, reference['lat'], reference['lon'], \
            reference['azim_motion'])
        for k, (lat, lon, a, c) in enumerate(zip(lat_arr, lon_arr, along, cross)):
            rows.append(dict(block=frag['block'], sep_ht_km=frag['sep_ht']/1000, fragment=frag['number'], \
                mass_kg=falls['mass'], \
                k_dark_m2_kg=frag['k_dark'], start_ht_km=run['state']['ht']/1000, draw=run['draw'], \
                kick=run['kick'], line=run['line'], \
                mc_run=k, lat_deg=lat, lon_deg=lon, along_km=float(a), cross_km=float(c)))

    # Save the falls
    if rows:
        with open(os.path.join(out_dir, 'fragment_falls.csv'), 'w') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    # Describe every breakup by its fragments, nominal fall and spread of the others
    out = ["Falls of the fragments, along and across the motion from the last observed point of the leading "
        "body at {:.5f}, {:.5f} (km):".format(reference['lat'], reference['lon'])]
    for i, sep_ht in enumerate(setup['settings']['sep_hts']):

        out.append("Breakup at {:.2f} km:".format(sep_ht/1000))
        for frag in [frag for frag in setup['fragments'] if frag['block'] == i]:

            frag_rows = [row for row in rows if (row['block'] == i) and (row['fragment'] == frag['number'])]
            if not frag_rows:
                continue

            nominal = ([row for row in frag_rows if _isNominal(row)] or frag_rows)[0]
            spread = spreadFalls(frag_rows)
            label = "K {:.3e} m^2/kg".format(frag['k_fixed']) if frag['k_fixed'] is not None \
                else "{:.4g} kg".format(frag['mass'])
            kick = ", kicks of {:.2f} deg".format(np.degrees(frag['kick_angle'])) \
                if setup['settings']['n_kicks'] > 0 else ""
            text = "    {:s}: dark flight from {:.2f} km{:s}, falls along {:+.2f}, across {:+.2f}".format( \
                label, nominal['start_ht_km'], kick, nominal['along_km'], nominal['cross_km'])

            # The range of the falls with the kicks and the errors of the wind
            if spread:
                a_arr = np.array([row['along_km'] for row in spread])
                c_arr = np.array([row['cross_km'] for row in spread])
                text += " (2.5-97.5% of {:d} falls: along {:+.2f} to {:+.2f}, across {:+.2f} to " \
                    "{:+.2f})".format(len(spread), np.percentile(a_arr, 2.5), np.percentile(a_arr, 97.5), \
                    np.percentile(c_arr, 2.5), np.percentile(c_arr, 97.5))
            out.append(text)

    out += setup['notes']
    if n_missing:
        out.append("{:d} runs without an output of OpenDarkflight".format(n_missing))
    if setup['settings'].get('n_lines', 0) > 0:
        out.append("The lines of the runs are drawn from {:s}".format( \
            LINES_FROM[setup['settings']['lines_from']]))
    if setup['settings']['n_kicks'] > 0:
        out.append("The sideways kicks follow Passey & Melosh (1980), Icarus 42, 211-233. Please cite it.")

    # Plot the falls by breakup
    if rows:
        _plotFragmentFalls(out_dir, setup, rows)


    return "\n".join(out)



def _plotFragmentFalls(out_dir, setup, rows):
    """ Plot the falls of the fragments by breakup, along and across the motion. """

    fig, ax = plt.subplots(figsize=(9, 6))
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, sep_ht in enumerate(setup['settings']['sep_hts']):

        colour = colours[i%len(colours)]
        block = [row for row in rows if row['block'] == i]
        spread = [row for row in block if not _isNominal(row)]
        ax.scatter([row['along_km'] for row in spread], [row['cross_km'] for row in spread], s=4, alpha=0.3, \
            color=colour, lw=0)

        # Join the nominal falls by mass
        nominal = sorted([row for row in block if _isNominal(row)], \
            key=lambda row: row['fragment'])
        ax.plot([row['along_km'] for row in nominal], [row['cross_km'] for row in nominal], 'o-', ms=4, \
            color=colour, label="Breakup at {:.2f} km".format(sep_ht/1000))
        for row in nominal:
            ax.annotate("{:.3g}".format(row['mass_kg']), (row['along_km'], row['cross_km']), fontsize=7, \
                textcoords='offset points', xytext=(3, 3))

    ax.set_xlabel("Along the motion from the last observed point (km)")
    ax.set_ylabel("Across the motion, + to the right (km)")
    ax.set_title("Falls of the fragments by breakup (labels: mass in kg)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'fragment_falls.png'), dpi=150)
    plt.close(fig)



def readMeteorites(file_path):
    """ Read the meteorites found from a CSV file.

    Arguments:
        file_path: [str] Path of the file, with the columns name, lat, lon (deg), mass_kg and density
            (kg/m^3), and optionally ground_ht_m, shape, drag_model, SHAPE_OPTIONS, or k_m2_kg for a fixed K.

    Return:
        [list] One dict per meteorite.
    """

    meteorites = []
    with open(file_path) as f:
        for row in csv.DictReader(f):

            row = {key.strip(): value.strip() for key, value in row.items() \
                if isinstance(key, str) and isinstance(value, str)}
            met = dict(name=row['name'], lat=float(row['lat']), lon=float(row['lon']), \
                mass=float(row['mass_kg']), density=float(row['density']))

            for key in ['shape', 'drag_model']:
                if row.get(key):
                    met[key] = row[key]
            for key in SHAPE_OPTIONS:
                if row.get(key):
                    met[key] = float(row[key])
            if row.get('ground_ht_m'):
                met['ground_ht'] = float(row['ground_ht_m'])
            if row.get('k_m2_kg'):
                met['k'] = float(row['k_m2_kg'])

            meteorites.append(met)


    return meteorites



def startHeights(curve, along):
    """ Find the starting heights whose falls are as far along the motion as the given one, by interpolating
        between the starts scanned.

    Arguments:
        curve: [ndarray] Starting heights (m), and the falls along and across the motion (km), one row per
            start.
        along: [float] Distance along the motion of the fall (km).

    Return:
        [list] For every start found, (height (m), distance across of the fall there (km), change of the
            starting height per km along (m/km)).
    """

    curve = np.asarray(curve)
    curve = curve[np.argsort(curve[:, 0])]
    hts, al, cr = curve[:, 0], curve[:, 1], curve[:, 2]

    found = []
    for i in range(len(hts) - 1):

        if ((al[i] - along)*(al[i + 1] - along) > 0) or (al[i] == al[i + 1]):
            continue

        # Skip the second side of a crossing at a scanned start
        if found and (al[i] == along):
            continue

        f = (along - al[i])/(al[i + 1] - al[i])
        found.append((hts[i] + f*(hts[i + 1] - hts[i]), cr[i] + f*(cr[i + 1] - cr[i]), \
            abs((hts[i + 1] - hts[i])/(al[i + 1] - al[i]))))


    return found



def _terrainHeights(meteorites, odf_python, odf_dir, cache_dir):
    """ Look up the terrain heights of the meteorites with the terrain lookup of OpenDarkflight (m). """

    code = "import json, sys; sys.path.insert(0, '.'); " \
        "from opendf.Routines.Terrain import terrainElevations; " \
        "print(json.dumps([float(h) for h in terrainElevations({:s}, {:s}, cache_dir={:s})]))".format( \
        json.dumps([met['lat'] for met in meteorites]), json.dumps([met['lon'] for met in meteorites]), \
        json.dumps(cache_dir))

    result = subprocess.run([odf_python, '-c', code], cwd=odf_dir, stdout=subprocess.PIPE, \
        stderr=subprocess.STDOUT, universal_newlines=True)
    if result.returncode != 0:
        raise RuntimeError("The terrain lookup failed:\n{:s}".format(result.stdout[-2000:]))

    return json.loads(result.stdout.strip().splitlines()[-1])



def _startGround(drag_fit, track, ht):
    """ Length along the ground track (m) where the leading body was at a height. """

    return groundLength(track, float(drag_fit.lengthAtHeight(np.array([ht]))[0]))



def _startState(drag_fit, track, ht, met, v_dark, kick=None, line_drawn=None):
    """ State of the dark flight of a meteorite starting at v_dark on the ground track at a height, or on a
        line drawn around it.
    """

    return trackState(drag_fit, track, _startGround(drag_fit, track, ht), met.get('k') or 1.0, v_dark, \
        line_drawn=line_drawn, kick=kick)



def _breakupState(drag_fit, track, ht, met, kick=None, line_drawn=None):
    """ State of a meteorite breaking off the leading body at a height, on its ground track and with its
        speed and direction there, or on a line drawn around it. The kick turns the motion there.
    """

    len_ground = _startGround(drag_fit, track, ht)
    vel = float(np.interp(len_ground, track['len'], track['vel_parent']))

    return trackState(drag_fit, track, len_ground, met.get('k') or 1.0, vel, line_drawn=line_drawn, kick=kick)



def _scanState(drag_fit, track, ht, met, from_breakup, v_dark, kick=None, line_drawn=None):
    """ State of a scanned height of the inverse: a breakup (see _breakupState()), or a start of the dark
        flight at v_dark (see _startState()).
    """

    if from_breakup:
        return _breakupState(drag_fit, track, ht, met, kick=kick, line_drawn=line_drawn)

    return _startState(drag_fit, track, ht, met, v_dark, kick=kick, line_drawn=line_drawn)



def _writeMeteoriteInput(run_dir, name, state, ref_time, atm_profile, met, mc_runs, seed, odf_options):
    """ Write the OpenDarkflight input of a meteorite from a state, ending at its terrain. """

    _cleanRunDir(run_dir)
    ini_path = os.path.join(run_dir, name + '_df_params.ini')
    writeOdfInput(ini_path, state, ref_time, atm_profile, met['density'], end_ht=met['ground_ht']/1000, \
        mc_runs=mc_runs, seed=seed, odf_options=odf_options, **_odfMeteorite(met, state, met['mass']))

    return ini_path



def _flyBack(work_dir, name, drag_fit, line, atm_profile, met, start_hts, v_dark, sigma, odf_python, odf_dir):
    """ Fly a meteorite back up from starts at v_dark to its breakup (see LuminousFlight.flyUp()). """

    flights = []
    for ht in start_hts:
        state = _startState(drag_fit, line, ht, met, v_dark)
        flights.append(dict(_fragmentFlight(met, met['mass'], sigma), direction='up', \
            len_start=_startGround(drag_fit, line, ht), vel_start=v_dark, zenith=90.0 - state['elev']))

    # Give the breakup its height, inertial length and speed relative to the ground
    results = _runLuminous(work_dir, name, line, flights, atm_profile, odf_python, odf_dir)
    for result in results:
        result['ht_end'] = float(np.interp(result['len_end'], line['len'], line['ht']))
        result['len_end_inertial'] = inertialLength(line, result['len_end'])
        result['vel_end_ground'] = result['vel_end']


    return results



def invertFalls(traj, out_dir, knot_hts, atm_profile, meteorites, odf_python, odf_dir, v_dark=V_DARK, \
    scan_top=None, scan_bottom=None, scan_step=SCAN_STEP, end_ht=None, mc_runs=30, n_kicks=0, pm_c=PM_C, \
    sigmas=(0.0,), n_lines=0, mc_lines=None, seed=0, odf_options=None, n_proc=2):
    """ Find where the meteorites found broke off: without ablation, scan the breakups on the line, flown
        from there with OpenDarkflight, take the one whose fall is as far along the motion as each meteorite,
        and get the spread of the falls there from the errors of the wind, the kicks and the line. With
        ablation, scan the starts of the dark flight on the line at v_dark instead, and fly back up from the
        one found to the breakup, along the ground track. Saves inverse.json, inverse.txt and inverse.png.

    Arguments:
        traj: [Trajectory object] Solved trajectory, with its orbit.
        out_dir: [str] Output directory.
        knot_hts: [list] Heights of the fragmentations of the drag fit (m), e.g. those of the flares.
        atm_profile: [AtmosphereProfile object] Atmosphere profile, with its winds, which has to cover the
            scanned starts.
        meteorites: [list] The meteorites (see readMeteorites()).
        odf_python: [str] Python of the environment of OpenDarkflight.
        odf_dir: [str] Directory of OpenDarkflight.

    Keyword arguments:
        v_dark: [float] Speed relative to the ground where the dark flight starts with ablation (m/s).
            V_DARK by default.
        scan_top: [float] Highest breakup (or start, with ablation) scanned (m). None by default, for
            BREAKUP_MARGIN above the highest fragmentation (the highest fragmentation, with ablation), or the
            first observed point, whichever is lower.
        scan_bottom: [float] Lowest breakup or start scanned (m). None by default, for the last observed
            point.
        scan_step: [float] Step of the heights scanned (m). SCAN_STEP by default.
        end_ht: [float or str] Terrain height of the meteorites without their own ground_ht (km), or 'auto' to
            look it up with OpenDarkflight. None by default, which needs ground_ht for all.
        mc_runs: [int] Monte Carlo runs of OpenDarkflight at the start found, for the errors of the wind. 30
            by default.
        n_kicks: [int] Kicks of Passey and Melosh drawn at the start found. 0 by default.
        pm_c: [float] Constant C of the kicks. PM_C by default.
        sigmas: [list] Ablation coefficients in the flight up to the breakup (s^2/m^2), each giving its own
            breakup. Only 0, without ablation, by default, and then the breakups are scanned.
        n_lines: [int] Lines drawn around the fitted one at the start found, for the error of the line. 0 by
            default.
        mc_lines: [list] Lines of the Monte Carlo of the solution to draw them from (see
            wmpl.Utils.DragFitDarkFlight.loadMcLines()). None by default, for the covariance of the
            line-of-sight fit.
        seed: [int] Seed of the kicks and of the Monte Carlo of OpenDarkflight. 0 by default.
        odf_options: [list] Further OpenDarkflight options as 'Section.option=value'. None by default.
        n_proc: [int] Runs of OpenDarkflight at the same time. 2 by default.

    Return:
        [dict] The results, as saved to inverse.json.
    """

    drag_fit, line, reference, ref_time = _setup(traj, knot_hts, atm_profile)
    os.makedirs(out_dir, exist_ok=True)

    # Without ablation OpenDarkflight flies the meteorites from the breakups scanned; with ablation the starts
    #   of the dark flight are scanned, and the breakups come from flying back up
    sigmas = [float(sigma) for sigma in sigmas]
    from_breakup = all(sigma == 0 for sigma in sigmas)

    # Heights scanned
    ht_first = float(drag_fit.heightAtLength(drag_fit.len_first)[0])
    if (scan_top is None) and from_breakup:
        scan_top = min(max(knot_hts) + BREAKUP_MARGIN, ht_first) if len(knot_hts) else ht_first
    elif scan_top is None:
        scan_top = max(knot_hts) if len(knot_hts) else ht_first
    if scan_bottom is None:
        scan_bottom = float(drag_fit.heightAtLength(drag_fit.len_last)[0])
    scan_hts = np.arange(scan_top, scan_bottom - 1e-6, -scan_step)
    if len(scan_hts) < 2:
        raise ValueError("Scan at least two heights, from {:.2f} down to {:.2f} km".format( \
            scan_top/1000, scan_bottom/1000))

    try:
        atm_profile.checkCoverage(scan_hts.min(), scan_hts.max())
    except ValueError as e:
        raise ValueError("The scanned starts are not covered by the atmosphere profile: {:s}".format(str(e)))

    # Terrain under every meteorite, and the shape and drag model of each
    meteorites = [resolveMeteorite(met) for met in meteorites]
    missing = [met for met in meteorites if met.get('ground_ht') is None]
    if missing and (end_ht == 'auto'):
        for met, ground_ht in zip(missing, _terrainHeights(missing, odf_python, odf_dir, out_dir)):
            met['ground_ht'] = ground_ht
    elif missing and (end_ht is not None):
        for met in missing:
            met['ground_ht'] = 1000*float(end_ht)
    elif missing:
        raise ValueError("Give the terrain height of every meteorite (ground_ht_m), end_ht or 'auto'")
    for met in meteorites:
        if not np.isfinite(met['ground_ht']):
            raise ValueError("No terrain height for the meteorite {:s}".format(met['name']))

    odf_options_mc = _odfOptions(odf_options, mc_runs)
    odf_options = list(odf_options or [])

    # Scan the starts of every meteorite
    ini_paths = []
    for m, met in enumerate(meteorites):
        for ht in scan_hts:
            name = 'm{:02d}_{:06.0f}'.format(m, ht)
            state = _scanState(drag_fit, line, ht, met, from_breakup, v_dark)
            ini_paths.append(_writeMeteoriteInput(os.path.join(out_dir, 'scan', name), name, state, \
                ref_time, atm_profile, met, 0, seed, odf_options))
    failed = runOdf(ini_paths, odf_python, odf_dir, n_proc=n_proc)

    # Find the starts whose falls are as far along as the meteorites
    for m, met in enumerate(meteorites):

        rows = []
        for ht in scan_hts:
            falls = readOdfRun(os.path.join(out_dir, 'scan', 'm{:02d}_{:06.0f}'.format(m, ht)))
            if falls is not None:
                a, c = alongCross([falls['lat']], [falls['lon']], reference['lat'], reference['lon'], \
                    reference['azim_motion'])
                rows.append((float(ht), float(a[0]), float(c[0])))

        a, c = alongCross([met['lat']], [met['lon']], reference['lat'], reference['lon'], \
            reference['azim_motion'])
        met.update(along=float(a[0]), cross=float(c[0]), scan=rows, starts=[dict(ht=s[0], cross_pred=s[1], \
            slope=s[2]) for s in startHeights(rows, float(a[0]))] if len(rows) > 1 else [])

    # With ablation, fly back up from every start found to the breakup, for every ablation coefficient
    for m, met in enumerate(meteorites):
        for k, sigma in enumerate(sigmas):
            if met['starts'] and (not from_breakup):
                backs = _flyBack(out_dir, 'back_m{:02d}_s{:d}'.format(m, k), drag_fit, line, atm_profile, \
                    met, [start['ht'] for start in met['starts']], v_dark, sigma, odf_python, odf_dir)
                for start, back in zip(met['starts'], backs):
                    start.setdefault('breakups', []).append(dict(sigma=sigma, back=back))
                    start['back'] = start['breakups'][0]['back']

    # Run the Monte Carlo of the wind and the kicks at the first start found
    rng = np.random.default_rng(seed)
    ini_paths = []
    for m, met in enumerate(meteorites):

        if not met['starts']:
            continue

        start = met['starts'][0]
        name = 'm{:02d}_mc'.format(m)
        state = _scanState(drag_fit, line, start['ht'], met, from_breakup, v_dark)
        ini_paths.append(_writeMeteoriteInput(os.path.join(out_dir, 'spread', name), name, state, \
            ref_time, atm_profile, met, mc_runs, seed + m, odf_options_mc))

        # The kick at the breakup, with the leading body there: the breakup found, or the one the flight back
        #   up reaches, if any
        back = start.get('back')
        if (n_kicks > 0) and from_breakup:

            len_ground = _startGround(drag_fit, line, start['ht'])
            vel = float(np.interp(len_ground, line['len'], line['vel_parent']))
            k_met = met.get('k')
            if k_met is None:
                flight = dict(_fragmentFlight(met, met['mass'], 0.0), direction='down', \
                    len_start=len_ground, vel_start=vel, vel_end=0.999*vel, zenith=90.0 - state['elev'])
                k_met = _runLuminous(out_dir, 'kick_m{:02d}'.format(m), line, [flight], atm_profile, \
                    odf_python, odf_dir)[0]['k_start']
            k_main = leadingBallisticCoeff(drag_fit, inertialLength(line, len_ground))
            rho_air = float(np.interp(len_ground, line['len'], line['rho']))
            start['kick_angle'] = pmKickAngle(vel, rho_air, met['density'], k_met/k_main, pm_c=pm_c)
            path = 0.0

        elif (n_kicks > 0) and (not back['ok']):
            start['kick_note'] = "no kicks, as no breakup was found"

        elif n_kicks > 0:

            k_main = leadingBallisticCoeff(drag_fit, back['len_end_inertial'])
            rho_air = float(np.interp(back['len_end'], line['len'], line['rho']))
            start['kick_angle'] = pmKickAngle(back['vel_end'], rho_air, met['density'], \
                back['k_end']/k_main, pm_c=pm_c)
            path = _startGround(drag_fit, line, start['ht']) - back['len_end']

        if 'kick_angle' in start:
            for kick in range(1, n_kicks + 1):
                name = 'm{:02d}_k{:03d}'.format(m, kick)
                state = _scanState(drag_fit, line, start['ht'], met, from_breakup, v_dark, \
                    kick=(start['kick_angle'], float(rng.uniform(0.0, 2*np.pi)), path))
                ini_paths.append(_writeMeteoriteInput(os.path.join(out_dir, 'spread', name), name, state, \
                    ref_time, atm_profile, met, 0, seed, odf_options))

        # The lines drawn around the fitted one, from the Monte Carlo of the solution or its covariance
        for k in range(1, n_lines + 1):

            line_drawn = drawLine(traj, rng, mc_lines=mc_lines)
            if line_drawn is None:
                break

            name = 'm{:02d}_l{:03d}'.format(m, k)
            state = _scanState(drag_fit, line, start['ht'], met, from_breakup, v_dark, line_drawn=line_drawn)
            ini_paths.append(_writeMeteoriteInput(os.path.join(out_dir, 'spread', name), name, state, \
                ref_time, atm_profile, met, 0, seed, odf_options))

    failed += runOdf(ini_paths, odf_python, odf_dir, n_proc=n_proc) if ini_paths else []

    # The spread of the falls, and the range of starts and breakups it gives
    for m, met in enumerate(meteorites):

        if not met['starts']:
            continue

        start = met['starts'][0]
        spread = readOdfRun(os.path.join(out_dir, 'spread', 'm{:02d}_mc'.format(m)))
        if spread is None:
            continue

        a0, c0 = alongCross([spread['lat']], [spread['lon']], reference['lat'], reference['lon'], \
            reference['azim_motion'])

        wind_a, wind_c = alongCross(spread['mc_lat'], spread['mc_lon'], reference['lat'], reference['lon'], \
            reference['azim_motion'])
        def spreadFalls(code, n):

            falls_list = [readOdfRun(os.path.join(out_dir, 'spread', 'm{:02d}_{:s}{:03d}'.format(m, code, \
                k))) for k in range(1, n + 1)]
            falls_list = [falls for falls in falls_list if falls is not None]
            along, cross = alongCross([falls['lat'] for falls in falls_list], [falls['lon'] for falls in \
                falls_list], reference['lat'], reference['lon'], reference['azim_motion'])

            # Their spread about the fall of the start found
            if not falls_list:
                return along, cross, (0.0, 0.0)
            return along, cross, (np.sqrt(np.mean((along - a0[0])**2)), np.sqrt(np.mean((cross - c0[0])**2)))

        kick_a, kick_c, sd_kick = spreadFalls('k', n_kicks)
        line_a, line_c, sd_line = spreadFalls('l', n_lines)

        # Add the spreads of the wind, the kicks and the line
        sd_wind = (np.std(wind_a), np.std(wind_c)) if len(wind_a) > 1 else (0.0, 0.0)
        sd_along = np.sqrt(sd_wind[0]**2 + sd_kick[0]**2 + sd_line[0]**2)
        sd_cross = np.sqrt(sd_wind[1]**2 + sd_kick[1]**2 + sd_line[1]**2)

        start.update(sd_wind=list(map(float, sd_wind)), sd_kick=list(map(float, sd_kick)), \
            sd_line=list(map(float, sd_line)), sd_along=float(sd_along), sd_cross=float(sd_cross), \
            n_wind=len(wind_a), n_kick=len(kick_a), n_line=len(line_a), \
            cross_sigma=float((met['cross'] - start['cross_pred'])/sd_cross) if sd_cross > 0 else None, \
            wind_along=wind_a.tolist(), wind_cross=wind_c.tolist(), kick_along=kick_a.tolist(), \
            kick_cross=kick_c.tolist(), line_along=line_a.tolist(), line_cross=line_c.tolist())

        # Fly back up from the ends of the range of starts
        if sd_along > 0:
            dht = sd_along*start['slope']
            start['ht_range'] = [max(start['ht'] - dht, float(scan_hts.min())), \
                min(start['ht'] + dht, float(scan_hts.max()))]
            start['ht_range_clipped'] = bool((start['ht'] - dht < scan_hts.min()) \
                or (start['ht'] + dht > scan_hts.max()))
            for k, breakup in enumerate(start.get('breakups', [])):
                breakup['back_range'] = _flyBack(out_dir, 'range_m{:02d}_s{:d}'.format(m, k), drag_fit, \
                    line, atm_profile, met, start['ht_range'], v_dark, breakup['sigma'], odf_python, odf_dir)
            if start.get('breakups'):
                start['back_range'] = start['breakups'][0]['back_range']

    results = dict(reference=reference, meteorites=meteorites, failed_runs=failed, \
        observed_hts=[float(drag_fit.heightAtLength(drag_fit.len_first)[0]), \
        float(drag_fit.heightAtLength(drag_fit.len_last)[0])], \
        atm_profile=os.path.abspath(atm_profile.path), atm_profile_sha256=atm_profile.sha256, \
        settings=dict(method='breakup' if from_breakup else 'start', knot_hts=list(map(float, knot_hts)), \
        v_dark=v_dark, scan_hts=scan_hts.tolist(), \
        end_ht=end_ht, mc_runs=mc_runs, n_kicks=n_kicks, pm_c=pm_c, sigmas=sigmas, n_lines=n_lines, \
        lines_from=('monte_carlo' if mc_lines else 'covariance') if n_lines > 0 else 'none', seed=seed, \
        odf_options=odf_options))

    with open(os.path.join(out_dir, 'inverse.json'), 'w') as f:
        json.dump(results, f, indent=1)
    with open(os.path.join(out_dir, 'inverse.txt'), 'w') as f:
        f.write(inverseReport(results) + "\n")
    _plotInverse(out_dir, results)


    return results



def inverseReport(results):
    """ Describe the results of invertFalls().

    Arguments:
        results: [dict] Results of invertFalls().

    Return:
        [str] The report.
    """

    settings = results['settings']
    ht_first, ht_last = results['observed_hts']
    from_breakup = (settings.get('method') == 'breakup')
    what = "Breakup" if from_breakup else "Start"

    def breakup(back):
        if back.get('status') == 'faster_at_start':
            return "none, as the leading body is not faster there"
        if not back['ok']:
            return "not met below the top of the drag fit"
        note = " (outside the observed part)" if not (ht_last <= back['ht_end'] <= ht_first) else ""
        return "{:.2f} km at {:.2f} km/s{:s}".format(back['ht_end']/1000, back['vel_end_ground']/1000, note)

    def ablation(sigma):
        if sigma > 0:
            return "with sigma {:.4f} s^2/km^2".format(1e6*sigma)
        return "without ablation"

    if from_breakup:
        out = ["Flown with OpenDarkflight from breakups on the line of the leading body, with its speed "
            "there, without ablation; distances along and across the motion from its last observed point "
            "(km)"]
    else:
        out = ["Dark flight started at {:.1f} km/s on the line of the leading body; distances along and "
            "across the motion from its last observed point (km)".format(settings['v_dark']/1000)]
    for met in results['meteorites']:

        out.append("{:s} ({:.4g} kg): along {:+.2f}, across {:+.2f}".format(met['name'], met['mass'], \
            met['along'], met['cross']))
        if not met['starts']:
            out.append("    No {:s} between {:.2f} and {:.2f} km makes it fall this far along".format( \
                what.lower(), settings['scan_hts'][-1]/1000, settings['scan_hts'][0]/1000))
            continue

        for i, start in enumerate(met['starts']):

            # The start, with its range from the spread of the fall along, and the breakups
            span = ""
            if 'ht_range' in start:
                span = " [{:.2f} to {:.2f}] (1 sigma of the fall along {:.2f} km: wind {:.2f}, kicks " \
                    "{:.2f}, line {:.2f})".format(start['ht_range'][0]/1000, start['ht_range'][1]/1000, \
                    start['sd_along'], start['sd_wind'][0], start['sd_kick'][0], start['sd_line'][0])
            out.append("    {:s}{:s} at {:.2f} km{:s}".format(what, "" if i == 0 else " (also)", \
                start['ht']/1000, span))

            for breakup_sigma in start.get('breakups', []):

                back = breakup_sigma['back']
                mass = ", {:.4g} kg then".format(back['mass_end']) if breakup_sigma['sigma'] > 0 else ""
                span = ""
                if 'back_range' in breakup_sigma:
                    span = " [{:s} to {:s}]".format(breakup(breakup_sigma['back_range'][0]), \
                        breakup(breakup_sigma['back_range'][1]))
                out.append("    Breakup {:s}: {:s}{:s}{:s}".format(ablation(breakup_sigma['sigma']), \
                    breakup(back), mass, span))

            if i > 0:
                continue
            if start.get('ht_range_clipped'):
                out.append("    The range of {:s}s is cut at the scanned heights".format(what.lower()))
            if start.get('kick_note'):
                out.append("    {:s}".format(start['kick_note'][0].upper() + start['kick_note'][1:]))
            if start.get('cross_sigma') is not None:
                out.append("    Across: predicted {:+.2f}, found {:+.2f} km, {:+.1f} sigma".format( \
                    start['cross_pred'], met['cross'], start['cross_sigma']))

    # The meteorites against each other, as the wind is the same for all
    found = [met for met in results['meteorites'] if met['starts']]
    for i in range(len(found)):
        for j in range(i + 1, len(found)):
            a, b = found[i]['starts'][0], found[j]['starts'][0]
            if from_breakup:
                out.append("{:s} - {:s}: breakup {:+.2f} km".format(found[i]['name'], found[j]['name'], \
                    (a['ht'] - b['ht'])/1000))
            else:
                out.append("{:s} - {:s}: start {:+.2f} km, breakup {:+.2f} km".format(found[i]['name'], \
                    found[j]['name'], (a['ht'] - b['ht'])/1000, \
                    (a['back']['ht_end'] - b['back']['ht_end'])/1000))

    if from_breakup:
        out.append("The breakup assumes the meteorite flew alone from the leading body, with its speed "
            "there.")
    else:
        out.append("Only the start of the dark flight follows from the fall; the breakup assumes the "
            "meteorite flew alone from the leading body, and the flight back up keeps to its ground track.")
    if settings.get('n_lines', 0) > 0:
        out.append("The lines are drawn from {:s}".format(LINES_FROM[settings['lines_from']]))
    if results['failed_runs']:
        out.append("{:d} runs of OpenDarkflight failed".format(len(results['failed_runs'])))
    if settings['n_kicks'] > 0:
        out.append("The sideways kicks follow Passey & Melosh (1980), Icarus 42, 211-233. Please cite it.")


    return "\n".join(out)



def _plotInverse(out_dir, results):
    """ Plot the falls of the scanned heights and the meteorites, and the height against the fall. """

    fig, (ax_map, ax_ht) = plt.subplots(1, 2, figsize=(13, 5.5))
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    what = "breakup" if results['settings'].get('method') == 'breakup' else "start of the dark flight"

    for m, met in enumerate(results['meteorites']):

        colour = colours[m%len(colours)]
        scan = np.array(met['scan'])
        if len(scan):
            ax_map.plot(scan[:, 1], scan[:, 2], '-', lw=1.2, color=colour)
            ax_ht.plot(scan[:, 1], scan[:, 0]/1000, '-', lw=1.2, color=colour, label=met['name'])

        if met['starts']:
            start = met['starts'][0]
            ax_map.scatter(start.get('wind_along', []), start.get('wind_cross', []), s=5, alpha=0.4, lw=0, \
                color=colour)
            ax_map.scatter(start.get('kick_along', []), start.get('kick_cross', []), s=5, alpha=0.4, lw=0, \
                marker='s', color=colour)
            ax_map.scatter(start.get('line_along', []), start.get('line_cross', []), s=6, alpha=0.4, lw=0, \
                marker='^', color=colour)
            ax_ht.plot(met['along'], start['ht']/1000, 'o', color=colour)

        ax_map.plot(met['along'], met['cross'], '*', ms=12, color=colour, mec='k', label=met['name'])
        ax_ht.axvline(met['along'], color=colour, lw=0.8, ls='--')

    ax_map.set_xlabel("Along the motion from the last observed point (km)")
    ax_map.set_ylabel("Across the motion, + to the right (km)")
    ax_map.set_title("Falls of the scanned heights, the spread at the {:s} found, and the "
        "meteorites".format(what))
    ax_ht.set_xlabel("Along the motion from the last observed point (km)")
    ax_ht.set_ylabel("{:s} (km)".format(what[0].upper() + what[1:]))
    ax_ht.set_title("{:s} against the fall".format(what[0].upper() + what[1:]))
    for ax in (ax_map, ax_ht):
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'inverse.png'), dpi=150)
    plt.close(fig)



if __name__ == "__main__":


    ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser( \
        description="Falls of the fragments of a fireball with OpenDarkflight: "
        "where the fragments which broke off at given heights fall (--direct), or where the meteorites found "
        "started their dark flight and broke off (--inverse).")

    arg_parser.add_argument('traj_path', type=str, \
        help="Path to the trajectory pickle file, solved with its orbit.")
    arg_parser.add_argument('out_dir', type=str, help="Output directory.")
    arg_parser.add_argument('-f', '--fragheights', nargs='+', type=float, default=[], metavar='HT', \
        help="Heights of the fragmentations of the drag fit (km), e.g. those of the flares.")

    arg_parser.add_argument('--direct', nargs='+', type=float, metavar='HT', \
        help="Direct: heights where the fragments broke off (km).")
    arg_parser.add_argument('--masses', nargs='+', type=float, metavar='KG', \
        help="Direct: masses of the fragments (kg).")
    arg_parser.add_argument('--ks', nargs='+', type=float, metavar='K', \
        help="Direct: fixed K of the fragments instead of a drag model (m^2/kg), one for every mass.")
    arg_parser.add_argument('--density', type=float, default=3500.0, \
        help="Direct: bulk density of the fragments (kg/m^3). 3500 by default.")
    arg_parser.add_argument('--shape', type=str, default='sphere', \
        help="Direct: shape of the fragments, as in OpenDarkflight (sphere, brick, cube, cylinder, cone or "
        "A=<value>). sphere by default.")
    arg_parser.add_argument('--drag_model', type=str, default='dfn', \
        help="Direct: drag model of OpenDarkflight. dfn by default.")
    for name in SHAPE_OPTIONS:
        arg_parser.add_argument('--' + name, type=float, \
            help="Direct: {:s} of the fragments, as in OpenDarkflight.".format(name))
    arg_parser.add_argument('--run', action='store_true', \
        help="Direct: run OpenDarkflight on the inputs and collect the falls.")
    arg_parser.add_argument('--collect', action='store_true', \
        help="Direct: only collect the falls of runs done before.")

    arg_parser.add_argument('--inverse', type=str, metavar='CSV', \
        help="Inverse: CSV file of the meteorites found (see the module documentation).")
    arg_parser.add_argument('--scan_top', type=float, \
        help="Inverse: highest breakup (or start of the dark flight, with ablation) scanned (km). {:.0f} km "
        "above the highest fragmentation (the highest fragmentation itself with ablation), or the first "
        "observed point, by default.".format(BREAKUP_MARGIN/1000))
    arg_parser.add_argument('--scan_bottom', type=float, \
        help="Inverse: lowest breakup or start scanned (km). The last observed point by default.")
    arg_parser.add_argument('--scan_step', type=float, default=SCAN_STEP/1000, \
        help="Inverse: step of the heights scanned (km). {:.2f} by default.".format(SCAN_STEP/1000))

    arg_parser.add_argument('--v_dark', type=float, default=V_DARK/1000, \
        help="Speed relative to the ground where the dark flight starts, with ablation (km/s). {:.1f} by "
        "default.".format(V_DARK/1000))
    arg_parser.add_argument('--kicks', type=int, default=0, \
        help="Sideways kicks of Passey and Melosh (1980) drawn for every fragment or meteorite. 0 by "
        "default.")
    arg_parser.add_argument('--pm_c', type=float, default=PM_C, \
        help="Constant C of the kicks. {:.1f} by default (they fit 0.02 to 1.52).".format(PM_C))
    arg_parser.add_argument('--sigma', nargs='+', type=float, default=[0.0], metavar='SIGMA', \
        help="Ablation coefficient before the dark flight (s^2/km^2), never in the dark flight. 0 by "
        "default, and then OpenDarkflight flies the fragments from the breakup. With ablation the flight "
        "down to --v_dark is along the ground track of the leading body. The inverse takes several, e.g. "
        "0 0.014, and gives the breakup for each.")
    arg_parser.add_argument('--lines', type=int, \
        help="Lines drawn around the fitted one for every fragment (direct) or at the start found (inverse), "
        "from the Monte Carlo of the solution, or the covariance of the line-of-sight fit without it. 0 for "
        "--direct and 30 for --inverse by default.")
    arg_parser.add_argument('--mc_uncertainties', type=str, metavar='PATH', \
        help="Monte Carlo uncertainties of the solution to draw the lines from. By default the "
        "<name>_mc_uncertainties.pickle in the 'Monte Carlo' directory next to the trajectory pickle, if "
        "any.")
    arg_parser.add_argument('--mc_runs', type=int, \
        help="Monte Carlo runs of OpenDarkflight for the errors of the wind. 0 for --direct and 30 for "
        "--inverse by default.")
    arg_parser.add_argument('--end_ht', type=str, \
        help="Height where the dark flight ends (km), or 'auto' for the terrain, by OpenDarkflight.")
    arg_parser.add_argument('--seed', type=int, default=0, help="Seed of the kicks and of the Monte Carlo.")
    arg_parser.add_argument('--odf_option', action='append', default=[], metavar='SECTION.OPTION=VALUE', \
        help="Further OpenDarkflight option, may be given several times.")
    arg_parser.add_argument('--odf_python', type=str, help="Python of the environment of OpenDarkflight.")
    arg_parser.add_argument('--odf_dir', type=str, help="Directory of OpenDarkflight.")
    arg_parser.add_argument('--nproc', type=int, default=2, \
        help="Runs of OpenDarkflight at the same time. 2 by default.")

    # Add the atmosphere profile
    addProfileArguments(arg_parser)

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    if (cml_args.direct is None) == (cml_args.inverse is None):
        arg_parser.error("Give either --direct or --inverse")

    traj = loadPickle(*os.path.split(os.path.abspath(cml_args.traj_path)))
    atm_profile = profileFromArguments(cml_args)
    if atm_profile is None:
        arg_parser.error("The dark flight needs the atmosphere profile (--atm_profile)")

    knot_hts = [1000*ht for ht in cml_args.fragheights]

    # Take the lines from the Monte Carlo of the solution, if it was solved with it
    mc_path = cml_args.mc_uncertainties or findMcUncertainties(os.path.abspath(cml_args.traj_path))
    mc_lines = loadMcLines(mc_path)
    end_ht = None
    if cml_args.end_ht is not None:
        end_ht = cml_args.end_ht if (cml_args.end_ht == 'auto') else float(cml_args.end_ht)

    if cml_args.direct is not None:

        if not cml_args.collect:

            if cml_args.masses is None:
                arg_parser.error("--direct needs the masses of the fragments (--masses)")
            if cml_args.run and ((cml_args.odf_python is None) or (cml_args.odf_dir is None)):
                arg_parser.error("--run runs OpenDarkflight, so give --odf_python and --odf_dir")
            if len(cml_args.sigma) > 1:
                arg_parser.error("--direct takes one ablation coefficient (--sigma)")

            meteorite = dict(density=cml_args.density, shape=cml_args.shape, drag_model=cml_args.drag_model, \
                k=cml_args.ks)
            for name in SHAPE_OPTIONS:
                if getattr(cml_args, name) is not None:
                    meteorite[name] = getattr(cml_args, name)

            setup, ini_paths = prepareFragmentFalls(traj, cml_args.out_dir, knot_hts, atm_profile, \
                [1000*ht for ht in cml_args.direct], cml_args.masses, meteorite, \
                v_dark=1000*cml_args.v_dark, \
                n_kicks=cml_args.kicks, pm_c=cml_args.pm_c, sigma=1e-6*cml_args.sigma[0], \
                n_lines=cml_args.lines or 0, mc_lines=mc_lines, end_ht=end_ht, \
                mc_runs=cml_args.mc_runs or 0, seed=cml_args.seed, odf_options=cml_args.odf_option, \
                odf_python=cml_args.odf_python, odf_dir=cml_args.odf_dir)
            print("Wrote {:d} OpenDarkflight inputs to {:s}".format(len(ini_paths), cml_args.out_dir))
            for note in setup['notes']:
                print(note)

            if cml_args.run:
                failed = runOdf(ini_paths, cml_args.odf_python, cml_args.odf_dir, n_proc=cml_args.nproc)
                if failed:
                    print("{:d} runs failed, see their run.log".format(len(failed)))

        if cml_args.run or cml_args.collect:
            print(summarizeFragmentFalls(cml_args.out_dir))

    else:

        if (cml_args.odf_python is None) or (cml_args.odf_dir is None):
            arg_parser.error("--inverse runs OpenDarkflight, so give --odf_python and --odf_dir")

        results = invertFalls(traj, cml_args.out_dir, knot_hts, atm_profile, \
            readMeteorites(cml_args.inverse), \
            cml_args.odf_python, cml_args.odf_dir, v_dark=1000*cml_args.v_dark, \
            scan_top=None if cml_args.scan_top is None else 1000*cml_args.scan_top, \
            scan_bottom=None if cml_args.scan_bottom is None else 1000*cml_args.scan_bottom, \
            scan_step=1000*cml_args.scan_step, end_ht=end_ht, \
            mc_runs=30 if cml_args.mc_runs is None else cml_args.mc_runs, n_kicks=cml_args.kicks, \
            pm_c=cml_args.pm_c, sigmas=[1e-6*sigma for sigma in cml_args.sigma], \
            n_lines=30 if cml_args.lines is None else cml_args.lines, mc_lines=mc_lines, seed=cml_args.seed, \
            odf_options=cml_args.odf_option, \
            n_proc=cml_args.nproc)
        print(inverseReport(results))
