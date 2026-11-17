""" Hand a drag fit over to OpenDarkflight, with joint clones of the ballistic coefficient and the motion at
the last observed point, and the fragments of the last fragmentation.

The dark flight starts at the last observed point, where the drag fit (wmpl.Utils.DragFit) already continues
as a drag-only flight with a constant ballistic coefficient K. OpenDarkflight takes K directly
(ballistic_coeff, in its pull request #45) and gives the meteorite the mass of that K for its shape and
density, so the motion stays the same across the hand-over whatever the shape. Both take the air density and
the winds from the same atmosphere profile (--atm_profile), which has to cover the starting heights, so the
dark flight starts in the air the drag fit ended in. The trajectory has to be solved with its orbit, which
gives the rotation of the Earth along it. The speed handed over is the speed relative to the ground the drag
of the fit acts on, so that K gives the same deceleration on both sides; the summary checks the speeds of both
0.25 to 1 km below the hand-over. The vertical wind of a profile is used by OpenDarkflight only. WMPL does not
import OpenDarkflight: it writes its inputs and, with --run, runs it from its own environment and reads its
outputs.

The uncertainty of K at the last point is dominated by how spread the last fragmentation was, which the
motion does not fix (DragFit.spreadScan). Every spread which fits as well as the best one (Delta chi2 <= 4) is
taken as equally likely, and every clone is one of these fits with its parameters drawn from their
covariance. K and the speed at the last point of a clone come from the same parameters, so they keep their
correlation: a larger K comes with a lower speed, both of which shorten the dark flight, while drawing them
separately would mix them and narrow the strewn field. The spread only makes K grow below the flares, so its
range is one-sided, which a log-normal K could not give either. The nominal state, from which the falls are
measured, is the fit of the best spread. The jumps, the K of DragFit, are the end of that range with the
smallest K and the longest flight, so most falls are short of theirs, and the summary also gives the median
fall of the clones, the centre of the strewn field.

The direction of motion relative to the ground is derotated with the initial velocity, as in
wmpl.Utils.SampleTrajectoryPositions, as drag keeps it. Every clone also gets its own line, which moves its
starting point across the line and turns its direction together: the line of one of the Monte Carlo
trajectories of the solution if it was solved with the Monte Carlo (the <name>_mc_uncertainties.pickle in its
'Monte Carlo' directory), or else a line drawn from the covariance of its line-of-sight fit. Only the line is
taken from the Monte Carlo, as the drag fit gives the motion along it. Every clone can get a shape drawn from
a list and a density drawn from a normal distribution, and its own OpenDarkflight Monte Carlo (--mc_runs) adds
the statistical wind errors of OpenDarkflight, unless --odf_option sets another wind_error_mode. The falls are
then those of the Monte Carlo runs.

The fragments of the last fragmentation start at its height with the state of the fit there, and a K between
that of the leading body after the fragmentation and --frag_k_max times it, log-uniform: smaller bodies of the
shape and density drawn for each, which the leading body they broke from is taken to share, down to
1/frag_k_max^3 of its mass. OpenDarkflight gives them the sideways kick of Passey and Melosh (1980) for that
ratio of radii in its Monte Carlo runs. It does not ablate them, while they start at the speed of the
fireball, so their K is the one at the breakup. If the fragmentation is above the top of the atmosphere
profile, the fragments are left out.

On the command line, with the heights of the flares, the atmosphere profile and where OpenDarkflight is:

    python -m wmpl.Utils.DragFitDarkFlight traj.pickle out_dir -f 35.0 30.0 --atm_profile profile.csv \\
        --drag_model dfn --shapes sphere brick --mc_runs 5 --fragments 100 --run \\
        --odf_python <environment of OpenDarkflight>/bin/python --odf_dir <OpenDarkflight directory>

In Python, the same three steps:

    mc_lines = loadMcLines(findMcUncertainties('traj.pickle'))
    handover, ini_paths = prepareHandover(traj, out_dir, [35000.0, 30000.0], AtmosphereProfile('profile.csv'),
        mc_lines=mc_lines, drag_model='dfn', mc_runs=5, n_fragments=100)
    runOdf(ini_paths, odf_python, odf_dir)
    print(summarizeFalls(out_dir, handover))

prepareHandover() writes the inputs and handover.json without running anything, so OpenDarkflight can also be
run by hand on them, and loadHandover() and summarizeFalls() collect the falls afterwards. Every
OpenDarkflight run with a Monte Carlo uses all the cores of the computer, so a few runs at the same time are
enough.

    Reference:
        Passey, Q. R., & Melosh, H. J. (1980). Effects of atmospheric breakup on crater field formation.
            Icarus, 42(2), 211-233.
"""

from __future__ import print_function, division, absolute_import

import argparse
import concurrent.futures
import configparser
import csv
import glob
import json
import os
import re
import subprocess

import matplotlib.pyplot as plt
import numpy as np

from wmpl.Trajectory.Trajectory import lineFromFrame
from wmpl.Utils.AtmosphereProfile import addProfileArguments, profileFromArguments
from wmpl.Utils.DragFit import DragFit, spreadScan
from wmpl.Utils.Math import vectNorm
from wmpl.Utils.Pickling import loadPickle
from wmpl.Utils.TrajConversions import cartesian2Geo, derotatedRadiantAltAz, jd2Date


# Largest Delta chi2 of an acceptable spread of the fragmentations, as in DragFit.spreadScan
DELTA_CHI2_MAX = 4.0

# Heights below the hand-over where the speeds of the drag fit and OpenDarkflight are compared (m)
CHECK_DEPTHS = [250.0, 500.0, 1000.0]

# Units of the masses in the names of the OpenDarkflight output files (kg)
MASS_UNITS = {'mg': 1e-6, 'g': 1e-3, 'kg': 1.0}

# Lowest bulk density drawn for a meteorite (kg/m^3)
DENSITY_MIN = 1000.0



def stateAtLength(drag_fit, params, len_val):
    """ Compute the ballistic coefficient, the inertial speed and the time at a length within the observed
        part for the given parameters of a drag fit.

    Arguments:
        drag_fit: [DragFit object] Drag fit of the motion along the trajectory.
        params: [ndarray] Parameters of the fit (see DragFit).
        len_val: [float] Length along the trajectory (m).

    Return:
        (k_val, vel, time):
            - k_val: [float] Ballistic coefficient (m^2/kg).
            - vel: [float] Inertial speed (m/s).
            - time: [float] Time relative to the reference time of the trajectory (s).
    """

    len_obs = drag_fit.len_grid[drag_fit.i_first:drag_fit.i_last + 1]
    vel_obs, time_obs = drag_fit._solveObserved(params)

    k_val = float(np.exp(drag_fit._logK(params, np.array([len_val]))[0]))
    vel = float(np.interp(len_val, len_obs, vel_obs))
    time = float(np.interp(len_val, len_obs, time_obs))


    return k_val, vel, time



def lineMove(traj, eci_point, line):
    """ Compute how a point of the trajectory and the radiant move from the line of the solution to another
        line.

    Arguments:
        traj: [Trajectory object] Solved trajectory.
        eci_point: [ndarray] ECI coordinates of the point (m).
        line: [tuple] ECI coordinates of a point on the other line (m) and its radiant vector.

    Return:
        (move, turn):
            - move: [ndarray] Move of the closest point on the line across it (m).
            - turn: [ndarray] Change of the unit radiant vector.
    """

    point, radiant = np.array(traj.state_vect_mini), vectNorm(traj.radiant_eci_mini)
    point_new, radiant_new = np.array(line[0]), vectNorm(line[1])

    # Take the closest points of both lines to the given point
    closest = point + np.dot(eci_point - point, radiant)*radiant
    closest_new = point_new + np.dot(eci_point - point_new, radiant_new)*radiant_new

    return closest_new - closest, radiant_new - radiant



def handoverState(drag_fit, len_val, k_val, vel, time, line=None, kick=None):
    """ Compute the state of the meteorite at a length along the trajectory for a dark flight.

    Arguments:
        drag_fit: [DragFit object] Drag fit of the motion along the trajectory.
        len_val: [float] Length along the trajectory (m).
        k_val: [float] Ballistic coefficient (m^2/kg).
        vel: [float] Inertial speed (m/s).
        time: [float] Time relative to the reference time of the trajectory (s).

    Keyword arguments:
        line: [tuple] Point (m) and radiant vector of a line drawn around the fitted one (see drawLine()).
            None by default, for the fitted line.
        kick: [tuple] Sideways kick of a breakup, as (angle, phase, path): the direction of motion is turned
            by the angle (radians) towards the phase (radians, 0 to the right of the motion and pi/2 up), and
            the point moves that way by path*sin(angle), for the path flown since the breakup (m). None by
            default.

    Return:
        [dict] Ballistic coefficient (m^2/kg), speed relative to the ground and inertial speed (m/s), time
            (s), latitude and longitude (deg), height above the mean sea level (m), and azimuth (+E of N) and
            elevation of the ground-fixed radiant (deg). The speed relative to the ground is the inertial
            speed less the rotation of the Earth along the trajectory, the speed the drag of the fit acts on
            in still air, so that K gives the same deceleration in the dark flight.
    """

    traj = drag_fit.traj
    len_arr = np.array([len_val])

    eci = drag_fit.eciAtLength(len_arr)[0]
    radiant_eci = drag_fit.radiantAtLength(len_arr)[0]
    ht = float(drag_fit.heightAtLength(len_arr)[0])

    # Move the point and turn the direction with the line drawn
    if line is not None:
        move, turn = lineMove(traj, eci, line)
        eci = eci + move
        radiant_eci = vectNorm(radiant_eci + turn)
        ht += np.dot(move, vectNorm(eci))

    # Turn the direction and move the point by the kick, which drag keeps
    if kick is not None:
        angle, phase, path = kick
        right = vectNorm(np.cross(vectNorm(eci), radiant_eci))
        side = np.cos(phase)*right + np.sin(phase)*vectNorm(np.cross(radiant_eci, right))
        move = path*np.sin(angle)*side
        eci = eci + move
        radiant_eci = vectNorm(np.cos(angle)*radiant_eci - np.sin(angle)*side)
        ht += np.dot(move, vectNorm(eci))

    jd = traj.jdt_ref + time/86400.0
    lat, lon, _ = cartesian2Geo(jd, *eci)

    # Derotate the direction with the initial velocity, as drag keeps the direction of motion above the
    #   ground, and take the speed relative to the ground the drag of the fit acts on
    azim, elev, _ = derotatedRadiantAltAz(traj.v_init*radiant_eci, eci, jd, lat, lon)
    vel_ground = vel - drag_fit.v_rotation


    return dict(k=k_val, vel=float(vel_ground), vel_inertial=vel, time=time, lat=float(np.degrees(lat)), \
        lon=float(np.degrees(lon)), ht=ht, azim=float(np.degrees(azim)), elev=float(np.degrees(elev)))



def findMcUncertainties(traj_pickle_path):
    """ Find the uncertainties the trajectory solver saved with its Monte Carlo, in the 'Monte Carlo'
        directory next to the trajectory pickle.

    Arguments:
        traj_pickle_path: [str] Path to the trajectory pickle, <name>_trajectory.pickle.

    Return:
        [str] Path to <name>_mc_uncertainties.pickle, None if there is none.
    """

    dir_path, file_name = os.path.split(os.path.abspath(traj_pickle_path))
    name = file_name.replace('_trajectory.pickle', '')
    mc_path = os.path.join(dir_path, 'Monte Carlo', name + '_mc_uncertainties.pickle')

    return mc_path if os.path.isfile(mc_path) else None



def loadMcLines(mc_path):
    """ Load the lines of the Monte Carlo trajectories of a solution.

    Arguments:
        mc_path: [str] Path to the <name>_mc_uncertainties.pickle of the solution, e.g. from
            findMcUncertainties(). None gives None.

    Return:
        [list] ECI coordinates of a point on every line (m) and its radiant vector, None if no path is given.
    """

    if mc_path is None:
        return None

    uncertainties = loadPickle(*os.path.split(mc_path))

    lines = []
    for mc_traj in (getattr(uncertainties, 'mc_traj_list', None) or []):

        point = getattr(mc_traj, 'state_vect_mini', None)
        radiant = getattr(mc_traj, 'radiant_eci_mini', None)
        if (point is not None) and (radiant is not None) and np.all(np.isfinite(np.r_[point, radiant])):
            lines.append((np.array(point), np.array(radiant)))


    return lines



def drawLine(traj, rng, mc_lines=None):
    """ Draw a line around the fitted one: one of the given Monte Carlo lines, or else one from the
        covariance of the line-of-sight fit.

    Arguments:
        traj: [Trajectory object] Solved trajectory.
        rng: [Generator object] Random number generator.

    Keyword arguments:
        mc_lines: [list] Lines of the Monte Carlo trajectories (see loadMcLines()). None by default.

    Return:
        [tuple] ECI coordinates of a point on the line (m) and its radiant vector, None if there are no Monte
            Carlo lines and the trajectory has no covariance of its line-of-sight fit.
    """

    if mc_lines:
        return mc_lines[rng.integers(len(mc_lines))]

    covariance = getattr(traj, 'los_covariance', None)
    frame = getattr(traj, 'los_frame', None)
    if (covariance is None) or (frame is None):
        return None

    # Move the line of the solution as the drawn coordinates move the line of the frame, about the middle
    #   of the meteor, as the frame may belong to the solution before its Monte Carlo
    point, radiant = np.array(traj.state_vect_mini), vectNorm(traj.radiant_eci_mini)
    middle = point + np.dot(frame[0] - point, radiant)*radiant
    zero = lineFromFrame(np.zeros(4), frame)
    line = lineFromFrame(rng.multivariate_normal(np.zeros(4), covariance, method='eigh'), frame)


    return middle + line[:3] - zero[:3], radiant + line[3:] - zero[3:]



def handoverClones(traj, knot_hts, n_clones=200, seed=0, geometry=True, mc_lines=None, n_fragments=0, \
    frag_k_max=10.0, **kwargs):
    """ Draw joint clones of the ballistic coefficient and the state at the last observed point, from the fits
        of every acceptable spread of the fragmentations, all equally likely, with their parameters drawn from
        their covariance, and the fragments of the last fragmentation.

    Arguments:
        traj: [Trajectory object] Solved trajectory.
        knot_hts: [list] Heights of the fragmentations (m), e.g. those of the flares.

    Keyword arguments:
        n_clones: [int] Number of clones. 200 by default.
        seed: [int] Seed of the random draws. 0 by default.
        geometry: [bool] Draw a line for every clone and fragment (see drawLine()). True by default.
        mc_lines: [list] Lines of the Monte Carlo trajectories of the solution (see loadMcLines()), drawn
            instead of the covariance of the line-of-sight fit. None by default.
        n_fragments: [int] Number of fragments of the last fragmentation within the observed part. 0 by
            default.
        frag_k_max: [float] Largest K of a fragment over that of the leading body after the last
            fragmentation. 10 by default.
        **kwargs: Keyword arguments for DragFit, e.g. atm_profile.

    The parameters of a clone are drawn from the covariance of its fit, a linear propagation of the noise of
    the points: the bounds of the fit, which keep K from decreasing, are not imposed on the draws, as clipping
    them would bias K, so K may decrease somewhere along a clone, while its K and speed at the last point keep
    the spread and the correlation the fit gives them.

    Return:
        (nominal, clones, fragments, scan):
            - nominal: [dict] State of the fit of the best spread (see handoverState()), the one with the
                smallest Delta chi2, with its spread (m), the speeds relative to
                the ground of its drag-only flight below the last point (m/s) at the heights CHECK_DEPTHS
                below it (m), where the lines were drawn from ('monte_carlo', 'covariance' or 'none'), and
                whether some parameters of the fits are not constrained by the data, along which the clones
                do not spread.
            - clones: [list] States of the clones, each with its spread (m) and number.
            - fragments: [list] States of the fragments, each with its number and its K over that of the
                leading body.
            - scan: [dict] Spreads (m), their Delta chi2 and K at the last point (m^2/kg).
    """

    if seed < 0:
        raise ValueError("The seed has to be zero or positive, as OpenDarkflight takes it, not " \
            "{:d}".format(seed))
    if frag_k_max < 1:
        raise ValueError("frag_k_max is the largest K of a fragment over that of the leading body, so at " \
            "least 1, not {:.3f}".format(frag_k_max))

    rng = np.random.default_rng(seed)

    # Find where the lines are drawn from
    lines_from = 'none'
    if geometry and mc_lines:
        lines_from = 'monte_carlo'
    elif geometry and (drawLine(traj, np.random.default_rng(0)) is not None):
        lines_from = 'covariance'
    geometry = lines_from != 'none'

    # Scan the spreads, and keep the fits of the acceptable ones
    _, spreads, delta_chi2, log_k_last, _, scan_fits = spreadScan(traj, knot_hts, return_fits=True, **kwargs)
    acceptable = np.where(delta_chi2 <= DELTA_CHI2_MAX)[0]
    fits = [scan_fits[i] for i in acceptable]

    # The nominal state is that of the best spread, with the speeds of its flight below the last point
    i_nominal = int(np.argmin(delta_chi2))
    nominal_fit = scan_fits[i_nominal]
    nominal = handoverState(nominal_fit, nominal_fit.len_last, *stateAtLength(nominal_fit, \
        nominal_fit.params, nominal_fit.len_last))
    check_hts = nominal['ht'] - np.array(CHECK_DEPTHS)
    check_vel = nominal_fit.velocityAtLength(nominal_fit.lengthAtHeight(check_hts)) - nominal_fit.v_rotation
    nominal.update(spread=float(spreads[i_nominal]), check_ht=check_hts.tolist(), \
        check_vel=check_vel.tolist(), lines_from=lines_from, \
        unconstrained=any(fit.params_null.size > 0 for fit in fits))

    # Draw every clone from one acceptable spread and the covariance of its fit
    clones = []
    for clone_num in range(1, n_clones + 1):

        i = rng.integers(len(fits))
        params = rng.multivariate_normal(fits[i].params, fits[i].params_cov, method='eigh')
        line = drawLine(traj, rng, mc_lines=mc_lines) if geometry else None

        state = handoverState(fits[i], fits[i].len_last, *stateAtLength(fits[i], params, fits[i].len_last), \
            line=line)
        state.update(clone=clone_num, spread=float(spreads[acceptable[i]]))
        clones.append(state)

    # Draw the fragments of the last fragmentation, from the state of the nominal fit there
    fragments = []
    knots = nominal_fit.pen_weight == 0
    if n_fragments and np.any(knots):

        i_knot = np.argmax(np.where(knots, nominal_fit.ramp_beg, -np.inf))
        len_knot = nominal_fit.ramp_beg[i_knot]
        len_after = min(len_knot + nominal_fit.ramp_len[i_knot], nominal_fit.len_last)
        k_after, _, _ = stateAtLength(nominal_fit, nominal_fit.params, len_after)
        _, vel, time = stateAtLength(nominal_fit, nominal_fit.params, len_knot)

        for frag_num in range(1, n_fragments + 1):

            k_factor = frag_k_max**rng.uniform()
            line = drawLine(traj, rng, mc_lines=mc_lines) if geometry else None

            state = handoverState(nominal_fit, len_knot, k_after*k_factor, vel, time, line=line)
            state.update(fragment=frag_num, k_factor=float(k_factor))
            fragments.append(state)


    scan = dict(spreads=spreads.tolist(), delta_chi2=delta_chi2.tolist(), k_last=np.exp(log_k_last).tolist())

    return nominal, clones, fragments, scan



def drawMeteorites(states, shapes, density, density_sigma=0.0, seed=0):
    """ Give every state a shape drawn from a list and a bulk density drawn from a normal distribution, kept
        above DENSITY_MIN. The leading body a fragment broke from is taken to have the shape and the density
        drawn for the fragment, so that the K of the fragment over the K measured for the leading body is
        the ratio of their radii, whatever the shape and the density.

    Arguments:
        states: [list] States (see handoverState()), changed in place.
        shapes: [list] Shapes to draw from, all equally likely. None or empty for the default of
            OpenDarkflight.
        density: [float] Mean bulk density (kg/m^3).

    Keyword arguments:
        density_sigma: [float] Standard deviation of the bulk density (kg/m^3). 0 by default.
        seed: [int] Seed of the random draws. 0 by default.
    """

    rng = np.random.default_rng([seed, 1])

    for state in states:

        state['shape'] = str(rng.choice(shapes)) if shapes else None
        state['density'] = float(max(rng.normal(density, density_sigma) if density_sigma > 0 else density, \
            DENSITY_MIN))



def writeOdfInput(file_path, state, ref_time, atm_profile, density, shape=None, drag_model=None, \
    end_ht=None, mc_runs=0, seed=0, odf_options=None, mass=None, shape_options=None):
    """ Write an OpenDarkflight input file which starts the dark flight from the given state with its
        ballistic coefficient, or with a given mass.

    Arguments:
        file_path: [str] Path of the input file.
        state: [dict] State at the hand-over (see handoverState()).
        ref_time: [str] Reference time of the trajectory, 'YYYY-MM-DD HH:MM:SS.ffffff'.
        atm_profile: [AtmosphereProfile object] Atmosphere profile, whose file gives the winds and the air.
        density: [float] Bulk density of the meteorite (kg/m^3).

    Keyword arguments:
        shape: [str] Shape of the meteorite for OpenDarkflight. None by default, for its default.
        drag_model: [str] Drag model of OpenDarkflight. None by default, for its default.
        end_ht: [float or str] Height where the dark flight ends (km), or 'auto' for the terrain below the
            falls. None by default, for its default.
        mc_runs: [int] Runs of the Monte Carlo of OpenDarkflight. 0 by default, for none.
        seed: [int] Seed of the Monte Carlo of OpenDarkflight. 0 by default.
        odf_options: [list] Further options as 'Section.option=value'. None by default.
        mass: [float] Mass of the meteorite (kg), run instead of the K of the state. It is given as the
            whole mass range too, as OpenDarkflight runs its default range otherwise. None by default.
        shape_options: [dict] Further options of the meteorite, such as its sphericity, flatness, elongation
            and volume_factor. None by default.
    """

    config = configparser.ConfigParser()
    config.optionxform = str

    config['Atmosphere'] = dict(winds_file=os.path.abspath(atm_profile.path), \
        winds_type=atm_profile.profile_type)

    config['Ejection'] = dict(lat='{:.7f}'.format(state['lat']), lon='{:.7f}'.format(state['lon']), \
        ht='{:.5f}'.format(state['ht']/1000), vel='{:.5f}'.format(state['vel']/1000), \
        az='{:.5f}'.format(state['azim']), alt='{:.5f}'.format(state['elev']), ref_time=ref_time, \
        ref_time_ds='{:.4f}'.format(state['time']))
    if isinstance(end_ht, str):
        config['Ejection']['end_ht'] = end_ht
    elif end_ht is not None:
        config['Ejection']['end_ht'] = '{:.4f}'.format(end_ht)

    # Start from the K of the state, or from the given mass alone
    if mass is None:
        config['Meteorite'] = dict(ballistic_coeff='{:.6e}'.format(state['k']), \
            density='{:.1f}'.format(density))
    else:
        mass_str = '{:.6e}'.format(mass)
        config['Meteorite'] = dict(mass_min=mass_str, mass_max=mass_str, mass_list=mass_str, \
            density='{:.1f}'.format(density))

    if shape is not None:
        config['Meteorite']['shape'] = shape
    if drag_model is not None:
        config['Meteorite']['drag_model'] = drag_model
    for name, value in (shape_options or {}).items():
        config['Meteorite'][name] = str(value)

    config['Monte Carlo'] = dict(mc_enable=str(mc_runs > 0), mc_runs=str(max(mc_runs, 1)), \
        random_seed=str(seed))

    # Add the further options
    for option in (odf_options or []):

        key, value = option.split('=', 1)
        section, name = key.rsplit('.', 1)
        if not config.has_section(section):
            config.add_section(section)
        config[section][name] = value


    with open(file_path, 'w') as f:
        config.write(f)



def _cleanRunDir(run_dir):
    """ Make the directory of an OpenDarkflight run, or empty it of the files of an earlier run, whose outputs
        are named after the mass and would be read again.
    """

    os.makedirs(run_dir, exist_ok=True)
    for file_path in glob.glob(os.path.join(run_dir, '*')):
        if os.path.isfile(file_path):
            os.remove(file_path)



def checkDarkFlightInputs(traj, atm_profile):
    """ Check that a trajectory and an atmosphere profile can be handed over to OpenDarkflight, and raise a
        ValueError if not.

    Arguments:
        traj: [Trajectory object] Solved trajectory.
        atm_profile: [AtmosphereProfile object] Atmosphere profile of the drag fit and the dark flight.
    """

    # Without an orbit the drag fit, like the solver, leaves out the rotation of the Earth, which the dark
    #   flight above the ground needs, unless the stations were fixed and the speeds are already relative to
    #   the ground
    no_orbit = getattr(getattr(traj, 'orbit', None), 'v_init_norot', None) is None
    if no_orbit and getattr(traj, 'orbit_los_success', True):
        raise ValueError("The trajectory has no orbit, so the drag fit leaves out the rotation of the Earth "
            "along it: solve it with the orbit for a dark flight")

    # The dark flight always takes the winds of the profile, so the drag fit has to as well
    if not atm_profile.use_winds:
        raise ValueError("OpenDarkflight takes the winds of the atmosphere profile, so the drag fit has to "
            "take them too for the dark flight to start in the same air")



def prepareHandover(traj, out_dir, knot_hts, atm_profile, n_clones=200, seed=0, geometry=True, \
    mc_lines=None, mc_source=None, n_fragments=0, frag_k_max=10.0, shapes=None, density=3500.0, \
    density_sigma=0.0, drag_model=None, end_ht=None, mc_runs=0, odf_options=None):
    """ Prepare the dark flight of a fireball: draw the clones and the fragments (see handoverClones()) and
        their meteorites (see drawMeteorites()), write an OpenDarkflight input for the nominal state and for
        every clone and fragment under out_dir, and save everything to handover.json. Run OpenDarkflight on
        the inputs with runOdf() and collect the falls with summarizeFalls().

    Arguments:
        traj: [Trajectory object] Solved trajectory.
        out_dir: [str] Output directory.
        knot_hts: [list] Heights of the fragmentations (m), e.g. those of the flares.
        atm_profile: [AtmosphereProfile object] Atmosphere profile of both the drag fit and the dark flight.

    Keyword arguments:
        n_clones: [int] Number of clones. 200 by default.
        seed: [int] Seed of the random draws and of the Monte Carlo of OpenDarkflight. 0 by default.
        geometry: [bool] Draw a line for every clone and fragment (see drawLine()). True by default.
        mc_lines: [list] Lines of the Monte Carlo trajectories of the solution (see loadMcLines()). None by
            default, for the covariance of the line-of-sight fit.
        mc_source: [str] Path of the Monte Carlo uncertainties the lines come from, kept in handover.json.
            None by default.
        n_fragments: [int] Number of fragments of the last fragmentation. 0 by default.
        frag_k_max: [float] Largest K of a fragment over that of the leading body. 10 by default.
        shapes: [list] Shapes of the meteorite for OpenDarkflight, one drawn for every clone and fragment.
            None by default, for its default.
        density: [float] Bulk density of the meteorite (kg/m^3). 3500 by default.
        density_sigma: [float] Standard deviation of the bulk density (kg/m^3). 0 by default.
        drag_model: [str] Drag model of OpenDarkflight. None by default, for its default.
        end_ht: [float] Height where the dark flight ends (km). None by default, for its default.
        mc_runs: [int] Monte Carlo runs of OpenDarkflight for every clone and fragment, with its statistical
            wind errors unless odf_options sets another 'Monte Carlo.wind_error_mode'. 0 by default.
        odf_options: [list] Further options of the OpenDarkflight inputs as 'Section.option=value'. None by
            default.

    Return:
        (handover, ini_paths):
            - handover: [dict] The nominal state, the clones, the fragments, the scan of the spreads and the
                settings, as saved to handover.json.
            - ini_paths: [list] Paths of the OpenDarkflight inputs, the nominal one first.
    """

    checkDarkFlightInputs(traj, atm_profile)

    # Give the Monte Carlo of OpenDarkflight its statistical wind errors, which the options can override
    odf_options = list(odf_options or [])
    if mc_runs > 0:
        odf_options = ['Monte Carlo.wind_error_mode=statistical'] + odf_options

    # Draw the clones and the fragments, and their meteorites
    nominal, clones, fragments, scan = handoverClones(traj, knot_hts, n_clones=n_clones, seed=seed, \
        geometry=geometry, mc_lines=mc_lines, n_fragments=n_fragments, frag_k_max=frag_k_max, \
        atm_profile=atm_profile)
    drawMeteorites([nominal], shapes[:1] if shapes else None, density)
    drawMeteorites(clones + fragments, shapes, density, density_sigma=density_sigma, seed=seed)

    # The profile has to cover the starting heights of the leading body, where OpenDarkflight matches K in its
    #   own air
    start_hts = [state['ht'] for state in [nominal] + clones]
    try:
        atm_profile.checkCoverage(min(start_hts), max(start_hts))
    except ValueError as e:
        raise ValueError("The dark flight would not start in the air of the drag fit: {:s}".format(str(e)))

    # Leave out the fragments if the profile does not cover the height of the fragmentation
    fragments_note = None
    if fragments and (max(state['ht'] for state in fragments) > atm_profile.ht_max):
        fragments_note = "the fragmentation at {:.2f} km is above the top of the atmosphere profile, " \
            "{:.2f} km".format(fragments[0]['ht']/1000, atm_profile.ht_max/1000)
        fragments = []

    # Remove the falls of an earlier hand-over
    for file_name in ('falls.csv', 'falls.png'):
        if os.path.isfile(os.path.join(out_dir, file_name)):
            os.remove(os.path.join(out_dir, file_name))

    ref_time = jd2Date(traj.jdt_ref, dt_obj=True).strftime('%Y-%m-%d %H:%M:%S.%f')

    # Write the input of the nominal state, without a Monte Carlo
    _cleanRunDir(os.path.join(out_dir, 'nominal'))
    ini_paths = [os.path.join(out_dir, 'nominal', 'nominal_df_params.ini')]
    writeOdfInput(ini_paths[0], nominal, ref_time, atm_profile, nominal['density'], shape=nominal['shape'], \
        drag_model=drag_model, end_ht=end_ht, odf_options=odf_options)

    # Write the inputs of the clones and the fragments, which get the sideways kick of their breakup
    for kind, states in (('clone', clones), ('fragment', fragments)):
        for state in states:

            run_dir = os.path.join(out_dir, kind + 's', '{:s}_{:04d}'.format(kind, state[kind]))
            _cleanRunDir(run_dir)
            ini_paths.append(os.path.join(run_dir, '{:s}_{:04d}_df_params.ini'.format(kind, state[kind])))

            run_options = list(odf_options)
            if kind == 'fragment':
                run_options += ['Deviating Fragments.v_perp_model=passey_melosh', \
                    'Deviating Fragments.pm_radius_ratio={:.5f}'.format(state['k_factor'])]

            writeOdfInput(ini_paths[-1], state, ref_time, atm_profile, state['density'], \
                shape=state['shape'], drag_model=drag_model, end_ht=end_ht, mc_runs=mc_runs, \
                seed=seed + n_clones*(kind == 'fragment') + state[kind], odf_options=run_options)

    # Save the hand-over and its settings
    handover = dict(nominal=nominal, clones=clones, fragments=fragments, scan=scan, \
        atm_profile=os.path.abspath(atm_profile.path), atm_profile_sha256=atm_profile.sha256, \
        atm_profile_winds=atm_profile.use_winds, mc_uncertainties=mc_source, fragments_note=fragments_note, \
        n_mc_lines=len(mc_lines or []), settings=dict(knot_hts=list(map(float, knot_hts)), \
        n_clones=n_clones, seed=seed, n_fragments=n_fragments, frag_k_max=frag_k_max, shapes=shapes, \
        density=density, density_sigma=density_sigma, drag_model=drag_model, end_ht=end_ht, mc_runs=mc_runs, \
        odf_options=odf_options))

    with open(os.path.join(out_dir, 'handover.json'), 'w') as f:
        json.dump(handover, f, indent=1)


    return handover, ini_paths



def loadHandover(out_dir):
    """ Load a hand-over saved by prepareHandover().

    Arguments:
        out_dir: [str] Output directory of the hand-over.

    Return:
        [dict] The hand-over (see prepareHandover()).
    """

    with open(os.path.join(out_dir, 'handover.json')) as f:
        return json.load(f)



def handoverReport(handover):
    """ Describe a hand-over prepared by prepareHandover(): the starting point, the nominal state, the
        acceptable spreads, the range of the clones, where the lines come from and the fragments.

    Arguments:
        handover: [dict] The hand-over (see prepareHandover()).

    Return:
        [str] The report.
    """

    nominal, clones, fragments, scan = [handover[key] for key in ('nominal', 'clones', 'fragments', 'scan')]
    settings = handover.get('settings', {})
    k_clones = np.array([state['k'] for state in clones])
    vel_clones = np.array([state['vel'] for state in clones])

    out = []
    out.append("Hand-over at the last point: {:.3f} km, {:.6f}, {:.6f}, radiant azim {:.3f}, elev {:.3f} " \
        "deg".format(nominal['ht']/1000, nominal['lat'], nominal['lon'], nominal['azim'], nominal['elev']))
    out.append("Nominal (best spread, {:.2f} km): K = {:.3e} m^2/kg, speed relative to the ground {:.1f} " \
        "m/s".format(nominal.get('spread', 0.0)/1000, nominal['k'], nominal['vel']))
    out.append("Jumps: K = {:.3e} m^2/kg, Delta chi2 {:.2f}".format(scan['k_last'][0], scan['delta_chi2'][0]))
    acceptable = [spread for spread, d_chi2 in zip(scan['spreads'], scan['delta_chi2']) \
        if d_chi2 <= DELTA_CHI2_MAX]
    out.append("Acceptable spreads: {:s} km".format(", ".join("{:.2f}".format(s/1000) for s in acceptable)))
    if clones:
        out.append("Clones: K {:.3e} to {:.3e} m^2/kg (median {:.3e}), speed {:.0f} to {:.0f} m/s, "
            "correlation of ln K and the speed {:+.2f}".format(np.min(k_clones), np.max(k_clones), \
            np.median(k_clones), np.min(vel_clones), np.max(vel_clones), \
            np.corrcoef(np.log(k_clones), vel_clones)[0, 1]))

    if nominal.get('unconstrained', False):
        out.append("Some parameters of the fits are not constrained by the data, and the clones do not "
            "spread along them")

    # Say where the lines come from
    lines_from = nominal.get('lines_from', 'none')
    if lines_from == 'monte_carlo':
        out.append("Lines drawn from the {:d} Monte Carlo trajectories of {:s}".format( \
            handover['n_mc_lines'], str(handover['mc_uncertainties'])))
    elif lines_from == 'covariance':
        out.append("Lines drawn from the covariance of the line-of-sight fit (no Monte Carlo uncertainties)")
    else:
        out.append("The lines are not drawn: no Monte Carlo uncertainties nor covariance of the "
            "line-of-sight fit, or the geometry was switched off")

    # Describe the fragments
    if fragments:
        out.append("Fragments: {:d} from {:.3f} km, K {:.3e} to {:.3e} m^2/kg".format(len(fragments), \
            fragments[0]['ht']/1000, min(state['k'] for state in fragments), \
            max(state['k'] for state in fragments)))
        if settings.get('mc_runs', 0) == 0:
            out.append("Without Monte Carlo runs of OpenDarkflight the fragments get no sideways kick")
    elif handover.get('fragments_note'):
        out.append("No fragments: {:s}".format(handover['fragments_note']))
    elif settings.get('n_fragments'):
        out.append("No fragments: no fragmentation within the observed part")


    return "\n".join(out)



def runOdf(ini_paths, odf_python, odf_dir, n_proc=2):
    """ Run OpenDarkflight on the given input files, in parallel. Its outputs and a run.log are written next
        to every input file. Every run with a Monte Carlo uses all the cores of the computer, so with a Monte
        Carlo a few runs at the same time are enough.

    Arguments:
        ini_paths: [list] Paths of the input files.
        odf_python: [str] Python interpreter of the environment where OpenDarkflight is installed.
        odf_dir: [str] Directory of OpenDarkflight, where it is run from.

    Keyword arguments:
        n_proc: [int] Number of runs at the same time. 2 by default.

    Return:
        [list] Paths of the input files whose run failed.
    """

    env = dict(os.environ, MPLBACKEND='Agg', OMP_NUM_THREADS='1')

    def _run(ini_path):

        with open(os.path.join(os.path.dirname(ini_path), 'run.log'), 'w') as log:

            try:
                result = subprocess.run([odf_python, '-m', 'opendf.Darkflight', os.path.abspath(ini_path)], \
                    cwd=odf_dir, stdout=log, stderr=subprocess.STDOUT, env=env)

            # A run which cannot start counts as failed
            except OSError as e:
                log.write("Could not run OpenDarkflight: {:s}\n".format(str(e)))
                return ini_path, -1

        return ini_path, result.returncode


    with concurrent.futures.ThreadPoolExecutor(max_workers=n_proc) as executor:
        results = list(executor.map(_run, ini_paths))

    return [ini_path for ini_path, code in results if code != 0]



def _readOdfRows(file_path):
    """ Read the data rows of an OpenDarkflight output file, as lists of floats. """

    rows = []
    with open(file_path) as f:
        for line in f:

            if line.startswith('#') or (not line.strip()):
                continue

            rows.append([float(value) for value in line.split('|')])

    return rows



def readOdfRun(run_dir):
    """ Read the fall points of an OpenDarkflight run.

    Arguments:
        run_dir: [str] Directory of the run.

    Return:
        [dict] Mass (kg), fall latitude and longitude (deg) and the heights (m) and speeds (m/s) along the
            nominal flight, and the fall latitudes and longitudes of the Monte Carlo runs (deg), None if there
            is no nominal output.
    """

    nominal_files = glob.glob(os.path.join(run_dir, 'nominal_*.txt'))
    if not nominal_files:
        return None

    if len(nominal_files) > 1:
        raise ValueError("{:d} nominal outputs in {:s}, of more than one run".format(len(nominal_files), \
            run_dir))

    # Take the mass from the name of the file and the fall from its last row
    match = re.search(r'nominal_([0-9.]+)_(mg|g|kg)\.txt$', nominal_files[0])
    rows = _readOdfRows(nominal_files[0])

    out = dict(mass=float(match.group(1))*MASS_UNITS[match.group(2)], lat=rows[-1][1], lon=rows[-1][2], \
        ht=[1000*row[3] for row in rows], vel=[row[4] for row in rows], mc_lat=[], mc_lon=[])

    # Add the falls of the Monte Carlo runs of the same mass
    mc_file = os.path.join(run_dir, 'mc_' + os.path.basename(nominal_files[0])[len('nominal_'):])
    if os.path.isfile(mc_file):
        for row in _readOdfRows(mc_file):
            out['mc_lat'].append(row[1])
            out['mc_lon'].append(row[2])


    return out



def alongCross(lat, lon, lat_ref, lon_ref, azim_motion):
    """ Compute the distances along and across the direction of motion of points from a reference point, on a
        local flat Earth.

    Arguments:
        lat: [ndarray] Latitudes (deg).
        lon: [ndarray] Longitudes (deg).
        lat_ref: [float] Latitude of the reference point (deg).
        lon_ref: [float] Longitude of the reference point (deg).
        azim_motion: [float] Azimuth of the direction of motion (deg, +E of N).

    Return:
        (along, cross): [tuple of ndarrays] Distances along the motion, positive downrange, and across it,
            positive to the right (km).
    """

    north = np.radians(np.asarray(lat) - lat_ref)*6371.0
    east = np.radians(np.asarray(lon) - lon_ref)*6371.0*np.cos(np.radians(lat_ref))
    azim = np.radians(azim_motion)

    return east*np.sin(azim) + north*np.cos(azim), east*np.cos(azim) - north*np.sin(azim)



def collectFalls(out_dir, kind, states, lat_ref, lon_ref, azim_motion):
    """ Collect the falls of the runs of the clones or the fragments, with their Monte Carlo runs.

    Arguments:
        out_dir: [str] Output directory of the hand-over (see the command line).
        kind: [str] 'clone' or 'fragment', the name of the runs.
        states: [list] States of the runs.
        lat_ref: [float] Latitude of the reference fall (deg).
        lon_ref: [float] Longitude of the reference fall (deg).
        azim_motion: [float] Azimuth of the direction of motion (deg, +E of N).

    Return:
        (falls, n_missing, n_unreadable):
            - falls: [list] One dict per fall, the nominal run of every run with mc_run 0.
            - n_missing: [int] Number of runs without an output.
            - n_unreadable: [int] Number of runs whose directory holds the outputs of several runs.
    """

    falls = []
    n_missing, n_unreadable = 0, 0
    for state in states:

        run_dir = os.path.join(out_dir, kind + 's', '{:s}_{:04d}'.format(kind, state[kind]))
        try:
            run = readOdfRun(run_dir)
        except ValueError:
            n_unreadable += 1
            continue

        if run is None:
            n_missing += 1
            continue

        lat_arr = [run['lat']] + run['mc_lat']
        lon_arr = [run['lon']] + run['mc_lon']
        along, cross = alongCross(lat_arr, lon_arr, lat_ref, lon_ref, azim_motion)
        for j, (lat, lon, a, c) in enumerate(zip(lat_arr, lon_arr, along, cross)):
            falls.append(dict(kind=kind, number=state[kind], mc_run=j, spread_m=state.get('spread', ''), \
                k_factor=state.get('k_factor', ''), k_m2_kg=state['k'], vel_ground_m_s=state['vel'], \
                shape=state.get('shape') or '', density=state.get('density', ''), mass_kg=run['mass'], \
                lat_deg=lat, lon_deg=lon, along_km=float(a), cross_km=float(c)))


    return falls, n_missing, n_unreadable



def perturbedFalls(falls):
    """ Keep the Monte Carlo falls of every run which has them, and the nominal fall of the others.

    Arguments:
        falls: [list] Falls (see collectFalls()).

    Return:
        (kept, n_nominal_only):
            - kept: [list] The falls kept.
            - n_nominal_only: [int] Number of runs without Monte Carlo falls.
    """

    with_mc = set(fall['number'] for fall in falls if fall['mc_run'] > 0)
    kept = [fall for fall in falls if (fall['mc_run'] > 0) or (fall['number'] not in with_mc)]


    return kept, len(set(fall['number'] for fall in falls) - with_mc)



def summarizeFalls(out_dir, handover):
    """ Collect the falls of the OpenDarkflight runs of a hand-over, save them to falls.csv and plot them to
        falls.png, and return a summary, with the speeds of the drag fit and of OpenDarkflight below the
        hand-over, which have to agree. The falls of a run are those of its Monte Carlo runs if it has any,
        and else its nominal fall; the nominal fall of every run, without wind errors nor kick, is kept in
        falls.csv with mc_run 0.

    Arguments:
        out_dir: [str] Output directory of the hand-over.
        handover: [dict] The hand-over (see prepareHandover() and loadHandover()).

    Return:
        [str] Summary of the falls.
    """

    nominal, clones, fragments = handover['nominal'], handover['clones'], handover['fragments']

    nominal_run = readOdfRun(os.path.join(out_dir, 'nominal'))
    if nominal_run is None:
        return "No nominal OpenDarkflight output in {:s}".format(out_dir)

    ref = (nominal_run['lat'], nominal_run['lon'], (nominal['azim'] + 180.0)%360.0)
    falls, n_missing, n_unreadable = collectFalls(out_dir, 'clone', clones, *ref)
    frag_falls, frag_missing, frag_unreadable = collectFalls(out_dir, 'fragment', fragments, *ref)
    if not (falls or frag_falls):
        return "No OpenDarkflight output of the clones nor the fragments in {:s}".format(out_dir)

    # Save the falls
    with open(os.path.join(out_dir, 'falls.csv'), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list((falls + frag_falls)[0].keys()))
        writer.writeheader()
        writer.writerows(falls + frag_falls)

    # Keep the Monte Carlo falls of the runs which have them
    falls, nominal_only = perturbedFalls(falls)
    frag_falls, frag_nominal_only = perturbedFalls(frag_falls)

    column = lambda rows, key: np.array([row[key] for row in rows], dtype=float)

    # Plot the falls on the ground around the nominal one, the clones coloured by K
    fig, ax = plt.subplots(figsize=(8, 6))
    if frag_falls:
        ax.scatter(column(frag_falls, 'cross_km'), column(frag_falls, 'along_km'), s=5, c='0.6', marker='x', \
            label='Fragments of the last fragmentation')
    if falls:
        points = ax.scatter(column(falls, 'cross_km'), column(falls, 'along_km'), \
            c=np.log10(column(falls, 'k_m2_kg')), s=6, cmap='viridis', label='Leading body')
        fig.colorbar(points, ax=ax, label='log10 K of the leading body (m$^2$/kg)')
        ax.plot(np.median(column(falls, 'cross_km')), np.median(column(falls, 'along_km')), 'kX', ms=11, \
            label='Median of the leading body')
    ax.plot(0, 0, 'r*', ms=14, label='Nominal (best spread)')
    ax.set_xlabel('Across the motion, right (km)')
    ax.set_ylabel('Along the motion, downrange (km)')
    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'falls.png'), dpi=120)
    plt.close(fig)

    # Summarize the falls
    pct = [2.5, 16, 50, 84, 97.5]
    pct_text = lambda values, fmt: " / ".join(fmt.format(v) for v in np.percentile(values, pct))
    out = []
    out.append("Nominal fall: {:.6f}, {:.6f} (mass for its K {:.3g} kg)".format(nominal_run['lat'], \
        nominal_run['lon'], nominal_run['mass']))
    vel_odf = np.interp(nominal['check_ht'], nominal_run['ht'][::-1], nominal_run['vel'][::-1])
    checks = ["{:.2f} km {:.0f}/{:.0f} m/s".format(ht/1000, v_fit, v_odf) for ht, v_fit, v_odf in \
        zip(nominal['check_ht'], nominal['check_vel'], vel_odf)]
    out.append("Speed relative to the ground below the hand-over at {:.3f} km, drag fit / OpenDarkflight: "
        "{:s}".format(nominal['ht']/1000, ", ".join(checks)))

    for label, rows, states, n_nominal_only, n_out, n_bad, k_label in ( \
        ("Leading body", falls, clones, nominal_only, n_missing, n_unreadable, "K (m^2/kg)"), \
        ("Fragments of the last fragmentation", frag_falls, fragments, frag_nominal_only, frag_missing, \
            frag_unreadable, None)):

        if not states:
            continue

        out.append("{:s}: {:d} falls of {:d} runs, percentiles {:s} of:".format(label, len(rows), \
            len(states), "/".join("{:g}".format(p) for p in pct)))
        if rows:
            out.append("    along the motion (km): " + pct_text(column(rows, 'along_km'), "{:+.2f}"))
            out.append("    across the motion (km): " + pct_text(column(rows, 'cross_km'), "{:+.2f}"))
            if k_label:
                out.append("    {:s}: ".format(k_label) + pct_text(column(rows, 'k_m2_kg'), "{:.3e}"))
            out.append("    mass for K (kg): " + pct_text(column(rows, 'mass_kg'), "{:.3g}"))
            out.append("    median fall: {:.6f}, {:.6f}".format(np.median(column(rows, 'lat_deg')), \
                np.median(column(rows, 'lon_deg'))))

        # Say which runs are missing or only have their nominal fall
        if n_out or n_bad:
            out.append("    {:d} runs without an output and {:d} with the outputs of several runs are left "
                "out".format(n_out, n_bad))
        if n_nominal_only and (len(rows) > n_nominal_only):
            out.append("    {:d} runs without Monte Carlo falls count with their nominal fall".format( \
                n_nominal_only))

    if frag_falls:
        out.append("The sideways kicks of the fragments follow Passey & Melosh (1980), Icarus 42, 211-233. "
            "Please cite it.")


    return "\n".join(out)



if __name__ == "__main__":


    ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description=""" Hand a drag fit with fragmentations over to
        OpenDarkflight: joint clones of K and the state at the last observed point, and the fragments of the
        last fragmentation, one OpenDarkflight input each, and with --run their falls. """)

    arg_parser.add_argument('traj_pickle_file', type=str, help='Path to the trajectory .pickle file.')

    arg_parser.add_argument('out_dir', type=str, help='Output directory.')

    arg_parser.add_argument('-f', '--fragheights', metavar='HT', type=float, nargs='+', \
        help='Heights of the fragmentations in km, e.g. those of the flares (see wmpl.Utils.DragFit). Needed '
        'unless --collect is given.')

    arg_parser.add_argument('-n', '--clones', metavar='N', type=int, default=200, \
        help='Number of clones. 200 by default.')

    arg_parser.add_argument('--seed', type=int, default=0, help='Seed of the random draws. 0 by default.')

    arg_parser.add_argument('--no_geometry', action='store_true', \
        help='Do not draw a line for every clone and fragment.')

    arg_parser.add_argument('--mc_uncertainties', metavar='PATH', type=str, default=None, \
        help="Monte Carlo uncertainties of the trajectory solution, whose lines are drawn for the clones and "
        "the fragments. By default the <name>_mc_uncertainties.pickle in the 'Monte Carlo' directory next to "
        "the trajectory pickle, if there is one; without it, the lines are drawn from the covariance of the "
        "line-of-sight fit.")

    arg_parser.add_argument('-d', '--dens', metavar='DENS', type=float, default=3500.0, \
        help='Bulk density of the meteorite in kg/m^3 for OpenDarkflight. 3500 by default.')

    arg_parser.add_argument('--dens_sigma', metavar='SIGMA', type=float, default=0.0, \
        help='Standard deviation of the bulk density in kg/m^3, drawn for every clone and fragment. 0 by '
        'default.')

    arg_parser.add_argument('--shapes', metavar='SHAPE', type=str, nargs='+', default=None, \
        help="Shapes of the meteorite for OpenDarkflight, e.g. 'sphere brick cylinder', one drawn for every "
        "clone and fragment. The default of OpenDarkflight if not given.")

    arg_parser.add_argument('--drag_model', type=str, default=None, \
        help="Drag model of OpenDarkflight, e.g. 'ceplecha' or 'dfn'. Its default if not given.")

    arg_parser.add_argument('--end_ht', metavar='HT', type=float, default=None, \
        help='Height where the dark flight ends, in km. The default of OpenDarkflight if not given.')

    arg_parser.add_argument('--fragments', metavar='N', type=int, default=0, \
        help='Number of fragments of the last fragmentation, which start at its height. 0 by default.')

    arg_parser.add_argument('--frag_k_max', metavar='FACTOR', type=float, default=10.0, \
        help='Largest K of a fragment over that of the leading body after the last fragmentation, drawn '
        'log-uniform from 1. 10 by default, for masses down to 1/1000.')

    arg_parser.add_argument('--mc_runs', metavar='N', type=int, default=0, \
        help='Monte Carlo runs of OpenDarkflight for every clone and fragment, with its statistical wind '
        'errors (15%% in speed, 5 deg in direction by default), the sideways kick of the fragments and the '
        'uncertainties given with --odf_option. 0 by default.')

    arg_parser.add_argument('--odf_option', metavar='SECTION.OPTION=VALUE', type=str, action='append', \
        default=[], help="Further option of the OpenDarkflight input, e.g. "
        "'Monte Carlo.wind_speed_sigma=0.15'. Can be repeated.")

    arg_parser.add_argument('--run', action='store_true', help='Run OpenDarkflight and collect the falls.')

    arg_parser.add_argument('--collect', action='store_true', \
        help='Only collect the falls of runs already done in the output directory. The trajectory pickle '
        'is not read then.')

    arg_parser.add_argument('--odf_python', type=str, default='python', \
        help='Python interpreter of the environment where OpenDarkflight is installed.')

    arg_parser.add_argument('--odf_dir', type=str, default='.', help='Directory of OpenDarkflight.')

    arg_parser.add_argument('--nproc', type=int, default=2, \
        help='OpenDarkflight runs at the same time, 2 by default. With --mc_runs every run uses all the '
        'cores.')

    # Add the options of an atmosphere profile
    addProfileArguments(arg_parser)

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    ############################


    # Only collect the falls of a hand-over already run
    if cml_args.collect:

        print(summarizeFalls(cml_args.out_dir, loadHandover(cml_args.out_dir)))
        raise SystemExit()

    if cml_args.fragheights is None:
        raise SystemExit("-f is needed: the heights of the fragmentations")

    atm_profile = profileFromArguments(cml_args)
    if atm_profile is None:
        raise SystemExit("--atm_profile is needed: the dark flight takes its winds and air from the same " \
            "file")


    # Load the trajectory, and the lines of its Monte Carlo trajectories if there are any
    traj = loadPickle(*os.path.split(cml_args.traj_pickle_file))
    mc_path = cml_args.mc_uncertainties or findMcUncertainties(cml_args.traj_pickle_file)
    mc_lines = None if cml_args.no_geometry else loadMcLines(mc_path)

    # Prepare the hand-over and the inputs of OpenDarkflight
    handover, ini_paths = prepareHandover(traj, cml_args.out_dir, 1000*np.array(cml_args.fragheights), \
        atm_profile, n_clones=cml_args.clones, seed=cml_args.seed, geometry=not cml_args.no_geometry, \
        mc_lines=mc_lines, mc_source=mc_path, n_fragments=cml_args.fragments, \
        frag_k_max=cml_args.frag_k_max, \
        shapes=cml_args.shapes, density=cml_args.dens, density_sigma=cml_args.dens_sigma, \
        drag_model=cml_args.drag_model, end_ht=cml_args.end_ht, mc_runs=cml_args.mc_runs, \
        odf_options=cml_args.odf_option)

    print(handoverReport(handover))
    print("OpenDarkflight inputs written to {:s}".format(cml_args.out_dir))

    # Run OpenDarkflight and summarize the falls
    if cml_args.run:

        failed = runOdf(ini_paths, cml_args.odf_python, cml_args.odf_dir, n_proc=cml_args.nproc)
        if failed:
            print("{:d} OpenDarkflight runs failed, see their run.log, e.g. {:s}".format(len(failed), \
                failed[0]))

        print(summarizeFalls(cml_args.out_dir, handover))
