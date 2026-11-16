""" Hand a drag fit over to OpenDarkflight, with joint clones of the ballistic coefficient and the motion at
the last observed point.

The dark flight starts at the last observed point, where the drag fit (wmpl.Utils.DragFit) already continues
as a drag-only flight with a constant ballistic coefficient K. OpenDarkflight takes K directly
(ballistic_coeff, in its pull request #45) and gives the meteorite the mass of that K for its shape and
density, so the motion stays the same across the hand-over whatever the shape. With the same atmosphere
profile in both (--atm_profile), the drag fit ends and the dark flight starts in the same air.

The uncertainty of K at the last point is dominated by how spread the last fragmentation was, which the motion
does not fix (DragFit.spreadScan). Every spread which fits as well as the best one (Delta chi2 <= 4) is taken
as equally likely, and every clone is one of these fits with its parameters drawn from their covariance. K and
the speed at the last point of a clone come from the same parameters, so they keep their correlation: a larger
K comes with a lower speed, both of which shorten the dark flight, while drawing them separately would mix
them and narrow the strewn field. The spread only makes K grow below the flares, so its range is one-sided,
which a log-normal K could not give either.

The direction of motion relative to the ground is derotated with the initial velocity, as in
wmpl.Utils.SampleTrajectoryPositions, and the speed relative to the ground with the speed at the point. Every
clone gets its own OpenDarkflight input, whose own Monte Carlo can then add the winds, the shape and the
density of the meteorite (--mc_runs and --odf_option).
"""

from __future__ import print_function, division, absolute_import

import argparse
import concurrent.futures
import configparser
import glob
import json
import os
import re
import subprocess

import matplotlib.pyplot as plt
import numpy as np

from wmpl.Utils.AtmosphereProfile import addProfileArguments, profileFromArguments
from wmpl.Utils.DragFit import DragFit, spreadScan
from wmpl.Utils.Pickling import loadPickle
from wmpl.Utils.TrajConversions import cartesian2Geo, derotatedRadiantAltAz, jd2Date


# Largest Delta chi2 of an acceptable spread of the fragmentations, as in DragFit.spreadScan
DELTA_CHI2_MAX = 4.0

# Heights below the hand-over where the speeds of the drag fit and OpenDarkflight are compared (m)
CHECK_DEPTHS = [250.0, 500.0, 1000.0]

# Units of the masses in the names of the OpenDarkflight output files (kg)
MASS_UNITS = {'mg': 1e-6, 'g': 1e-3, 'kg': 1.0}



def stateAtLastPoint(drag_fit, params):
    """ Compute the ballistic coefficient, the inertial speed and the time at the last observed point for the
        given parameters of a drag fit.

    Arguments:
        drag_fit: [DragFit object] Drag fit of the motion along the trajectory.
        params: [ndarray] Parameters of the fit (see DragFit).

    Return:
        (k_last, vel_last, time_last):
            - k_last: [float] Ballistic coefficient at the last point (m^2/kg).
            - vel_last: [float] Inertial speed at the last point (m/s).
            - time_last: [float] Time of the last point relative to the reference time of the trajectory (s).
    """

    len_obs = drag_fit.len_grid[drag_fit.i_first:drag_fit.i_last + 1]
    vel, time = drag_fit._solveObserved(params)

    k_last = float(np.exp(drag_fit._logK(params, np.array([drag_fit.len_last]))[0]))
    vel_last = float(np.interp(drag_fit.len_last, len_obs, vel))
    time_last = float(np.interp(drag_fit.len_last, len_obs, time))


    return k_last, vel_last, time_last



def handoverState(drag_fit, k_last, vel_last, time_last):
    """ Compute the state of the meteorite at the last observed point for a dark flight.

    Arguments:
        drag_fit: [DragFit object] Drag fit of the motion along the trajectory.
        k_last: [float] Ballistic coefficient at the last point (m^2/kg).
        vel_last: [float] Inertial speed at the last point (m/s).
        time_last: [float] Time of the last point relative to the reference time of the trajectory (s).

    Return:
        [dict] Ballistic coefficient (m^2/kg), speed relative to the ground and inertial speed (m/s), time
            (s), latitude and longitude (deg), height above the mean sea level (m), and azimuth (+E of N) and
            elevation of the ground-fixed radiant (deg).
    """

    traj = drag_fit.traj
    len_last = np.array([drag_fit.len_last])

    eci = drag_fit.eciAtLength(len_last)[0]
    radiant_eci = drag_fit.radiantAtLength(len_last)[0]
    jd = traj.jdt_ref + time_last/86400.0
    lat, lon, _ = cartesian2Geo(jd, *eci)

    # Derotate the direction with the initial velocity and the speed with the speed at the point
    azim, elev, _ = derotatedRadiantAltAz(traj.v_init*radiant_eci, eci, jd, lat, lon)
    _, _, vel_ground = derotatedRadiantAltAz(vel_last*radiant_eci, eci, jd, lat, lon)

    ht = float(drag_fit.heightAtLength(len_last)[0])

    return dict(k=k_last, vel=float(vel_ground), vel_inertial=vel_last, time=time_last, \
        lat=float(np.degrees(lat)), lon=float(np.degrees(lon)), ht=ht, azim=float(np.degrees(azim)), \
        elev=float(np.degrees(elev)))



def handoverClones(traj, knot_hts, n_clones=200, seed=0, **kwargs):
    """ Draw joint clones of the ballistic coefficient and the state at the last observed point, from the fits
        of every acceptable spread of the fragmentations, all equally likely, with their parameters drawn from
        their covariance.

    Arguments:
        traj: [Trajectory object] Solved trajectory.
        knot_hts: [list] Heights of the fragmentations (m), e.g. those of the flares.

    Keyword arguments:
        n_clones: [int] Number of clones. 200 by default.
        seed: [int] Seed of the random draws. 0 by default.
        **kwargs: Keyword arguments for DragFit, e.g. atm_profile.

    Return:
        (nominal, clones, scan):
            - nominal: [dict] State of the fit with the fragmentations as jumps (see handoverState()), with
                the speeds relative to the ground of its drag-only flight below the last point (m/s) at the
                heights CHECK_DEPTHS below it (m).
            - clones: [list] States of the clones, each with its spread (m) and number.
            - scan: [dict] Spreads (m), their Delta chi2 and K at the last point (m^2/kg).
    """

    # Scan the spreads, and fit every acceptable one
    jump_fit, spreads, delta_chi2, log_k_last, _ = spreadScan(traj, knot_hts, **kwargs)
    acceptable = np.where(delta_chi2 <= DELTA_CHI2_MAX)[0]
    fits = [jump_fit if (spreads[i] == 0) else DragFit(traj, knot_hts=knot_hts, knot_spread=spreads[i], \
        **kwargs) for i in acceptable]

    # The nominal state is that of the fit with jumps, with the speeds of its flight below the last point
    nominal = handoverState(jump_fit, *stateAtLastPoint(jump_fit, jump_fit.params))
    check_hts = nominal['ht'] - np.array(CHECK_DEPTHS)
    nominal['check_ht'] = check_hts.tolist()
    nominal['check_vel'] = jump_fit.groundVelocityAtLength(jump_fit.lengthAtHeight(check_hts)).tolist()

    # Draw every clone from one acceptable spread and the covariance of its fit
    rng = np.random.default_rng(seed)
    clones = []
    for clone_num in range(1, n_clones + 1):

        i = rng.integers(len(fits))
        params = rng.multivariate_normal(fits[i].params, fits[i].params_cov, method='eigh')

        state = handoverState(fits[i], *stateAtLastPoint(fits[i], params))
        state.update(clone=clone_num, spread=float(spreads[acceptable[i]]))
        clones.append(state)


    scan = dict(spreads=spreads.tolist(), delta_chi2=delta_chi2.tolist(), k_last=np.exp(log_k_last).tolist())

    return nominal, clones, scan



def writeOdfInput(file_path, state, ref_time, atm_profile, density, shape=None, drag_model=None, \
    end_ht=None, mc_runs=0, seed=0, odf_options=None):
    """ Write an OpenDarkflight input file which starts the dark flight from the given state with its
        ballistic coefficient.

    Arguments:
        file_path: [str] Path of the input file.
        state: [dict] State at the hand-over (see handoverState()).
        ref_time: [str] Reference time of the trajectory, 'YYYY-MM-DD HH:MM:SS.ffffff'.
        atm_profile: [AtmosphereProfile object] Atmosphere profile, whose file gives the winds and the air.
        density: [float] Bulk density of the meteorite (kg/m^3).

    Keyword arguments:
        shape: [str] Shape of the meteorite for OpenDarkflight. None by default, for its default.
        drag_model: [str] Drag model of OpenDarkflight. None by default, for its default.
        end_ht: [float] Height where the dark flight ends (km). None by default, for its default.
        mc_runs: [int] Runs of the Monte Carlo of OpenDarkflight. 0 by default, for none.
        seed: [int] Seed of the Monte Carlo of OpenDarkflight. 0 by default.
        odf_options: [list] Further options as 'Section.option=value'. None by default.
    """

    config = configparser.ConfigParser()
    config.optionxform = str

    config['Atmosphere'] = dict(winds_file=os.path.abspath(atm_profile.path), \
        winds_type=atm_profile.profile_type)

    config['Ejection'] = dict(lat='{:.7f}'.format(state['lat']), lon='{:.7f}'.format(state['lon']), \
        ht='{:.5f}'.format(state['ht']/1000), vel='{:.5f}'.format(state['vel']/1000), \
        az='{:.5f}'.format(state['azim']), alt='{:.5f}'.format(state['elev']), ref_time=ref_time, \
        ref_time_ds='{:.4f}'.format(state['time']))
    if end_ht is not None:
        config['Ejection']['end_ht'] = '{:.4f}'.format(end_ht)

    config['Meteorite'] = dict(ballistic_coeff='{:.6e}'.format(state['k']), density='{:.1f}'.format(density))
    if shape is not None:
        config['Meteorite']['shape'] = shape
    if drag_model is not None:
        config['Meteorite']['drag_model'] = drag_model

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



def runOdf(ini_paths, odf_python, odf_dir, n_proc=4):
    """ Run OpenDarkflight on the given input files, in parallel. Its outputs and a run.log are written next
        to every input file.

    Arguments:
        ini_paths: [list] Paths of the input files.
        odf_python: [str] Python interpreter of the environment where OpenDarkflight is installed.
        odf_dir: [str] Directory of OpenDarkflight, where it is run from.

    Keyword arguments:
        n_proc: [int] Number of runs at the same time. 4 by default.

    Return:
        [list] Paths of the input files whose run failed.
    """

    env = dict(os.environ, MPLBACKEND='Agg', OMP_NUM_THREADS='1')

    def _run(ini_path):

        with open(os.path.join(os.path.dirname(ini_path), 'run.log'), 'w') as log:
            result = subprocess.run([odf_python, '-m', 'opendf.Darkflight', os.path.abspath(ini_path)], \
                cwd=odf_dir, stdout=log, stderr=subprocess.STDOUT, env=env)

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

    # Take the mass from the name of the file and the fall from its last row
    match = re.search(r'nominal_([0-9.]+)_(mg|g|kg)\.txt$', nominal_files[0])
    rows = _readOdfRows(nominal_files[0])

    out = dict(mass=float(match.group(1))*MASS_UNITS[match.group(2)], lat=rows[-1][1], lon=rows[-1][2], \
        ht=[1000*row[3] for row in rows], vel=[row[4] for row in rows], mc_lat=[], mc_lon=[])

    # Add the falls of the Monte Carlo runs
    for mc_file in glob.glob(os.path.join(run_dir, 'mc_*.txt')):
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



def summarizeFalls(out_dir, nominal, clones):
    """ Collect the falls of the runs, save them to falls.csv and plot them to falls.png, and return a
        summary.

    Arguments:
        out_dir: [str] Output directory of the hand-over (see the command line).
        nominal: [dict] Nominal state (see handoverClones()).
        clones: [list] States of the clones.

    Return:
        [str] Summary of the falls.
    """

    nominal_run = readOdfRun(os.path.join(out_dir, 'nominal'))
    if nominal_run is None:
        return "No nominal OpenDarkflight output in {:s}".format(out_dir)

    azim_motion = (nominal['azim'] + 180.0)%360.0

    # Collect the falls of the clones, with the Monte Carlo runs of every clone
    rows = []
    for state in clones:

        run = readOdfRun(os.path.join(out_dir, 'clones', 'clone_{:04d}'.format(state['clone'])))
        if run is None:
            continue

        lat_arr = [run['lat']] + run['mc_lat']
        lon_arr = [run['lon']] + run['mc_lon']
        along, cross = alongCross(lat_arr, lon_arr, nominal_run['lat'], nominal_run['lon'], azim_motion)
        for j, (lat, lon, a, c) in enumerate(zip(lat_arr, lon_arr, along, cross)):
            rows.append((state['clone'], j, state['spread'], state['k'], state['vel'], run['mass'], lat, \
                lon, a, c))

    if not rows:
        return "No OpenDarkflight output of the clones in {:s}".format(out_dir)

    rows = np.array(rows)
    header = "clone,mc_run,spread_m,k_m2_kg,vel_ground_m_s,mass_kg,lat_deg,lon_deg,along_km,cross_km"
    np.savetxt(os.path.join(out_dir, 'falls.csv'), rows, delimiter=',', header=header, comments='', \
        fmt=['%d', '%d', '%.1f', '%.6e', '%.2f', '%.6g', '%.7f', '%.7f', '%.4f', '%.4f'])

    along, cross, k_arr = rows[:, 8], rows[:, 9], rows[:, 3]

    # Plot the falls on the ground around the nominal one, coloured by K
    fig, ax = plt.subplots(figsize=(8, 6))
    points = ax.scatter(cross, along, c=np.log10(k_arr), s=6, cmap='viridis')
    ax.plot(0, 0, 'r*', ms=14, label='Nominal (jumps)')
    fig.colorbar(points, ax=ax, label='log10 K (m$^2$/kg)')
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
    out = []
    out.append("Nominal fall: {:.6f}, {:.6f} (mass for its K {:.3g} kg)".format(nominal_run['lat'], \
        nominal_run['lon'], nominal_run['mass']))
    vel_odf = np.interp(nominal['check_ht'], nominal_run['ht'][::-1], nominal_run['vel'][::-1])
    checks = ["{:.2f} km {:.0f}/{:.0f} m/s".format(ht/1000, v_fit, v_odf) for ht, v_fit, v_odf in \
        zip(nominal['check_ht'], nominal['check_vel'], vel_odf)]
    out.append("Speed relative to the ground below the hand-over at {:.3f} km, drag fit / OpenDarkflight: "
        "{:s}".format(nominal['ht']/1000, ", ".join(checks)))
    out.append("{:d} falls of {:d} clones, percentiles {:s} of:".format(len(rows), \
        len(np.unique(rows[:, 0])), "/".join("{:g}".format(p) for p in pct)))
    for label, values, fmt in (("along the motion (km)", along, "{:+.2f}"), \
        ("across the motion (km)", cross, "{:+.2f}"), ("K (m^2/kg)", k_arr, "{:.3e}"), \
        ("mass for K (kg)", rows[:, 5], "{:.3g}")):

        out.append("    {:s}: ".format(label) + " / ".join(fmt.format(v) for v in np.percentile(values, pct)))


    return "\n".join(out)




if __name__ == "__main__":


    ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description=""" Hand a drag fit with fragmentations over to
        OpenDarkflight: joint clones of K and the state at the last observed point, one OpenDarkflight input
        each, and with --run their falls. """)

    arg_parser.add_argument('traj_pickle_file', type=str, help='Path to the trajectory .pickle file.')

    arg_parser.add_argument('out_dir', type=str, help='Output directory.')

    arg_parser.add_argument('-f', '--fragheights', metavar='HT', type=float, nargs='+', required=True, \
        help='Heights of the fragmentations in km, e.g. those of the flares (see wmpl.Utils.DragFit).')

    arg_parser.add_argument('-n', '--clones', metavar='N', type=int, default=200, \
        help='Number of clones. 200 by default.')

    arg_parser.add_argument('--seed', type=int, default=0, help='Seed of the random draws. 0 by default.')

    arg_parser.add_argument('-d', '--dens', metavar='DENS', type=float, default=3500.0, \
        help='Bulk density of the meteorite in kg/m^3 for OpenDarkflight. 3500 by default.')

    arg_parser.add_argument('--shape', type=str, default=None, \
        help="Shape of the meteorite for OpenDarkflight, e.g. 'sphere' or 'brick'. Its default if not given.")

    arg_parser.add_argument('--drag_model', type=str, default=None, \
        help="Drag model of OpenDarkflight, e.g. 'ceplecha' or 'dfn'. Its default if not given.")

    arg_parser.add_argument('--end_ht', metavar='HT', type=float, default=None, \
        help='Height where the dark flight ends, in km. The default of OpenDarkflight if not given.')

    arg_parser.add_argument('--mc_runs', metavar='N', type=int, default=0, \
        help='Monte Carlo runs of OpenDarkflight for every clone, for the uncertainties it is given with '
        '--odf_option (e.g. of the wind). 0 by default.')

    arg_parser.add_argument('--odf_option', metavar='SECTION.OPTION=VALUE', type=str, action='append', \
        default=[], help="Further option of the OpenDarkflight input, e.g. "
        "'Monte Carlo.wind_speed_sigma=0.15'. Can be repeated.")

    arg_parser.add_argument('--run', action='store_true', help='Run OpenDarkflight and collect the falls.')

    arg_parser.add_argument('--collect', action='store_true', \
        help='Only collect the falls of runs already done in the output directory.')

    arg_parser.add_argument('--odf_python', type=str, default='python', \
        help='Python interpreter of the environment where OpenDarkflight is installed.')

    arg_parser.add_argument('--odf_dir', type=str, default='.', help='Directory of OpenDarkflight.')

    arg_parser.add_argument('--nproc', type=int, default=4, help='OpenDarkflight runs at the same time.')

    # Add the options of an atmosphere profile
    addProfileArguments(arg_parser)

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    ############################


    atm_profile = profileFromArguments(cml_args)
    if atm_profile is None:
        raise SystemExit("--atm_profile is needed: the dark flight takes its winds and air from the same " \
            "file")

    state_file = os.path.join(cml_args.out_dir, 'handover.json')

    if cml_args.collect:

        with open(state_file) as f:
            saved = json.load(f)

        print(summarizeFalls(cml_args.out_dir, saved['nominal'], saved['clones']))
        raise SystemExit()


    # Load the trajectory and draw the clones
    traj = loadPickle(*os.path.split(cml_args.traj_pickle_file))
    nominal, clones, scan = handoverClones(traj, 1000*np.array(cml_args.fragheights), \
        n_clones=cml_args.clones, seed=cml_args.seed, atm_profile=atm_profile)

    ref_time = jd2Date(traj.jdt_ref, dt_obj=True).strftime('%Y-%m-%d %H:%M:%S.%f')
    odf_kwargs = dict(shape=cml_args.shape, drag_model=cml_args.drag_model, end_ht=cml_args.end_ht, \
        odf_options=cml_args.odf_option)

    # Write the inputs of the nominal state and of the clones
    os.makedirs(os.path.join(cml_args.out_dir, 'nominal'), exist_ok=True)
    ini_paths = [os.path.join(cml_args.out_dir, 'nominal', 'nominal_df_params.ini')]
    writeOdfInput(ini_paths[0], nominal, ref_time, atm_profile, cml_args.dens, **odf_kwargs)

    for state in clones:

        clone_dir = os.path.join(cml_args.out_dir, 'clones', 'clone_{:04d}'.format(state['clone']))
        os.makedirs(clone_dir, exist_ok=True)
        ini_paths.append(os.path.join(clone_dir, 'clone_{:04d}_df_params.ini'.format(state['clone'])))
        writeOdfInput(ini_paths[-1], state, ref_time, atm_profile, cml_args.dens, mc_runs=cml_args.mc_runs, \
            seed=cml_args.seed + state['clone'], **odf_kwargs)

    with open(state_file, 'w') as f:
        json.dump(dict(nominal=nominal, clones=clones, scan=scan, atm_profile=atm_profile.path, \
            atm_profile_sha256=atm_profile.sha256), f, indent=1)

    # Report the hand-over
    k_clones = np.array([state['k'] for state in clones])
    vel_clones = np.array([state['vel'] for state in clones])
    print("Hand-over at the last point: {:.3f} km, {:.6f}, {:.6f}, radiant azim {:.3f}, elev {:.3f} " \
        "deg".format(nominal['ht']/1000, nominal['lat'], nominal['lon'], nominal['azim'], nominal['elev']))
    print("Nominal (jumps): K = {:.3e} m^2/kg, speed relative to the ground {:.1f} m/s".format(nominal['k'], \
        nominal['vel']))
    print("Acceptable spreads: {:s} km".format(", ".join("{:.2f}".format(s/1000) for s, d in \
        zip(scan['spreads'], scan['delta_chi2']) if d <= DELTA_CHI2_MAX)))
    print("Clones: K {:.3e} to {:.3e} m^2/kg (median {:.3e}), speed {:.0f} to {:.0f} m/s, correlation of "
        "ln K and the speed {:+.2f}".format(np.min(k_clones), np.max(k_clones), np.median(k_clones), \
        np.min(vel_clones), np.max(vel_clones), np.corrcoef(np.log(k_clones), vel_clones)[0, 1]))
    print("OpenDarkflight inputs written to {:s}".format(cml_args.out_dir))

    # Run OpenDarkflight and summarize the falls
    if cml_args.run:

        failed = runOdf(ini_paths, cml_args.odf_python, cml_args.odf_dir, n_proc=cml_args.nproc)
        if failed:
            print("{:d} OpenDarkflight runs failed, see their run.log, e.g. {:s}".format(len(failed), \
                failed[0]))

        print(summarizeFalls(cml_args.out_dir, nominal, clones))
