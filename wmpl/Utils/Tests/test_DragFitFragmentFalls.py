""" Tests of the falls of the fragments in wmpl.Utils.DragFitFragmentFalls and of their flight down to the
    dark flight in wmpl.Utils.LuminousFlight: the flight against the exact solution in an exponential
    atmosphere, the flight back up to the breakup, the ablation, the sideways kick, the starting heights of
    the inverse, and the direct and the inverse runs, from the breakup without ablation and from where the
    fragment has slowed down with it, with the outputs of OpenDarkflight written in place of its runs. The
    trajectory is the synthetic fireball of test_DragFitFragmentations.py.

Run with pytest:
    python -m pytest wmpl/Utils/Tests/test_DragFitFragmentFalls.py -v
"""

from __future__ import print_function, division, absolute_import

import configparser
import contextlib
import csv
import io
import os
import sys

import numpy as np
import pytest

from wmpl.Utils.DragFitDarkFlight import alongCross, handoverState, writeOdfInput
from wmpl.Utils.DragFitFragmentFalls import BREAKUP_MARGIN, V_DARK, _breakupState, _startGround, \
    _startState, fragmentLine, groundLength, inertialLength, invertFalls, leadingBallisticCoeff, \
    pmKickAngle, prepareFragmentFalls, readMeteorites, startHeights, summarizeFragmentFalls, trackState
from wmpl.Utils.LuminousFlight import ablatedMass, fixedBallisticCoeff, flyDown, flyUp, runJobs
from wmpl.Utils.Tests.test_DragFitAtmosphereProfile import loadProfile
from wmpl.Utils.Tests.test_DragFitDarkFlight import writeFakeOdfRun
from wmpl.Utils.Tests.test_DragFitFragmentations import HT_FRAG, fitDrag, fragmentedTrajectory


# Elevation (deg), scale height (m), sea level density (kg/m^3) and top of the line of the exact solution
LINE_ELEV = 40.0
SCALE_HT = 7200.0
RHO_SEA = 1.225
HT_TOP = 60000.0



def exponentialLine(k_parent, vel_top, wind_along=0.0):
    """ A straight line down from HT_TOP through an exponential atmosphere, with the exact speed of a leading
        body of a fixed K without wind, in a frame which does not rotate, without gravity.

    Arguments:
        k_parent: [float] Ballistic coefficient of the leading body (m^2/kg).
        vel_top: [float] Its speed at the top (m/s).

    Keyword arguments:
        wind_along: [float] Wind along the motion (m/s). 0 by default.

    Return:
        [dict] The line (see wmpl.Utils.LuminousFlight.lineAt()).
    """

    sin_elev = np.sin(np.radians(LINE_ELEV))
    len_arr = np.linspace(0.0, HT_TOP/sin_elev, 20001)
    ht_arr = HT_TOP - len_arr*sin_elev
    rho_arr = RHO_SEA*np.exp(-ht_arr/SCALE_HT)

    # The integral of the density along the line from the top
    column = SCALE_HT/sin_elev*(rho_arr - rho_arr[0])

    return dict(len=len_arr, ht=ht_arr, rho=rho_arr, air_offset=np.full_like(len_arr, wind_along), \
        vel_parent=vel_top*np.exp(-k_parent*column), g_along=np.zeros_like(len_arr))



def exactLength(line, k_val, len_start, vel_start, vel_end):
    """ Length where a body of fixed K slows down from vel_start to vel_end without wind, exactly. """

    sin_elev = np.sin(np.radians(LINE_ELEV))
    rho_start = float(np.interp(len_start, line['len'], line['rho']))
    rho_end = rho_start + np.log(vel_start/vel_end)*sin_elev/(k_val*SCALE_HT)

    return (HT_TOP + SCALE_HT*np.log(rho_end/RHO_SEA))/sin_elev



def testTheFlightDownFollowsTheExactSolution():
    """ A fragment of fixed K slows down along the line as v = v0*exp(-K*integral of rho), which the flight
        down follows to within one step.
    """

    line = exponentialLine(1e-3, 17000.0)

    for k_val in [3e-3, 3e-2]:

        len_end, vel, time, mass, reached = flyDown(line, fixedBallisticCoeff(k_val, 1.0), 20000.0, 15000.0, \
            V_DARK, 1.0)

        assert reached
        assert vel <= V_DARK
        assert abs(len_end - exactLength(line, k_val, 20000.0, 15000.0, V_DARK)) < 2.5
        assert mass == 1.0
        assert time > (len_end - 20000.0)/15000.0

    # A tail wind lowers the speed through the air, so the fragment slows down further on
    len_still, _, _, _, _ = flyDown(line, fixedBallisticCoeff(3e-3, 1.0), 20000.0, 15000.0, V_DARK, 1.0)
    len_tail, _, _, _, _ = flyDown(exponentialLine(1e-3, 17000.0, wind_along=100.0), \
        fixedBallisticCoeff(3e-3, 1.0), 20000.0, 15000.0, V_DARK, 1.0)
    assert len_tail > len_still + 20.0



def testAFragmentWithTheKOfTheLeadingBodyFliesWithIt(tmp_path):
    """ A fragment with the K of the leading body, from its speed at a point, flies on as the leading body
        does in the drag fit: along its ground track with the speed relative to the ground, through the air of
        the profile with its wind, with the gravity along the line.
    """

    traj, _ = fragmentedTrajectory()
    profile = loadProfile(tmp_path/'profile.csv', wind_speed=30.0, wind_dir=45.0)
    drag_fit = fitDrag(traj, knot_hts=[HT_FRAG], atm_profile=profile)

    # Below the last point, where the drag fit flies on with a constant K
    len_a, len_b = drag_fit.lengthAtHeight(np.array([46000.0, 40000.0]))
    k_a = drag_fit.ballisticCoeffAtLength(np.array([len_a]))[0][0]
    vel_a, vel_b = drag_fit.velocityAtLength(np.array([len_a, len_b]))

    line = fragmentLine(drag_fit)
    ground_a, ground_b = groundLength(line, len_a), groundLength(line, len_b)
    len_end, _, time, _, reached = flyDown(line, fixedBallisticCoeff(k_a, 1.0), ground_a, \
        vel_a - drag_fit.v_rotation, vel_b - drag_fit.v_rotation, 1.0)
    time_ab = np.diff(drag_fit.timeAtLength(np.array([len_a, len_b])))[0]

    assert reached
    assert abs(len_end - ground_b) < 5.0
    assert abs(time/time_ab - 1) < 2e-3



def testTheFlightUpReturnsToTheBreakup():
    """ Flying a fragment down from a breakup to its dark flight and back up from there finds the breakup
        again, with the speed of the leading body and its mass there, also with a wind and ablation.
    """

    line = exponentialLine(1e-3, 17000.0, wind_along=-30.0)
    len_sep = 25000.0
    vel_sep = float(np.interp(len_sep, line['len'], line['vel_parent']))

    for sigma in [0.0, 0.014e-6]:

        k_func = fixedBallisticCoeff(1e-2, 0.05)
        len_dark, vel_dark, time_down, mass_dark, _ = flyDown(line, k_func, len_sep, vel_sep, V_DARK, 0.05, \
            sigma=sigma)
        len_back, vel_back, time_up, mass_back, met = flyUp(line, k_func, len_dark, vel_dark, mass_dark, \
            sigma=sigma)

        assert met
        assert abs(len_back - len_sep) < 4.0
        assert abs(vel_back/vel_sep - 1) < 2e-3
        assert abs(mass_back/0.05 - 1) < 1e-3
        assert abs(time_up/time_down - 1) < 1e-2

    # A fragment slower than the leading body everywhere above never meets it
    _, _, _, _, met = flyUp(line, fixedBallisticCoeff(1e-2, 0.05), 25000.0, 1000.0, 0.05)
    assert not met



def testAFlightFromABadStartGetsNowhere():
    """ A flight up which starts as fast as the leading body has no breakup, and a flight from a value which
        is not finite stops at once, both telling why.
    """

    line = exponentialLine(1e-3, 17000.0)
    k_func = fixedBallisticCoeff(1e-2, 0.05)
    vel_parent = float(np.interp(25000.0, line['len'], line['vel_parent']))

    len_end, _, time, _, met = flyUp(line, k_func, 25000.0, vel_parent + 10.0, 0.05)
    assert (not met) and (len_end == 25000.0) and (time == 0.0)
    _, _, _, _, reached = flyDown(line, k_func, np.nan, 15000.0, V_DARK, 0.05)
    assert not reached

    jobs = dict(line={key: value.tolist() for key, value in line.items()}, profile=None, flights=[ \
        dict(direction='up', len_start=25000.0, vel_start=vel_parent + 10.0, mass=0.05, k=1e-2, \
        zenith=50.0), \
        dict(direction='down', len_start=np.nan, vel_start=15000.0, vel_end=V_DARK, mass=0.05, k=1e-2, \
        zenith=50.0)])
    assert [result['status'] for result in runJobs(jobs)] == ['faster_at_start', 'not_finite']



def testTheMeteoritesAreRead(tmp_path):
    """ The meteorites are read from their CSV file, also from rows which end with a comma. """

    met_path = tmp_path/'meteorites.csv'
    with open(met_path, 'w') as f:
        f.write("name,lat,lon,mass_kg,density,ground_ht_m,sphericity\n")
        f.write("A,-30.1,140.2,0.12,3400,180,0.78,\n")
        f.write("B,-30.2,140.3,0.05,3300,,\n")

    meteorites = readMeteorites(str(met_path))
    assert [met['name'] for met in meteorites] == ['A', 'B']
    assert (meteorites[0]['ground_ht'], meteorites[0]['sphericity']) == (180.0, 0.78)
    assert 'ground_ht' not in meteorites[1]



def testTheAblatedMassFollowsTheSpeed():
    """ An ablating fragment loses its mass as m = m0*exp(sigma/2*(v^2 - v0^2)), whatever its path. """

    line = exponentialLine(1e-3, 17000.0)
    sigma = 0.014e-6

    _, vel, _, mass, _ = flyDown(line, fixedBallisticCoeff(1e-2, 0.05), 20000.0, 15000.0, V_DARK, 0.05, \
        sigma=sigma)

    assert abs(mass/(0.05*np.exp(0.5*sigma*(vel**2 - 15000.0**2))) - 1) < 1e-9
    assert abs(ablatedMass(0.05, 15000.0, vel, sigma)/mass - 1) < 1e-12
    assert mass < 0.05



def testTheKickTurnsTheMotionByThePasseyMeloshAngle():
    """ The kick turns the direction by atan(U/V) of Passey and Melosh (1980), and moves the start sideways by
        the path flown since the breakup times its sine, to the right of the motion for a phase of 0 and up
        for pi/2.
    """

    rho_ratio = 1e-5
    assert abs(np.tan(pmKickAngle(15000.0, 0.035, 3500.0, 8.0)) - np.sqrt(1.5*8.0*rho_ratio)) < 1e-12
    assert pmKickAngle(15000.0, 0.035, 3500.0, 8.0, pm_c=0.0) == 0.0
    assert pmKickAngle(15000.0, 0.035, 3500.0, 0.5) == pmKickAngle(15000.0, 0.035, 3500.0, 1.0)

    traj, _ = fragmentedTrajectory()
    drag_fit = fitDrag(traj, knot_hts=[HT_FRAG])
    len_val = drag_fit.lengthAtHeight(np.array([50000.0]))[0]

    # Above the first point, where the drag fit has no drag, the K of the leading body is that of the first
    k_first = drag_fit.ballisticCoeffAtLength(np.array([drag_fit.len_first]))[0][0]
    assert drag_fit.ballisticCoeffAtLength(np.array([drag_fit.len_first - 3000.0]))[0][0] == 0
    assert leadingBallisticCoeff(drag_fit, drag_fit.len_first - 3000.0) == k_first > 0
    vel = drag_fit.velocityAtLength(np.array([len_val]))[0]
    angle, path = 0.01, 2000.0

    plain = handoverState(drag_fit, len_val, 1e-3, vel, 1.0)
    right = handoverState(drag_fit, len_val, 1e-3, vel, 1.0, kick=(angle, 0.0, path))
    up = handoverState(drag_fit, len_val, 1e-3, vel, 1.0, kick=(angle, np.pi/2, path))

    # Turn of the radiant, from the ground-fixed azimuth and elevation
    def unit(state):
        azim, elev = np.radians(state['azim']), np.radians(state['elev'])
        return np.array([np.cos(elev)*np.sin(azim), np.cos(elev)*np.cos(azim), np.sin(elev)])

    for kicked in (right, up):
        turn = np.arccos(np.clip(np.dot(unit(plain), unit(kicked)), -1.0, 1.0))
        assert abs(turn/angle - 1) < 0.03

    # Moving right turns the radiant to a larger azimuth, and moving up lowers it
    assert right['azim'] > plain['azim']
    assert up['elev'] < plain['elev']
    assert abs((up['ht'] - plain['ht'])/(path*np.sin(angle)*np.cos(np.radians(plain['elev']))) - 1) < 0.05

    along, cross = alongCross([right['lat']], [right['lon']], plain['lat'], plain['lon'], \
        (plain['azim'] + 180.0)%360)
    assert abs(1000*cross[0]/(path*np.sin(angle)) - 1) < 0.05
    assert abs(1000*along[0]) < 1.0



def testTheStartHeightsAreFoundBetweenTheScannedStarts():
    """ The starting heights whose falls are as far along as a meteorite are interpolated between the scanned
        starts, every crossing once, and none outside the scanned falls.
    """

    hts = np.arange(30000.0, 15000.0 - 1, -500.0)
    along = -(hts - 15000.0)/1500.0
    curve = np.column_stack([hts, along, 0.001*hts])

    found = startHeights(curve, -4.3)
    assert len(found) == 1
    assert abs(found[0][0] - (15000.0 + 4.3*1500.0)) < 1e-6
    assert abs(found[0][1] - 0.001*found[0][0]) < 1e-9
    assert abs(found[0][2] - 1500.0) < 1e-9

    assert len(startHeights(curve, along[5])) == 1
    assert startHeights(curve, 5.0) == []

    # A fall reached from two starts gives both
    bent = np.column_stack([hts, -np.abs(hts - 22500.0)/1000.0, np.zeros_like(hts)])
    assert len(startHeights(bent, -3.0)) == 2



def testTheInputRunsTheGivenMass(tmp_path):
    """ An OpenDarkflight input of a given mass runs that mass alone, with its shape options, and no K. """

    traj, _ = fragmentedTrajectory()
    drag_fit = fitDrag(traj, knot_hts=[HT_FRAG])
    state = handoverState(drag_fit, drag_fit.len_last, 1e-3, drag_fit.velocityAtLength( \
        np.array([drag_fit.len_last]))[0], 1.0)

    ini_path = str(tmp_path/'run.ini')
    writeOdfInput(ini_path, state, '2020-01-01 00:00:00.000000', loadProfile(tmp_path/'profile.csv'), \
        3500.0, shape='sphere', drag_model='dfn', end_ht='auto', mass=0.12, \
        shape_options=dict(sphericity=0.78))

    config = configparser.ConfigParser()
    config.read(ini_path)

    meteorite = config['Meteorite']
    assert float(meteorite['mass_min']) == float(meteorite['mass_max']) == 0.12
    assert float(meteorite['mass_list']) == 0.12
    assert 'ballistic_coeff' not in meteorite
    assert float(meteorite['sphericity']) == 0.78
    assert config['Ejection']['end_ht'] == 'auto'



def fragmentKs(drag_fit, sep_ht, factors):
    """ Fixed K of fragments as factors of the K of the leading body at a breakup. """

    len_sep = drag_fit.lengthAtHeight(np.array([sep_ht]))[0]

    return [factor*drag_fit.ballisticCoeffAtLength(np.array([len_sep]))[0][0] for factor in factors]



def unitVector(azim, elev):
    """ Unit vector of an azimuth and an elevation (deg) in the local east, north, up frame. """

    azim, elev = np.radians(azim), np.radians(elev)

    return np.array([np.cos(elev)*np.sin(azim), np.cos(elev)*np.cos(azim), np.sin(elev)])



def testWithoutAblationTheFragmentsFlyFromTheBreakup(tmp_path):
    """ Without ablation, OpenDarkflight flies every fragment from the breakup, with the place, the speed,
        the direction and the time of the leading body there and the K of the fragment; the kick turns the
        motion there by the angle of Passey and Melosh without moving it.
    """

    traj, _ = fragmentedTrajectory()
    drag_fit = fitDrag(traj, knot_hts=[HT_FRAG])
    profile = loadProfile(tmp_path/'profile.csv', wind_speed=20.0, wind_dir=90.0)
    ks = fragmentKs(drag_fit, 60000.0, [5.0, 50.0])

    with contextlib.redirect_stdout(io.StringIO()):
        setup, ini_paths = prepareFragmentFalls(traj, str(tmp_path/'direct'), [HT_FRAG], profile, \
            [60000.0, 50000.0], [0.5, 0.005], dict(density=3500.0, shape='sphere', k=ks), n_kicks=2, \
            n_lines=1, mc_runs=3, seed=5)

    assert setup['settings']['start'] == 'breakup'
    assert len(ini_paths) == len(setup['runs']) == 4*3

    profile_fit = fitDrag(traj, knot_hts=[HT_FRAG], atm_profile=profile)
    states = {(run['block'], run['number'], run['draw']): run['state'] for run in setup['runs']}
    for frag in setup['fragments']:

        len_sep = np.array([frag['len_sep']])
        sep = handoverState(profile_fit, frag['len_sep'], 1.0, profile_fit.velocityAtLength(len_sep)[0], \
            float(profile_fit.timeAtLength(len_sep)[0]))

        # The run without a kick nor a line starts at the breakup, as the leading body
        plain = states[(frag['block'], frag['number'], 0)]
        for key in ('lat', 'lon', 'azim', 'elev'):
            assert abs(plain[key] - sep[key]) < 1e-9
        assert abs(plain['ht'] - frag['sep_ht']) < 1.0
        assert abs(plain['vel'] - frag['vel_sep']) < 1e-6
        assert plain['time'] == frag['time_sep']

        # The kick turns the motion by its angle at the breakup, without moving the start. It turns the inertial
        #   motion, so the motion relative to the ground, which the state gives, turns by as much within the
        #   ratio of the speed of the rotation of the Earth to that of the fragment, here ~2 %
        kicked = states[(frag['block'], frag['number'], 2)]
        assert abs(kicked['lat'] - plain['lat']) < 1e-9
        assert abs(kicked['lon'] - plain['lon']) < 1e-9
        turn = np.arccos(np.clip(np.dot(unitVector(kicked['azim'], kicked['elev']), \
            unitVector(plain['azim'], plain['elev'])), -1.0, 1.0))
        assert abs(turn/frag['kick_angle'] - 1) < 0.03

        # The line drawn moves it
        lined = states[(frag['block'], frag['number'], 1)]
        assert abs(lined['lat'] - plain['lat']) + abs(lined['lon'] - plain['lon']) > 1e-7

    # OpenDarkflight gets the K of the fragment at the breakup
    config = configparser.ConfigParser()
    config.read(ini_paths[0])
    assert abs(float(config['Meteorite']['ballistic_coeff'])/setup['fragments'][0]['k_sep'] - 1) < 1e-5
    assert abs(setup['fragments'][0]['k_sep']/ks[0] - 1) < 1e-9



def testWithAblationTheFragmentFallsStartWhereTheyHaveSlowedDown(tmp_path):
    """ With ablation, every fragment starts its dark flight at v_dark where it has slowed down to it, the
        smaller one sooner after the breakup, with a run without a kick and its kicks, and the falls are
        collected by breakup.
    """

    traj, _ = fragmentedTrajectory()
    drag_fit = fitDrag(traj, knot_hts=[HT_FRAG])
    profile = loadProfile(tmp_path/'profile.csv', wind_speed=20.0, wind_dir=90.0)
    ks = fragmentKs(drag_fit, 60000.0, [5.0, 50.0])
    out_dir = str(tmp_path/'direct')

    with contextlib.redirect_stdout(io.StringIO()):
        setup, ini_paths = prepareFragmentFalls(traj, out_dir, [HT_FRAG], profile, [60000.0, 50000.0], \
            [0.5, 0.005], dict(density=3500.0, shape='sphere', k=ks), n_kicks=2, n_lines=1, mc_runs=3, \
            seed=5, sigma=0.014e-6)

    # On the ground track the air moves with the wind of the profile along the motion only
    profile_fit = fitDrag(traj, knot_hts=[HT_FRAG], atm_profile=profile)
    wind_along = np.array(fragmentLine(profile_fit)['air_offset'])
    assert np.max(np.abs(wind_along)) <= 20.0
    assert np.ptp(wind_along) < 1e-6
    assert np.max(np.abs(wind_along)) > 1.0

    fragments = setup['fragments']
    assert len(fragments) == 4
    assert len(ini_paths) == len(setup['runs']) == 4*3
    assert setup['settings']['lines_from'] == 'covariance'
    assert setup['settings']['start'] == 'v_dark'

    # The dark flight starts at v_dark below the breakup, sooner for the larger K
    for frag in fragments:
        assert frag['vel_dark'] <= V_DARK
        assert frag['vel_dark'] > V_DARK - 20.0
        assert frag['len_dark'] > frag['len_sep']
    paths = [frag['len_dark'] - frag['len_sep'] for frag in fragments]
    assert paths[1] < paths[0]

    # The kick is that of Passey and Melosh for the ratio of the radii of the leading body and the fragment
    line = fragmentLine(profile_fit)
    for frag in fragments:
        k_main = leadingBallisticCoeff(profile_fit, frag['len_sep'])
        rho_air = np.interp(frag['ground_sep'], line['len'], line['rho'])
        assert abs(frag['kick_angle'] - pmKickAngle(frag['vel_sep'], rho_air, 3500.0, frag['k_sep']/k_main)) \
            < 1e-12
        assert frag['k_sep'] > k_main

    # The run without a kick nor a line starts on the ground track, the kicks turn it and the lines move it
    states = {(run['block'], run['number'], run['draw']): run for run in setup['runs']}
    plain = trackState(profile_fit, line, fragments[0]['ground_dark'], 1.0, fragments[0]['vel_dark'])
    assert abs(states[(0, 0, 0)]['state']['azim'] - plain['azim']) < 1e-9
    assert abs(states[(0, 0, 0)]['state']['lat'] - plain['lat']) < 1e-9
    assert abs(states[(0, 0, 0)]['state']['vel'] - fragments[0]['vel_dark']) < 1e-6
    assert states[(0, 0, 0)]['state']['time'] == fragments[0]['time_dark']
    assert (states[(0, 0, 1)]['kick'], states[(0, 0, 1)]['line']) == (1, 1)
    assert (states[(0, 0, 2)]['kick'], states[(0, 0, 2)]['line']) == (2, 0)
    assert abs(states[(0, 0, 2)]['state']['azim'] - plain['azim']) > 1e-6
    lined = states[(0, 0, 1)]['state']
    assert abs(lined['lat'] - plain['lat']) + abs(lined['lon'] - plain['lon']) > 1e-7

    # The inverse starts the same dark flight at the same place
    start = _startState(profile_fit, line, plain['ht'], dict(k=1.0), plain['vel'])
    assert abs(start['lat'] - plain['lat']) < 1e-7
    assert abs(start['lon'] - plain['lon']) < 1e-7

    # Without a shape, the input of OpenDarkflight gets the sphere and the drag model of the flight down
    config = configparser.ConfigParser()
    config.read(ini_paths[0])
    assert (config['Meteorite']['shape'], config['Meteorite']['drag_model']) == ('sphere', 'dfn')

    # The inputs start from the state with the K of the fragment there, and the statistical wind errors
    config = configparser.ConfigParser()
    config.read(ini_paths[0])
    assert abs(float(config['Meteorite']['ballistic_coeff'])/fragments[0]['k_dark'] - 1) < 1e-5
    assert config['Monte Carlo']['wind_error_mode'] == 'statistical'

    # Collect falls written in place of the runs
    for ini_path in ini_paths:
        writeFakeOdfRun(os.path.dirname(ini_path), -30.0, 140.0, \
            mc_falls=[(-30.001, 140.0), (-30.0, 140.001)])

    summary = summarizeFragmentFalls(out_dir)
    with open(os.path.join(out_dir, 'fragment_falls.csv')) as f:
        rows = list(csv.DictReader(f))

    assert len(rows) == 12*3
    assert "of 4 falls" in summary
    assert set(row['sep_ht_km'] for row in rows) == {'60.0', '50.0'}
    assert os.path.isfile(os.path.join(out_dir, 'fragment_falls.png'))
    assert "Breakup at 60.00 km" in summary
    assert "kicks of" in summary
    assert "Passey & Melosh" in summary



def fakeRunOdf(ini_paths, odf_python, odf_dir, n_proc=2):
    """ Write a fall for every OpenDarkflight input, straight below its start, 0.001 deg apart for the Monte
        Carlo runs, in place of running it.
    """

    for ini_path in ini_paths:

        config = configparser.ConfigParser()
        config.read(ini_path)
        lat, lon = float(config['Ejection']['lat']), float(config['Ejection']['lon'])

        mc_falls = []
        if config['Monte Carlo'].getboolean('mc_enable'):
            mc_falls = [(lat + 0.001, lon), (lat - 0.001, lon), (lat, lon + 0.001), (lat, lon - 0.001)]

        writeFakeOdfRun(os.path.dirname(ini_path), lat, lon, mc_falls=mc_falls)

    return []



def testTheInverseFindsTheStartAndTheBreakup(tmp_path, monkeypatch):
    """ A meteorite which fell where the dark flight starting at a given height falls is traced back to that
        start, and from there up to the breakup where it had the speed of the leading body, from which it
        flies down to the same start again.
    """

    # Run the fake OpenDarkflight in the module as loaded
    monkeypatch.setattr(sys.modules['wmpl.Utils.DragFitFragmentFalls'], 'runOdf', fakeRunOdf)

    traj, _ = fragmentedTrajectory()
    profile = loadProfile(tmp_path/'profile.csv')
    drag_fit = fitDrag(traj, knot_hts=[HT_FRAG], atm_profile=profile)
    k_met = fragmentKs(drag_fit, 55000.0, [20.0])[0]

    # Put the meteorite below the start of its dark flight at 52.3 km
    ht_true = 52300.0
    len_true = drag_fit.lengthAtHeight(np.array([ht_true]))[0]
    start = handoverState(drag_fit, len_true, k_met, V_DARK + drag_fit.v_rotation, \
        float(drag_fit.timeAtLength(np.array([len_true]))[0]))
    met_path = tmp_path/'meteorites.csv'
    with open(met_path, 'w') as f:
        f.write("name,lat,lon,mass_kg,density,ground_ht_m,k_m2_kg\n")
        f.write("A,{:.9f},{:.9f},0.05,3500,0,{:.6e}\n".format(start['lat'], start['lon'], k_met))

    with contextlib.redirect_stdout(io.StringIO()):
        results = invertFalls(traj, str(tmp_path/'inverse'), [HT_FRAG], profile, \
            readMeteorites(str(met_path)), None, None, scan_top=56000.0, scan_bottom=48000.0, mc_runs=4, \
            sigmas=[0.0, 0.014e-6], n_lines=3, seed=3)

    met = results['meteorites'][0]
    assert len(met['starts']) == 1
    found = met['starts'][0]
    assert abs(found['ht'] - ht_true) < 20.0

    # The breakup is above the start, where the meteorite had the speed of the leading body
    back = found['back']
    assert back['ok']
    assert back['ht_end'] > ht_true
    vel_parent = drag_fit.velocityAtLength(np.array([back['len_end_inertial']]))[0] - drag_fit.v_rotation
    assert abs(back['vel_end']/vel_parent - 1) < 1e-3

    line = fragmentLine(drag_fit)
    len_down, _, _, _, _ = flyDown(line, fixedBallisticCoeff(k_met, 0.05), back['len_end'], back['vel_end'], \
        V_DARK, 0.05)
    assert abs(len_down - groundLength(line, drag_fit.lengthAtHeight(np.array([found['ht']]))[0])) < 10.0

    # Ablating on the way down, the meteorite was heavier at the breakup and broke off higher up
    ablated = found['breakups'][1]['back']
    assert found['breakups'][0]['back'] == back
    assert ablated['ht_end'] > back['ht_end']
    assert ablated['mass_end'] > 0.05

    # The spread of the falls gives a range of starts and breakups around them, with the lines drawn
    assert found['n_line'] == 3
    assert found['sd_line'][0] > 0
    assert abs(found['sd_along']**2 - found['sd_wind'][0]**2 - found['sd_line'][0]**2) < 1e-9
    assert found['sd_along'] > 0
    assert abs(found['ht_range'][1] - found['ht'] - found['sd_along']*found['slope']) < 1e-6
    assert abs(found['ht'] - found['ht_range'][0] - found['sd_along']*found['slope']) < 1e-6
    assert found['back_range'][0]['ht_end'] < back['ht_end'] < found['back_range'][1]['ht_end']
    assert os.path.isfile(str(tmp_path/'inverse'/'inverse.png'))
    with open(str(tmp_path/'inverse'/'inverse.txt')) as f:
        report = f.read()
    assert "A (0.05 kg)" in report
    assert "Breakup without ablation" in report
    assert "Breakup with sigma 0.0140 s^2/km^2" in report

    # The terrain of every meteorite is needed
    with pytest.raises(ValueError, match="terrain"):
        invertFalls(traj, str(tmp_path/'inverse2'), [HT_FRAG], profile, [dict(name='B', lat=0.0, lon=0.0, \
            mass=0.05, density=3500.0, k=k_met)], None, None)



def testWithoutAblationTheInverseFindsTheBreakup(tmp_path, monkeypatch):
    """ Without ablation, a meteorite which fell where the flight from a breakup at a given height falls is
        traced back to that breakup, scanned from above the highest fragmentation, with the kick of Passey
        and Melosh for its K and that of the leading body there, and the spread of the falls around it.
    """

    # Run the fake OpenDarkflight in the module as loaded
    monkeypatch.setattr(sys.modules['wmpl.Utils.DragFitFragmentFalls'], 'runOdf', fakeRunOdf)

    traj, _ = fragmentedTrajectory()
    profile = loadProfile(tmp_path/'profile.csv')
    drag_fit = fitDrag(traj, knot_hts=[HT_FRAG], atm_profile=profile)
    line = fragmentLine(drag_fit)
    k_met = fragmentKs(drag_fit, 55000.0, [20.0])[0]

    # Put the meteorite below the breakup at 52.3 km
    ht_true = 52300.0
    breakup = _breakupState(drag_fit, line, ht_true, dict(k=k_met))
    met_path = tmp_path/'meteorites.csv'
    with open(met_path, 'w') as f:
        f.write("name,lat,lon,mass_kg,density,ground_ht_m,k_m2_kg\n")
        f.write("A,{:.9f},{:.9f},0.05,3500,0,{:.6e}\n".format(breakup['lat'], breakup['lon'], k_met))

    with contextlib.redirect_stdout(io.StringIO()):
        results = invertFalls(traj, str(tmp_path/'inverse'), [HT_FRAG], profile, \
            readMeteorites(str(met_path)), None, None, scan_bottom=48000.0, mc_runs=4, n_kicks=2, n_lines=3, \
            seed=3)

    assert results['settings']['method'] == 'breakup'
    assert results['settings']['scan_hts'][0] == HT_FRAG + BREAKUP_MARGIN

    met = results['meteorites'][0]
    assert len(met['starts']) == 1
    found = met['starts'][0]
    assert abs(found['ht'] - ht_true) < 20.0
    assert ('back' not in found) and ('breakups' not in found)

    # The kick is that of Passey and Melosh for the K of the meteorite and of the leading body there
    len_ground = _startGround(drag_fit, line, found['ht'])
    vel = np.interp(len_ground, line['len'], line['vel_parent'])
    rho_air = np.interp(len_ground, line['len'], line['rho'])
    k_main = leadingBallisticCoeff(drag_fit, inertialLength(line, len_ground))
    assert abs(found['kick_angle']/pmKickAngle(vel, rho_air, 3500.0, k_met/k_main) - 1) < 1e-6
    assert found['n_kick'] == 2

    # The spread of the falls gives a range of breakups around the one found, with the lines drawn
    assert found['n_line'] == 3
    assert abs(found['sd_along']**2 - found['sd_wind'][0]**2 - found['sd_kick'][0]**2 \
        - found['sd_line'][0]**2) < 1e-9
    assert abs(found['ht_range'][1] - found['ht'] - found['sd_along']*found['slope']) < 1e-6
    with open(str(tmp_path/'inverse'/'inverse.txt')) as f:
        report = f.read()
    assert "Breakup at" in report
    assert "with its speed there" in report
    assert "Start" not in report
