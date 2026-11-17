""" Tests of the hand-over of a drag fit to OpenDarkflight in wmpl.Utils.DragFitDarkFlight: the joint clones
    of K and the state at the last observed point, the lines drawn around the fitted one, the fragments of the
    last fragmentation, the meteorites, the OpenDarkflight inputs and the reading of its outputs. The
    trajectory is the synthetic fireball of test_DragFitFragmentations.py, with a fragmentation spread over
    1.5 km.

Run with pytest:
    python -m pytest wmpl/Utils/Tests/test_DragFitDarkFlight.py -v
"""

from __future__ import print_function, division, absolute_import

import configparser
import contextlib
import csv
import functools
import glob
import io
import os
import sys
from types import SimpleNamespace

import copy

import numpy as np
import pytest

from wmpl.Trajectory.Trajectory import lineFromFrame
from wmpl.Utils.AtmosphereProfile import AtmosphereProfile
from wmpl.Utils.DragFitDarkFlight import DENSITY_MIN, alongCross, drawLine, drawMeteorites, \
    findMcUncertainties, handoverClones, handoverReport, handoverState, lineMove, loadHandover, loadMcLines, \
    prepareHandover, readOdfRun, stateAtLength, summarizeFalls, writeOdfInput
from wmpl.Utils.Pickling import savePickle
from wmpl.Utils.Tests.test_DragFitAtmosphereProfile import loadProfile
from wmpl.Utils.SampleTrajectoryPositions import sampleTrajectory
from wmpl.Utils.Tests.test_DragFitFragmentations import HT_FRAG, fitDrag, fragmentedTrajectory, \
    lastBallisticCoeff


# Spread of the fragmentation (m), number of clones and of fragments, and largest K factor of a fragment
SPREAD = 1500.0
N_CLONES = 100
N_FRAGMENTS = 20
FRAG_K_MAX = 10.0



@functools.lru_cache(maxsize=None)
def clonedHandover(geometry=True):
    """ The hand-over of the fireball with the spread fragmentation, with or without the lines drawn
        (computed once for each).
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    with contextlib.redirect_stdout(io.StringIO()):
        return handoverClones(traj, [HT_FRAG], n_clones=N_CLONES, seed=1, geometry=geometry, \
            n_fragments=N_FRAGMENTS, frag_k_max=FRAG_K_MAX)



def testTheClonesComeFromTheAcceptableSpreads():
    """ Every clone comes from a spread within Delta chi2 <= 4, a larger K comes with a lower speed, and as
        the jump does not fit this spread fragmentation, the nominal state is that of the best spread.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    nominal, clones, _, scan = clonedHandover()

    spreads = np.array(scan['spreads'])
    acceptable = spreads[np.array(scan['delta_chi2']) <= 4.0]
    k_arr = np.array([state['k'] for state in clones])
    vel_arr = np.array([state['vel'] for state in clones])

    assert len(clones) == N_CLONES
    assert set(state['spread'] for state in clones) <= set(acceptable)
    assert len(set(state['spread'] for state in clones)) > 1
    assert np.corrcoef(np.log(k_arr), vel_arr)[0, 1] < -0.5

    delta_chi2 = np.array(scan['delta_chi2'])
    assert delta_chi2[0] > 4.0
    assert nominal['spread'] == spreads[np.argmin(delta_chi2)]
    best_fit = fitDrag(traj, knot_hts=[HT_FRAG], knot_spread=nominal['spread'])
    assert abs(nominal['k']/lastBallisticCoeff(best_fit) - 1) < 1e-9



def testTheNominalStateIsThatOfTheBestSpreadAlsoWhenTheJumpsFit(monkeypatch):
    """ The nominal state is that of the best spread also when the jumps fit the data within Delta chi2 <= 4,
        as they are the end of the range of the spreads, with the smallest K.
    """

    module = sys.modules['wmpl.Utils.DragFitDarkFlight']
    scan = module.spreadScan


    def jumpsAcceptable(*args, **kwargs):
        """ Scan the spreads, with the jumps acceptable but not the best. """

        result = list(scan(*args, **kwargs))
        result[2] = np.concatenate([[2.0], result[2][1:]])

        return tuple(result)


    monkeypatch.setattr(module, 'spreadScan', jumpsAcceptable)

    traj, _ = fragmentedTrajectory(SPREAD)
    with contextlib.redirect_stdout(io.StringIO()):
        nominal, _, _, scan_out = module.handoverClones(traj, [HT_FRAG], n_clones=1, geometry=False)

    assert scan_out['delta_chi2'][0] == 2.0
    assert nominal['spread'] == scan_out['spreads'][int(np.argmin(scan_out['delta_chi2']))]
    assert nominal['spread'] > 0



def testTheNominalStateIsThatOfTheSamples():
    """ For a fragmentation which is a jump, the nominal state at the last point is that of the fit with
        jumps, at the position and in the ground-fixed direction which SampleTrajectoryPositions gives there
        for the same fit, with the speed relative to the ground the drag of the fit acts on, the inertial
        speed less the rotation of the Earth along the trajectory.
    """

    traj, _ = fragmentedTrajectory()
    with contextlib.redirect_stdout(io.StringIO()):
        nominal, _, _, _ = handoverClones(traj, [HT_FRAG], n_clones=1)
        samples = sampleTrajectory(traj, nominal['time'], nominal['time'], 1.0, by_time=True, \
            knot_hts=[HT_FRAG])
    best_fit = fitDrag(traj, knot_hts=[HT_FRAG])

    assert nominal['spread'] == 0.0

    assert abs(np.degrees(samples.lat[0]) - nominal['lat']) < 1e-6
    assert abs(np.degrees(samples.lon[0]) - nominal['lon']) < 1e-6
    assert abs(samples.ht[0] - nominal['ht']) < 1.0
    assert abs(np.degrees(samples.azim_norot[0]) - nominal['azim']) < 1e-6
    assert abs(np.degrees(samples.elev_norot[0]) - nominal['elev']) < 1e-6
    vel_ground = best_fit.velocityAtLength(np.array([best_fit.len_last]))[0] - best_fit.v_rotation
    assert abs(nominal['vel'] - vel_ground) < 0.1
    assert abs(samples.vel_norot[0] - nominal['vel']) < 0.1



def testTheInputCarriesTheState(tmp_path):
    """ The OpenDarkflight input starts from the state with its K, in the winds of the profile, with the
        further options added.
    """

    nominal, _, _, _ = clonedHandover()

    profile_path = tmp_path/'profile.csv'
    with open(str(profile_path), 'w') as f:
        f.write("height,temperature,pressure,relative_humidity,wind_horizontal,wind_direction,wind_east," \
            + "wind_north,wind_up,density\n")
        for ht in (0.0, 50000.0, 100000.0):
            f.write("{:.1f},250.0,1000.0,0.0,5.0,90.0,0.0,0.0,0.0,{:.6e}\n".format(ht, 1.2*np.exp(-ht/7200)))

    ini_path = str(tmp_path/'clone_df_params.ini')
    writeOdfInput(ini_path, nominal, '2020-01-01 00:00:00.000000', AtmosphereProfile(str(profile_path)), \
        3500.0, shape='sphere', mc_runs=5, seed=7, odf_options=['Monte Carlo.wind_speed_sigma=0.2'])

    config = configparser.ConfigParser()
    config.read(ini_path)

    assert config['Atmosphere']['winds_file'] == str(profile_path)
    assert abs(float(config['Meteorite']['ballistic_coeff'])/nominal['k'] - 1) < 1e-5
    assert abs(1000*float(config['Ejection']['vel']) - nominal['vel']) < 0.01
    assert abs(1000*float(config['Ejection']['ht']) - nominal['ht']) < 0.01
    assert abs(float(config['Ejection']['az']) - nominal['azim']) < 1e-4
    assert abs(float(config['Ejection']['alt']) - nominal['elev']) < 1e-4
    assert config['Meteorite']['shape'] == 'sphere'
    assert config['Monte Carlo']['mc_enable'] == 'True'
    assert config['Monte Carlo']['mc_runs'] == '5'
    assert config['Monte Carlo']['wind_speed_sigma'] == '0.2'



def testTheFallsAreRead(tmp_path):
    """ The mass comes from the name of the nominal output, its fall from the last row and its heights and
        speeds from all rows, and the Monte Carlo falls from their table.
    """

    with open(str(tmp_path/'nominal_250.5_g.txt'), 'w') as f:
        f.write("# header\n")
        f.write("0.000 | 10.0 | 20.0 | 15.0 | 4000.0 | 83.0 | 27.0\n")
        f.write("0.002 | 10.0 | 20.0 | 14.9 | 3994.0 | 83.0 | 27.0\n")
        f.write("80.000 | 10.1 | 19.9 | 0.2 | 120.0 | 250.0 | 85.0\n")

    with open(str(tmp_path/'mc_250.5_g.txt'), 'w') as f:
        f.write("# header\n")
        f.write("81.0 | 10.11 | 19.91 | 0.2 | 121.0\n")
        f.write("82.0 | 10.12 | 19.92 | 0.2 | 122.0\n")

    run = readOdfRun(str(tmp_path))

    assert abs(run['mass'] - 0.2505) < 1e-12
    assert (run['lat'], run['lon']) == (10.1, 19.9)
    assert np.allclose(run['ht'], [15000.0, 14900.0, 200.0])
    assert run['vel'] == [4000.0, 3994.0, 120.0]
    assert run['mc_lat'] == [10.11, 10.12]
    assert run['mc_lon'] == [19.91, 19.92]
    assert readOdfRun(str(tmp_path/'empty')) is None



def testAlongAndAcrossTheMotion():
    """ For a motion to the north, a point to the north is downrange and a point to the east is to the
        right.
    """

    along, cross = alongCross(np.array([0.01, 0.0]), np.array([0.0, 0.01]), 0.0, 0.0, 0.0)

    assert np.allclose(along, [1.112, 0.0], atol=1e-3)
    assert np.allclose(cross, [0.0, 1.112], atol=1e-3)



def testALineMovesThePointAndTurnsTheDirection():
    """ A rotation of the line moves a point by its distance from the middle over the half length and turns
        the radiant by as much, and a shift moves every point by the shift without turning the radiant.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    middle, half_length, radiant, axes = traj.los_frame
    point = middle - 0.8*half_length*radiant

    line = lineFromFrame(np.array([10.0, 0.0, 0.0, 0.0]), traj.los_frame)
    move, turn = lineMove(traj, point, (line[:3], line[3:]))
    assert np.allclose(move, -0.8*10.0*axes[0], atol=1e-2)
    assert np.allclose(turn, 10.0/half_length*axes[0], atol=1e-6)

    line = lineFromFrame(np.array([0.0, 0.0, 0.0, 5.0]), traj.los_frame)
    move, turn = lineMove(traj, point, (line[:3], line[3:]))
    assert np.allclose(move, 5.0*axes[1], atol=1e-9)
    assert np.allclose(turn, 0.0, atol=1e-12)



def testTheLinesOfTheClonesAreDrawn():
    """ The lines drawn from the covariance of the line-of-sight fit move the starting points of the clones
        by about its size at the end of the meteor. Without them the points of the clones of one spread only
        move with the rotation of the Earth between their times, by much less.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    _, clones, _, _ = clonedHandover()
    _, clones_fixed, _, _ = clonedHandover(geometry=False)

    def _moves(states):
        lat = np.radians(np.array([state['lat'] for state in states]))
        lon = np.radians(np.array([state['lon'] for state in states]))
        north = (lat - lat[0])*6371e3
        east = (lon - lon[0])*6371e3*np.cos(lat[0])
        return np.hypot(north, east)

    # The clones of one spread share the point of the fitted line, up to the moves of their lines
    same_spread = [state for state in clones if state['spread'] == clones[0]['spread']]
    same_spread_fixed = [state for state in clones_fixed if state['spread'] == clones_fixed[0]['spread']]
    size = np.sqrt(np.trace(traj.los_covariance))

    assert 0.2*size < np.median(_moves(same_spread)[1:]) < 5*size
    assert np.median(_moves(same_spread_fixed)[1:]) < 0.1*size



def testTheFragmentsStartAtTheLastFragmentation():
    """ The fragments start at the height of the last fragmentation with the speed of the nominal fit there,
        and a K between that of the leading body after the fragmentation and FRAG_K_MAX times it.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    nominal, _, fragments, _ = clonedHandover(geometry=False)
    nominal_fit = fitDrag(traj, knot_hts=[HT_FRAG], knot_spread=nominal['spread'])

    knot = np.where(nominal_fit.pen_weight == 0)[0][0]
    len_knot = nominal_fit.ramp_beg[knot]
    len_after = len_knot + nominal_fit.ramp_len[knot]
    k_after = np.exp(nominal_fit._logK(nominal_fit.params, np.array([len_after]))[0])
    k_factors = np.array([state['k']/k_after for state in fragments])

    assert len(fragments) == N_FRAGMENTS
    assert np.allclose([state['ht'] for state in fragments], HT_FRAG, atol=1.0)
    assert np.allclose([state['vel_inertial'] for state in fragments], \
        nominal_fit.velocityAtLength(np.array([len_knot]))[0], atol=1.0)
    assert np.all(k_factors > 1 - 1e-6) and np.all(k_factors < FRAG_K_MAX)
    assert np.allclose(k_factors, [state['k_factor'] for state in fragments], rtol=1e-6)



def testTheMeteoritesAreDrawn():
    """ Every state gets one of the shapes and a density around the mean, never below DENSITY_MIN. """

    states = [dict() for _ in range(400)]
    drawMeteorites(states, ['sphere', 'brick'], 3500.0, density_sigma=500.0, seed=3)
    densities = np.array([state['density'] for state in states])

    assert set(state['shape'] for state in states) == {'sphere', 'brick'}
    assert abs(np.mean(densities) - 3500.0) < 100.0
    assert abs(np.std(densities) - 500.0) < 60.0
    assert np.min(densities) >= DENSITY_MIN

    drawMeteorites(states, None, 3000.0)
    assert all((state['shape'] is None) and (state['density'] == 3000.0) for state in states)



def testTheMonteCarloLinesAreFoundAndDrawn(tmp_path):
    """ The uncertainties the solver saves with its Monte Carlo are found next to the trajectory pickle, their
        lines are loaded and drawn, and a line shifted across moves the starting point by the shift without
        turning the direction, other than the local north turns over the shift.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    radiant = traj.radiant_eci_mini/np.linalg.norm(traj.radiant_eci_mini)
    across = np.cross(radiant, [0.0, 0.0, 1.0])
    across /= np.linalg.norm(across)
    shifted = (np.array(traj.state_vect_mini) + 7.0*across, np.array(traj.radiant_eci_mini))

    # Save the uncertainties as the solver does, with one trajectory without a line
    mc_traj_list = [SimpleNamespace(state_vect_mini=shifted[0], radiant_eci_mini=shifted[1]), \
        SimpleNamespace(state_vect_mini=None, radiant_eci_mini=None)]
    savePickle(SimpleNamespace(mc_traj_list=mc_traj_list), str(tmp_path/'Monte Carlo'), \
        'event_mc_uncertainties.pickle')
    traj_path = str(tmp_path/'event_trajectory.pickle')

    mc_path = findMcUncertainties(traj_path)
    mc_lines = loadMcLines(mc_path)

    assert mc_path == os.path.join(str(tmp_path), 'Monte Carlo', 'event_mc_uncertainties.pickle')
    assert findMcUncertainties(str(tmp_path/'other_trajectory.pickle')) is None
    assert len(mc_lines) == 1
    assert np.array_equal(drawLine(traj, np.random.default_rng(0), mc_lines=mc_lines)[0], shifted[0])

    # The starting point moves by the shift and the direction stays
    jump_fit = fitDrag(traj, knot_hts=[HT_FRAG])
    args = (jump_fit, jump_fit.len_last) + stateAtLength(jump_fit, jump_fit.params, jump_fit.len_last)
    fixed, moved = handoverState(*args), handoverState(*args, line=mc_lines[0])

    north = np.radians(moved['lat'] - fixed['lat'])*6371e3
    east = np.radians(moved['lon'] - fixed['lon'])*6371e3*np.cos(np.radians(fixed['lat']))
    assert abs(np.sqrt(north**2 + east**2 + (moved['ht'] - fixed['ht'])**2) - 7.0) < 0.1
    # The direction stays, up to the turn of the local north over the shift
    assert abs(moved['azim'] - fixed['azim']) < 1e-3
    assert abs(moved['elev'] - fixed['elev']) < 1e-3



def writeFakeOdfRun(run_dir, lat, lon, ht_arr=None, vel_arr=None, mc_falls=()):
    """ Write the outputs of an OpenDarkflight run of a 1.5 kg meteorite which falls at the given point, with
        the given heights (m) and speeds (m/s) before, and the given Monte Carlo falls (lat, lon).
    """

    ht_arr = [15000.0, 200.0] if ht_arr is None else ht_arr
    vel_arr = [4000.0, 100.0] if vel_arr is None else vel_arr

    with open(os.path.join(run_dir, 'nominal_1.5_kg.txt'), 'w') as f:
        f.write("# header\n")
        for i, (ht, vel) in enumerate(zip(ht_arr, vel_arr)):
            f.write("{:.3f} | {:.6f} | {:.6f} | {:.4f} | {:.2f}\n".format(float(i), lat, lon, ht/1000, vel))

    with open(os.path.join(run_dir, 'mc_1.5_kg.txt'), 'w') as f:
        f.write("# header\n")
        for lat_mc, lon_mc in mc_falls:
            f.write("90.0 | {:.6f} | {:.6f} | 0.2 | 100.0\n".format(lat_mc, lon_mc))



def testTheHandoverIsPreparedAndItsFallsCollected(tmp_path):
    """ prepareHandover() writes an input for the nominal state and every clone and fragment and saves the
        hand-over, which loadHandover() reads back and handoverReport() describes. With the outputs of
        OpenDarkflight in place, summarizeFalls() collects every fall to falls.csv and compares the speeds
        below the hand-over.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    profile = loadProfile(tmp_path/'profile.csv')
    out_dir = str(tmp_path/'handover')

    # Leave the outputs of an earlier run of another mass in a run directory
    stale_dir = os.path.join(out_dir, 'clones', 'clone_0001')
    os.makedirs(stale_dir)
    writeFakeOdfRun(stale_dir, 9.0, 9.0)
    os.rename(os.path.join(stale_dir, 'nominal_1.5_kg.txt'), os.path.join(stale_dir, 'nominal_850.0_g.txt'))

    with contextlib.redirect_stdout(io.StringIO()):
        handover, ini_paths = prepareHandover(traj, out_dir, [HT_FRAG], profile, n_clones=6, seed=2, \
            n_fragments=3, shapes=['sphere', 'brick'], density_sigma=100.0, mc_runs=2, drag_model='dfn', \
            odf_options=['Monte Carlo.wind_speed_sigma=0.2'])

    assert os.listdir(stale_dir) == ['clone_0001_df_params.ini']

    # The inputs: the nominal one without a Monte Carlo, the fragments with their kicks
    assert len(ini_paths) == 1 + 6 + 3
    assert ini_paths[0].endswith(os.path.join('nominal', 'nominal_df_params.ini'))
    configs = []
    for ini_path in ini_paths:
        config = configparser.ConfigParser()
        config.read(ini_path)
        configs.append(config)

    assert configs[0]['Monte Carlo']['mc_enable'] == 'False'
    assert all(config['Monte Carlo']['mc_enable'] == 'True' for config in configs[1:])
    assert len(set(config['Monte Carlo']['random_seed'] for config in configs[1:])) == 9
    assert all(config['Meteorite']['drag_model'] == 'dfn' for config in configs)
    assert all(config['Monte Carlo']['wind_error_mode'] == 'statistical' for config in configs[1:])
    assert all(config['Monte Carlo']['wind_speed_sigma'] == '0.2' for config in configs[1:])
    for config, state in zip(configs[7:], handover['fragments']):
        assert config['Deviating Fragments']['v_perp_model'] == 'passey_melosh'
        assert abs(float(config['Deviating Fragments']['pm_radius_ratio']) - state['k_factor']) < 1e-4

    # The hand-over saved and described
    saved = loadHandover(out_dir)
    assert saved['clones'] == handover['clones']
    assert saved['settings']['n_fragments'] == 3
    report = handoverReport(handover)
    assert "Hand-over at the last point" in report
    assert "Jumps: K = {:.3e} m^2/kg".format(handover['scan']['k_last'][0]) in report
    assert "Fragments: 3 from" in report

    # Put the outputs of OpenDarkflight in place, the nominal one with the speeds of the drag fit below
    nominal = handover['nominal']
    ht_track = [nominal['ht']] + nominal['check_ht'] + [200.0]
    vel_track = [nominal['vel']] + nominal['check_vel'] + [100.0]
    writeFakeOdfRun(os.path.dirname(ini_paths[0]), 0.0, 0.0, ht_arr=ht_track, vel_arr=vel_track)
    for ini_path in ini_paths[1:]:
        writeFakeOdfRun(os.path.dirname(ini_path), 0.001, 0.0, mc_falls=[(0.002, 0.0), (0.003, 0.0)])

    summary = summarizeFalls(out_dir, handover)

    with open(os.path.join(out_dir, 'falls.csv')) as f:
        falls = list(csv.DictReader(f))

    assert len(falls) == (6 + 3)*3
    assert sum(row['kind'] == 'fragment' for row in falls) == 3*3
    assert os.path.isfile(os.path.join(out_dir, 'falls.png'))
    for v_fit in nominal['check_vel']:
        assert "{:.0f}/{:.0f} m/s".format(v_fit, v_fit) in summary
    assert "Leading body: 12 falls of 6 runs" in summary
    assert "median fall: 0.002500, 0.000000" in summary
    assert "Fragments of the last fragmentation: 6 falls of 3 runs" in summary
    assert len(glob.glob(os.path.join(out_dir, 'clones', 'clone_*'))) == 6

    # A run without Monte Carlo falls counts with its nominal fall
    os.remove(os.path.join(out_dir, 'clones', 'clone_0002', 'mc_1.5_kg.txt'))
    summary = summarizeFalls(out_dir, handover)
    assert "Leading body: 11 falls of 6 runs" in summary
    assert "1 runs without Monte Carlo falls count with their nominal fall" in summary





def testTheHandoverRefusesAnIncoherentAtmosphereOrSettings(tmp_path):
    """ The profile has to cover the starting heights and keep its winds, as OpenDarkflight takes them, the
        trajectory needs an orbit for the rotation of the Earth, and the seed and the largest K factor of the
        fragments have to make sense.
    """

    traj, _ = fragmentedTrajectory(SPREAD)

    with contextlib.redirect_stdout(io.StringIO()):

        with pytest.raises(ValueError, match="would not start in the air"):
            prepareHandover(traj, str(tmp_path/'low'), [HT_FRAG], loadProfile(tmp_path/'low.csv', \
                ht_max=40000.0), n_clones=2)

        with pytest.raises(ValueError, match="winds"):
            prepareHandover(traj, str(tmp_path/'still'), [HT_FRAG], loadProfile(tmp_path/'still.csv', \
                use_winds=False), n_clones=2)

        no_orbit = copy.copy(traj)
        no_orbit.orbit = None
        with pytest.raises(ValueError, match="no orbit"):
            prepareHandover(no_orbit, str(tmp_path/'no_orbit'), [HT_FRAG], loadProfile(tmp_path/'p.csv'), \
                n_clones=2)

        with pytest.raises(ValueError, match="seed"):
            handoverClones(traj, [HT_FRAG], n_clones=2, seed=-1)

        with pytest.raises(ValueError, match="frag_k_max"):
            handoverClones(traj, [HT_FRAG], n_clones=2, frag_k_max=0.5)



def testTheCovarianceLinesAreCentredOnTheSolution():
    """ The lines drawn from the covariance are centred on the line of the solution, also when it is not
        the line of the frame, as after the Monte Carlo of the solver.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    middle, _, radiant, axes = traj.los_frame
    point = middle - 0.8*traj.los_frame[1]*radiant

    # Move the line of the solution 50 m away from that of the frame
    moved = copy.copy(traj)
    moved.state_vect_mini = np.array(traj.state_vect_mini) + 50.0*axes[0]

    rng = np.random.default_rng(4)
    moves = np.array([lineMove(moved, point, drawLine(moved, rng))[0] for _ in range(400)])

    assert np.linalg.norm(np.mean(moves, axis=0)) < 1.0
    assert 0.5*np.sqrt(np.trace(traj.los_covariance)) < np.median(np.linalg.norm(moves, axis=1))



def testSeveralNominalOutputsAreRefused(tmp_path):
    """ Two nominal outputs in one run directory belong to different runs, and are not read. """

    writeFakeOdfRun(str(tmp_path), 1.0, 2.0)
    os.rename(str(tmp_path/'nominal_1.5_kg.txt'), str(tmp_path/'nominal_850.0_g.txt'))
    writeFakeOdfRun(str(tmp_path), 1.0, 2.0)

    with pytest.raises(ValueError, match="nominal outputs"):
        readOdfRun(str(tmp_path))



def testFragmentsAboveTheProfileAreLeftOut(tmp_path):
    """ A profile which covers the starting point of the leading body but not the fragmentation above it
        leaves the fragments out, and the report says why.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    profile = loadProfile(tmp_path/'profile.csv', ht_max=HT_FRAG - 1000.0)

    with contextlib.redirect_stdout(io.StringIO()):
        handover, ini_paths = prepareHandover(traj, str(tmp_path/'handover'), [HT_FRAG], profile, \
            n_clones=2, n_fragments=3)

    assert handover['fragments'] == []
    assert len(ini_paths) == 1 + 2
    assert "No fragments: the fragmentation at" in handoverReport(handover)
