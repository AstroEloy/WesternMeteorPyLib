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
import functools
import io

import numpy as np

from wmpl.Utils.AtmosphereProfile import AtmosphereProfile
from wmpl.Utils.DragFitDarkFlight import DENSITY_MIN, alongCross, drawMeteorites, handoverClones, lineMove, \
    readOdfRun, writeOdfInput
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
    """ Every clone comes from a spread within Delta chi2 <= 4, the nominal K is that of the fit with jumps,
        and a larger K comes with a lower speed.
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
    assert abs(nominal['k']/lastBallisticCoeff(fitDrag(traj, knot_hts=[HT_FRAG])) - 1) < 1e-9
    assert np.corrcoef(np.log(k_arr), vel_arr)[0, 1] < -0.5



def testTheNominalStateIsThatOfTheSamples():
    """ The nominal state at the last point is the one SampleTrajectoryPositions gives there for the same fit:
        position, ground-fixed radiant and speed relative to the ground.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    nominal, _, _, _ = clonedHandover()

    with contextlib.redirect_stdout(io.StringIO()):
        samples = sampleTrajectory(traj, nominal['time'], nominal['time'], 1.0, by_time=True, \
            knot_hts=[HT_FRAG])

    assert abs(np.degrees(samples.lat[0]) - nominal['lat']) < 1e-6
    assert abs(np.degrees(samples.lon[0]) - nominal['lon']) < 1e-6
    assert abs(samples.ht[0] - nominal['ht']) < 1.0
    assert abs(np.degrees(samples.azim_norot[0]) - nominal['azim']) < 1e-6
    assert abs(np.degrees(samples.elev_norot[0]) - nominal['elev']) < 1e-6
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

    move, turn = lineMove(traj, point, np.array([10.0, 0.0, 0.0, 0.0]))
    assert np.allclose(move, -0.8*10.0*axes[0], atol=1e-2)
    assert np.allclose(turn, 10.0/half_length*axes[0], atol=1e-6)

    move, turn = lineMove(traj, point, np.array([0.0, 0.0, 0.0, 5.0]))
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
    """ The fragments start at the height of the last fragmentation with the speed of the fit with jumps
        there, and a K between that of the leading body after the jump and FRAG_K_MAX times it.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    _, _, fragments, _ = clonedHandover(geometry=False)
    jump_fit = fitDrag(traj, knot_hts=[HT_FRAG])

    len_knot = jump_fit.ramp_beg[jump_fit.pen_weight == 0][0]
    k_after = jump_fit.ballisticCoeffAtLength(np.array([len_knot + 2.0]))[0][0]
    k_factors = np.array([state['k']/k_after for state in fragments])

    assert len(fragments) == N_FRAGMENTS
    assert np.allclose([state['ht'] for state in fragments], HT_FRAG, atol=1.0)
    assert np.allclose([state['vel_inertial'] for state in fragments], \
        jump_fit.velocityAtLength(np.array([len_knot]))[0], atol=1.0)
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

