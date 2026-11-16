""" Tests of the hand-over of a drag fit to OpenDarkflight in wmpl.Utils.DragFitDarkFlight: the joint clones
    of K and the state at the last observed point, the OpenDarkflight inputs and the reading of its outputs.
    The trajectory is the synthetic fireball of test_DragFitFragmentations.py, with a fragmentation spread
    over 1.5 km.

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
from wmpl.Utils.DragFitDarkFlight import alongCross, handoverClones, readOdfRun, writeOdfInput
from wmpl.Utils.SampleTrajectoryPositions import sampleTrajectory
from wmpl.Utils.Tests.test_DragFitFragmentations import HT_FRAG, fitDrag, fragmentedTrajectory, \
    lastBallisticCoeff


# Spread of the fragmentation (m) and number of clones
SPREAD = 1500.0
N_CLONES = 100



@functools.lru_cache(maxsize=None)
def clonedHandover():
    """ The hand-over of the fireball with the spread fragmentation (computed once). """

    traj, _ = fragmentedTrajectory(SPREAD)
    with contextlib.redirect_stdout(io.StringIO()):
        return handoverClones(traj, [HT_FRAG], n_clones=N_CLONES, seed=1)



def testTheClonesComeFromTheAcceptableSpreads():
    """ Every clone comes from a spread within Delta chi2 <= 4, the nominal K is that of the fit with jumps,
        and a larger K comes with a lower speed.
    """

    traj, _ = fragmentedTrajectory(SPREAD)
    nominal, clones, scan = clonedHandover()

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
    nominal, _, _ = clonedHandover()

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

    nominal, _, _ = clonedHandover()

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
