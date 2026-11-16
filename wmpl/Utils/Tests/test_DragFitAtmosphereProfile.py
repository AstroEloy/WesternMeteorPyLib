""" Tests of the atmosphere profile in the drag model of wmpl.Utils.DragFit: the air density of the profile,
    with the MSIS model joined to it outside its heights, and its winds along the direction of motion. The
    profiles are written as OpenDarkflight's wrf CSV files, of an exponential atmosphere.

Run with pytest:
    python -m pytest wmpl/Utils/Tests/test_DragFitAtmosphereProfile.py -v
"""

from __future__ import print_function, division, absolute_import

import contextlib
import functools
import io

import numpy as np

from wmpl.Trajectory.Tests.SyntheticTrajectory import DragMeteor, observeMeteor
from wmpl.Utils.AtmosphereProfile import AtmosphereProfile
from wmpl.Utils.DragFit import DragFit



@functools.lru_cache(maxsize=None)
def solvedTrajectory():
    """ Simulate, observe and solve the synthetic fireball, seen from 86 to 48 km (computed once). """

    traj = observeMeteor(DragMeteor(18000.0, 120.0, 45.0, 0.01, 0.05, 2.7, 6.0))
    with contextlib.redirect_stdout(io.StringIO()):
        traj = traj.run()

    return traj



def fitDrag(traj, **kwargs):
    """ Fit the drag model without printing. """

    with contextlib.redirect_stdout(io.StringIO()):
        return DragFit(traj, **kwargs)



def loadProfile(path, density_scale=1.0, wind_speed=0.0, wind_dir=0.0, ht_max=120000.0, use_winds=True):
    """ Write a wrf profile of an exponential atmosphere with a constant wind, and load it.

    Arguments:
        path: [str] Path of the file to write.

    Keyword arguments:
        density_scale: [float] Factor of the density. 1 by default.
        wind_speed: [float] Wind speed (m/s). 0 by default.
        wind_dir: [float] Direction the wind blows from (deg, +E of N). 0 by default.
        ht_max: [float] Top of the profile (m). 120 km by default.
        use_winds: [bool] Use the winds of the profile. True by default.

    Return:
        [AtmosphereProfile object] The profile.
    """

    heights = np.arange(0.0, ht_max + 1.0, 250.0)
    densities = density_scale*1.225*np.exp(-heights/7200.0)

    with open(str(path), 'w') as f:
        f.write("height,temperature,pressure,relative_humidity,wind_horizontal,wind_direction,wind_east," \
            + "wind_north,wind_up,density\n")
        for ht, rho in zip(heights, densities):
            f.write("{:.1f},250.0,1000.0,0.0,{:.3f},{:.3f},0.0,0.0,0.0,{:.6e}\n".format(ht, wind_speed, \
                wind_dir, rho))

    atm_profile = AtmosphereProfile(str(path))
    atm_profile.use_winds = use_winds

    return atm_profile



def lastBallisticCoeff(drag_fit):
    """ Ballistic coefficient at the last point (m^2/kg). """

    return drag_fit.ballisticCoeffAtLength(np.array([drag_fit.len_last]))[0][0]



def testTheDensityOfTheProfileIsUsed(tmp_path):
    """ The air density of the model is that of the profile, and in a profile 20% denser the ballistic
        coefficient comes out 20% smaller, as only their product sets the drag.
    """

    traj = solvedTrajectory()
    profile = loadProfile(tmp_path/'profile.csv')
    denser = loadProfile(tmp_path/'denser.csv', density_scale=1.2)

    drag_fit = fitDrag(traj, atm_profile=profile)
    drag_fit_denser = fitDrag(traj, atm_profile=denser)

    len_arr = np.linspace(drag_fit.len_first, drag_fit.len_last, 20)
    assert np.allclose(drag_fit.airDensityAtLength(len_arr), profile.density(drag_fit.heightAtLength(len_arr)), \
        rtol=1e-6)
    assert abs(1.2*lastBallisticCoeff(drag_fit_denser)/lastBallisticCoeff(drag_fit) - 1) < 0.01



def testTheMsisModelJoinsTheProfileAboveItsTop(tmp_path):
    """ Above the top of a profile which ends within the observed part, the density is that of the MSIS model
        scaled to join the profile at its top, without a step.
    """

    traj = solvedTrajectory()
    profile = loadProfile(tmp_path/'profile.csv', ht_max=60000.0)

    drag_fit = fitDrag(traj, atm_profile=profile)
    msis_fit = fitDrag(traj)

    len_top = drag_fit.lengthAtHeight(np.array([60000.0 + 1.0, 60000.0 - 1.0]))
    rho_top = drag_fit.airDensityAtLength(len_top)
    assert abs(rho_top[0]/rho_top[1] - 1) < 1e-3

    len_above = drag_fit.lengthAtHeight(np.array([70000.0, 80000.0]))
    ratio = drag_fit.airDensityAtLength(len_above)/msis_fit.airDensityAtLength(len_above)
    assert abs(ratio[1]/ratio[0] - 1) < 1e-6
    assert abs(ratio[0]*msis_fit.airDensityAtLength(len_top[1])[0]/profile.density(60000.0) - 1) < 1e-2



def testATailWindLowersTheAirSpeed(tmp_path):
    """ A wind blowing along the direction of motion lowers the speed relative to the air by its speed, so the
        ballistic coefficient at the last point grows by the squared ratio of the air speeds, and switching
        the winds off gives back the fit in still air.
    """

    traj = solvedTrajectory()
    still = fitDrag(traj, atm_profile=loadProfile(tmp_path/'still.csv'))

    # Blow the wind from behind the meteor, along its horizontal direction of motion
    wind_speed = 100.0
    motion_azim = np.degrees(np.arctan2(still.motion_en[0], still.motion_en[1]))
    profile = loadProfile(tmp_path/'tail.csv', wind_speed=wind_speed, wind_dir=(motion_azim + 180.0)%360.0)
    drag_fit = fitDrag(traj, atm_profile=profile)

    wind_along = wind_speed*np.linalg.norm(still.motion_en)
    assert np.allclose(drag_fit.v_air_offset, drag_fit.v_rotation + wind_along)

    vel_air = still.velocityAtLength(np.array([still.len_last]))[0] - still.v_rotation
    expected = (vel_air/(vel_air - wind_along))**2
    assert abs(lastBallisticCoeff(drag_fit)/lastBallisticCoeff(still)/expected - 1) < 0.02

    profile.use_winds = False
    assert np.array_equal(fitDrag(traj, atm_profile=profile).params, still.params)
