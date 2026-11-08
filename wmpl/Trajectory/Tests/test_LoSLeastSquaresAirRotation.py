""" Test of the least squares line-of-sight fit (los_fit='least_squares', wmpg#203) with the displacement of a
    decelerating meteor by drag against the rotating air (air_rotation, wmpg#192): the residuals move the
    points as the angles of the cost do.

Run with pytest:
    python -m pytest wmpl/Trajectory/Tests/test_LoSLeastSquaresAirRotation.py -v
"""

from __future__ import print_function, division, absolute_import

import numpy as np

import matplotlib
matplotlib.use('Agg')

from wmpl.Trajectory.Tests.SyntheticTrajectory import JDT_REF, BentDragMeteor, observeMeteor


# Fireball of the air rotation tests: 18 km/s at 45 deg, slowing down to about 10 km/s in 3.3 s
FIREBALL = (18000.0, 120.0, 45.0, 0.01, 0.05, 2.7, 6.0)



def testLeastSquaresRecoversAFireballInTheRotatingAir():
    """ Check that the least squares fit recovers the radiant of the noise-free fireball which drag against
        the rotating air pulls across its path, as the fit of the angles does. Without the displacement in
        the residuals, the least squares radiant is far off.
    """

    meteor = BentDragMeteor(*FIREBALL, air_rotation=True)

    traj = observeMeteor(meteor, noise=0.0, gravity_correction=True, air_rotation=True, calc_orbit=False, \
        los_fit='least_squares')
    traj.run()

    assert traj.air_rotation_shift is not None

    # The radiant is the direction of motion at the first point of the solution
    t_beg = (traj.jdt_ref - JDT_REF)*86400.0
    chord = np.linalg.norm(traj.radiant_eci_mini - meteor.radiantAtTime(t_beg))

    assert 3600*np.degrees(2*np.arcsin(chord/2)) < 2.0
