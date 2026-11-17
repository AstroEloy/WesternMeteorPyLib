""" Drag-only flight of a fragment along the line of a fireball while it is still fast: down from where it
broke off the leading body until it has slowed down to the speed where its dark flight starts, or up from
there until it has the speed of the leading body again, which is where it broke off.

The fragment keeps to the line of the leading body, as drag does not turn it, and slows down by its ballistic
coefficient K = Gamma*S/m, dv/dt = -K*rho_air*v_air^2 + g_along, as in wmpl.Utils.DragFit. WMPL gives the line
as the ground track of the leading body, with the lengths along it and the speeds relative to the ground, as
drag keeps the direction of motion relative to the ground, and the drag acts on the speed relative to the air,
which is less by the wind along the motion (the air offset of the line). The air density, the air offset, the
speed of the leading body and the gravity along the line, g_along, come with the line. The bending of the path
by gravity is left out, as it is small while the fragment is fast. The fragment can ablate, with the mass
following its speed through the air as m = m_ref*exp(sigma/2*(v_air^2 - v_air_ref^2)), and its K growing as
m^(-1/3), but the dark flight never ablates. A flight up which starts as fast as the leading body, or any
flight which leaves the line or starts from a value which is not finite, does not get there.

K comes either from a fixed value at a reference mass, or from the drag models and shapes of OpenDarkflight
(odfBallisticCoeff), at the speed, height and mass of every step. For the latter, this file is run as a script
by the Python of OpenDarkflight, from its directory, so it imports nothing from WMPL:

    python LuminousFlight.py jobs.json results.json

The jobs file holds the line, the atmosphere profile and the flights (see runJobs()). WMPL writes it and reads
the results in wmpl.Utils.DragFitFragmentFalls.
"""

from __future__ import print_function, division, absolute_import

import json
import os
import sys

import numpy as np


# Step of the integration along the line (m)
LEN_STEP = 2.0



def lineAt(line, len_val):
    """ Interpolate the line of the leading body at a length along it.

    Arguments:
        line: [dict] Arrays along the line: 'len' (m, increasing), 'ht' (m), 'rho' (kg/m^3), 'air_offset',
            the speed less the speed relative to the air (m/s), 'vel_parent', the speed of the leading body
            (m/s), and 'g_along', the gravity along the motion (m/s^2).
        len_val: [float] Length along the line (m).

    Return:
        (ht, rho, air_offset, vel_parent, g_along): [tuple of floats] The values at that length.
    """

    return tuple(float(np.interp(len_val, line['len'], line[key])) for key in ('ht', 'rho', 'air_offset', \
        'vel_parent', 'g_along'))



def ablatedMass(mass_ref, vel_ref, vel, sigma):
    """ Compute the mass of an ablating fragment at a speed through the air, from its mass at a reference
        speed.

    Arguments:
        mass_ref: [float] Mass at the reference speed (kg).
        vel_ref: [float] Reference speed through the air (m/s).
        vel: [float] Speed through the air (m/s).
        sigma: [float] Ablation coefficient (s^2/m^2).

    Return:
        [float] Mass (kg).
    """

    return mass_ref*np.exp(0.5*sigma*(vel**2 - vel_ref**2))



def _fly(line, kFunc, len_start, vel_start, mass, sigma, direction, vel_end=None, len_step=LEN_STEP):
    """ Integrate the flight along the line, down (direction 1) until the speed falls to vel_end, or up
        (direction -1) until it reaches that of the leading body, by midpoint steps in length.
    """

    len_min, len_max = line['len'][0], line['len'][-1]
    vel_air_start = vel_start - lineAt(line, len_start)[2]

    def massAt(len_val, vel):
        return ablatedMass(mass, vel_air_start, vel - lineAt(line, len_val)[2], sigma)

    def deceleration(len_val, vel):

        ht, rho, air_offset, _, g_along = lineAt(line, len_val)
        vel_air = vel - air_offset

        return (kFunc(massAt(len_val, vel), ht, vel_air)*rho*vel_air**2 - g_along)/vel

    len_val, vel, time = float(len_start), float(vel_start), 0.0

    # A start which is not finite, off the line, or as fast as the leading body on the way up gets nowhere
    if not (np.isfinite(len_val) and np.isfinite(vel) and (len_min <= len_val <= len_max)):
        return len_val, vel, time, mass, False
    if (direction < 0) and (vel >= lineAt(line, len_val)[3]):
        return len_val, vel, time, mass, False

    while True:

        # Stop where the flight ends
        if (direction > 0) and (vel <= vel_end):
            return len_val, vel, time, massAt(len_val, vel), True
        if (direction < 0) and (vel >= lineAt(line, len_val)[3]):
            return len_val, vel, time, massAt(len_val, vel), True

        # Stop at the ends of the line, or where the speed is lost
        len_next = len_val + direction*len_step
        if (not (len_min <= len_next <= len_max)) or (not np.isfinite(vel)) or (vel <= 0):
            return len_val, vel, time, massAt(len_val, vel), False

        # Take a midpoint step, as dv/dL = (-K rho v_air^2 + g_along)/v going down
        vel_mid = vel - direction*0.5*len_step*deceleration(len_val, vel)
        vel_next = vel - direction*len_step*deceleration(len_val + direction*0.5*len_step, vel_mid)

        time += direction*len_step/vel_mid
        len_val, vel = len_next, vel_next



def flyDown(line, kFunc, len_start, vel_start, vel_end, mass, sigma=0.0, len_step=LEN_STEP):
    """ Fly a fragment down the line from where it broke off until its speed falls to the speed where its dark
        flight starts.

    Arguments:
        line: [dict] The line of the leading body (see lineAt()).
        kFunc: [function] Ballistic coefficient (m^2/kg) of the fragment for its mass (kg), the height (m)
            and the speed relative to the air (m/s).
        len_start: [float] Length along the line where it broke off (m).
        vel_start: [float] Inertial speed there (m/s).
        vel_end: [float] Inertial speed where the dark flight starts (m/s).
        mass: [float] Mass where it broke off (kg).

    Keyword arguments:
        sigma: [float] Ablation coefficient (s^2/m^2). 0 by default.
        len_step: [float] Step along the line (m). LEN_STEP by default.

    Return:
        (len_end, vel, time, mass_end, reached):
            - len_end: [float] Length where the dark flight starts (m).
            - vel: [float] Inertial speed there (m/s).
            - time: [float] Time of flight from the breakup (s).
            - mass_end: [float] Mass there (kg).
            - reached: [bool] False if the line ended before the speed fell to vel_end.
    """

    return _fly(line, kFunc, len_start, vel_start, mass, sigma, 1, vel_end=vel_end, len_step=len_step)



def flyUp(line, kFunc, len_start, vel_start, mass, sigma=0.0, len_step=LEN_STEP):
    """ Fly a fragment up the line, back in time, from where its dark flight starts until its speed is that of
        the leading body, which is where it broke off.

    Arguments:
        line: [dict] The line of the leading body (see lineAt()).
        kFunc: [function] Ballistic coefficient (m^2/kg) of the fragment for its mass (kg), the height (m)
            and the speed relative to the air (m/s).
        len_start: [float] Length along the line where the dark flight starts (m).
        vel_start: [float] Inertial speed there (m/s).
        mass: [float] Mass there (kg).

    Keyword arguments:
        sigma: [float] Ablation coefficient (s^2/m^2). 0 by default.
        len_step: [float] Step along the line (m). LEN_STEP by default.

    Return:
        (len_sep, vel, time, mass_sep, met):
            - len_sep: [float] Length where it broke off (m).
            - vel: [float] Inertial speed there (m/s).
            - time: [float] Time of flight from the breakup to the start of the dark flight (s).
            - mass_sep: [float] Mass where it broke off (kg).
            - met: [bool] False if the line ended before the fragment had the speed of the leading body, or
                if it was as fast as the leading body where it started.
    """

    len_sep, vel, time, mass_sep, met = _fly(line, kFunc, len_start, vel_start, mass, sigma, -1, \
        len_step=len_step)

    return len_sep, vel, -time, mass_sep, met



def fixedBallisticCoeff(k_ref, mass_ref):
    """ Make the function of a ballistic coefficient fixed at a reference mass, which grows as m^(-1/3) as the
        fragment ablates.

    Arguments:
        k_ref: [float] Ballistic coefficient at the reference mass (m^2/kg).
        mass_ref: [float] Reference mass (kg).

    Return:
        [function] K for the mass (kg), the height (m) and the speed relative to the air (m/s).
    """

    def kFunc(mass, ht, vel_air):
        return k_ref*(mass/mass_ref)**(-1.0/3)

    return kFunc



def odfBallisticCoeff(meteorite, atm, zenith):
    """ Make the function of the ballistic coefficient of a meteorite as OpenDarkflight integrates it, with
        its drag model and shape at the Mach and Reynolds numbers of the flow. Only runs in the environment of
        OpenDarkflight.

    Arguments:
        meteorite: [dict] Its density (kg/m^3), and its shape, drag_model, sphericity, flatness, elongation
            and volume_factor as in an OpenDarkflight input, where given.
        atm: [AtmosphereModel object] Atmosphere of OpenDarkflight, for the speed of sound and the
            temperature.
        zenith: [float] Zenith angle of the motion (radians), for the area of a shape which is not a sphere.

    Return:
        [function] K for the mass (kg), the height (m) and the speed relative to the air (m/s).
    """

    from opendf.Routines.Atmosphere import DragCoefficients, knudsenNumber, reynoldsNumber
    from opendf.Routines.MeteoriteShapes import computeShapeFactor, darkSurfaceToMassRatio, \
        volumeEquivalentDiameter

    shape = meteorite.get('shape', 'sphere')
    drag_options = {key: meteorite[key] for key in ('sphericity', 'flatness', 'elongation', 'volume_factor') \
        if meteorite.get(key) is not None}
    drag = DragCoefficients(model_type=meteorite.get('drag_model', 'ceplecha'), shape=shape, **drag_options)

    def kFunc(mass, ht, vel_air):

        radius, _, _ = computeShapeFactor(meteorite['density'], mass, shape=shape)
        sx, sl, sh = darkSurfaceToMassRatio(zenith, radius, mass, shape)

        # Take the Mach and Reynolds numbers in the air of OpenDarkflight
        rho_air, sound_speed, temp = atm.model(ht)
        mach = vel_air/sound_speed
        reynolds = reynoldsNumber(rho_air, vel_air, volumeEquivalentDiameter(radius, shape), temp)
        gamma = drag.computeDragCoeff(mach, reynolds=reynolds, knudsen=knudsenNumber(mach, reynolds), \
            vel=vel_air)

        # The area facing the flow along the line, as OpenDarkflight combines its areas
        return gamma*np.hypot(sl*np.sin(zenith), sh*np.cos(zenith))

    return kFunc



def runJobs(jobs):
    """ Fly the fragments of a jobs file.

    Arguments:
        jobs: [dict] The 'line' (see lineAt()), the 'profile' ({'path', 'type'}) of the atmosphere, needed
            for the drag models of OpenDarkflight, and the 'flights', each with 'direction' ('down' or 'up'),
            'len_start' (m), 'vel_start' and 'vel_end' (inertial, m/s, the latter down only), 'mass' (kg),
            'sigma' (s^2/m^2),
            'zenith' (deg), and either 'k' (m^2/kg, fixed at that mass) or 'meteorite' (see
            odfBallisticCoeff()).

    Return:
        [list] For every flight, the length, speed, time and mass where it ends, whether it got there ('ok')
            and why not ('status': 'line_end', 'faster_at_start' or 'not_finite'), and K where it starts
            ('k_start') and ends ('k_end').
    """

    line = {key: np.array(value) for key, value in jobs['line'].items()}

    atm = None
    results = []
    for flight in jobs['flights']:

        zenith = np.radians(flight['zenith'])
        if flight.get('k') is not None:
            kFunc = fixedBallisticCoeff(flight['k'], flight['mass'])

        else:

            # Load the atmosphere of OpenDarkflight once, from its own directory
            if atm is None:
                sys.path.insert(0, os.getcwd())
                from opendf.Routines.Atmosphere import AtmosphereModel, WindMeasurements
                atm = AtmosphereModel(WindMeasurements(jobs['profile']['path'], \
                    wind_data_file_type=jobs['profile']['type']))

            kFunc = odfBallisticCoeff(flight['meteorite'], atm, zenith)

        if flight['direction'] == 'down':
            len_end, vel, time, mass_end, ok = flyDown(line, kFunc, flight['len_start'], \
                flight['vel_start'], flight['vel_end'], flight['mass'], sigma=flight.get('sigma', 0.0))
        else:
            len_end, vel, time, mass_end, ok = flyUp(line, kFunc, flight['len_start'], flight['vel_start'], \
                flight['mass'], sigma=flight.get('sigma', 0.0))

        # K where the flight starts and ends, with the air speed there
        k_start, k_end = np.nan, np.nan
        if np.isfinite(flight['len_start']) and np.isfinite(flight['vel_start']):
            ht_start, _, offset_start, _, _ = lineAt(line, flight['len_start'])
            k_start = kFunc(flight['mass'], ht_start, flight['vel_start'] - offset_start)
        if np.isfinite(len_end) and np.isfinite(vel) and (vel > 0):
            ht_end, _, offset_end, _, _ = lineAt(line, len_end)
            k_end = kFunc(mass_end, ht_end, vel - offset_end)

        # Tell why a flight did not get there
        status = 'ok'
        if not ok:
            status = 'line_end'
            if not (np.isfinite(flight['len_start']) and np.isfinite(flight['vel_start'])):
                status = 'not_finite'
            elif (flight['direction'] != 'down') and (flight['vel_start'] >= lineAt(line, \
                flight['len_start'])[3]):
                status = 'faster_at_start'

        results.append(dict(len_end=len_end, vel_end=vel, time=time, mass_end=mass_end, ok=bool(ok), \
            status=status, k_start=float(k_start), k_end=float(k_end)))


    return results



if __name__ == "__main__":


    ### COMMAND LINE ARGUMENTS

    if len(sys.argv) != 3:
        print("Usage: python LuminousFlight.py jobs.json results.json")
        sys.exit(1)

    with open(sys.argv[1]) as f:
        jobs = json.load(f)

    with open(sys.argv[2], 'w') as f:
        json.dump(runJobs(jobs), f)
