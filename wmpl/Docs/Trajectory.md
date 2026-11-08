# Solving a trajectory

**The trajectory solver** turns multi-station meteor observations into a trajectory and a heliocentric
orbit, and saves the result as a pickle file. That pickle is the input to most of the rest of the
library, including [the REBOUND integrator](REBOUND.md), so this page is mainly about how to produce
one and what is inside it.

## Table of Contents

- [Running the solver](#running-the-solver)
- [The trajectory pickle](#the-trajectory-pickle)
- [Loading a pickle in Python](#loading-a-pickle-in-python)
- [The gravity drop](#the-gravity-drop)
- [Drag against the rotating air](#drag-against-the-rotating-air)

---

## Running the solver

The solver is driven through the input format of your data, not through a single generic entry point.
Each format module under `wmpl/Formats/` has its own command line interface that reads the
observations, runs the solver and writes the results:

```
python -m wmpl.Formats.ECSV /path/to/event/*.ecsv
python -m wmpl.Formats.RMSJSON /path/to/event/*.json
python -m wmpl.Formats.CAMS /path/to/event/
python -m wmpl.Formats.Met /path/to/event.met
```

Run any of them with `--help` for the full list of options. The ones that matter most often:

| Flag | Description |
| :--- | :--- |
| `-s`, `--solver` | Which solver to use. `original` (default) is the Monte Carlo solver; `gural0`–`gural3` are the Gural multi-parameter fits with constant, linear, quadratic and exponential deceleration. |
| `-r`, `--mcruns` | Number of Monte Carlo runs used to estimate the uncertainties. |
| `-d`, `--disablemc` | Only compute the geometric solution, skipping the Monte Carlo step. Fast, but gives no uncertainties. |
| `-t`, `--maxtoffset` | Maximum timing offset allowed between stations, in seconds. |
| `-v`, `--vinitht` | Estimate the initial velocity as the average above this height, in km. |
| `--vinitdrag` | Estimate the initial velocity from a single-body drag and ablation fit to the points selected by `--vinitdragtime` and `--vinitdraght`, for meteors that already decelerate in their first part (e.g. fireballs first seen below 60-70 km), where the default straight-line fit underestimates it. Erosion of the body is absorbed by the fitted ablation coefficient (sigma + eta), but a significant wake of eroded grains shifts the measured centroids behind the body and can bias the velocity by hundreds of m/s; see `wmpl.Utils.DragInitialVelocity`. The straight line is kept if the fit fails, does not fit better or measures no deceleration (then it is only a straight line over its points). The fit estimates the time offsets of the stations again, as the solver's come from its own lag model, except those given as fixed (`--fixedtimes`), which it keeps. The report gives the dynamic pressure and the energy received over the fitted part, with a note when they cross typical fragmentation pressures (0.04-0.12 and 0.5-5 MPa, Borovicka et al. 2020) or erosion onset energies (1-2 MJ/m^2, Buccongello et al. 2024), a cue to check the light curve. Off by default; the solver then warns, in the output and the report, when a parabola over the straight line's points puts the initial velocity more than 2 sigma above it. |
| `--vinitdragtime` | Only points within this time from the first point, in seconds, are used in the `--vinitdrag` fit, 1 s by default, or no time limit if only `--vinitdraght` is given; with both, the fit ends at whichever is reached first. The fit does not model fragmentation, so the fitted part should end before the first one; but the longer the part that can be trusted to be free of fragmentation, the better the velocity (median uncertainty 97, 33 and 9 m/s fitting 0.5 s, 1 s and all points in synthetic tests). Extend it as far as the light curve shows no flare. |
| `--vinitdraght` | Only points above this height, in km, are used in the `--vinitdrag` fit, e.g. above a flare in the light curve or where the report notes typical fragmentation pressures. Given alone, it is the only limit; with `--vinitdragtime`, the fit ends at whichever is reached first. No limit by default. |
| `-l`, `--plotallspatial` | Save the full set of diagnostic plots. |

To see the solver working on a synthetic example without any data of your own:

```
python -m wmpl.Trajectory.Trajectory
```

### Additional fragments in ECSV files

ECSV files can describe more fragments than the main one, as in Appendix 4 of the GFE standard: one row
per frame, with the columns of each additional fragment carrying its number as a suffix (`datetime1`,
`ra1`, `dec1`, `azimuth1`, `altitude1`, ...), empty on the frames where it was not measured. The main
fragment, fragment 0, is always solved from the columns without a suffix (or with the optional suffix 0),
leaving out its points with `trajectory_use` = False (e.g. saturated flares kept in the file).

Points of the main fragment flagged with `flare` = True are used in the trajectory like the others,
unless they are also flagged with `trajectory_use` = False. Either way, the solver reports the time and
height of every flare point in the terminal and at the end of the saved report. A point used in the
trajectory is projected on it. A point not used is only a time, as where it was clicked is not a
measurement (e.g. a saturated flare in which the meteor cannot be seen), and takes the height of the
trajectory at that time. The report gives the flares of every station, with the uncertainties of their
times and heights, and how the flares of different stations compare, as they need not be the same flare
(two cameras can see a fireball flaring at different times). The flares of a station are its flare points
in consecutive rows of its file, among the points of the main fragment: to separate two flares, leave a
point between them that is not flagged. This works in the same way for any frame rate, and when only every
second frame or fewer were picked. Two flares on consecutive points cannot be told apart, but two cameras
could not tell them apart either, as the comparison below allows a frame between the same flare. The
uncertainties come from the Monte Carlo runs and from the duration of the frame the flare was seen in, as
it happened at any time within it. The duration of a frame comes from the frame numbers of the file, if it
gives them, so it is that of one frame also when not every frame was picked. The times are relative to the
clock of the reference station, whose time offset is not estimated, so the Monte Carlo runs give it no
uncertainty in time, but the frame does.

Every flare also gets the speed of the meteor relative to the ground at its time, from the drag fit of the
trajectory (`wmpl.Utils.DragFit`, with the Earth's rotation along the path taken off its inertial speeds),
and the dynamic pressure p = rho*v^2 with the NRLMSISE air density at its height, with uncertainties from
the Monte Carlo runs and the frame as for the height. The error of the air density model is not included.
On fireballs simulated with MetSim, the drag fit gave the dynamic pressure at the fragmentations within
0.7 % in 90 % of the cases (0.2 % in the median), and better than a robust alpha-beta fit late in the
flight, where that one gave it about 4 % too high.

The noise of the Monte Carlo runs only moves the radiant, so the runs alone do not change the speed: refitting
the drag model on every run gave 2 to 3 m/s of scatter in these fireballs, against a real error of about
30 m/s. As `Trajectory.sampleMCVelocities()` does with the initial velocity, every run therefore shifts the
speed within its uncertainty from the covariance of the drag fit (`DragFit.motionStdAtTime()`), and the
pressure of the run follows. Near the fragmentations of the simulated fireballs that uncertainty, about
30 m/s, held the real error within 1 sigma in 50 to 55 % of the cases and within 2 sigma in 77 to 82 %,
so it is somewhat small: the bias of the fit where the meteor fragments is not included. Its share of the
uncertainty of the pressure, about 0.3 %, is small next to that of the frame.

When the solver estimated the initial velocity from a drag fit of the first part of the trajectory
(`--vinitdrag`, `wmpl.Utils.DragInitialVelocity`), which models a single body without fragmentation, the
report warns if a flare point lies at or above the lowest fitted point, giving the highest flare point of
every station, and suggests a `--vinitdraght` that keeps only the points above the highest flare.

Two cameras can mark the same flare a frame apart, e.g. as it lasts longer than a frame or peaks between
two, so the comparison takes one frame of the slower camera off the time between two flares before it
weighs it against its uncertainty. In the simulated fireballs, with the frames of each camera at a random
phase and each marking the frame with the most light, 4.9 % of the same flares of a 25 and a 13 fps camera
were further apart than 3 sigma of the duration of their frames, and 2.7 % with that frame taken off (2.1 %
and 1.1 % with two 25 fps cameras). When the solver estimated the time offsets, the uncertainty of the time
between two stations also includes their real error, 7 ms (68 % of the errors in the same simulations, with
a median of 4.4 ms), against about 0.24 ms from the fit and the Monte Carlo runs, which do not predict it.
It also enters the heights compared at the same time. It changed the rates above by 0.2 to 0.3 points, with
cameras of 25 fps; with faster cameras the frames weigh less and it matters more.

The report then proposes the flares of the event. Two flares of different stations are taken as the same
flare when they are consistent and each one is the only consistent flare of its station for the other, and
the joined flares make one flare of the event. A flare consistent with two or more flares of another station
is left alone as ambiguous, and so are the flares of a group that would hold two flares of one station. Every
flare of the event gives the values of its first point in time among its stations, how many frames apart its
stations saw it begin and how far apart their pressures are, which says how much joining them matters.
When two of the flares an ambiguous flare could be are flares of one station with a single point not flagged
between them, the report says so: a frame left unflagged inside a flare splits it in two.
Without the Monte Carlo runs the flares are not joined. In the simulated fireballs, with every camera flagging
the frames in which the light of the dust of a fragmentation exceeds the rest by 0.25 mag, this joined no
flares of different fragmentations. Of the fragmentations seen by two or more stations whose flares held only
that fragmentation, 92 % were joined and 8 % left ambiguous; windows of 3 and 5 frames instead of 3 sigma left
14 and 26 % ambiguous. A window of one frame did better (99 %) when every frame of a flare is flagged, but
missed more of the same flares marked only on their brightest frame. A third of the fragmentations seen by two
or more stations had a flare that already held another fragmentation, a median 85 ms apart, which no
comparison of stations can split.

Every flare of the event takes the values of its first point in time among its stations. A camera flags a
flare once the light it integrates over a frame stands out, on average after the fragmentation, and the
earliest camera is the least late. In the simulated fireballs (threshold 0.25 mag, 54 957 joined flares), with
the times compared to that of the fragmentation and the pressures to the true pressure of the main at the time
taken:

| Time of the joined flare | Bias in time | Bias in pressure | Error of the pressure, median / 90th pct. |
|---|---|---|---|
| first point of any camera (used) | -2 ms | -0.3 % | 1.5 % / 4.3 % |
| mean of the beginnings of the stations | +10 ms | +1.7 % | 1.8 % / 5.1 % |
| one station | +10 ms | +1.7 % | 2.0 % / 5.8 % |
| mean of the centres of the flares | +35 ms | +6.4 % | 6.4 % / 31 % |

The centres of the flares come after the fragmentation, as its dust shines for several frames. With a
threshold of 0.5 mag the first point does better still (1.5 % / 4.6 % against 2.6 % / 6.6 % for the mean of
the beginnings). With 0.1 mag, which flags almost any brightening, the mean of the beginnings does slightly
better (1.4 % / 5.9 % against 1.7 % / 7.0 %), as the first point is then early. The simulated frames are timed
at the middle of their exposure, and RMS times the frames of FF files at their beginning, but the flares and
the trajectory share those times and the solver fits the time offsets of the stations, so the heights and
pressures follow the position of the meteor in the frame either way.

```
python -m wmpl.Formats.ECSV /path/to/event/*.ecsv --fragments
```

| Flag | Description |
| :--- | :--- |
| `--fragments` | After the main fragment, solve every additional fragment seen from at least 2 stations with the same options (Monte Carlo included), save each one in a `fragment_k` folder of the output directory, and print how it differs from the main one: convergence angle of its stations, radiant, velocities, heights, distance from the main path and how far ahead of the main fragment it is at the same time. |
| `--fragtimefit` | Estimate the station time offsets for every fragment, instead of reusing those of the main trajectory, which its larger number of points constrains better. Without the timing estimation (`-e`), the fragments take the same fixed offsets as the main fragment. |
| `--fragnosys` | Leave out of the Monte Carlo uncertainties of the fragments the systematic errors of the stations, which are otherwise included when the main trajectory has 3 stations or more (see below). |

Two solutions differ slightly even for identical data, as each goes through its own iterations of the
timing and velocity estimation: in tests, by about 0.001 deg and tens of metres with well spread stations,
and up to about 0.1 deg and a few hundred metres with nearly parallel ones. With the Monte Carlo runs, the
solutions compared are still those with the original picks (not the best Monte Carlo runs, which differ by
about 1 sigma even for identical data), and the differences are also given in units of the combined
uncertainty. The time offsets reused by the fragments are also those of the main solution with the original
picks, the ones in its report.

The comparison notes a fragment more than 1.5 km off the main path at its first or last point. No fragment
of the main went further in fireballs simulated with MetSim (2 to 7 fragments, with the lateral speeds of
real fragmentations, stations converging by 20 deg), while fragments of another body 2 km away did in 97 %
of the cases. Such a fragment may come from another body, the stations may have measured different fragments
under the same number, or its trajectory may be poorly constrained: with stations converging by 8 deg, 8 to
14 % of the fragments of the main went further. Two stations measuring different fragments mostly give a line
close to the main path anyway, and often no trajectory at all.

With the time offsets of the main trajectory, the comparison also gives how far apart along the path of the
fragment every two of its stations see it at the same times, from their points on its line, as a constant
difference in a joint fit of their positions in time, and notes more than 1 km. Two stations measuring
different fragments see them where the line passes closest to their lines of sight, at different places
along it. In the simulated fireballs, 0.2 to 0.3 % of the fragments of the main went further, against 70 %
of the fragments whose 2 stations measured different fragments, and 56 % with 3 stations, one of which did.
Time offsets estimated for the fragment (`--fragtimefit`) would absorb the difference, so it is then not
given.

The comparison also gives the convergence angle of the stations of every fragment, against that of the main
trajectory, and notes one below 3 deg, the least the trajectory correlator of WMPL takes. Two stations that
hardly converge for a fragment, e.g. close to each other and looking the same way, still give a line, but one
which the small errors of their picks can turn far from the true path.

The pointing of a station can be off by the same amount at all its points (its calibration, or the centroid
of a saturated meteor), which does not average out with their number, so the Monte Carlo runs do not
include it. With 3 or more stations in the main trajectory, the others constrain its line, and the residuals
of each station across the line show its own error, as a constant and a slope along the trajectory. Their
typical size is added to the Monte Carlo runs of every fragment, as a random constant and slope of each
station in every run. The solution of the fragment is the one with the original picks, as without them, and
only its uncertainties change. On fireballs simulated with MetSim, with 2 to 7 fragments, noise and
saturation, the true direction of a fragment was within the 95 % region of its uncertainties in 93 % of the
cases with them, and in 75 % without. With 2 stations in the main trajectory there is nothing to measure,
and the solver says so.

A fragment seen with enough points (4 or more) from only one station cannot be solved, as one station gives a
plane and not a line, and the trajectory of a fragment seen from more can still fail, or say little: its
stations can converge by less than 3 deg, or it can be more than 1.5 km off the main path (see above). The
points of such a fragment are projected on the main trajectory, each station on its own, at the times of their
frames corrected by the time offset of the station, as the flares are, and the report says why. The terminal
and the end of the saved report give the height and time of every point, how far off the main path it is and
how fast it moves away from it, and how far ahead of the main fragment it is at the same time. This assumes
that the fragment is on the main path, which one station cannot test. Only the part of a separation across the
line of sight is seen, so the distances off the path and the speeds are lower limits. The rest moves the point
along the line of sight, which the report gives as the error of the height, and of how far ahead it is, per
100 m off the path. The uncertainties come from the Monte Carlo runs of the main trajectory, with the scatter
of the points of the fragment as their noise, and do not include a separation from the main path. A station
which is not in the main trajectory is skipped, as its time offset is unknown.

When several fragments are close together, two stations can measure them under crossed numbers, which their
trajectory hardly shows (see above). So the fragments of every two stations are also matched by where they
are along the main path at the same times, whatever their numbers in the files. The points of every station
are projected on the main trajectory on their own, as above, and the distance of two fragments is how far
apart along the path their projections are, as a constant difference in a joint fit of their positions in
time. Pairs more than 1.5 km apart are not matched, and the pairs are chosen by the least total distance
(Hungarian algorithm). The end of the saved report lists the matches of every two stations with the distance
of the closest other candidate of either fragment (which can be closer than the match, as the pairs are
chosen together), calls a match ambiguous when it is less than 0.5 km closer than that one, and notes a
match of different numbers, and a fragment more than 1.5 km from the one with its number at the other
station. The numbers in the files are not changed. On fireballs simulated with MetSim (2 or 3 stations, 2
to 7 fragments, noise and saturation, 100 events per case), 96 to 99 % of the true pairs were matched. The
1 to 3 % of the matches which joined different fragments were all ambiguous, and 0.2 % of the true pairs
were more than 1.5 km apart. As the matching does not use the numbers, crossed numbers are noted without
ambiguity about two times out of three, and as ambiguous in most of the rest.

## The trajectory pickle

A successful run writes `<event>_trajectory.pickle` (plus plots and a text report) into the output
directory. The pickled object is a `wmpl.Trajectory.Trajectory.Trajectory`. The attributes used by the
other tools are:

| Attribute | Meaning |
| :--- | :--- |
| `traj.jdt_ref` | Reference epoch of the trajectory, as a Julian date in TDB. |
| `traj.traj_id` | Identifier of the event, e.g. `20191023_091225`. Used as the object name downstream. |
| `traj.orbit` | The computed heliocentric orbit, including the geocentric radiant and velocity. |
| `traj.state_vect_mini`, `traj.v_init`, `traj.radiant_eci_mini` | Position, speed and radiant direction at the reference point. Together these are the state vector REBOUND starts from. |
| `traj.uncertainties` | One-sigma uncertainties from the Monte Carlo runs, used to draw clones. |
| `traj.observations` | The per-station observations that went into the fit. |

## Loading a pickle in Python

```python
import os
from wmpl.Utils.Pickling import loadPickle

pickle_path = "/path/to/20191023_091225_trajectory.pickle"
traj = loadPickle(*os.path.split(pickle_path))

print(traj.traj_id, traj.jdt_ref)
print("a = {:.4f} AU, e = {:.4f}".format(traj.orbit.a, traj.orbit.e))
```

Once you have the pickle, see [the REBOUND manual](REBOUND.md) for integrating the orbit backwards or
forwards in time.

## The gravity drop

With the gravity correction on (the default, `-g` turns it off), the solver fits the line of the trajectory
together with the drop of the meteor below it caused by gravity, from the first point. The line is the
direction of motion at the first point.

Drag also slows down the velocity which gravity adds across the path, like the rest of the velocity, so a
decelerating meteor drops less than one at a constant velocity. Its path bends at `dtheta/dt = g_perp/v`, as
in [the drag fit](DragFit.md). This depends only on the velocity: the size, shape and density of the meteoroid
and the density of the air act on the drop only through the deceleration they cause. The solver scales the
drop at a constant velocity by the ratio of the two. In the pass with the drop, it first fits the line with
the drop at a constant velocity, fits the velocity to the lengths of the points on that line, and fits the
line again with the scaled drop. The lengths on the straight line of the first pass are not used, as they also
hold the speed gravity adds along the path, which does not bend it: without drag, the drop across the first
direction is `g_perp*t^2/2` whatever the speed along it.

The velocity is fitted without a model of the deceleration, of the atmosphere or of the heights:

- It is constant, or linear in time over up to 8 segments joined end to end, and the Bayesian information
  criterion (Schwarz 1978) chooses the number of segments. Points which show no change of velocity keep a
  constant velocity, the drop is not scaled, and the solution is the one with the drop at a constant velocity.
  A deceleration of any shape, such as that of a fireball which fragments, gets the segments it needs.
- The velocity does not increase faster than gravity can speed a meteor up (9.83 m/s^2), as drag only slows
  a meteor down. This keeps the fit from following the scatter of few or noisy points, while it follows a
  meteor which gravity alone speeds up or slows down.
- It is held above 3 km/s, where ablation stops.

For Winchcombe (8.2 s, 13.8 to 3.6 km/s), the drop at the last point is 189 m instead of 245 m, the same as
from the bending of the drag fit, and along the whole path the two drops differ by at most 0.9 m. The end
height moves up by 41 m and the radiant by 27 arcsec. On fireballs integrated in 3D with drag and gravity,
which slow down to 4-10 km/s, the radiant error of a noise-free solution falls from 5-65 arcsec to 0.2-0.8
arcsec, and the largest error of the positions across the path from 4-71 m to 0.1-1.3 m.

The meteor falls from the line along the whole path, while the vertical turns with the curvature of the
Earth, so dropping every point by `g*t^2/2` along its own vertical tilts the line, by about `0.08*t^2`
arcsec after `t` seconds. The solver drops the points by the free fall of the meteor from the line under the
gravity of the Earth instead, a point mass with its oblateness (J2 of WGS84, around the axis of the ECI frame
of date):

    displacement(t) = integral from 0 to t of g(L(q)) (l(t) - l(q))/v(q) dq

where `L(q)` is the point of the line at the time `q`, `l` the length travelled and `v` the velocity of drag
alone, the one fitted above. At a constant velocity this is the free fall from the line, and in a uniform
field it is `g` times the scaled drop. Gravity is interpolated by a cubic in the distance along the line, and
the velocity enters through integrals of the motion tabulated from the fit, so the displacement is exact where
the fitted velocity changes its slope. It follows fireballs integrated in 3D with drag and gravity to 0.5-2.4
cm, against 0.2-2.1 m for the scaled drop along the vertical. On the same fireballs, without noise, the
radiant error falls from 0.6-3.6 arcsec to 0.00-0.18 arcsec and the largest error across the path from
0.08-1.3 m to 0.01-0.26 m. The radiant of Winchcombe moves by 7.7 arcsec, where the two drops differ by up to
2.8 m across the 97 km of the path. The oblateness accounts for 1 arcsec of it: on a sphere the drop across
its path differs by up to 0.45 m.

- `Trajectory(..., gravity_decel=False)` keeps the drop at a constant velocity. Trajectories pickled before
  the option was added keep it too.
- `Trajectory(..., gravity_free_fall=False)` drops every point along its own vertical, as before.
  Trajectories pickled before the option was added keep it too. `traj.gravityFreeFall()` gives the table of
  the motion along the line for the free fall, and `traj.gravityDropPoint()` drops a point of the line as the
  solver does.
- `traj.gravity_drop_scale` holds the times from the first point and the scale at those times, and
  `traj.gravityDropFactor(t)` gives the factor of the drop at the time `t` from the first point. The drag
  fit and the dynamic mass fit apply the same drop as the solver with it.
- With few points, gaps between the stations and a large scatter of the lengths, the points may only support
  a constant deceleration, which then stands for a deceleration concentrated at the end. On such a fireball
  (3.4 s, 57 points, 400 m of scatter), the drop differs from that of the drag fit by up to 2.7 m out of 53 m.
- As in the drag fit, the velocity is that of the solver, not relative to the air.

## Drag against the rotating air

The air rotates with the Earth, eastwards at up to 465 m/s at the equator and about 330 m/s at 45 deg, in the
ECI frame the solver works in. Drag acts against the velocity of the meteor relative to the air, not against
its ECI velocity, so it pulls the velocity of a decelerating meteor across its path towards the velocity of
the air. In ECI coordinates a decelerating meteor then does not move along a straight line, even without
gravity: the fireballs below, which slow down to 4-10 km/s, are pulled by 80-600 m. A meteor at a constant
velocity is not pulled.

With the gravity correction on, the solver also moves the points of its line by this displacement. The
velocity across the line, `p`, follows `dp/dt = (a/w)*(u_perp - p)`, where `a` is the deceleration by drag,
`w` the speed relative to the air and `u_perp` the velocity of the air across the line. It is 0 at the first
point, so the line stays the direction of motion at the first point in ECI coordinates, which the orbit is
computed from. Like the drop, it depends only on the velocity, which comes from the same fit to the lengths
of the points, without a model of the deceleration. These are the lengths on the line of the pass with the
drop, which do not hold the speed gravity adds along the path, so a velocity changed by gravity alone does not
move the points. The velocity of the air is taken at the points along the line, which changes the
displacement by 0.9-13 m in the fireballs checked.

In the pass with the drop, after fitting the line with the scaled drop, the solver computes the displacement
from the lengths of the points on that line and fits the line again with it. The lengths on a line depend on
the displacement, as moving the line across also moves the points of the lines of sight closest to it along
the line, so this is repeated until the displacement changes by less than 0.1 m, at most 3 times. On the
fireballs below, the displacement from the first fit only gives radiant errors of up to 15 arcsec and errors
of the positions of up to 9.4 m, against 7.4 arcsec and 3.8 m. The median radiant error stays at 2.4 arcsec,
and the median error of the positions falls from 0.9 to 0.3 m.

- On 20 fireballs integrated in 3D with gravity and drag against the rotating air (5 fireballs of 14-25 km/s
  slowing down to 4-10 km/s in 3.3-9 s, each coming from 4 azimuths), the radiant error of a noise-free
  solution falls from 110-735 arcsec (median 360) to 0.5-7.4 arcsec (median 2.4), and the largest error of
  the positions across the path from 56-390 m to 0.1-3.8 m. The same fireballs integrated with drag against
  the ECI velocity are solved without the displacement to 0.5-3.5 arcsec and 0.05-1.2 m. The displacement
  matches the deviation of the integrated paths to 0.02-0.26 m out of 100-600 m, and the rest comes from the
  velocity fit and the drop.
- With 10 and 30 arcsec of noise (4 of the fireballs, 5 seeds each), the radiant errors are the same as those
  of the fireballs integrated without the rotating air and solved without the displacement, with medians of
  5-10 and 15-38 arcsec, against 170-540 and 220-570 arcsec without the displacement.
- For Winchcombe (8.2 s, 13.8 to 3.6 km/s), the points move by up to 237 m, the radiant by 154 arcsec and the
  end height up by 147 m, and the RMS of the angular residuals falls from 153 to 75 arcsec.

- `Trajectory(..., air_rotation=False)` keeps the dropped line. Trajectories pickled before the option was
  added keep it too, and so do solutions without the gravity correction.
- `traj.air_rotation_shift` holds the times from the first point and the displacements at those times, and
  `traj.airRotationShiftAt(t, radiant_eci)` gives the displacement at the time `t` from the first point,
  across the given line. The drag fit and the dynamic mass fit apply the same displacement as the solver.
- With the timing offsets estimated, the displacement also moves the lengths, so the offsets and the initial
  velocity can change too, when the stations see different parts of a fireball.
- A decelerating meteor needs two more fits of the line in each pass with the drop, and a better fit can
  reject more picks and solve the trajectory again: Winchcombe takes 4.2 s instead of 1.4 s.
- Winds, of tens of m/s at the heights of meteors, are not included. They would add to the velocity of the
  air.
