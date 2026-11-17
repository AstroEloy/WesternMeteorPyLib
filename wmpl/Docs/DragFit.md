# Drag fit of the motion along a trajectory, and sampling it

**`wmpl.Utils.DragFit`** fits a physical drag model to the lengths and times of a solved
trajectory. It gives smooth positions, times and velocities along the whole trajectory, extrapolates
them above the beginning and below the end, and shows how the ballistic coefficient of the meteor
changed during the flight. **`wmpl.Utils.SampleTrajectoryPositions`** uses it to sample the trajectory
by height or by time, which is the usual first step of a dark flight model.

The ballistic coefficient is fitted from the dynamics alone, unlike in
[`wmpl.Utils.DynamicMassFit`](../Utils/DynamicMassFit.py), which computes a mass from the drag equation
for an assumed shape and density.

## Table of Contents

- [Quick start](#quick-start)
- [The model](#the-model)
- [Command line reference](#command-line-reference)
- [Reading the output](#reading-the-output)
- [Using it as a library](#using-it-as-a-library)
- [Interpreting the ballistic coefficient](#interpreting-the-ballistic-coefficient)
- [The dynamic mass](#the-dynamic-mass)
- [Limits of the model](#limits-of-the-model)
- [Troubleshooting](#troubleshooting)

---

## Quick start

Fit the model to a trajectory, print the table of the ballistic coefficient and show the plots:

```
python -m wmpl.Utils.DragFit /path/to/20191023_091225_trajectory.pickle
```

Sample the trajectory every 0.5 km between the observed beginning and end heights (`-1` means the
observed height), or between 30 and 15 km, below the end, for a dark flight model:

```
python -m wmpl.Utils.SampleTrajectoryPositions traj.pickle -1 -1 0.5
python -m wmpl.Utils.SampleTrajectoryPositions traj.pickle 30 15 0.1
```

Sample by time instead, every 0.05 s from 0 to 2 s after the reference time of the trajectory:

```
python -m wmpl.Utils.SampleTrajectoryPositions traj.pickle 0 2 0.05 --time
```

Save the samples to a CSV file, with all the sampled quantities:

```
python -m wmpl.Utils.SampleTrajectoryPositions traj.pickle 30 15 0.1 --csv samples.csv
```

An example trajectory is included in the repository:
`wmpl/Dynesty/examples/20191023_091225/20191023_091225_trajectory.pickle`.

---

## The model

The meteor moves along the solved trajectory and is decelerated by atmospheric drag:

```
dv/dt = -K*rho*v_air^2 + g_along
```

where `rho` is the NRLMSISE-00 air density along the path, `g_along` the component of gravity along
the path, `v_air` the speed relative to the air, which moves with the ground (the inertial speed `v`
minus the speed of the Earth's rotation along the trajectory), and `K` the ballistic coefficient,

```
K = Gamma*A / (m^(1/3) * rho_m^(2/3))
```

with the drag coefficient `Gamma`, the shape factor `A`, the mass `m` and the bulk density `rho_m`.
`K` has units of m²/kg. It is the drag area per unit mass, `K = C_D*S/(2m)`.

- **The fit.** `ln K` is a sum of ramps, one per segment of the observed part, and is only allowed to
  grow, as ablation and fragmentation both increase it. The equation is solved exactly in the specific
  kinetic energy, and the model is fitted to the lengths along the track at the observed times of all
  used points (ignored points and stations are skipped), with a robust loss and a small penalty on
  the growth of `K`.
- **Above the beginning** the meteor moves without drag, only accelerated by gravity, as it barely
  decelerates in the thin air. With `--dragabove` the ballistic coefficient at the first point is
  kept there instead. As `K` can only have been smaller above, this is the largest possible
  deceleration, and the model then stops where the velocity would have grown by 10%.
- **Below the end** `K` is kept at its last value, which is drag-only flight with gravity, as at the
  beginning of the dark flight. Ablation essentially stops below about 3 km/s, so below that speed
  this is exact.
- **Gravity** bends the path by the angle `theta`, `dtheta/dt = g_perp/v`. Drag damps the velocity
  gained from gravity like the rest of the velocity, so the bending grows quickly as the meteor slows
  down. The model is computed until the path has bent by 5 deg from the first point. Further on, a
  dark flight model has to be used.
- **Gravity drop.** If the trajectory was solved with the gravity correction, the positions in the
  observed part are dropped from the line in the same way as in the trajectory solver, which scales the
  drop for the deceleration ([the gravity drop](Trajectory.md#the-gravity-drop)). Below the end
  the path continues along the direction bent by `theta`, and the radiant of the samples is bent by
  `theta`.
- **Heights** are above the mean sea level, as everywhere in WMPL.

On 338 trajectories of real fireballs, the fit residuals are 173 m (median), and hiding the first or
last 20% of the points and predicting them from the rest gives errors of 233 m and 463 m.

---

## Command line reference

`python -m wmpl.Utils.DragFit traj.pickle [options]`

| Flag | Description |
| :--- | :--- |
| `-k`, `--knots` | Number of segments of the ballistic coefficient. One for every 25 points by default, between 3 and 15. |
| `-a`, `--dragabove` | Keep the ballistic coefficient of the first point above the beginning, instead of moving without drag. |
| `-g`, `--ga` | The product of the drag coefficient Gamma and the shape factor A, assumed for the dynamic mass. 0.78 by default. |
| `-d`, `--dens` | The bulk density assumed for the dynamic mass (kg/m³). 3500 by default. |
| `-f`, `--fragheights` | Heights of the fragmentations (km), e.g. those of the flares, where `K` can jump without penalty. The spread of the fragmentations is then scanned, and the uncertainty of `K` and of the dynamic mass in the table and the plot includes it. |
| `--atm_profile` | Atmosphere profile file to take the air density and the winds from (experimental), read as OpenDarkflight does. Above and below the heights it covers, the MSIS model is scaled to join it. The wind along the direction of motion changes the speed relative to the air. |
| `--atm_profile_type` | Format of the profile: `wrf` (default, the CSV of every model profile), `wyoming` or `supracenter`. |
| `--no_winds` | Ignore the winds of the profile. |

`python -m wmpl.Utils.SampleTrajectoryPositions traj.pickle beg end step [options]`

| Argument | Description |
| :--- | :--- |
| `beg`, `end`, `step` | Heights of the first and last sample and the step (km). `-1` uses the observed beginning or end height. With `--time` they are times relative to the reference time of the trajectory (s), and `-1` has no special meaning. |
| `-t`, `--time` | Sample by time instead of height. |
| `-f`, `--fragheights` | As above: the samples follow the fit with the jumps at the fragmentations. |
| `-a`, `--dragabove` | As above. |
| `-c`, `--csv` | Save the samples to the given CSV file. |
| `--atm_profile`, `--atm_profile_type`, `--no_winds` | As above. |

### Handing the fit over to OpenDarkflight (experimental)

Which tool to use, from the solved trajectory to the strewn field:

| Step | Tool | What it gives |
| :--- | :--- | :--- |
| Motion and `K` along the trajectory | `wmpl.Utils.DragFit` (`-f` with the heights of the flares, `--atm_profile`) | `K` and its uncertainty, the final dynamic mass for an assumed shape and density |
| Positions, speeds and `K` at given heights or times | `wmpl.Utils.SampleTrajectoryPositions` (same options) | The samples, e.g. to start a dark flight by hand |
| Dark flight of the leading body with its uncertainties | `wmpl.Utils.DragFitDarkFlight` | OpenDarkflight inputs from the last observed point, the falls of the clones and of the fragments of the last fragmentation |
| Falls of the fragments which broke off at given heights | `wmpl.Utils.DragFitFragmentFalls --direct` | The falls by breakup and mass, with the kicks and the wind errors |
| Where meteorites found broke off | `wmpl.Utils.DragFitFragmentFalls --inverse` | The breakup of each, with its range (with ablation, also the start of the dark flight) |

The dark flight does not need a mass, nor `wmpl.Utils.DynamicMassFit`: OpenDarkflight takes `K` and gives the
meteorite the mass of that `K` for its own shape and density.

`python -m wmpl.Utils.DragFitDarkFlight traj.pickle out_dir -f HT [HT ...] --atm_profile profile.csv [options]`

The dark flight starts at the last observed point, where the drag fit already flies on with a constant `K`.
OpenDarkflight takes `K` directly (`ballistic_coeff`, its pull request #45) and gives the meteorite the mass
of that `K` for its shape and density. Both take the density and the winds of the same profile, which has to
cover the starting height of the leading body (the hand-over stops otherwise, and `--no_winds` is refused), so
the dark flight starts in the air the drag fit ended in. The trajectory has to be solved with its orbit, which
gives the rotation of the Earth along it. The speed handed over is the speed relative to the ground the drag
of the fit acts on, so that `K` gives the same deceleration on both sides, and the summary compares the speeds
of the two models 0.25 to 1 km below the hand-over, which have to agree. The vertical wind of a profile is
only used by OpenDarkflight. WMPL does not import OpenDarkflight: it writes its inputs and, with `--run`, runs
it from its own environment and reads its outputs.

The uncertainty of `K` at the end is dominated by how spread the last fragmentation was. Every spread within
`Delta chi2 <= 4` is equally likely, and every clone is one of those fits with its parameters drawn from their
covariance, so `K` and the speed keep their correlation (a larger `K` comes with a lower speed) and the
one-sided range the spread gives. On the synthetic sets of the knots, this mixture covered the true `K` in
82-85% of the events within its central 68% and in 90-95% within its central 95%. The parameters are drawn as a
linear propagation of the noise of the points, without the bounds of the fit, which clipping would bias. The
falls are measured from the nominal state, the fit with jumps (the `K` of DragFit), or the best spread if the
jumps do not fit the data.

Every clone also gets its own line, which moves its starting point across the line and turns its direction
together: the line of one of the Monte Carlo trajectories of the solution if it was solved with the Monte
Carlo (the `<name>_mc_uncertainties.pickle` the solver saves in the `Monte Carlo` directory next to the
trajectory pickle, found by default), or else a line drawn from the covariance of the line-of-sight fit. Only
the line is taken from the Monte Carlo, as the drag fit gives the motion along it. Every clone can also get a
shape and a density drawn for its meteorite, and with `--mc_runs` its own Monte Carlo of OpenDarkflight adds
its statistical wind errors (15% in speed and 5 deg in direction by default, set with `--odf_option`). The
falls are then those of the Monte Carlo runs; the unperturbed run of every clone stays in `falls.csv` with
`mc_run` 0. Every OpenDarkflight run with a Monte Carlo uses all the cores, so `--nproc` is 2 by default.

The fragments of the last fragmentation (`--fragments`) start at its height with the state of the fit there
and a `K` from that of the leading body after the fragmentation to `--frag_k_max` times it, log-uniform:
smaller bodies of the shape and density drawn for each, which the leading body they broke from is taken to
share, so that the ratio of their `K` is the ratio of their radii. Their Monte Carlo runs get the sideways kick
of Passey and Melosh (1980) for that ratio of radii. OpenDarkflight does not ablate them, so their `K` stays
that of the breakup. If the fragmentation is above the top of the profile, the fragments are left out, while
the starting point of the leading body has to be covered. For the fragments of any breakup, by mass and
shape, `wmpl.Utils.DragFitFragmentFalls --direct` (below) flies each one from the breakup as well, with its
kick across the motion, while OpenDarkflight draws the direction of its kick over all directions.

| Flag | Description |
| :--- | :--- |
| `-f`, `--fragheights` | Heights of the fragmentations (km), required. |
| `-n`, `--clones` | Number of clones. 200 by default. |
| `--seed` | Seed of the draws. 0 by default. |
| `--mc_uncertainties` | Monte Carlo uncertainties of the solution to draw the lines from. By default the one next to the trajectory pickle, if any. |
| `--no_geometry` | Do not draw a line for every clone and fragment. |
| `-d`, `--dens`, `--dens_sigma` | Bulk density of the meteorite (3500 kg/m³ by default) and its standard deviation, drawn for every clone and fragment (0 by default). |
| `--shapes`, `--drag_model` | Shapes of the meteorite for OpenDarkflight, one drawn for every clone and fragment, and its drag model. Those of OpenDarkflight by default. |
| `--fragments`, `--frag_k_max` | Number of fragments of the last fragmentation (0 by default) and the largest `K` of a fragment over that of the leading body (10 by default). |
| `--end_ht` | Height where the dark flight ends (km). |
| `--mc_runs` | Monte Carlo runs of OpenDarkflight for every clone and fragment, with its statistical wind errors, the kicks of the fragments and anything set with `--odf_option`. |
| `--odf_option` | Further OpenDarkflight option, e.g. `'Monte Carlo.wind_speed_sigma=0.15'`. Can be repeated. |
| `--run` | Run OpenDarkflight (`--odf_python` of its environment, from `--odf_dir`, `--nproc` at a time, 2 by default) and collect the falls. |
| `--collect` | Only collect the falls of runs already done; the trajectory pickle is not read. |

The output directory gets `handover.json` (the nominal state, the clones, the fragments, the scan and the
settings), one input per clone under `clones/` and per fragment under `fragments/`, the nominal one under
`nominal/`, and after the runs `falls.csv` and `falls.png`, with the falls along and across the direction of
motion from the nominal fall.

In Python, the same three steps:

```python
import os

from wmpl.Utils.AtmosphereProfile import AtmosphereProfile
from wmpl.Utils.DragFitDarkFlight import findMcUncertainties, loadMcLines, prepareHandover, runOdf, \
    summarizeFalls
from wmpl.Utils.Pickling import loadPickle

traj = loadPickle(*os.path.split("/path/to/traj.pickle"))
mc_lines = loadMcLines(findMcUncertainties("/path/to/traj.pickle"))

# Draw the clones and the fragments and write the inputs of OpenDarkflight, without running anything
handover, ini_paths = prepareHandover(traj, "out_dir", [35000.0, 30000.0], AtmosphereProfile("profile.csv"),
    mc_lines=mc_lines, drag_model='dfn', shapes=['sphere', 'brick'], mc_runs=5, n_fragments=100)

# Run OpenDarkflight from its own environment, or by hand on the inputs
runOdf(ini_paths, "/path/to/envs/opendf/bin/python", "/path/to/OpenDarkflight", n_proc=2)

# Collect the falls, also later with loadHandover("out_dir")
print(summarizeFalls("out_dir", handover))
```

`prepareHandover()` empties the run directories first, as OpenDarkflight names its outputs after the mass and
those of an earlier run would be read again.

### Falls of the fragments, direct and inverse (experimental)

`python -m wmpl.Utils.DragFitFragmentFalls traj.pickle out_dir -f HT [HT ...] --atm_profile profile.csv --odf_python PY --odf_dir DIR (--direct HT [HT ...] --masses KG [KG ...] | --inverse meteorites.csv) [options]`

A fragment breaks off the leading body with its speed and slows down by its own `K`, faster than the leading
body as it is smaller, with the `K` of the drag model and shape of OpenDarkflight at every step (or a fixed
`K`).

**Without ablation** (`--sigma 0`, the default), OpenDarkflight flies the fragment from the breakup, with the
speed and the direction of the leading body there: the whole flight in three dimensions, its path bent by
gravity and pushed by the whole wind, exact within the drag model. A run takes some tens of seconds, as
OpenDarkflight steps 0.1 m in height above 5.8 km/s.

**With ablation** (`--sigma` above 0), which OpenDarkflight does not do, the flight down to where the fragment
has slowed down to a few km/s (`--v_dark`, 3.5 km/s by default) is computed in WMPL (`wmpl.Utils.LuminousFlight`),
in the density and the winds of the same profile, and OpenDarkflight takes it over from there. A small fragment
gets there soon after the breakup, a large one far below it. That flight keeps to the ground track of the
leading body: it takes the gravity and the wind along the motion, but not the bending of the path by gravity
nor the wind across it, so it moves the falls of large fragments (see the validation below). The dark flight
never ablates.

Two shortcuts give wrong falls, and are not offered:

- **Starting the dark flight at the breakup at a few km/s.** The fragment does not slow down on the spot: it
  flies on along the line, kilometres for a large one, before it is that slow. Starting it at the breakup at a
  few km/s moves its fall back along the trajectory by up to several km. The breakup is started at the speed
  of the leading body.
- **Scanning the speed at which the dark flight starts, with ablation.** It does not change the fall, as long
  as the dark flight starts where the fragment has that speed: at 4.5 km/s it starts higher up than at
  2.5 km/s, and flies the rest of the way the same. One speed is enough.

**Validation** on a synthetic fireball (18 km/s, 45 deg, `K` = 1e-3 m²/kg, seen from 82 to 19 km, with a jump of
`K` x3 at 30 km), in an exponential atmosphere with a 20 m/s wind, against OpenDarkflight from the breakup,
for spheres of 3.8 and 60 g breaking off at 30 and 25 km:

| Breakup, fragment | Without ablation, from the breakup | With ablation, along the ground track, dark flight from 2.5 / 3.5 / 4.5 km/s |
| :--- | :--- | :--- |
| 30 km, 3.8 g | 0 / 0 m | +0.7 / +9.9, +0.4 / +6.2, +0.4 / +4.3 m |
| 25 km, 3.8 g | 0 / 0 m | +16 / -8, +12 / -7, +9 / -5 m |
| 30 km, 60 g | 0 / 0 m | +73 / -29, +52 / -21, +39 / -15 m |
| 25 km, 60 g | 0 / 0 m | +85 / -48, +57 / -32, +39 / -22 m |

The differences are along / across the motion, the forced ablation path run with a negligible `sigma`. From the
breakup the run is the reference itself. Along the ground track the large fragments land up to 85 m along and
48 m across off, less for a higher `--v_dark`, while two fragments of similar `K` are moved alike. The inverse
on the reference falls finds the breakup within 5 m without ablation (-4.5 and +0.6 m for breakups between
the heights scanned every 250 m), and 3 to 29 m (3.8 g) and 303 to 306 m (60 g) too high along the ground track
with the flight back up.

**The sideways kick** of the breakup (`--kicks N`, N draws for every fragment or meteorite) follows Passey and
Melosh (1980), Eq. 24: `U = V*sqrt(3/2*C*(R1/R2)*rho_air/rho_met)`, with `C = 1` (`--pm_c`; they fit 0.02 to
1.52), `V` and `rho_air` at the breakup and `R1/R2` the ratio of the radii of the leading body and the
fragment, taken from their `K` as if they had the same shape and density. It turns the direction of motion by
`atan(U/V)` in a random direction around it at the breakup, which drag keeps. With ablation, it also moves
the start of the dark flight sideways by the path flown since the breakup times its sine. The kick turns the
inertial motion, so the motion relative to the ground turns by as much within the ratio of the speed of the
rotation of the Earth to that of the fragment (~2%).

**Ablation before the dark flight** (`--sigma`, in s²/km², 0 by default). The dark flight never ablates: below
3-4 km/s the ablation has all but stopped. Before it, a fragment which broke off fast does, and its mass follows
its speed through the air as `m = m0*exp(sigma/2*(v^2 - v0^2))`, its `K` growing as `m^(-1/3)`. So `sigma` is
the fraction of its mass a body loses for every unit of kinetic energy per kg it dissipates,
`sigma = -d(ln m)/d(v^2/2)`: 0.01 s²/km² is 1% of the mass for every MJ/kg, and slowing down from 15 to
3.5 km/s, which dissipates 106 MJ/kg, leaves 35% of it. It is `sigma = Lambda/(2*Gamma*Q)`, with `Lambda` the
fraction of the energy of the air that heats and ablates the surface, `Gamma` the drag coefficient and `Q` the
energy to ablate 1 kg. The fireball types of Ceplecha et al. (1998, Table XVII) give apparent values of 0.014,
0.042, 0.10 and 0.21 s²/km² for the types I (ordinary chondrites), II, IIIA and IIIB, which hold mostly the
fragmentation. A value for the fireball itself follows from the alpha-beta fit of its leading body
(`python -m wmpl.Utils.AlphaBeta traj.pickle -v -e --fitvinit`), as `beta = (1 - mu)*sigma*V_e^2/2`
(Gritsevich 2009): `sigma = 2*beta/((1 - mu)*V_e^2)`, with `mu` from 0 (the area does not change) to 2/3 (the body
shrinks keeping its shape). It holds all the mass loss of the leading body, its fragmentation included, so for a
fragment it is an upper value; take `mu = 0` for it, as `mu = 2/3` triples it. With ablation the inverse scans
the starts of the dark flight, which do not depend on `sigma`, and flies back up for the breakup and the mass
there, so `--sigma 0 0.014` gives both, the breakup without ablation then along the ground track as well; the
exact one is that of `--sigma 0` alone. In the direct mode the masses are those at the breakup.

**The lines** (`--lines N`) are drawn as for the clones of the hand-over: from the Monte Carlo of the solution
if it was solved with it (the `<name>_mc_uncertainties.pickle` in the `Monte Carlo` directory next to the
trajectory pickle, found by default, or `--mc_uncertainties`), or else from the covariance of the line-of-sight
fit. Each moves the breakup (or the start of the dark flight) across the line and turns its direction with
it.

**Direct** (`--direct HT [HT ...]`, the breakup heights in km): every fragment of `--masses` (kg, at the
breakup) with `--density`, `--shape` (`sphere`, `brick`, `cube`, `cylinder`, `cone` or `A=<value>`),
`--drag_model` and, where measured, `--sphericity` or `--flatness`, `--elongation` and `--volume_factor` (for
`--drag_model bagheri`), as in OpenDarkflight. `--ks` gives fixed values of `K` instead, one for every mass.
Every fragment gets a run without a kick nor a line drawn, from the breakup without ablation and from where
it has slowed down with it, and max(`--kicks`, `--lines`) runs which each draw a
kick and a line around the fitted one (below), as many as asked of each, every run with `--mc_runs` Monte Carlo
runs of OpenDarkflight for its statistical wind errors. The range of the falls in the summary is taken from
the falls which hold all these errors. The output directory gets `fragment_falls.json`, one input per run under
`fragments/`, and after `--run` (or `--collect` later) `fragment_falls.csv` and `fragment_falls.png`, with the
falls by breakup, along and across the motion from the last observed point of the leading body.

**Inverse** (`--inverse meteorites.csv`): the CSV has the columns `name`, `lat`, `lon` (deg), `mass_kg` and
`density` (kg/m³), and optionally `ground_ht_m` (or `--end_ht` in km, or `auto` for the terrain lookup of
OpenDarkflight), `shape`, `drag_model`, `sphericity`, `flatness`, `elongation` and `volume_factor`, or `k_m2_kg`
for a fixed `K`. Without a shape or a drag model, a meteorite is a sphere with the `dfn` drag, in the flight
down and in OpenDarkflight alike (whose own default is a brick). For every meteorite:

1. OpenDarkflight flies the meteorite from breakups on the line, with the speed of the leading body there, at
   heights every `--scan_step` (0.25 km), from `--scan_top` (5 km above the highest fragmentation, or the
   first observed point, by default) down to `--scan_bottom` (the last observed point), and the breakup whose
   fall is as far along the motion as the meteorite is interpolated. What is left over across the motion tells
   how well the wind and the line agree with it.
2. `--mc_runs` (30 by default) Monte Carlo runs of OpenDarkflight at that breakup, the `--kicks` and the
   `--lines` (30 by default) give the spread of the fall along the motion, added in quadrature, which turns
   into a range of breakups.

With ablation (`--sigma` above 0), the dark flight is started at `--v_dark` on the line at the scanned
heights instead, from the highest fragmentation down; the start whose fall is as far along is found, and from
it the meteorite flies back up the ground track until it has the speed of the leading body, for the breakup
and the mass it had there, for every `sigma` given.

The breakup assumes that the meteorite flew alone from the leading body, with its speed: if it broke off a
larger fragment later, it broke off the leading body higher up and was born as a piece lower down. The shape matters most: the same meteorite traced
back as a sphere or as a brick can give breakups kilometres apart, so its measured shape (`sphericity`, or
better the form for `--drag_model bagheri`) and density are worth having. As the wind is the same for all the
meteorites of a fall, their differences in breakup height are better known than where each broke off, and
the report gives them for every pair. The output directory gets `inverse.json`, `inverse.txt` and `inverse.png`.

| Flag | Description |
| :--- | :--- |
| `-f`, `--fragheights` | Heights of the fragmentations of the drag fit (km), e.g. those of the flares. |
| `--direct`, `--masses`, `--ks` | Breakup heights (km), masses of the fragments (kg), and optionally their fixed `K` (m²/kg). |
| `--density`, `--shape`, `--drag_model` | Bulk density (3500 kg/m³), shape (`sphere`) and drag model (`dfn`) of the fragments. |
| `--sphericity`, `--flatness`, `--elongation`, `--volume_factor` | Measured form of the fragments, as in OpenDarkflight. |
| `--inverse`, `--scan_top`, `--scan_bottom`, `--scan_step` | CSV of the meteorites found, and the breakups (with ablation, the starts) scanned (km). |
| `--v_dark` | Speed relative to the ground where the dark flight starts with ablation (km/s). 3.5 by default. |
| `--kicks`, `--pm_c` | Kicks of Passey and Melosh drawn for every fragment or meteorite (0 by default), and their constant `C` (1). |
| `--lines`, `--mc_uncertainties` | Lines drawn around the fitted one (0 for `--direct`, 30 for `--inverse`), and the Monte Carlo uncertainties to draw them from (by default the one next to the trajectory pickle, if any). |
| `--sigma` | Ablation coefficient before the dark flight (s²/km²), 0 by default, and then OpenDarkflight flies from the breakup. The inverse takes several and gives the breakup for each, along the ground track. |
| `--mc_runs` | Monte Carlo runs of OpenDarkflight for the wind errors. 0 for `--direct`, 30 for `--inverse`. |
| `--end_ht` | Height where the dark flight ends (km), or `auto`. |
| `--run`, `--collect` | Direct: run OpenDarkflight and collect the falls, or only collect them. The inverse always runs it. |

In Python:

```python
import os

from wmpl.Utils.AtmosphereProfile import AtmosphereProfile
from wmpl.Utils.DragFitDarkFlight import runOdf
from wmpl.Utils.DragFitFragmentFalls import invertFalls, prepareFragmentFalls, readMeteorites, \
    summarizeFragmentFalls
from wmpl.Utils.Pickling import loadPickle

traj = loadPickle(*os.path.split("/path/to/traj.pickle"))
odf = ("/path/to/envs/opendf/bin/python", "/path/to/OpenDarkflight")
profile = AtmosphereProfile("profile.csv")

# Direct: fragments of 1 g to 1 kg from two breakups, with 10 kicks each
setup, ini_paths = prepareFragmentFalls(traj, "direct", [35000.0, 30000.0], profile, [35000.0, 30000.0],
    [0.001, 0.01, 0.1, 1.0], dict(density=3500.0, shape='sphere', drag_model='dfn'), n_kicks=10,
    mc_runs=5, end_ht=0.25, odf_python=odf[0], odf_dir=odf[1])
runOdf(ini_paths, *odf)
print(summarizeFragmentFalls("direct"))

# Inverse: the meteorites found
results = invertFalls(traj, "inverse", [35000.0, 30000.0], profile, readMeteorites("meteorites.csv"), *odf,
    n_kicks=10)
```

---

## Reading the output

**The table** printed by `DragFit` (an example from a fireball over Ontario, 2022):

```
Drag model of the motion: dv/dt = -K*rho*v² + g
    Points: 267, segments of K: 10, residuals RMS = 62.4 m, robust scale = 50.3 m
    First point: t = 0.002 s, v = 14.821 km/s, h = 89.749 km
    Gravity drop: yes, above the first point: no drag
    K end/beg = 1.64, single-body m end/beg = 0.226
    Final dynamic mass = 0.126 kg (1 sigma 0.0803 - 0.197 kg), for Gamma*A = 0.78 and bulk density = 3500 kg/m³

 Seg, Ht beg, Ht end, Time beg,  Vel beg,  Dyn pres,     K beg,     K end, K std, Growth, Signif
    ,   (km),   (km),      (s),   (km/s), beg (kPa),   (m²/kg),   (m²/kg),   (%),    (x), (sigma)
   1,  89.75,  83.80,    0.002,   14.821,       0.8,  4.11e-03,  4.12e-03,    81,   1.00,    0.0
 ...
   9,  42.18,  36.23,    3.749,   13.345,     511.9,  4.12e-03,  4.12e-03,     8,   1.00,    0.0
  10,  36.23,  30.28,    4.289,   11.734,     958.7,  4.12e-03,  6.75e-03,    15,   1.64,    2.2 *
```

- **Residuals RMS** is the scatter of the observed lengths around the model, along the track.
- **Dyn pres** is the dynamic pressure `rho*v^2` at the top of the segment, with the velocity
  relative to the ground, without the winds of an atmosphere profile.
- **K std** is the relative uncertainty of `K` at the bottom of the segment (1 sigma). It is large
  high up, where the meteor barely decelerates, and `inf` where the data do not constrain `K`.
- **Growth** is the factor by which `K` grows in the segment, and **Signif** its significance in
  sigma. A star marks the segments where `K` grows by more than 2 sigma.
- **K end/beg** and **single-body m end/beg** are the total growth of `K` and the mass ratio a single
  body of a constant shape would need to explain it (`m ~ K^-3`).
- **Final dynamic mass** is the dynamic mass of the leading body at the last point, for the assumed
  `Gamma*A` and density (see [the dynamic mass](#the-dynamic-mass)).

**The plots** show the height against time, the velocity against height and the ballistic
coefficient against height with its 1 sigma band. The data are the solver's points, and the
velocities of the data are point-to-point differences, which are noisy. The model is solid in the
observed part and dashed where it is extrapolated.

![Drag fit of the 2022-04-18 fireball over Ontario](figures/DragFit_example.png)

The plots for the fireball of the table above:

- **Height against time.** The model goes through the points of all six cameras. Above the beginning
  it continues without drag, and below the end at 30 km it continues as a drag flight, which goes
  beyond the plotted time range down to 20 km, where the path has bent by 5 deg.
- **Velocity against height.** The meteor barely decelerates above 50 km, and slows from 14.3 km/s at
  50 km to 7.6 km/s at the end. The drag flight below the end slows it down to 0.25 km/s at 20 km.
- **Ballistic coefficient against height.** `K` is constant down to 36 km and grows by a factor of
  1.64 in the last segment (2.2 sigma). The MetSim model of this fireball has a fragmentation of 70% of
  the mass of the main fragment at 33.8 km, which grows `K` by 1.56. The band widens above 60 km,
  where the deceleration is too small to constrain `K`. The title of the panel gives the final dynamic mass of the
  leading body, 0.13 kg (0.08 to 0.20 kg) for `Gamma*A = 0.78` and `rho_m = 3500 kg/m³`.

**The samples** printed by `SampleTrajectoryPositions` give, for every sample, the time, the
geographical coordinates, the ground-fixed azimuth and elevation of the radiant, the velocity in
the ECI frame (`Vel`) with its 1 sigma uncertainty (`Vstd`) and relative to the ground (`Vgnd`), the
deceleration of the model (`Decel`) and the dynamic pressure (`Dyn pres`). Samples outside the
observed part are marked with a star. Samples outside the modelled part are skipped, and a message
gives the modelled range of heights and times. The largest dynamic pressure among the samples is
printed at the end. For the example trajectory, sampled every 2 km from 118 to 94 km:

```
 Time(s), Sample ht (m),  Lat (deg),   Lon (deg), Height (m), Azim (deg), Elev (deg), Vel (km/s), Vstd (m/s), Vgnd (km/s), Decel (m/s^2), Dyn pres (kPa)
* -0.032,      118000.0,  44.121908,  -81.316691,   118000.0, 162.213031,  61.012481,     67.426,       30.2,      67.476,         -8.28,          0.124
   0.002,      116000.0,  44.131234,  -81.320845,   116000.0, 162.210702,  61.002759,     67.426,       30.1,      67.476,        141.91,          0.161
 ...
   0.308,       98000.0,  44.215458,  -81.358420,    98000.0, 162.189576,  60.916085,     66.584,       64.5,      66.633,      19239.82,          3.658
*  0.342,       96000.0,  44.224849,  -81.362618,    96000.0, 162.187211,  60.906425,     65.566,      379.3,      65.616,      40396.81,          5.101
```

- **Azim** and **Elev** are the radiant above a fixed ground. Drag acts against the velocity
  relative to the air, so it does not turn the direction of motion above the ground, and the radiant
  is derotated with the initial velocity, as in `Orbit.calcOrbit` and `DynamicMassFit`. Only the
  bending by gravity changes it. Derotating it with the instantaneous velocity would turn it by about
  `v_rot*(1/v - 1/v_init)` radians, where `v_rot` is the component of the Earth's rotation velocity
  across the path: 0.5 deg at the end and 2 deg at 5 km/s for the synthetic fireball of the tests.
- **Decel** is `K*rho*v^2 - g_along`, positive when the meteor slows down. Above the beginning, without
  drag, it is the small negative acceleration by gravity.
- **Dyn pres** is `rho*v^2` with the velocity relative to the ground, `Vgnd`, as in the `DragFit`
  table, without the winds of an atmosphere profile. `Vgnd` is the modelled inertial speed less the
  rotation of the Earth along the trajectory, the speed the drag acts on, also handed to a dark flight
  by `DragFitDarkFlight`. The ECI
  velocity differs from it by up to the speed of the Earth's rotation (465 m/s at the equator). For
  the synthetic fireball of the tests, sampled from 18 to 7 km/s at 45 deg latitude, the dynamic
  pressure with the ECI velocity would be 2 to 6% lower.
- **Vstd** is propagated from the covariance of the fitted parameters. When sampling by height, the
  uncertainty of the time is also computed. When sampling by time, it is the uncertainty of the
  length along the trajectory. The sampled coordinate itself is exact.

With `--csv`, the samples are saved with the Julian date, the length along the trajectory, the ECI
coordinates, the uncertainties of the time and the length, the ballistic coefficient and its
uncertainty, and a flag of the extrapolated samples.

---

## Using it as a library

```python
import os

import numpy as np

from wmpl.Utils.DragFit import DragFit, plotDragFit
from wmpl.Utils.Pickling import loadPickle
from wmpl.Utils.SampleTrajectoryPositions import sampleTrajectory

traj = loadPickle(*os.path.split("/path/to/20191023_091225_trajectory.pickle"))

# Fit the model
drag_fit = DragFit(traj)
print(drag_fit.summaryTable())

# Smooth velocities and heights at any time (s, relative to traj.jdt_ref) or height (m)
time_arr = np.linspace(0.0, 1.0, 11)
vel = drag_fit.velocityAtTime(time_arr)
ht = drag_fit.heightAtTime(time_arr)
vel_30km = drag_fit.velocityAtHeight(30000.0)

# The ballistic coefficient and the uncertainty of its logarithm at lengths along the trajectory
k_arr, log_k_std = drag_fit.ballisticCoeffAtLength(drag_fit.lengthAtTime(time_arr))

# The final dynamic mass with its 1 sigma range, for an assumed Gamma*A and bulk density
mass, mass_lo, mass_hi = drag_fit.dynamicMassAtLength(drag_fit.len_last, gamma_a=0.78, bulk_density=3500.0)
print(drag_fit.summaryTable(gamma_a=0.78, bulk_density=3500.0))

# Positions and the direction of the radiant, bent by gravity
eci = drag_fit.eciAtLength(drag_fit.lengthAtTime(time_arr))
radiant = drag_fit.radiantAtLength(drag_fit.lengthAtTime(time_arr))

# The deceleration (m/s^2) and the air density (kg/m^3) of the model, and the velocity relative to the
#   ground (m/s)
len_arr = drag_fit.lengthAtTime(time_arr)
decel = drag_fit.decelerationAtLength(len_arr)
rho = drag_fit.airDensityAtLength(len_arr)
vel_gnd = drag_fit.groundVelocityAtLength(len_arr)

# The 1 sigma uncertainties of the velocity and the time at given lengths, or of the velocity and the
#   length at given times
vel_std, time_std = drag_fit.motionStdAtLength(len_arr)
vel_std, len_std = drag_fit.motionStdAtTime(time_arr)

# The plots, with the final dynamic mass
plotDragFit(drag_fit, gamma_a=0.78, bulk_density=3500.0)

# Samples every 100 m of height between the observed beginning and end (-1)
samples = sampleTrajectory(traj, -1, -1, 100.0)
print(samples.t_est[0], samples.ht[0], samples.vel[0], samples.azim_norot[0], samples.elev_norot[0])

# All sampled quantities as numpy arrays: t_est, ht, lat, lon, ele_geo, azim_norot, elev_norot, vel,
#   vel_norot, length, eci, radiant_eci, vel_std, time_std, len_std, decel, dyn_pressure, k, log_k_std
#   and extrapolated
arrays = samples.toArrays()
print(arrays['vel_std'], arrays['decel'], arrays['dyn_pressure'])

# The dynamic mass of every sample, from the fitted ballistic coefficient
mass, mass_lo, mass_hi = samples.drag_fit.dynamicMassAtLength(arrays['length'])

# Save the samples to a CSV file
samples.saveCSV("samples.csv")
```

All functions return `NaN` outside the modelled part of the trajectory. The lengths are measured
along the trajectory from the state vector, as `obs.state_vect_dist` in the trajectory solution.

---

## Interpreting the ballistic coefficient

`K` is fitted only from the deceleration, so it is independent of the photometry. How to read it:

- **Constant K** means the deceleration is driven only by the air density and the velocity, with no
  measurable change of the meteor.
- **Slowly growing K** is expected from ablation. For a single body of a constant shape `K ~ m^(-1/3)`,
  so losing 90% of the mass doubles `K`.
- **A fast growth over one or two segments**, marked with a star, points to a fragmentation, as the
  leading fragment has less mass for its area. A jump of `K` by a factor `f` means the leading
  fragment kept `1/f^3` of the mass. The dynamic pressure at that height is an estimate of the
  strength of the meteoroid.
- **The value of K** gives a dynamic mass for an assumed shape and density (see
  [the dynamic mass](#the-dynamic-mass)). Comparing it with the photometric mass constrains the shape
  and the density.

Compared with nine MetSim fragmentation models of fireballs, the fit found the large fragmentations
where the deceleration is measurable. For example, the 70% fragmentation of the Ontario 2022
fireball at 33.8 km shows up as a growth by 1.71 against 1.56 in MetSim. In the well constrained part,
the fitted `K` agreed with the `K` of the main fragment in MetSim within a factor of 0.9 to 1.6.

Keep in mind:

- **High up, K is poorly constrained.** Above about 60 km the deceleration of a fireball is too
  small to measure, so fragmentations there are not seen. Check the 1 sigma band and the K std
  column.
- **The fit smooths the steps.** The penalty on the growth of `K` turns a sharp fragmentation into a
  ramp, often a few km below the real break, so the dynamic pressure at the start of the growth
  tends to overestimate the strength. Steps closer than a segment cannot be separated, and growths
  by less than about 20% (a mass loss below about 40%) are not detectable.
- **Known fragmentations can be given.** If the heights of the fragmentations are known, e.g. from the
  flares, `-f` (or `knot_hts`) lets `K` jump there without penalty. This matters most when the last
  fragmentation is close to the end: the smoothed ramp then overestimates `K` at the last point and
  underestimates the final dynamic mass. On 165 MetSim fireballs with 2 or 3 fragmentations, the last
  one 1.5 to 9 km above the last point, the median error of `K` at the last point went from +8-12% to
  -3 to +0.5%, and `K` was within 10% of the truth in 61-76% of them instead of 31-44%.
- **How spread the fragmentation was is not known.** The motion alone does not tell a jump from a
  fragmentation spread over a few km, so with `-f` the fit is repeated with the growth spread below
  the heights, in steps of 250 m down to the last point, and the farthest `K` at the last point among
  the spreads which fit within `Delta chi2 <= 4` is added to its uncertainty. On the same 165
  fireballs, the uncertainty of the fit alone covered the true `K` in 61% (1 sigma) and 81% (2 sigma)
  of them, and with the spread in 86% and 96%.
- **K cannot decrease.** A change of shape that lowers the drag cannot be followed.
- **The centroid is not the leading fragment.** After a large fragmentation the measured positions
  follow the centroid of the luminous cloud, which can lag the leading fragment and inflate `K`.
- **The air density model** (NRLMSISE-00) can be off by 10-20%, which goes directly into `K`.

---

## The dynamic mass

For an assumed product of the drag coefficient and the shape factor `Gamma*A` and bulk density
`rho_m`, the ballistic coefficient gives the mass of the leading body,

```
m = (Gamma*A / (K*rho_m^(2/3)))^3
```

This is the same dynamic mass as in `wmpl.Utils.Physics.dynamicMass` and
`wmpl.Utils.DynamicMassFit`, as `K = decel/(rho*v^2)`. Here it comes from the fit of the whole
trajectory instead of the deceleration at one point. By default `Gamma*A = 0.78` and
`rho_m = 3500 kg/m³`, which can be changed with `--ga` and `--dens` (note that `DynamicMassFit` uses
`Gamma*A = 0.55` by default).

The **final dynamic mass** at the last point is the best estimate of the mass that leaves the
luminous part of the trajectory. Below the end `K` is kept constant, so this is also the mass of the
drag flight. The terminal velocity of the body at the ground is `sqrt(g/(K*rho_0))`, about 34 m/s for
the 0.13 kg of the example.

**The uncertainty.** The mass goes as `K^-3`, so the uncertainty of `ln K` from the fit is tripled,
`sigma(ln m) = 3*sigma(ln K)`, and the 1 sigma range is asymmetric. In the example, `K` at the end is
known to 15%, so the mass is known to a factor of 1.56 (0.08 to 0.20 kg). Only the uncertainty of the
fit is propagated. The assumptions usually matter more: the mass scales as `(Gamma*A)^3/rho_m^2`, so
changing `Gamma*A` by 20% changes the mass by a factor of 1.7, and an iron (`rho_m = 7800 kg/m³`)
has 0.2 of the mass of a stone with the same `K`.

The dynamic mass is meant for fireballs. For faint meteors, which fragment and erode into small
grains, the deceleration follows the grains and the dynamic mass comes out orders of magnitude below
the photometric mass (e.g. 1e-10 kg for the faint example trajectory in the repository).

---

## Limits of the model

- **The tail below the end is a drag flight of a surviving body.** For a meteor which ends by
  ablating completely, it is not a real trajectory.
- **The modelled part ends** at 200 km height (or 100 km above the beginning, if higher), where the
  path has bent by 5 deg below the end, or at the ground. With `--dragabove` it also ends where the
  velocity above the beginning would have grown by 10%.
- **Without the gravity correction** in the trajectory solution, the positions and the radiant stay
  on the straight line.
- **The uncertainties of the velocity, the time and the length** only include the noise of the
  points. They do not include the uncertainty of the trajectory solution or the bias of the model,
  e.g. from the smoothing of `K` by the penalty on its growth. On the synthetic fireball of the tests
  (18 to 10 km/s, `K` growing five times), over 20 realisations of the noise the errors of the
  velocity were 0.8 to 0.9 sigma (RMS) in the first quarter of the observed part and 0.6 to 0.7 sigma
  below the end, but 5 sigma at 70 to 90% of the observed part, where `K` grows fastest and the model
  is biased by 70 to 100 m/s. They are linearised, which is less accurate far below the end.

---

## Troubleshooting

**`ValueError: The drag model could not be fitted, the lengths and times of the points of the
trajectory are not consistent`.** The lengths along the track and the times of the solution do not
describe a physical motion, usually because of a bad solution or a wrong timing of a station. Check
the lag and the residual plots of the trajectory solution.

**Samples are skipped.** They are outside the modelled part of the trajectory; the message gives the
modelled range of heights and times. Below that, use a dark flight model.

**The K std is `inf` and the band covers the whole plot.** The data do not constrain `K` there,
usually because the meteor did not decelerate measurably.
