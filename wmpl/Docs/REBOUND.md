# REBOUND orbital integration

**`wmpl.Rebound.REBOUND`** integrates a meteoroid orbit backwards or forwards in time through the
real solar system, starting from a solved trajectory. It reports the orbital elements at the end of
the integration, finds close encounters and impacts, propagates the measurement uncertainty with
Monte Carlo clones, and can optionally include radiation forces and measure how chaotic the orbit is.

## Table of Contents

- [Quick start](#quick-start)
- [Requirements](#requirements)
- [What it reads](#what-it-reads)
- [Command line reference](#command-line-reference)
- [Choosing an integrator](#choosing-an-integrator)
- [Radiation forces](#radiation-forces)
- [Close encounters, impacts and ejections](#close-encounters-impacts-and-ejections)
- [Monte Carlo clones](#monte-carlo-clones)
- [Chaos: MEGNO](#chaos-megno)
- [What it writes](#what-it-writes)
- [The orbit to use at the end](#the-orbit-to-use-at-the-end)
- [Worked examples](#worked-examples)
- [Using it as a library](#using-it-as-a-library)
- [Troubleshooting](#troubleshooting)

---

## Quick start

```
python -m wmpl.Rebound.REBOUND /path/to/20191023_091225_trajectory.pickle
```

That integrates the nominal orbit 60 days backwards with IAS15 and writes a report, a plot and a JSON
file next to the pickle. A more typical run, 100 years backwards with 100 Monte Carlo clones on 8
cores and a fixed seed:

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 36525 --outputs 2000 --mc 100 --cores 8 --seed 42
```

---

## Requirements

REBOUND is an optional dependency and is not installed by the standard wmpl setup:

```
conda install -y -c conda-forge astropy
conda install -y -c conda-forge rebound
pip install --no-build-isolation reboundx
```

Check that both imported **and that REBOUNDx attached**:

```
python -c "import rebound, reboundx; s = rebound.Simulation(); reboundx.Extras(s); print('ok' if s.extras else 'reboundx not attached: reinstall it')"
```

**Both are needed.** `rebound` alone is not enough — the general relativity correction, the Earth's
gravitational harmonics and the radiation forces all come from `reboundx`.

**Use REBOUNDx ≥ 4.4.0.** The simulations run in the ecliptic J2000 frame, and the Earth's J2/J4 are
made symmetric about its rotation pole of date (precessed from J2000, about 20″ per year) through the
`Omega` spin-axis parameter that REBOUNDx added in 4.4.0. Older versions ignore it and make them
symmetric about the ecliptic pole, 23.4° off; the run then warns, and the summary printed at the end
says so. The effect is a few m/s on the departure of a slow meteoroid (5.2 m/s for Winchcombe, at 14
km/s).

**`--no-build-isolation` is not optional.** `reboundx` publishes no wheels, so `pip` always compiles
it, by default against the newest `rebound` on PyPI rather than the one you installed. When the two
differ, their C simulation structures differ, `reboundx` writes its pointer outside the structure,
and every parameter assignment fails with `Need to attach reboundx.Extras instance to simulation
before setting params` — or `import reboundx` fails outright with a missing symbol. The same
mismatch appears if `rebound` is upgraded after `reboundx` was compiled. To repair an existing
install:

```
pip install --force-reinstall --no-deps --no-build-isolation --no-binary reboundx reboundx
```

That pulls the latest `reboundx`, which needs `rebound` 5 or newer, so update `rebound` first. To
keep an older `rebound`, pin the `reboundx` released alongside it instead (e.g. `reboundx==4.3.0`
for `rebound==4.3.0`); a newer `reboundx` does not compile against an older `rebound`.

Both are listed unpinned in `requirements.txt`, so `pip install -r requirements.txt` also works:
with no pin, `pip` builds `reboundx` against the same newest `rebound` it installs. Pinning
`rebound` there does not work, since `reboundx` builds against the newest one regardless.

**Platform note:** this works on Linux and macOS. `reboundx` has no Windows wheel and does not compile
with MSVC (it uses C99 variable-length arrays, which MSVC does not support), so on Windows use WSL2:
install Ubuntu, set up the wmpl conda environment inside it, and run from that shell.

**Version note:** the code works with both REBOUND 4 and REBOUND 5. The TRACE integrator was added
in REBOUND 4.4.0; on an older REBOUND, `--integrator trace` fails immediately with a message naming
the installed version, before any work is done. A REBOUNDx that cannot attach is also reported
before any work is done, with the command that repairs it.

---

## What it reads

**The trajectory pickle** (the positional argument). See [the trajectory manual](Trajectory.md) for how
to produce one. The integration uses `traj.jdt_ref` as the starting epoch, `traj.state_vect_mini`,
`traj.v_init` and `traj.radiant_eci_mini` as the starting state, `traj.traj_id` as the object name,
and `traj.uncertainties` to draw the Monte Carlo clones.

**The planetary ephemeris.** By default the local DE430 kernel at `wmpl/share/de430.bsp`, which the
wmpl setup downloads. With `--horizons` the planet states come from the JPL Horizons web service
instead, which is slower and needs network access, but does not need the kernel.

The Sun, all eight planets and the Moon are included as massive bodies (the Moon separately, since a meteoroid can pass through the Earth-Moon system). The meteoroid is a
massless test particle, so it does not perturb them, which is what makes the Monte Carlo clones
independent and cheap to run in parallel.

---

## Command line reference

```
python -m wmpl.Rebound.REBOUND --help
```

| Flag | Description |
| :--- | :--- |
| `pickle_path` | **Required.** Path to the trajectory pickle. All output is written next to it. |
| `--days DAYS` | Length of the integration in days. Default: `60`. |
| `--forward [DAYS]` | Integrate forward in time instead of backward. With a value, that many days; without one, `--days` is used. |
| `--outputs N` | Number of times along the integration at which the state is saved. Default: `500`. Increase it for long runs — 500 samples over a century is one sample per 10 weeks. |
| `--mc N` | Number of Monte Carlo clones drawn from the trajectory covariance. Default: `0`. |
| `--seed SEED` | Random seed for the Monte Carlo sampling, making a run exactly reproducible. If not given, a seed is drawn and reported. |
| `--cores N` | Number of parallel processes used for the clones. Default: all but one core. |
| `--geocentric` | Report in the geocentric frame instead of the heliocentric one. Switches the semi-major axis and perihelion units from AU to km. |
| `--horizons` | Use the JPL Horizons web service for planet positions instead of the local DE430 kernel. |
| `--integrator {ias15,whfast,trace}` | Which integrator to use. Default: `ias15`. See [Choosing an integrator](#choosing-an-integrator). |
| `--dt DAYS` | Timestep in days for `whfast` and `trace`. Default: `0.5`. Must be positive. Ignored by `ias15`. |
| `--beta BETA` | Include solar radiation pressure and Poynting–Robertson drag with this beta. |
| `--radius METRES` | Object radius in metres. With `--density`, beta is computed from it. Purely gravitational if none of `--beta`, `--radius` and `--mass` is given. |
| `--mass [KG]` | Object mass in kg. With `--density`, beta is computed from it. Without a value, the photometric mass of the light curve is used. |
| `--density KG_M3` | Bulk density in kg/m³, used with `--radius` or `--mass`. Default: `3000`. |
| `--compute_megno` | After the integration, measure the MEGNO chaos indicator of the nominal orbit. See [Chaos: MEGNO](#chaos-megno). |
| `--verbose` | Print the progress of the simulation. |

`--beta`, `--radius` and `--mass` are mutually exclusive: give the beta directly, or give a size or a
mass (with the density) and let the code compute it, but only one of them.

---

## Choosing an integrator

| Integrator | Timestep | Resolves close encounters | Use it when |
| :--- | :--- | :--- | :--- |
| `ias15` (default) | Adaptive | Yes | Almost always. Accurate to machine precision, and the only safe choice if you are unsure. |
| `whfast` | Fixed (`--dt`) | **No** | Long integrations of an orbit that does not pass near a planet, where speed matters. |
| `trace` | Fixed (`--dt`) | Yes | Long integrations that do have close encounters. Needs REBOUND ≥ 4.4. |

Two things happen automatically with the fixed-step integrators, and it is worth knowing about both:

**The departure from the Earth is always integrated with IAS15.** The object starts on the Earth's
surface, which is as deep a close encounter as it gets. WHFast cannot resolve it at all, and TRACE
resolves the encounter but applies the Earth's J2/J4 harmonics without sub-stepping them, which leaves
the orbit about 1–2 % off in the semi-major axis. Both therefore start with IAS15 and take over at the
first output after the object has left the Earth's neighbourhood (3 Earth Hill radii), where the
harmonics are negligible. The report states when the handover happened, and the JSON records it as
`run.fixed_step_from_days`.

**A backward TRACE run is integrated in reversed time.** TRACE mishandles close encounters when the
timestep is negative — a REBOUND bug, present in 4.6.0 and 5.1.1, where a flyby integrated backward
comes out about 10 % off in the semi-major axis while the same flyby forward matches IAS15 to 1e-7.
The workaround is to integrate forward with every velocity reversed, which is exact for gravity, for
general relativity and for the Earth's harmonics. Poynting–Robertson drag is not time-reversible, so
it is handled by flipping the sign of `c` on the radiation force, which gives exactly `a(r, −v)`.
This is all internal; the times reported stay on the signed axis you asked for.

**On accuracy.** On the example trajectory `20191023_091225` integrated 100 years backward, WHFast and
TRACE at `--dt 0.5` both ended within about 1e-8 of IAS15 in the semi-major axis, at roughly 4× and 3×
the speed. **That number is not a guarantee.** On a synthetic Earth-departure test with a lunar
encounter, the same 0.5-day step gave differences of 2e-4 (WHFast) and 7e-5 (TRACE) over the same
span. An orbit with a small perihelion needs a shorter step, because the fixed step has to resolve the
fastest part of the orbit. **Before trusting a long fixed-step run, do one shorter run both ways and
compare** — the run is cheap relative to being wrong:

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 3652 --integrator ias15
python -m wmpl.Rebound.REBOUND traj.pickle --days 3652 --integrator trace --dt 0.5
```

WHFast in particular does not resolve close encounters at all: measured through flybys inside 0.1 Hill
radii, its error in the semi-major axis afterwards was 0.3 % to 200 %, where TRACE stayed within 4e-6.
The report therefore prints a warning whenever a close encounter happened while WHFast was
integrating, for the nominal solution or for any clone. **Take that warning seriously** — rerun with
`--integrator trace` or `ias15`.

**A small perihelion is a close encounter for this purpose too.** A fixed step that is fine along most
of an eccentric orbit is far too long through perihelion, and the Sun passages inside 0.1 AU
(`WHFAST_SUN_WARNING_AU`) are counted in the warning for that reason. Sun passages are listed out to
0.2 AU, but beyond 0.1 AU the default step is enough: on orbits with a = 1.27 and 2.6 AU (i = 25°,
100 years, the Sun and five planets) WHFast at 0.5 d was 9e-4 to 1e-2 off in a at q = 0.07 AU, 5e-6 to
6e-6 at 0.10 AU, 1e-6 to 4e-6 at 0.14 AU and below 7e-7 from 0.17 AU. Measured on a δ-Aquariid-like
orbit (a = 2.6 AU, q = 0.07 AU, i = 27°), 100 years, with the Sun and five planets:

| Integrator | `--dt` | Final a against IAS15 |
| :--- | :--- | :--- |
| IAS15 | adaptive | reference (a = 2.5587 AU) |
| WHFast | 0.5 d (default) | **ejected** (a = −8.9 AU) |
| WHFast | 0.1 d | 1.7e-4 |
| TRACE | 0.5 d (default) | 7e-4 |
| TRACE | 0.1 d | 2e-9 |

For an orbit with perihelion inside ~0.1 AU use IAS15, or TRACE with `--dt 0.1`.

REBOUNDx will warn that the `gr_full` general relativity correction is velocity-dependent when used
with WHFast. Its measured effect on the orbit was below 1e-8 in the semi-major axis over 100 years, so
the warning can be ignored here.

---

## Radiation forces

By default the integration is purely gravitational, since a trajectory solution does not constrain the
object's size or density. To include solar radiation pressure and Poynting–Robertson drag, give the
dimensionless **beta**, the ratio of radiation pressure to solar gravity:

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 36525 --beta 0.01
```

or let it be computed from a size and a density:

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 36525 --radius 1e-5 --density 3000
```

A size is not what a meteor solution gives, so a mass works too, either explicitly or — with no value
— as the photometric mass computed from the trajectory's own light curve:

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 36525 --mass 1e-6 --density 3000
python -m wmpl.Rebound.REBOUND traj.pickle --days 36525 --mass --density 3000
```

The mass is turned into the radius of the sphere of that mass and density, `s = (3m/(4πρ))^(1/3)`.
At a fixed mass, beta goes as `m^(-1/3)*rho^(-2/3)`. The mass, and through it the luminous efficiency
of a photometric mass, enters as a cube root, so a factor of 2 in it changes beta by a factor of 1.26.
The density enters more strongly: halving it raises beta by a factor of 1.59. `--mass 0` is the same
as `--mass` without a value.

For a spherical grain, `beta = 5.7425e-4*Q_pr/(rho*s)` with the density in kg/m³ and the radius in m.
A 1 µm grain at 3000 kg/m³ gives beta = 0.19; a 1 cm meteoroid of the same density gives 1.9e-5. In
other words, radiation forces matter for small grains and are negligible for fireball-producing
bodies — but they are not negligible over long integrations, where they changed the semi-major axis by
25 % to 190 % in the cases measured here.

---

## Close encounters, impacts and ejections

Closest approaches are tracked at **every internal integrator timestep**, not at the output samples,
and then refined *within* the step. Both steps matter:

- Near the Earth and the Moon the object can cross a whole detection sphere between two output
  samples, so a scan of the sampled output would miss the encounter entirely or badly overestimate
  its distance.
- Away from the Earth, IAS15 takes steps of one to two days, so the object can pass a planet
  completely inside a single step. The distance at the step ends alone can then be far above the
  true minimum. The minimum is therefore found on the cubic Hermite interpolant of the relative
  position and velocity across the step, which reproduces straight-line motion exactly and so is
  accurate even when one step spans the whole flyby.

- **Close encounters** are reported for any body approached within **3 Hill radii**, and for the Sun
  within **0.2 AU** (it has no Hill sphere, so its `n_hill` is `None`), listed in the order they
  happened. **A body can appear more than once**, if the object passed it more than once: each local
  minimum of the distance is refined and listed separately. A meteoroid in resonance with a planet
  commonly meets it many times.
  - *Why 0.2 AU for the Sun.* It is the limit of the Sun-approaching orbital class of
    [Borovička et al. 2005](https://doi.org/10.1016/j.icarus.2004.09.011), the one meteor studies use:
    almost all meteoroids inside it are Na-free or Na-poor
    ([Vojáček et al. 2019](https://doi.org/10.1051/0004-6361/201833289), section 3.3). Inside 0.1 AU all
    of them are, the sodium being thermally desorbed
    ([Kasuga et al. 2006](https://doi.org/10.1051/0004-6361:20065330)), asteroids disrupt at
    q ≈ 0.076 AU ([Granvik et al. 2016](https://doi.org/10.1038/nature16934)), and the sunskirter zone
    (q < 33 R☉ = 0.153 AU) of [Jones et al. 2018](https://doi.org/10.1007/s11214-017-0446-5) lies inside
    it. Kasuga et al. found no sodium loss down to q = 0.14 AU among larger shower meteoroids, so the
    limit is that of the small meteoroids. A Sun passage is a thermal event, not a gravitational
    perturbation, and is listed for that reason.
  - *Repeated Sun passages.* An orbit with q < 0.2 AU passes the Sun every revolution, so a
    low-perihelion stream lists one Sun passage every few years (measured: 24 in 100 yr for a
    δ-Aquariid-like orbit, against 6 planetary flybys). The terminal and the text report collapse more
    than three Sun passages into one line with the count, the perihelion range and the time span,
    placed where the first passage was. The JSON lists every passage.
- **Impacts** are detected against the physical radii of the Sun, the Earth and the Moon, using
  REBOUND's line-of-travel collision mode so a fast mover cannot tunnel through. Detection is armed
  only once the object has left the Earth's neighbourhood, so the trivial fact that it starts on the
  Earth is not reported as an impact. An impact ends the integration and is printed as a banner at
  the top of the report.
- **Ejections** are recorded if the object leaves the simulation volume (1000 AU by default), after
  which it stops being integrated.

On the plot every encounter is a star in the colour of the body it passed (see
[The figure](#the-figure)).

The first close encounter going backwards is almost always the Earth itself — that is the meteoroid
arriving. What matters is what came before it.

---

## Monte Carlo clones

`--mc N` draws N clones from the trajectory's measurement covariance and integrates each one, in
parallel, through an independent copy of the simulation. This is what turns a single number into a
number with an uncertainty:

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 36525 --mc 200 --cores 8 --seed 42
```

The report gives the final elements with 95 % confidence intervals from the clone ensemble, how many
clones impacted which body, how many escaped, and how many had a close encounter with each body. If
some clones did not survive, the confidence intervals are computed from the survivors only and the
report says so (`ci_uses_survivors_only` in the JSON). Under the final elements come the Tisserand
parameter and the orbital class at the end of the run: the meteoroid class of Borovička et al. (2005),
from `wmpl.Utils.OrbitClassification.classifyBorovicka` (`final_elements.orbital_class` in the JSON, see
[Classifying a single orbit](OrbitClassification.md)). Like the elements, they are the nominal's, and
with clones T_J has their 1σ and 95% interval and the class the fraction of the clones in each class
(`final_elements.clones_tisserand_and_class`). For Winchcombe, 1000 years back with 200 clones:

```
Tisserand parameter w.r.t. Jupiter and orbital class (Borovicka et al. 2005) at the end of the run,
nominal +/- 1 sigma [95% CI] of the 200 clones:
  T_J   = 3.083673 +/- 0.029853 [  2.995651,   3.112674]
  class = asteroidal-chondritic (nominal); clones: asteroidal-chondritic 97.0%, jupiter-family 3.0%
```

In a run longer than the Lyapunov time the nominal's class is that of one clone among many, and the
clones' fractions are the ones to quote. The class of the observed orbit is that of the trajectory
solution, before any integration. The ecliptic-shower class needs a shower association, which the
integration does not have, so it is never given.

The ensemble also gives a divergence estimate for free: the clones are trajectories started from
slightly different state vectors, so their spread over time is exactly what a Lyapunov analysis
measures. The report fits their root-mean-square separation both as exponential growth (the chaotic
case) and as linear growth (the regular case) and says which describes the data better. This is a
finite-time estimate driven by the real measurement uncertainty, so it is only meaningful while the
ensemble stays compact.

**Always pass `--seed`** if you intend to quote the result. Without it a seed is drawn and reported,
which is recoverable, but passing it is less work than reading it back out of the report.

---

## Chaos: MEGNO

`--compute_megno` runs a second, separate integration of the nominal solution over the same span, with
Newtonian gravity only and always with IAS15, and computes the MEGNO chaos indicator. In a run with clones
it does the same for 32 of them, which give the Lyapunov time to quote (see [below](#over-the-clones)).

MEGNO (Mean Exponential Growth factor of Nearby Orbits) follows how fast a small deviation from the
orbit grows. Its running mean ⟨Y⟩ tends to **2** for a regular, quasi-periodic orbit, tends to **0**
for a stable periodic one (for example an orbit librating in a resonance), and grows without bound as
(λ/2)·t for a chaotic one, where λ is the Lyapunov exponent.

Only the object's own deviation is followed. REBOUND's built-in `init_megno()` perturbs every body in
the simulation, and the linearly growing deviations of the Moon and the planets then swamp the
object's exponential growth for centuries — an orbit with a 64-year Lyapunov time still read
⟨Y⟩ = 2.01 after 300 years that way, against 4.5 when only the object is followed.

The verdict is based on the final ⟨Y⟩ and its mean over the last quarter of the run:

| Verdict | Meaning |
| :--- | :--- |
| `regular` | Both within 0.1 of 2. Quasi-periodic motion over the integrated span. |
| `chaotic` | Both above 2.5. A Lyapunov time is fitted to the growth of the deviation (see below). |
| `periodic` | Both below 1. The deviation stayed bounded over the span, e.g. libration in a mean-motion resonance. An Earth-crossing orbit can also do this for a while before its first deep encounter, so this describes the nominal orbit over the span, not the clones (see [Over the clones](#over-the-clones)). |
| `not_converged` | Anything else. Integrate for longer to decide. |
| `too_short` | The run covers fewer than 20 orbital periods. No verdict; the report says how long the run would need to be. |
| `unbound` | The orbit is hyperbolic, so MEGNO — which describes bounded motion — does not apply. |

**Lyapunov time.** It is fitted to ln|δ|, the logarithm of the deviation, which grows as λ·t on a
chaotic orbit, rather than to ⟨Y⟩, a double running average that lags the growth. The fit starts
when ⟨Y⟩ rises above 2.5 for the last time: before that the deviation can stay flat or even shrink for
a long time, e.g. on an Earth-crossing orbit until its first deep encounter, and including it biases
the time high. Chaos driven by close encounters grows in steps, so the growth phase is also split
into four consecutive windows of equal length, its quarters, and each is fitted on its own. The
report lists every quarter with its time; a quarter in which ln|δ| does not grow (it can shrink for a
while between encounters) is listed as "no growth".

- **1σ.** The scatter of the quarters' growth rates gives the uncertainty of the time: their standard
  deviation over √4 is the standard error of λ, and since the time is 1/λ its interval is asymmetric,
  e.g. 69 (+175/−29) years, with no upper limit when λ is within 1σ of 0.
- **Has it settled?** The time fitted from the onset to the end of each quarter (cumulative) shows
  whether it has: the report warns when it still changes by more than 20% over the last half of the
  growth phase, and then the run is too short to time the chaos well. The last quarter alone misses a
  slow drift: over 4000 years, 21 of 21 Earth-crossing orbits passed that test while the median of their
  times still rose by 21% from 2000 years.

In the JSON the quarters are the verdict's `lyapunov_windows` (`start_years`, `end_years` and
`lyapunov_time_years`, `null` where ln|δ| does not grow), the 1σ interval is
`lyapunov_time_1sigma_years` (its upper end `null` when unbounded) and the cumulative fits are
`lyapunov_cumulative`. `rebound_megno_lyapunov.png` has two panels: ⟨Y⟩ with the lines at 2 and
2.5, and ln|δ| with the straight fit over the growth phase, labelled with the Lyapunov time and its
1σ. The start of the fit, where ⟨Y⟩ last crosses 2.5, is marked in both. No time is given if the growth
phase covers fewer than 10 orbital periods (or a tenth of the run), since the fit is not reliable there.
This Lyapunov time is of the nominal orbit alone, in its Newtonian re-integration; for the clones, see
below.

<a id="over-the-clones"></a>
**Over the clones.** The Lyapunov time of one orbit is one draw. In a run with clones and with enough
orbits for a verdict, `--compute_megno` also integrates the first 32 clones in the same way, sharing the
cores. The report gives their verdicts and the median of their times with its 16–84% range, which is the
value to quote. Measured on Winchcombe over 1000 years with 64 clones, the median is 42 years (16–84%:
30–66), and the nominal orbit's 69 years is at the 86th percentile. With 32 clones the median is good to
9%, against 12% with 16 and 18% with 8. The 32 clones are integrated with Newtonian gravity only, so
they cost less than as many clones of the main integration.

**A rate over a finite time.** The Lyapunov exponent is the limit of this rate for t → ∞, and the
report gives its value over the span of the run. The median over the first half of the run says how much
it changes when the span doubles: within 10% it has settled, and otherwise it is a value for that span,
not yet a limit (that change has a scatter of 4–9% with 32 clones). On Winchcombe it is 39 years from
500 years and 46 from 1000. On a similar orbit with synthetic planets it keeps rising, at 20, 24, 28 and
35 years from 500, 1000, 2000 and 4000 years, which is more than the polynomial term of the deviation
explains (Tancredi et al. 2001, eq. 20).

In the JSON this is `clones` inside `megno`: `n_clones`, `status_counts`, `lyapunov_times_years`,
`lyapunov_time_median_years`, `lyapunov_time_p16_years`, `lyapunov_time_p84_years`,
`half_run_median_years` and `change_on_doubling`. The orbit figure marks this median as the Lyapunov
time (see [The figure](#the-figure)). The number of clones is the constant `MEGNO_CLONES`, which can be
set to 0 from a script to skip them.

**Why the deviation and not a second orbit.** The time comes from the variational equations, not from
the distance to a shadow orbit. That two-particle method gives a false positive exponent on regular and
quasi-regular orbits (Tancredi, Sánchez & Roig 2001, AJ 121, 1171). On Ceres around the Sun alone, a shadow
at 10⁻⁴ of the state, renormalised every 1000 years, levels off at 5·10⁻⁴ per year, i.e. 2000 years, while
`computeMegno` gives ⟨Y⟩ = 2.000 and `regular`. For Lyapunov times of tens of years both agree, to 0–2% as measured with a shadow particle.

Measured on an orbit with its aphelion near Jupiter (a = 3.6 AU, e = 0.45), whose growth phase starts
171 years back:

| Span | This fit (1σ) | Each quarter, from the start of the growth phase | Fitted from the start to each quarter's end | ⟨Y⟩ over the second half (before) |
| ---: | ---: | :--- | :--- | ---: |
| 200 years | not given (4 periods of growth) | | | 17 years |
| 400 years | 54 years (−29/+∞) | no growth, 19, 74, no growth | no growth, 26, 44, 54 | 34 years |
| 1000 years | 66 years (+30/−16) | 44, 42, 125, 174 | 44, 58, 55, 66 | 71 years |

The deviation is rescaled to unit length at the outputs once it exceeds 10⁸, so that it never overflows.
The variational equations are linear, so this does not change its growth, except for the round-off that
IAS15 keeps at the old scale in its compensated summation. That is why the limit is low: the jump of
ln|δ| at a rescale is 5·10⁻⁹ at 10⁸, but 0.6 at 10¹⁶ and 8 at 10²⁰, which shortened the Lyapunov time by
20–25% over 2000–4000 years.

Round-off also limits how far the orbit itself can be followed. Past about 36 Lyapunov times (ln 10¹⁶),
the integration follows an orbit like the one that started there, not that one. Even the seed of the
initial deviation changes it, since IAS15's convergence test includes the variational particles. On an
orbit with a Lyapunov time of about 20 years, five seeds gave the same time at 500 years, while their
semi-major axes differed by 10⁻³ AU by 520–790 years and their times by 20–33 years at 1000 years. The
same seed reproduces a run exactly.

**MEGNO needs tens of orbital periods to settle**, so set `--days` accordingly. For an orbit with
a = 3 AU (a period of about 5.2 years) that is at least 100 years, i.e. `--days 36525`. A run that is
too short gets no verdict rather than a wrong one.

Two caveats. Every verdict holds **over the integrated span only** — a sticky chaotic orbit can look
regular for a long time before it does not. And this is a gravity-only measurement: radiation forces
are dissipative, and MEGNO is not defined for them.

The extra integration costs roughly what a gravity-only run of the same span costs, which is small
next to a full run with the Monte Carlo clones. It is skipped, with a warning, if you call
`reboundSimulate` with `return_diagnostics=False`, since that is the only way the result comes back.

For example, Winchcombe integrated 1000 years back with 200 clones (from a covariance of 100 m and
1 m/s per axis):

```
  MEGNO of the nominal orbit (Newtonian gravity only, IAS15):
    <Y> at the end: 10.526   mean over the last quarter: 10.598
    span: 211 orbital periods
    Converges to 2: NO. <Y> does not converge to 2 and keeps growing: the orbit is chaotic.
    Lyapunov time ~ 69 years (+175/-29, 1 sigma), of the nominal orbit alone (for the clones, see
    below). Fitted to the growth of ln|delta| from 32 to 1000 years back,
    i.e. since <Y> last exceeded 2.5.
    Encounter-driven chaos grows in steps, so the same fit in each quarter of that span gives:
        32-274 years back: 24 years
       274-516 years back: 104 years
       516-758 years back: 48 years
      758-1000 years back: no growth
    The 1 sigma is the scatter of the quarters' growth rates over the square root of their number.
    Fitted from the onset to the end of each quarter: 24, 53, 52, 69 years.
    It still changes by 24% over the last half of the growth phase, so it has not settled: integrate for longer.
    The same over 32 clones of the run:
      Verdicts: 32 chaotic.
      Lyapunov time: median 46 years (16-84% of the clones: 30-68), from the 32 with a time.
      The nominal orbit is one draw among them, so this is the value to quote.
      From the first half of the run the median is 39 years: it changes by 15% when the span doubles,
      so this is a value for this span, not yet a limit.
```

Most of the growth comes in the first quarter, which holds the nominal orbit's deepest encounter with
the Earth (0.13 Hill radii, 43 years back); over the last quarter its deviation even shrinks. The fit
from the onset holds at 52–53 years over the middle quarters and ends at 69: close to settled, but not
yet. The clones say that this orbit is on the slow side: their median is 46 years.

---


---

## What it writes

All output goes into the directory containing the input pickle, under **fixed names**. Two runs on
the same pickle overwrite each other, so copy the results elsewhere before rerunning with different
options.

| File | Contents |
| :--- | :--- |
| `rebound_simulation_results.txt` | The human-readable report: run settings, the integrator that actually ran, final elements with confidence intervals, close encounters, impacts, clone outcomes, energy drift, divergence, and the MEGNO verdict. |
| `rebound_simulation.png` | The evolution of the orbit, 3×3 panels: a, e, q; i, node, argument of perihelion; longitude of perihelion, aphelion distance Q, and the close approaches to every body in Hill radii (log scale). See [The figure](#the-figure). |
| `rebound_megno_lyapunov.png` | Only with `--compute_megno`: the evolution of ⟨Y⟩ against time in years, with the lines at 2 and 2.5, and below it ln\|δ\| with the fit of the Lyapunov time over its growth phase (see [Chaos: MEGNO](#chaos-megno)). |
| `rebound_simulation_results.json` | The same results in machine-readable form, so downstream analysis does not have to parse the text. |

### The figure

- **Clones as bands.** With 20 clones or more, they are drawn as their 2σ and 1σ bands (the
  2.275–97.725% and 15.87–84.13% percentiles, those of a normal distribution) and median at each time, with 10 of them drawn individually to show how single trajectories behave; with
  fewer, every clone is drawn. The nominal solution is the black line. Clones that ended early (e.g.
  impacted) count only over their own span.
- **Angles within 0–360°, without jumps.** The longitude of the node Ω, the argument of perihelion ω
  and the longitude of perihelion ϖ are labelled within 0–360°, but drawn continuous, so crossing
  0/360° draws no vertical line. The axis covers the data alone when they span less than a turn, and
  its labels then run on through 0 (…, 330, 0, 30, …); otherwise it covers one full turn, 0–360°, and
  each line is drawn at every whole turn that falls on it. The bands are taken about the clones'
  circular mean at each time, so a cloud that straddles 0/360°, or lies away from the nominal
  solution, is not cut through. For a nearly ecliptic orbit Ω and ω are poorly defined, and
  ϖ = Ω + ω is the angle to read: clones whose inclination comes near 0 scatter their Ω and ω and
  widen those bands, while ϖ stays narrow. The report's 95% intervals of the final Ω, ω and f are
  taken about the clones' circular mean as well, and an interval that runs through 0 is given with
  low > high and marked "(across 0)".
- **From the departure.** The elements before the object leaves 3 Earth Hill radii (dotted vertical
  line) carry the Earth's pull and are not drawn, and the vertical axes are scaled to the rest of the
  run, so they cover the clones as well as the nominal solution.
- **Time** in years for spans of two years or more, in days otherwise.
- **Lyapunov time.** With `--compute_megno`, the Lyapunov time is the dark red dashed vertical line on
  every panel, if the run is longer, and the legend says where it comes from: the median over the MEGNO
  clones when they give one, otherwise the nominal orbit's (`lyapunovTimeToQuote`).
- **Aphelion distance.** The bottom middle panel gives Q = a(1 + e), with the clones' bands, next to the
  perihelion distance q of the top row. It is left out for an open orbit.
- **Close approaches in Hill radii.** The bottom right panel gives every distance in units of the
  Hill radius of the body, so that one dotted line is the encounter threshold (3 Hill radii by default)
  for all of them, and a Jupiter encounter at 0.9 AU (2.5 Hill radii) reads as deep as an Earth one at
  0.025 AU. It shows:
  - the nominal's **encounters** as stars in the body's colour;
  - the **clones' encounters** as faint dots, which show when and with which planets the ensemble meets,
    as opposed to the single nominal history;
  - for a body with none of these, its **closest approach** of the run as a ring ("never closer": the
    object never came closer to that body), if it is within 100 Hill radii; further out the body never
    mattered.

  The Sun has no Hill sphere: its distance is scaled so that its distance gate falls on the same line.
  In runs shorter than two years the distance to the Earth at the outputs is drawn too, as it shows the
  departure; in longer ones the output sampling misses the approaches, which the encounters give
  exactly.

The JSON keys, at the top level:

| Key | Contents |
| :--- | :--- |
| `traj_id` | Event identifier. |
| `run` | Run settings: `integration_days`, `direction`, `reference_frame`, `ephemeris`, `n_outputs`, `beta`, `start_epoch_jd_tdb` and `start_epoch_jd_utc`, `final_epoch_jd_tdb`, `final_epoch_jd_utc` and `final_epoch_utc` (ISO), `mc_runs`, `random_seed`, `runtime_s`, `integrator`, `dt_days`, `fixed_step_from_days`, and the close-encounter gates `n_hill_threshold` and `sun_encounter_au`. |
| `final_elements` | `a`, `q`, `e`, `incl_deg`, `peri_deg`, `node_deg`, `f_deg`, `tisserand_jupiter`, with `a_units`/`q_units` (AU heliocentric, km geocentric). |
| `encounters` | Every close encounter of the nominal solution, ordered by time, each with `body`, `min_dist_au`, `time_days`, `hill_radius_au` and `n_hill`. A body can appear more than once. For the Sun, `hill_radius_au` and `n_hill` are `null`. |
| `closest_approaches_au`, `closest_approach_times_days` | Closest approach to every tracked body, whether or not it counted as an encounter. This is the deepest approach only, one value per body. |
| `impact`, `escaped` | The impact or ejection of the nominal solution, or `null`. |
| `clone_outcomes` | Clone statistics: counts and fractions of impacts and encounters per body, `ci_uses_survivors_only`, `n_hill_threshold`, `sun_encounter_au`. Under `close_encounters`, `count` is how many clones met the body and `n_encounters` how many passages they made in total. `null` without `--mc`. |
| `clone_closest_approaches_au` | Per-clone closest approach to every body. |
| `clone_encounters` | Each clone's full list of encounters, in the same format as `encounters`. |
| `energy_rel_drift` | Relative energy drift of the massive subsystem, as an integrator-quality check. |
| `divergence` | The exponential and linear fits to the clone spread, with `r2_exponential` and `r2_linear`. |
| `whfast_encounter_warning` | The WHFast close-encounter warning string, or `null`. |
| `megno` | With `--compute_megno`: `times_days`, `megno`, `ln_delta` (logarithm of the deviation), `a_au`, `departure_days`, `escaped` and the `verdict` dict, with `lyapunov_time_years`, `lyapunov_time_range_years`, `lyapunov_onset_years` and `lyapunov_method`. `null` otherwise. |
| `representative` | With 10 clones or more over the whole span: the representative clone at the end of the run (see [The orbit to use at the end](#the-orbit-to-use-at-the-end)). `null` otherwise. |
| `nominal`, `monte_carlo` | Full element series: `time_days`, `a_au`, `e`, `incl_deg`, `peri_deg`, `node_deg`, `f_deg`, `body_distances_au`. |

A quick sanity check on any run is `energy_rel_drift`. It measures the massive subsystem the object
moves through, not the object itself (which is massless), but a large drift means the integration is
not to be trusted.

---

## The orbit to use at the end

The final elements of the nominal solution mean something on their own only while the run is shorter
than the Lyapunov time (see [Chaos: MEGNO](#chaos-megno)). Beyond it the nominal orbit is one
realization among the clones, no better than any of them: in Winchcombe integrated 1000 years back
with 200 clones, 86% of the clones end closer to the centre of the clones than the nominal does, about
what any one clone would give. Nor do the medians of the elements make an orbit: each comes from a
different clone, so that a(1 - e) is not the median q, and near i = 0 the median i is not even the
tilt of the central orbit plane (see below). In a run longer than the Lyapunov time the report therefore
flags the final elements with a warning: their central values are the nominal's (not a median), one
realization among the clones, while their intervals, the clones' spread, still hold. The Lyapunov time
comes from `--compute_megno`, so without it there is no warning.

With 10 clones or more over the whole span, the report then gives the representative orbit at the end,
next to the clones' median and 95% CI of each element, how far its initial state was from the nominal's,
and the share of the clones more central than the nominal at the end.

The representative is the clones' median in the space of the orbits, taken among the clones themselves.
The orbits are compared by their angular momentum and eccentricity vectors (the Milankovitch elements,
`orbitVectors`): six numbers that fix the size, the shape and the orientation of an orbit, but not where
the body is on it, with no wrap at 360 deg and no singularity at i = 0 or e = 0. Each direction is
measured in units of the clones' spread along it, so that the clones have a unit covariance (a
Mahalanobis distance; the direction normal to the surface h.e = 0, on which every orbit lies, is left
out). The representative is the clone with the smallest sum of distances to all the others, the medoid:
in one dimension that is the median, and here it is the median corrected for the spread in each direction
and for the correlations between the elements. Where the body is on its orbit (f) is not used: the clones
spread along the whole orbit long before their orbits do.

Being a clone, it is a real trajectory: the equations of motion integrated from a state drawn from the
measurement's uncertainty. The report gives the Mahalanobis distance of that state from the nominal's in
the spread of the clones' initial states, and the share of the draws that are closer (a chi-square in
six dimensions, in which a typical draw is 2.3 sigma away, not 1).

It is chosen at the end of the run only, and the report says so. Once the dynamics is chaotic, no
trajectory stays the most central at every time. The clone most central on average over the run is not a
better choice: with the clones
resampled (100 bootstrap draws), the final q of the representative varies by 0.008 AU, a sixth of the
clones' spread (0.050 AU), but that of the clone most central on average varies by 0.050 AU, as much as
the whole cloud. For a representative orbit at another epoch, end the run there.

Read the representative's elements together, as one orbit. Taken one at a time they can fall outside the
clones' 95% CI, and the report then says so. This happens at low inclination: i is the tilt of the orbit
plane whatever its direction, so for tilts spread around the pole of the frame the median i is larger than
the tilt of the central plane, and the node and the argument of perihelion are not defined on their own,
only their sum, the longitude of perihelion varpi, which the report adds. For example, Winchcombe 150
years back with 40 clones and `--compute_megno` (`--seed 20210228`; the pickle has no trajectory ID,
hence the clone names `None_MC_...`), as the text report gives it:

```
a    =   2.664466 +/- 0.012507 [  2.630404,   2.679430] AU
q    =   1.004697 +/- 0.029462 [  0.922710,   1.038202] AU
e    =   0.622928 +/- 0.009898 [  0.611316,   0.650118]
i    =   1.031627 +/- 0.250442 [  0.051990,   1.033724] deg
peri = 176.789616 +/- 22.587088 [176.737973, 243.430466] deg
node = -27.164289 +/- 22.197804 [266.203819, 332.840025] deg
f    = 152.999549 +/- 24.206251 [134.622213, 261.818563] deg
WARNING: the run is longer than the Lyapunov time (~36 years). The values above are the nominal orbit's
(not a median), one realization among the clones: they do not represent the orbit at the end. The intervals,
the clones' spread, do. For a real orbit that represents the clones, see the representative orbit below.

Representative orbit: clone None_MC_24, chosen at the end of the run only (150 years back).
How: the orbits of the 40 clones (their size, shape and orientation, not where the body is on them)
are compared two by two, each element in units of the clones' spread, and the representative is the
clone with the smallest sum of distances to all the others: the median of the clones, but a real orbit.
                 clone  clones: median [95% CI]
  a    =      2.661729        2.661108 [  2.630404,   2.679430]  AU
  q    =      1.002524        0.991919 [  0.922710,   1.038202]  AU
  e    =      0.623356        0.625776 [  0.611316,   0.650118]
  i    =      0.070011        0.123819 [  0.051990,   1.033724]  deg
  peri =    215.428699      190.230168 [176.737973, 243.430466]  deg
  node =    293.557626      318.609375 [266.203819, 332.840025]  deg
  varpi=    148.986325      148.916754 [147.506457, 149.670754]  deg
  f    =    142.234220      150.356284 [134.622213, 261.818563]  deg
It is a real trajectory, from an initial state 0.8 sigma from the nominal's (6-D: 0% of the draws
from the measurement are closer to the nominal).
At the end 95% of the clones are more central than the nominal orbit (about 50% for a typical clone).
```

There the representative's node and peri are 25 deg from their medians, in opposite directions, but its
varpi is 0.07 deg from the median, within an interval 2.2 deg wide.

`representativeOrbit` is saved under `representative` in the JSON: its `clone` name, its `elements`
(lengths in AU, angles in degrees, keyed `a_au`, `q_au`, `e`, `incl_deg`, `peri_deg`, `node_deg`,
`varpi_deg`, `f_deg`), the clones' `median`, `sigma` and `ci95` of each under `clones`, `nominal_rank`,
`start_sigma`, `start_dof`, `start_closer` and `outside_ci95`. A rank is the share of the clones more
central than the orbit. The figure draws the representative's path in fuchsia with a white edge: no other
line of the figure is pink, and the edge keeps it visible over the blue bands with protanopia or
deuteranopia, which turn the fuchsia into a dark grey-blue.

---

## Worked examples

**1. The default run — where did it come from in the last two months?**

```
python -m wmpl.Rebound.REBOUND traj.pickle
```

60 days backwards, IAS15, nominal solution only. Fast, and enough to see the geocentric encounter and
the pre-encounter heliocentric orbit.

**2. A century backwards, with enough sampling to see it.**

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 36525 --outputs 5000
```

The default 500 outputs over a century is one sample every ten weeks, which is too coarse to see a
planetary encounter in the plots. The encounter *detection* is unaffected — that runs at every
internal timestep — but the plots and the JSON series are.

**3. Uncertainty propagation.**

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 36525 --outputs 5000 --mc 200 --cores 8 --seed 42
```

200 clones on 8 cores. The report now gives 95 % confidence intervals on every element, and says how
many clones impacted something, escaped, or passed near a planet.

**4. A fast long run — and checking that it was safe.**

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 365250 --outputs 5000 --integrator trace --dt 0.5
```

A thousand years with TRACE, which resolves close encounters. Use `--integrator whfast` instead only
if you know the orbit stays away from the planets; if it does not, the report will tell you so, and
the run has to be redone. Before quoting a fixed-step result, compare a shorter run against IAS15 as
shown in [Choosing an integrator](#choosing-an-integrator).

**5. A small grain, where radiation forces matter.**

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 36525 --radius 1e-5 --density 3000
```

A 10 µm grain. The report prints the beta that was computed. Compare against the same run without
`--radius` to see how much the radiation forces actually moved the orbit.

**6. Is the orbit chaotic?**

```
python -m wmpl.Rebound.REBOUND traj.pickle --days 365250 --outputs 5000 --compute_megno
```

A thousand years, which covers enough orbital periods for a verdict for most meteoroid orbits. The
report ends with a MEGNO block:

```
  MEGNO of the nominal orbit (Newtonian gravity only, IAS15):
    <Y> at the end: 0.512   mean over the last quarter: 0.498
    span: 91.2 orbital periods
    Converges to 2: NO (tends to 0). <Y> tends to 0 rather than 2: the deviation from the orbit
    stayed bounded over this span, as around a stable periodic orbit (e.g. libration in a
    mean-motion resonance). An Earth-crossing orbit can also do this for a while before its first
    deep encounter: this is about the nominal orbit over this span, not the clones (see their
    Lyapunov time).
```

A chaotic orbit also gets a Lyapunov time, with the time fitted in each quarter of its growth phase
(see [Chaos: MEGNO](#chaos-megno) for an example).

If it says `too_short`, the message tells you how many years the run needs.

**7. Reading the results back in Python.**

```python
import json

with open("rebound_simulation_results.json") as f:
    res = json.load(f)

print("Integrator:", res["run"]["integrator"], "over", res["run"]["integration_days"], "days")
print("a = {:.4f} {:s}".format(res["final_elements"]["a"], res["final_elements"]["a_units"]))

for enc in res["encounters"]:
    # The Sun has no Hill sphere, so its "n_hill" is None
    hill_str = "" if enc["n_hill"] is None else " ({:.2f} Hill radii)".format(enc["n_hill"])
    print("{:s}: {:.5f} AU{:s} at t = {:.1f} d".format(
        enc["body"], enc["min_dist_au"], hill_str, enc["time_days"]))

if res["megno"] is not None:
    print("MEGNO verdict:", res["megno"]["verdict"]["status"])

# The full element series of the nominal solution
t_days = res["nominal"]["time_days"]
a_au = res["nominal"]["a_au"]
```

---

## Using it as a library

The command line is a wrapper around `reboundSimulate`, which can be called directly:

```python
import os
from wmpl.Utils.Pickling import loadPickle
from wmpl.Rebound.REBOUND import reboundSimulate

traj = loadPickle(*os.path.split("/path/to/traj.pickle"))

outputs, outputs_mc, diagnostics = reboundSimulate(
    None, None, traj=traj,
    direction="backward", sim_days=36525, n_outputs=5000,
    mc_runs=100, random_seed=42,
    integrator="trace", dt_days=0.5,
    compute_megno=True, return_diagnostics=True
    )

# Orbital elements at the end of the nominal integration
final = outputs[traj.traj_id][-1][2]
print("a = {:.4f} AU, e = {:.4f}".format(final.a, final.e))

# The MEGNO verdict lives in the nominal solution's diagnostics
print(diagnostics[traj.traj_id]["megno"]["verdict"]["status"])
```

`compute_megno=True` needs `return_diagnostics=True`, since the diagnostics are the only place the
result is returned.

**Units.** Each output row is `[t, state_vect_hel, orbit, planet_dists]`, in REBOUND's units: G = 1,
lengths in AU, masses in solar masses, so the time unit is sqrt(AU³/GM☉) = 58.13244 days
(`DAYS_PER_SIM_TIME`; REBOUND calls it `yr2pi`, a Gaussian year over 2π). Convert times with
`simTimeToDays(t)` and `daysToSimTime(days)`. Do not use 365.25/(2π), the Julian year over 2π: it is
1.9e-5 shorter, so with G = 1 it makes the Sun 3.8e-5 too massive and moves the Earth ~700 000 km off
the ephemeris in 10 years.

Other functions worth knowing about, all in `wmpl.Rebound.REBOUND`:

| Function | Purpose |
| :--- | :--- |
| `radiationPressureBeta(radius_m, density_kgm3)` | Beta for a spherical grain. |
| `equivalentSphereRadius(mass_kg, density_kgm3)` | The radius that beta needs, from a mass and a bulk density. |
| `encountersFromMinDistances(min_dist_au, min_time_days, n_hill=3.0)` | The close-encounter list from the closest-approach dictionaries, one entry per body, with the same gates as the integration (Hill radii, and `SUN_ENCOUNTER_AU` for the Sun). The per-particle `encounters` diagnostic lists every passage instead. |
| `sunPassageSummary(encounters, max_listed=3)` / `formatSunPassageSummary(summary)` | Collapse more than `max_listed` Sun passages into one summary, and its report line. |
| `cloneEncounterSummary(clone_diag)` | Per body, how many clones met it, how many passages they made, and the closest approach over all of them. |
| `hermiteClosestApproach(t0, r0, v0, t1, r1, v1)` | Closest approach of a relative trajectory inside one integrator step. |
| `whfastEncounterWarning(diagnostics, n_hill=3.0)` | The WHFast warning string, or `None`. |
| `computeMegno(task, seed=1)` / `classifyMegno(times_days, megno, a_au)` | MEGNO series and verdict. |
| `estimateLyapunovFromMC(sim_outputs, sim_outputs_mc)` | The divergence fit over the clone ensemble. |
| `checkIntegratorAvailable(integrator)` | Raises if the installed REBOUND does not have it. |
| `findEarthDepartureIndex(sim_outputs, n_hill=3.0)` | First output at which the object left the Earth's neighbourhood. |

---

## Troubleshooting

**`reboundx` will not install on Windows.** It has no Windows wheel and does not compile with MSVC.
Use WSL2. Installing only `rebound` is not enough.

**`Need to attach reboundx.Extras instance to simulation before setting params`,** or
**`REBOUNDx did not attach to the REBOUND simulation`.** `reboundx` was compiled against a different
`rebound` than the installed one. Recompile it with the `--force-reinstall` command in
[Requirements](#requirements).

**`import reboundx` fails with a missing symbol.** Same cause, same repair.

**`The installed REBOUND (4.3.0) does not provide the TRACE integrator.`** TRACE was added in REBOUND
4.4.0. Upgrade, or use `--integrator ias15`.

**`REBOUNDx ... cannot orient the Earth's J2/J4 along its rotation axis`.** REBOUNDx is older than
4.4.0. Upgrade it with the `--force-reinstall` command in [Requirements](#requirements).

**`--dt must be positive.`** A zero or negative step would silently collapse to a single step per
output interval, so it is refused rather than accepted.

**`WARNING: N close encounter(s) ... happened while WHFast was integrating.`** WHFast does not resolve
close encounters, so the orbit after them may be badly wrong. Rerun with `--integrator trace` or
`ias15`.

**`REBOUNDx: Passing a velocity-dependent force to WHFAST.`** Expected, and harmless here: the
measured effect of the general relativity correction on the orbit was below 1e-8 in the semi-major
axis over 100 years.

**The run hangs at the end.** It is not hung — the script ends with a blocking `plt.show()`. Close the
figure windows. The files are already written by that point.

**Results changed between two runs with the same arguments.** Pass `--seed`. Without it the Monte
Carlo seed is drawn at random (and reported, so an old run can be reproduced from its report).

**Two runs overwrote each other.** Output names are fixed and go next to the input pickle. Copy the
results out, or keep the pickles in separate directories.
