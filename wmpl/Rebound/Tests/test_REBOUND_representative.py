""" Tests for the orbit to use at the end of a REBOUND run: the representative clone, the clones' median in
the space of the orbits (the angular momentum and eccentricity vectors in units of the clones' spread), and
the warning on the final elements of the nominal orbit past the Lyapunov time.

Synthetic outputs only, so they run without REBOUND or REBOUNDx installed.
"""

import os
import json
import importlib.util
from types import SimpleNamespace

import numpy as np
import pytest
import scipy.stats


REBOUND_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "REBOUND.py")


@pytest.fixture(scope="module")
def reb():
    """ Load REBOUND.py under an isolated module name. """

    spec = importlib.util.spec_from_file_location("_wmpl_test_rebound_representative", REBOUND_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module



@pytest.fixture
def agg():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.use("Agg")
    yield
    plt.close("all")


# REBOUND time unit in days (DAYS_PER_SIM_TIME)
DAYS_PER_UNIT = np.sqrt(149597870700.0**3/1.3271244004193938e20)/86400.0

# Spread of the clones in a, e, i, node and peri (degrees for the angles), and their centre
SCALES = np.array([0.05, 0.01, 0.5, 2.0, 2.0])
CENTRE = np.array([2.5, 0.6, 5.0, 100.0, 200.0])


def _orbit(x, frac):
    """ A REBOUND-like orbit from (a, e, i, node, peri), angles in degrees. """

    return SimpleNamespace(a=x[0], e=x[1], inc=np.radians(x[2]), Omega=np.radians(x[3] % 360.0), \
                           omega=np.radians(x[4] % 360.0), f=np.radians(360.0*frac % 360.0))



def _ensemble(n_clones=60, n_out=120, years=100.0, k_scramble=40, seed=4, nominal_after=None, \
              states=True):
    """ Backward outputs of a nominal orbit and its clones, whose elements keep the order of their
    initial conditions until output k_scramble and are drawn anew from then on, so that they lose
    that order there. The nominal sits at the centre of the clones until
    then, and at nominal_after (in units of the spread, default the centre) after. The object leaves
    the Earth at the third output. With states, each output carries a state vector, the clones' drawn
    about the nominal's from a known covariance. """

    t_days = -np.linspace(0, years*365.25, n_out)
    rng = np.random.default_rng(seed)
    before = rng.normal(size=(n_clones, 5))
    after = rng.normal(size=(n_clones, 5))
    start = rng.normal(size=(n_clones, 6))
    chol = np.diag([1e-9, 2e-9, 1e-9, 1e-4, 3e-4, 2e-4])
    chol[4, 3] = 1e-4

    def rows(z_before, z_after, w):
        out = []
        for k in range(n_out):
            frac = k/(n_out - 1)
            z = z_before*0.1 if k < k_scramble else z_after
            state = list(np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0]) + chol @ w) if states else None
            out.append([t_days[k]/DAYS_PER_UNIT, state, _orbit(CENTRE + SCALES*z, frac), \
                        {"Earth": 1e-5 if k < 2 else 1.0, "Sun": 1.0}])
        return out

    nominal = rows(np.zeros(5), np.zeros(5) if nominal_after is None else np.asarray(nominal_after, float), \
                   np.zeros(6))
    clones = {"obj_MC_{:d}".format(j): rows(before[j], after[j], start[j]) for j in range(n_clones)}

    return nominal, clones



def _kepler(a, e, inc, Omega, omega, f):
    """ Position and velocity with GM = 1 from the elements (radians). """

    p = a*(1 - e**2)
    r = p/(1 + e*np.cos(f))
    r_pf = r*np.array([np.cos(f), np.sin(f), 0.0])
    v_pf = np.sqrt(1/p)*np.array([-np.sin(f), e + np.cos(f), 0.0])
    cO, sO = np.cos(Omega), np.sin(Omega)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(omega), np.sin(omega)
    rot = np.array([[cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci, sO*si], \
                    [sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si], \
                    [sw*si, cw*si, ci]])

    return rot @ r_pf, rot @ v_pf



@pytest.mark.parametrize("elements", [(2.5, 0.6, 0.3, 1.0, 4.0, 0.7), (-3.0, 1.4, 2.8, 5.5, 0.4, 0.3)])
def testOrbitVectorsAreTheAngularMomentumAndEccentricityVectors(reb, elements):
    """ sqrt(p) times the orbit normal is r x v with GM = 1, and the eccentricity vector is v x h - r/|r|,
    for an ellipse and for a retrograde hyperbola. """

    r, v = _kepler(*elements)
    h = np.cross(r, v)
    e_vec = np.cross(v, h) - r/np.linalg.norm(r)

    vec = reb.orbitVectors(*elements[:5])

    assert np.allclose(vec[:3], h, atol=1e-12)
    assert np.allclose(vec[3:], e_vec, atol=1e-12)
    assert abs(np.dot(vec[:3], vec[3:])) < 1e-12



def testOrbitVectorsHaveNoWrapAndNoSingularityInThePlane(reb):
    """ A whole turn in the node changes nothing, and in the plane of the frame only the longitude of
    the perihelion (node + peri) counts. """

    a, e = 2.5, 0.6
    assert np.allclose(reb.orbitVectors(a, e, 0.2, 0.1, 1.0), reb.orbitVectors(a, e, 0.2, 0.1 + 2*np.pi, 1.0))
    assert np.allclose(reb.orbitVectors(a, e, 0.0, 0.3, 1.0), reb.orbitVectors(a, e, 0.0, 1.1, 0.2), \
        atol=1e-14)



def testSpreadUnitsGiveTheClonesAUnitCovarianceInFiveDimensions(reb):
    """ The direction normal to h.e = 0 is left out, and the other five have unit covariance. """

    rng = np.random.default_rng(1)
    x = CENTRE + SCALES*rng.normal(size=(300, 5))
    vec = reb.orbitVectors(x[:, 0], x[:, 1], np.radians(x[:, 2]), np.radians(x[:, 3]), np.radians(x[:, 4]))

    centre, matrix = reb._spreadUnits(vec)
    z = (vec - centre) @ matrix

    assert z.shape == (300, 5)
    assert np.allclose(np.cov(z.T), np.eye(5), atol=1e-9)



def testRepresentativeIsTheMedoidClone(reb):
    """ The representative is one of the clones, with its own elements at the end, and the most central
    of them: no clone has a smaller sum of distances in units of the clones' spread. """

    nominal, clones = _ensemble()
    rep = reb.representativeOrbit(nominal, clones)

    assert rep["clone"] in clones
    assert rep["n_clones"] == len(clones)
    end = clones[rep["clone"]][-1][2]
    assert rep["elements"]["a_au"] == pytest.approx(end.a)
    assert rep["elements"]["q_au"] == pytest.approx(end.a*(1 - end.e))
    assert rep["elements"]["varpi_deg"] == pytest.approx(np.degrees(end.Omega + end.omega) % 360.0)

    names = list(clones)
    a, e, inc, Omega, omega, _ = reb._orbitArrays([clones[n][-1][2] for n in names])
    _, ranks, _ = reb._centralityRanks(reb.orbitVectors(a, e, inc, Omega, omega))
    assert ranks[names.index(rep["clone"])] == 0
    assert np.sum(ranks == 0) == 1



def testRepresentativeIsTheMedianCorrectedForTheSpread(reb):
    """ With the inclination the narrowest element, the medoid of the unscaled orbit vectors moved four
    sigma off in i stays their medoid, since i barely moves the vectors, but it is no longer the
    representative, which measures each direction in units of the clones' spread. """

    nominal, clones = _ensemble(n_clones=40, seed=7)
    for rows in clones.values():
        o = rows[-1][2]
        o.inc = np.radians(CENTRE[2] + (np.degrees(o.inc) - CENTRE[2])*0.02)
    names = list(clones)

    def unscaledMedoid():
        a, e, inc, Omega, omega, _ = reb._orbitArrays([clones[n][-1][2] for n in names])
        vec = reb.orbitVectors(a, e, inc, Omega, omega)
        return names[int(np.argmin(reb._distanceSums(vec, vec)))]

    medoid = unscaledMedoid()
    clones[medoid][-1][2].inc = np.radians(CENTRE[2] + 4*SCALES[2]*0.02)

    assert unscaledMedoid() == medoid
    assert reb.representativeOrbit(nominal, clones)["clone"] != medoid



def testNominalRankIsLowAtTheCentreAndHighWhenOffIt(reb):
    """ The share of clones more central than the nominal: small while the nominal is the centre, large
    when it ends three sigma off in every element. """

    nominal, clones = _ensemble(k_scramble=200)
    assert reb.representativeOrbit(nominal, clones)["nominal_rank"] < 0.1

    nominal, clones = _ensemble(nominal_after=[3, 3, 3, 3, 3])
    assert reb.representativeOrbit(nominal, clones)["nominal_rank"] > 0.9



def testStartSigmaIsTheMahalanobisDistanceOfTheInitialState(reb):
    """ The distance of the representative's initial state from the nominal's in the spread of the
    clones' initial states, and the chi-square share of draws closer than it. """

    nominal, clones = _ensemble()
    rep = reb.representativeOrbit(nominal, clones)

    offsets = np.array([rows[0][1] for rows in clones.values()]) - np.array(nominal[0][1])
    d = np.array(clones[rep["clone"]][0][1]) - np.array(nominal[0][1])
    expected = np.sqrt(d @ np.linalg.solve(np.cov(offsets.T), d))

    assert rep["start_dof"] == 6
    assert rep["start_sigma"] == pytest.approx(expected, rel=1e-4)
    assert rep["start_closer"] == pytest.approx(scipy.stats.chi2.cdf(expected**2, 6), rel=1e-4)



def testStartSigmaNeedsTheStateVectors(reb):
    """ Outputs read back from a results file have no state vectors. """

    nominal, clones = _ensemble(states=False)
    rep = reb.representativeOrbit(nominal, clones)

    assert rep["start_sigma"] is None
    assert rep["clone"] in clones



def testRepresentativeNeedsTenClonesOverTheWholeSpan(reb):
    """ Clones that ended early do not count. """

    nominal, clones = _ensemble(n_clones=12)
    for name in list(clones)[:3]:
        clones[name] = clones[name][:50]

    assert reb.representativeOrbit(nominal, clones) is None



def testAnglesSpreadAcrossZeroKeepTheirInterval(reb):
    """ A node spread across 0/360 deg has its median near 0 and an interval through 0. """

    rng = np.random.default_rng(3)
    orbits = [_orbit(CENTRE + SCALES*z - np.array([0, 0, 0, 100.0, 0]), 0.5) \
        for z in rng.normal(size=(200, 5))]
    spread = reb._cloneSpread(orbits)

    lo, hi = spread["node_deg"]["ci95"]
    assert lo > hi
    assert min(spread["node_deg"]["median"], 360.0 - spread["node_deg"]["median"]) < 1.0
    assert spread["node_deg"]["sigma"] == pytest.approx(2.0, rel=0.15)



def testResultsAreJsonSerialisable(reb):
    """ The representative goes into the results file. """

    nominal, clones = _ensemble()
    rep = reb.representativeOrbit(nominal, clones)

    back = json.loads(json.dumps(rep))
    assert back["clone"] == rep["clone"]
    assert back["elements"]["e"] == pytest.approx(rep["elements"]["e"])



def testReportWithoutRepresentative(reb):
    """ Without enough clones there is no representative and nothing to report. """

    assert reb._representativeReportLines(None) == []



def testReportGivesTheRepresentative(reb):
    """ The report gives the representative clone, how it was chosen, its distance from the nominal at the
    start, and how central the nominal is at the end. """

    nominal, clones = _ensemble(nominal_after=[3, 3, 3, 3, 3])
    rep = reb.representativeOrbit(nominal, clones)
    text = "\n".join(reb._representativeReportLines(rep, dist_unit_multiplier=1.0, length_units="AU"))

    assert "Representative orbit: clone {:s}, chosen at the end of the run only".format(rep["clone"]) in text
    assert "How: the orbits of the {:d} clones".format(len(clones)) in text
    assert "smallest sum of distances to all the others" in text
    assert "varpi=" in text
    assert "sigma from the nominal's (6-D" in text
    assert "more central than the nominal" in text



def testReportExplainsAnElementOutsideItsInterval(reb):
    """ Near i = 0 the representative's i can fall outside the clones' interval of i: the report says
    the orbit is central as a whole and why. """

    nominal, clones = _ensemble()
    rep = reb.representativeOrbit(nominal, clones)
    rep["outside_ci95"] = ["incl_deg"]
    rep["clones"]["incl_deg"]["ci95"] = [0.05, 0.9]
    text = " ".join(reb._representativeReportLines(rep))

    assert "Its i is outside the clones' 95% CI taken alone" in text
    assert "only varpi" in text



def testFigureDrawsTheRepresentative(reb, agg):
    """ The representative's path is drawn on the element panels when there is one, and not without. """

    nominal, clones = _ensemble()
    rep = reb.representativeOrbit(nominal, clones)
    _, axes = reb.plotReboundOrbits(nominal, clones, [], representative=rep)

    lines = [l for l in axes["a"].get_lines() if l.get_label() == "representative clone"]
    assert len(lines) == 1
    assert lines[0].get_color() == reb.PLOT_REPRESENTATIVE_COLOR
    assert lines[0].get_path_effects() == reb.PLOT_REPRESENTATIVE_EDGE
    a_rep = np.array([row[2].a for row in clones[rep["clone"]]])
    assert np.allclose(lines[0].get_ydata()[2:], a_rep[2:])

    _, axes = reb.plotReboundOrbits(nominal, clones, [])
    assert not [l for l in axes["a"].get_lines() if l.get_label() == "representative clone"]



def testFinalElementsAreFlaggedPastTheLyapunovTime(reb):
    """ In a run longer than the Lyapunov time the final elements of the nominal are flagged as one
    realization, not a median nor a representative orbit; in a shorter run, or without a time, they are
    not. """

    span = -150*365.25
    text = " ".join(reb._finalElementsWarningLines(36.0, span))
    assert text.startswith("WARNING: the run is longer than the Lyapunov time (~36 years)")
    assert "The values below are the nominal orbit's (not a median)" in text
    assert "do not represent the orbit at the end" in text
    assert "above are" in " ".join(reb._finalElementsWarningLines(36.0, span, above=True))

    assert reb._finalElementsWarningLines(200.0, span) == []
    assert reb._finalElementsWarningLines(None, span) == []
