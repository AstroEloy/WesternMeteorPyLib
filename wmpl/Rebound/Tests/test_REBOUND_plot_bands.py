""" Tests for the bands of the REBOUND figure: the clones' 1 and 2 sigma bands, the angles drawn
continuous and labelled within 0-360 deg with their bands about the clones' circular mean, the aphelion
panel, and the Lyapunov time on every panel.

Synthetic outputs only, so they run without REBOUND or REBOUNDx installed.
"""

import os
import importlib.util
from types import SimpleNamespace

import numpy as np
import pytest


REBOUND_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "REBOUND.py")


@pytest.fixture(scope="module")
def reb():
    """ Load REBOUND.py under an isolated module name. """

    spec = importlib.util.spec_from_file_location("_wmpl_test_rebound_plot_bands", REBOUND_PATH)
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


def _outputs(node_func, n_clones=40, n_out=200, years=100.0, seed=2):
    """ Backward outputs whose node is node_func(fraction of the run, clone index or None for the
    nominal) in degrees; the object leaves the Earth at the third output. """

    t_days = -np.linspace(0, years*365.25, n_out)
    rng = np.random.default_rng(seed)

    def rows(j):
        da = 0.0 if j is None else 0.05*rng.normal()
        out = []
        for k in range(n_out):
            frac = k/(n_out - 1)
            orbit = SimpleNamespace(a=2.0 + da*frac, e=0.6, inc=np.radians(1.0), \
                Omega=np.radians(node_func(frac, j) % 360.0), omega=np.radians(100.0), f=0.0)
            out.append([t_days[k]/DAYS_PER_UNIT, None, orbit, {"Earth": 1e-5 if k < 2 else 1.0, "Sun": 1.0}])
        return out

    return rows(None), {"obj_MC_{:d}".format(j): rows(j) for j in range(n_clones)}



def testCircularPercentilesOfACloudAcross360(reb):
    """ A cloud centred on 355 deg with a 10 deg spread: its percentiles are taken about the mean, so
    the 1 sigma band is 20 deg wide, where percentiles taken in 0-360 deg would span the circle. """

    rng = np.random.default_rng(0)
    stack = ((355.0 + 10.0*rng.normal(size=(4000, 3))) % 360.0)
    p = reb._circularPercentiles(stack, reb.PLOT_SIGMA_PERCENTILES)

    assert np.all(np.abs(((p[2] - 355.0 + 180.0) % 360.0) - 180.0) < 1.0)
    assert np.all(np.abs((p[3] - p[1]) - 20.0) < 1.5)
    assert np.all(np.abs((p[4] - p[0]) - 40.0) < 3.0)

    naive = np.percentile(stack, [15.866, 84.134], axis=0)
    assert np.all(naive[1] - naive[0] > 300.0)



def testCircularPercentilesAreAboutTheCloudNotOneClone(reb):
    """ The cut at 180 deg from the centre falls opposite the cloud, wherever a single clone is: a
    clone 170 deg away does not split the cloud. """

    rng = np.random.default_rng(1)
    stack = (100.0 + 5.0*rng.normal(size=(2000, 1))) % 360.0
    stack[0] = 270.0
    p = reb._circularPercentiles(stack, reb.PLOT_SIGMA_PERCENTILES)

    assert p[2, 0] == pytest.approx(100.0, abs=0.5)
    assert p[3, 0] - p[1, 0] == pytest.approx(10.0, abs=1.0)



def testCircularPercentilesFollowAMeanThatCircles(reb):
    """ A cloud whose mean goes once round the circle gives continuous percentiles. """

    rng = np.random.default_rng(2)
    centre = np.linspace(0.0, 720.0, 300)
    stack = (centre + 3.0*rng.normal(size=(500, 300))) % 360.0
    p = reb._circularPercentiles(stack, reb.PLOT_SIGMA_PERCENTILES)

    assert np.max(np.abs(np.diff(p[2]))) < 5.0
    assert p[2, -1] - p[2, 0] == pytest.approx(720.0, abs=2.0)



def testTurnShiftsCoverTheWindow(reb):

    assert reb._turnShifts([np.array([355.0, 365.0])], 350.0, 370.0) == [0.0]
    assert reb._turnShifts([np.array([355.0, 365.0])], -10.0, 10.0) == [360.0]
    assert reb._turnShifts([np.array([10.0, 700.0])], 0.0, 360.0) == [0.0, 360.0]
    assert reb._turnShifts([np.array([np.nan])], 0.0, 360.0) == [0.0]



def testBandsAreOneAndTwoSigma(reb, agg):

    nominal, clones = _outputs(lambda frac, j: 30.0)
    _, axes = reb.plotReboundOrbits(nominal, clones, [])

    labels = axes["a"].get_legend_handles_labels()[1]
    assert {"clones 2$\\sigma$", "clones 1$\\sigma$", "clones median"} <= set(labels)
    assert reb.PLOT_SIGMA_PERCENTILES == pytest.approx((2.275, 15.866, 50.0, 84.134, 97.725), abs=1e-3)



def testElementsAreLabelledByNameAndSymbol(reb, agg):

    nominal, clones = _outputs(lambda frac, j: 30.0)
    _, axes = reb.plotReboundOrbits(nominal, clones, [])

    assert axes["a"].get_ylabel() == "semi-major axis $a$ [AU]"
    assert axes["e"].get_ylabel() == "eccentricity $e$"
    assert axes["node"].get_ylabel() == "longitude of the node $\\Omega$ [deg]"
    assert axes["peri"].get_ylabel() == "argument of perihelion $\\omega$ [deg]"
    assert axes["varpi"].get_ylabel() == "longitude of perihelion $\\varpi$ [deg]"



def testAWideAngleCloudTakesOneFullTurn(reb, agg):
    """ Clones whose node spreads over the whole circle: the axis is one turn, 0-360 deg, and every
    line on it stays within it. """

    nominal, clones = _outputs(lambda frac, j: 40.0 + (0.0 if j is None else 600.0*frac*np.sin(j)))
    fig, axes = reb.plotReboundOrbits(nominal, clones, [])
    fig.canvas.draw()

    lo, hi = axes["node"].get_ylim()
    assert (lo, hi) == pytest.approx((0.0, 360.0))
    ticks = [y for y in axes["node"].get_yticks() if lo <= y <= hi]
    assert ticks == pytest.approx([0.0, 90.0, 180.0, 270.0, 360.0])



def testFewClonesAcross360AreContinuous(reb, agg):
    """ Without bands each clone is drawn: a clone whose node crosses 360 deg draws no vertical line
    (no step larger than its own motion between outputs). """

    nominal, clones = _outputs(lambda frac, j: 350.0 + 20.0*frac + (0.0 if j is None else 2.0*j), n_clones=5)
    _, axes = reb.plotReboundOrbits(nominal, clones, [])

    lo, hi = axes["node"].get_ylim()
    series = [line for line in axes["node"].get_lines() if len(line.get_ydata()) == len(nominal)]
    assert len(series) >= 6
    for line in series:
        y = np.array(line.get_ydata(), dtype=float)
        y = y[np.isfinite(y)]
        if np.any((y > lo) & (y < hi)):
            assert np.max(np.abs(np.diff(y))) < 1.0



def testAphelionPanelIsQ(reb, agg):
    """ The aphelion distance Q = a(1 + e) has its own panel, with the clones' bands, and is left out for an
    open orbit, which has none. """

    nominal, clones = _outputs(lambda frac, j: 30.0)
    _, axes = reb.plotReboundOrbits(nominal, clones, [])

    nominal_line = [line for line in axes["Q"].get_lines() if line.get_linewidth() > 1][0]
    q_big = np.array(nominal_line.get_ydata(), dtype=float)
    assert q_big[2:] == pytest.approx(2.0*(1 + 0.6))
    assert len(axes["Q"].collections) >= 2

    for row in nominal:
        row[2].e = 1.2
    _, axes = reb.plotReboundOrbits(nominal, {}, [])
    nominal_line = [line for line in axes["Q"].get_lines() if line.get_linewidth() > 1][0]
    assert np.all(np.isnan(np.array(nominal_line.get_ydata(), dtype=float)))



def testCircularIntervalOfTheReport(reb):
    """ The report's 95% interval of an angle: about the circular mean, with low > high across 0, where
    percentiles in 0-2 pi would cover almost the whole circle; away from 0 it is the usual one. """

    rng = np.random.default_rng(4)
    across = np.radians((350.0 + 10.0*rng.normal(size=20000)) % 360.0)
    lo, hi = reb._circularInterval(across)
    assert np.degrees(lo) == pytest.approx(350.0 - 19.6, abs=0.6)
    assert np.degrees(hi) == pytest.approx(350.0 + 19.6 - 360.0, abs=0.6)
    assert lo > hi

    naive = np.degrees(np.percentile(across, [2.5, 97.5]))
    assert naive[1] - naive[0] > 300.0

    away = np.radians(120.0 + 10.0*rng.normal(size=20000))
    assert reb._circularInterval(away) == pytest.approx(tuple(np.percentile(away, [2.5, 97.5])), abs=1e-6)



def testTheLyapunovTimeIsDarkRedOnEveryPanel(reb, agg):
    """ The Lyapunov time is drawn in dark red on every panel, the encounters' included. """

    from matplotlib.colors import to_hex

    nominal, clones = _outputs(lambda frac, j: 30.0, n_clones=25)
    _, axes = reb.plotReboundOrbits(nominal, clones, [], lyapunov_years=31.0, \
                                    lyapunov_text="median of the clones")

    for name in ("a", "node", "Q", "encounters"):
        dashed = [x for x in axes[name].get_lines() if x.get_linestyle() == "--" and len(x.get_xdata()) == 2]
        assert any(to_hex(x.get_color()) == to_hex(reb.PLOT_LYAPUNOV_COLOR) for x in dashed)
