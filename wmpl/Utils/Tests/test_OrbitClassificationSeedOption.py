""" Tests that --seed of the OrbitClassification command line serves both --clones and --dynamics. """

import os
import subprocess
import sys

import pytest


EXAMPLE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "Dynesty", "examples",
    "20191023_091225", "20191023_091225_trajectory.pickle")


def _run(*args):

    return subprocess.run([sys.executable, "-W", "ignore", "-m", "wmpl.Utils.OrbitClassification", EXAMPLE]
        + list(args), capture_output=True, text=True)


@pytest.mark.skipif(not os.path.isfile(EXAMPLE), reason="example trajectory not available")
def test_seed_goes_with_dynamics():

    # The span is refused before integrating, so this reaches the dynamics without running them
    result = _run("--dynamics", "shober2024", "--years", "2", "--seed", "1")

    assert result.returncode == 1
    assert "no output after the start" in result.stderr + result.stdout
    assert "--seed goes with" not in result.stderr + result.stdout


@pytest.mark.skipif(not os.path.isfile(EXAMPLE), reason="example trajectory not available")
def test_seed_alone_is_still_refused():

    result = _run("--seed", "1")

    assert result.returncode == 1
    assert "--seed goes with --clones or --dynamics" in result.stderr + result.stdout
