import math

import pytest

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import (
    ConstantLineDensityProfile,
    GaussianLineDensityProfile,
    InterpolatedLineDensityProfile,
    SpaceChargeCalcAnalyticGaussian,
)
from orbit.lattice import AccLattice
from orbit.space_charge.analytical import setSCanalyticalAccNodes
from orbit.teapot import DriftTEAPOT


def _make_bunch():
    bunch = Bunch()
    bunch.mass(0.93827231)
    bunch.charge(1.0)
    bunch.macroSize(1.0)
    bunch.getSyncParticle().kinEnergy(1.0)
    bunch.addParticle(1.0e-3, 0.0, 2.0e-3, 0.0, 0.0, 0.0)
    bunch.compress()
    return bunch


def test_line_density_profiles_expose_expected_values():
    gaussian = GaussianLineDensityProfile(0.2)
    assert gaussian.getLocalLineDensityFactor(0.0) > gaussian.getLocalLineDensityFactor(0.2)

    constant = ConstantLineDensityProfile(2.0)
    assert constant.getLocalLineDensityFactor(-10.0) == pytest.approx(0.5)
    assert constant.getLocalLineDensityFactor(10.0) == pytest.approx(0.5)

    interpolated = InterpolatedLineDensityProfile(-1.0, 1.0, [0.0, 2.0, 4.0])
    assert interpolated.getLocalLineDensityFactor(-2.0) == pytest.approx(0.0)
    assert interpolated.getLocalLineDensityFactor(-0.5) == pytest.approx(1.0)
    assert interpolated.getLocalLineDensityFactor(0.0) == pytest.approx(2.0)


def test_analytical_gaussian_calculator_tracks_bunch_with_finite_kick():
    profile = GaussianLineDensityProfile(0.2)
    calc = SpaceChargeCalcAnalyticGaussian(1.0e11, 1.0e-6, 2.0e-6, 1.0e-4, profile)
    calc.setLatticeParameters(10.0, 12.0, 0.0, 0.0, 0.0, 0.0)
    bunch = _make_bunch()

    calc.trackBunch(bunch, 0.25)

    assert math.isfinite(bunch.xp(0))
    assert math.isfinite(bunch.yp(0))
    assert bunch.xp(0) != pytest.approx(0.0)
    assert bunch.yp(0) != pytest.approx(0.0)


def test_set_analytical_space_charge_nodes_inserts_frozen_sc_nodes():
    lattice = AccLattice("analytical-sc")
    drift = DriftTEAPOT("drift")
    drift.setLength(2.0)
    drift.addParam("betax", 10.0)
    drift.addParam("betay", 12.0)
    lattice.addNode(drift)
    lattice.initialize()
    calc = SpaceChargeCalcAnalyticGaussian(
        1.0e11,
        1.0e-6,
        2.0e-6,
        1.0e-4,
        GaussianLineDensityProfile(0.2),
    )

    nodes = setSCanalyticalAccNodes(lattice, 0.5, calc)

    assert len(nodes) == 1
    assert nodes[0].getType() == "FrozenSC"
    assert nodes[0].getLength() == 0.0
    assert nodes[0].sc_length == pytest.approx(2.0)
    assert nodes[0].lattice_functions["betax"] == pytest.approx(10.0)
    assert nodes[0].lattice_functions["etay"] == pytest.approx(0.0)
