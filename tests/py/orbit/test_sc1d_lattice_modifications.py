from orbit.core.spacecharge import LSpaceChargeCalc
from orbit.space_charge.sc1d import addLongitudinalSpaceChargeNode
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice


def test_lspacechargecalc_exposes_phase_length():
    calc = LSpaceChargeCalc(1.0, 12.5, 100, 1, 16)

    assert calc.getLength() == 12.5


def test_add_longitudinal_space_charge_accepts_raw_calculator():
    lattice = AccLattice("sc1d")
    drift = DriftTEAPOT("drift")
    drift.setLength(2.0)
    lattice.addNode(drift)
    lattice.initialize()
    calc = LSpaceChargeCalc(1.0, lattice.getLength(), 100, 1, 16)

    nodes = addLongitudinalSpaceChargeNode(lattice, 0.5, calc)

    assert len(nodes) == 3
    assert nodes[1].getType() == "long sc node"
    assert nodes[1].getLength() == 0.0
    assert nodes[1].lspacecharge is calc
    assert lattice.getLength() == 2.0
