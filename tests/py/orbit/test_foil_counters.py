from orbit.core.bunch import Bunch
from orbit.core.foil import Foil


def _make_bunch():
    bunch = Bunch()
    bunch.mass(0.93827231)
    bunch.charge(1.0)
    bunch.getSyncParticle().kinEnergy(1.0)
    bunch.addParticle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    bunch.addParticle(10.0, 0.0, 10.0, 0.0, 0.0, 0.0)
    bunch.compress()
    return bunch


def test_foil_counter_api_exposes_initial_values_and_local_setters():
    foil = Foil(-1.0, 1.0, -1.0, 1.0, 1.0e-6)

    assert foil.getFoilHitsLocal() == 0
    assert foil.getFoilHitsGlobal() == 0
    assert foil.getFoilLossesLocal() == 0
    assert foil.getFoilLossesGlobal() == 0

    foil.setFoilHitsLocal(7)
    foil.setFoilLossesLocal(3)

    assert foil.getFoilHitsLocal() == 7
    assert foil.getFoilLossesLocal() == 3


def test_simple_scatter_updates_foil_hit_counters():
    foil = Foil(-1.0, 1.0, -1.0, 1.0, 1.0e-6)
    bunch = _make_bunch()

    foil.traverseFoilSimpleScatter(bunch)

    assert foil.getFoilHitsLocal() == 1
    assert foil.getFoilHitsGlobal() == 1
    assert foil.getFoilLossesLocal() == 0
    assert foil.getFoilLossesGlobal() == 0
