import math

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTuneAnalysis
from orbit.diagnostics import TeapotTuneAnalysisNode


def _particle_phase_attributes(bunch):
    return lambda index: bunch.partAttrValue("ParticlePhaseAttributes", 0, index)


def test_bunch_tune_analysis_accepts_vertical_dispersion_and_closed_orbit():
    bunch = Bunch()
    bunch.mass(0.93827231)
    bunch.getSyncParticle().kinEnergy(1.0)
    bunch.addParticle(1.2, 0.1, 2.3, 0.2, 0.0, 0.0)
    bunch.compress()
    analysis = BunchTuneAnalysis()

    analysis.assignTwiss(4.0, 0.5, 0.0, 0.0, 9.0, -0.25, 0.0, 0.0)
    analysis.assignClosedOrbit(1.0, 0.1, 2.0, 0.2)
    analysis.analyzeBunch(bunch)

    phase = _particle_phase_attributes(bunch)
    assert phase(4) > 0.0
    assert phase(5) > 0.0


def test_bunch_tune_analysis_vertical_dispersion_affects_y_action():
    bunch = Bunch()
    bunch.mass(0.93827231)
    bunch.getSyncParticle().kinEnergy(1.0)
    bunch.addParticle(0.0, 0.0, 1.5, 0.0, 0.0, 0.1)
    bunch.compress()
    analysis = BunchTuneAnalysis()

    analysis.assignTwiss(1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    analysis.analyzeBunch(bunch)
    action_without_etay = _particle_phase_attributes(bunch)(5)

    analysis.assignTwiss(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.5 / _dpp(bunch), 0.0)
    analysis.analyzeBunch(bunch)
    action_with_etay = _particle_phase_attributes(bunch)(5)

    assert action_without_etay > 0.0
    assert math.isclose(action_with_etay, 0.0, abs_tol=1.0e-12)


def test_teapot_tune_analysis_node_exposes_full_twiss_and_closed_orbit_api():
    bunch = Bunch()
    bunch.mass(0.93827231)
    bunch.getSyncParticle().kinEnergy(1.0)
    bunch.addParticle(1.2, 0.1, 2.3, 0.2, 0.0, 0.0)
    bunch.compress()
    node = TeapotTuneAnalysisNode("tune_analysis")

    node.assignTwiss(4.0, 0.5, 0.0, 0.0, 9.0, -0.25, 0.0, 0.0)
    node.assignClosedOrbit(1.0, 0.1, 2.0, 0.2)
    node.track({"bunch": bunch})

    phase = _particle_phase_attributes(bunch)
    assert phase(4) > 0.0
    assert phase(5) > 0.0


def _dpp(bunch):
    sync_part = bunch.getSyncParticle()
    beta = sync_part.beta()
    etot = sync_part.kinEnergy() + sync_part.mass()
    return bunch.dE(0) / (beta * beta * etot)
