from pathlib import Path

import pytest


def _ptc_flat_file():
    root = Path(__file__).resolve().parents[3]
    candidates = [
        root.parent / "PTC" / "ptc_standalone_readiness" / "inputs" / "PTC-PyORBIT_flat_file.madx.flt",
        root.parent / "ptc_pyorbit" / "PTC-PyORBIT_flat_file.madx.flt",
    ]
    for path in candidates:
        if path.exists():
            return path
    pytest.skip("PTC flat file fixture is unavailable")


def _load_ptc_lattice():
    ptc_orbit = pytest.importorskip("ext.ptc_orbit.ptc_orbit")
    lattice = ptc_orbit.PTC_Lattice()
    lattice.readPTC(_ptc_flat_file())
    return ptc_orbit, lattice


def test_all_fibre_apertures_returns_native_rows_after_setting_rectangular_aperture():
    ptc_orbit, _lattice = _load_ptc_lattice()

    ptc_orbit.set_fibre_aperture(
        1,
        kind=2,
        half_x=0.012,
        half_y=0.034,
        x_offset=0.001,
        y_offset=-0.002,
    )

    rows = ptc_orbit.all_fibre_apertures()
    row = next(row for row in rows if row["index"] == 1)

    assert row["name"]
    assert row["kind"] == 2
    assert row["half_x"] == pytest.approx(0.012)
    assert row["half_y"] == pytest.approx(0.034)
    assert row["x_offset"] == pytest.approx(0.001)
    assert row["y_offset"] == pytest.approx(-0.002)
    assert row["s"] >= 0.0
    assert row["source"] == "PTC fibre aperture"


def test_all_fibre_apertures_is_stable_across_twiss_update():
    ptc_orbit, _lattice = _load_ptc_lattice()

    ptc_orbit.set_fibre_aperture(1, kind=2, half_x=0.021, half_y=0.022)
    before = ptc_orbit.all_fibre_apertures()
    ptc_orbit.ptc_update_twiss_()
    after = ptc_orbit.all_fibre_apertures()

    assert after == before
