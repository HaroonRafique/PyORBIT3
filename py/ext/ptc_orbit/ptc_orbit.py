"""PTC lattice helpers for the integrated PyORBIT3 + PTC build."""

import math
from os import fspath

from orbit.teapot import BaseTEAPOT, TEAPOT_Lattice
from orbit.utils import consts

try:
    from pylibptc_orbit import *  # noqa: F401,F403
except ImportError as exc:
    raise ImportError(
        "ext.ptc_orbit requires PyORBIT3 to be built with PTC enabled "
        "and the pylibptc_orbit extension installed."
    ) from exc


def _ptc_path(path):
    return fspath(path)


def _ptc_path_len(path):
    return len(path) - 1


class PTC_Lattice(TEAPOT_Lattice):
    """TEAPOT-compatible lattice backed by PTC flat-file data."""

    def __init__(self, name="no name"):
        super().__init__(name)

    def readPTC(self, PTC_File):
        """Read a PTC flat file and initialize lattice nodes."""
        ptc_file = _ptc_path(PTC_File)
        self.setName(ptc_file)
        ptc_init_(ptc_file, _ptc_path_len(ptc_file))

        (
            self.betax0,
            self.betay0,
            self.alphax0,
            self.alphay0,
            self.etax0,
            self.etapx0,
            self.etay0,
            self.etapy0,
            self.orbitx0,
            self.orbitpx0,
            self.orbity0,
            self.orbitpy0,
        ) = ptc_get_twiss_init_()

        self.nNodes, self.nHarm, self.lRing, self.gammaT = ptc_get_ini_params_()

        for node_index in range(self.nNodes):
            (
                length,
                betax,
                betay,
                alphax,
                alphay,
                etax,
                etapx,
                etay,
                etapy,
                orbitx,
                orbitpx,
                orbity,
                orbitpy,
            ) = ptc_get_twiss_for_node_(node_index)

            elem = PTC_Node("PTC_Node")
            elem.setparams(
                node_index,
                length,
                betax,
                betay,
                alphax,
                alphay,
                etax,
                etapx,
                etay,
                etapy,
                orbitx,
                orbitpx,
                orbity,
                orbitpy,
            )
            self.addNode(elem)

        self.initialize()


class PTC_Node(BaseTEAPOT):
    """PTC-backed lattice node."""

    def __init__(self, name="ptc_node"):
        super().__init__(name)
        self.setType("ptc_node")

    def setparams(
        self,
        orbit_ptc_node_index,
        length,
        betax,
        betay,
        alphax,
        alphay,
        etax,
        etapx,
        etay,
        etapy,
        orbitx,
        orbitpx,
        orbity,
        orbitpy,
    ):
        self.addParam("node_index", orbit_ptc_node_index)
        self.setLength(length)
        self.addParam("betax", betax)
        self.addParam("betay", betay)
        self.addParam("alphax", alphax)
        self.addParam("alphay", alphay)
        self.addParam("etax", etax)
        self.addParam("etapx", etapx)
        self.addParam("etay", etay)
        self.addParam("etapy", etapy)
        self.addParam("orbitx", orbitx)
        self.addParam("orbitpx", orbitpx)
        self.addParam("orbity", orbity)
        self.addParam("orbitpy", orbitpy)

    def track(self, paramsDict):
        bunch = paramsDict["bunch"]
        phase_length = paramsDict["length"]
        orbit_ptc_node_index = self.getParam("node_index")
        ptc_trackBunch(bunch, phase_length, orbit_ptc_node_index)


def setBunchParamsPTC(bunch):
    """Set synchronous particle parameters from PTC."""
    mass, charge, kin_energy = ptc_get_syncpart_()
    mass *= consts.mass_proton
    sync_part = bunch.getSyncParticle()
    sync_part.kinEnergy(kin_energy)
    bunch.charge(charge)
    bunch.mass(mass)


def readAccelTablePTC(acc_file):
    acc_file = _ptc_path(acc_file)
    ptc_read_accel_table_(acc_file, _ptc_path_len(acc_file))


def readScriptPTC(script_file):
    script_file = _ptc_path(script_file)
    ptc_script_(script_file, _ptc_path_len(script_file))


def updateParamsPTC(lattice, bunch):
    (
        lattice.betax0,
        lattice.betay0,
        lattice.alphax0,
        lattice.alphay0,
        lattice.etax0,
        lattice.etapx0,
        lattice.etay0,
        lattice.etapy0,
        lattice.orbitx0,
        lattice.orbitpx0,
        lattice.orbity0,
        lattice.orbitpy0,
    ) = ptc_get_twiss_init_()

    lattice.nNodes, lattice.nHarm, lattice.lRing, lattice.gammaT = ptc_get_ini_params_()

    for node in lattice.getNodes():
        node_index = node.getParam("node_index")
        (
            _,
            betax,
            betay,
            alphax,
            alphay,
            etax,
            etapx,
            etay,
            etapy,
            orbitx,
            orbitpx,
            orbity,
            orbitpy,
        ) = ptc_get_twiss_for_node_(node_index)

        node.setParam("betax", betax)
        node.setParam("betay", betay)
        node.setParam("alphax", alphax)
        node.setParam("alphay", alphay)
        node.setParam("etax", etax)
        node.setParam("etapx", etapx)
        node.setParam("etay", etay)
        node.setParam("etapy", etapy)
        node.setParam("orbitx", orbitx)
        node.setParam("orbitpx", orbitpx)
        node.setParam("orbity", orbity)
        node.setParam("orbitpy", orbitpy)

    setBunchParamsPTC(bunch)


def _require_ptc_status_ok(operation, status):
    if status != 0:
        raise RuntimeError(f"{operation} failed with PTC status {status}")


def fibre_count():
    """Return the number of native PTC fibres in the active loaded layout."""
    status, count = ptc_get_fibre_count_()
    _require_ptc_status_ok("ptc_get_fibre_count_", status)
    return count


def fibre_name(index):
    """Return the native PTC fibre name for a 1-based fibre index."""
    status, name = ptc_get_fibre_name_(index)
    _require_ptc_status_ok("ptc_get_fibre_name_", status)
    return name


def set_fibre_aperture(index, kind, half_x, half_y, x_offset=0.0, y_offset=0.0, r1=0.0, r2=0.0):
    """Set a native PTC fibre aperture.

    This helper is intentionally thin; it is mainly useful for tests and for
    explicit user-side aperture application before querying readback state.
    """
    status = ptc_set_fibre_aperture_(
        index,
        kind,
        r1,
        r2,
        half_x,
        half_y,
        x_offset,
        y_offset,
    )
    _require_ptc_status_ok("ptc_set_fibre_aperture_", status)


def fibre_aperture(index):
    """Return native PTC aperture data for one fibre, or None if inactive."""
    status, kind, r1, r2, half_x, half_y, x_offset, y_offset, s = ptc_get_fibre_aperture_(index)
    if status in (2, 3):
        return None
    _require_ptc_status_ok("ptc_get_fibre_aperture_", status)
    if kind <= 0 or not math.isfinite(half_x) or not math.isfinite(half_y):
        return None
    return {
        "index": index,
        "name": fibre_name(index),
        "s": s,
        "kind": kind,
        "half_x": half_x,
        "half_y": half_y,
        "x_offset": x_offset,
        "y_offset": y_offset,
        "r1": r1,
        "r2": r2,
        "source": "PTC fibre aperture",
    }


def all_fibre_apertures():
    """Return active native PTC fibre aperture rows from the loaded layout."""
    rows = []
    for index in range(1, fibre_count() + 1):
        row = fibre_aperture(index)
        if row is not None:
            rows.append(row)
    return rows


def synchronousSetPTC(ival):
    if ival >= 0:
        raise ValueError("synchronousSetPTC requires ival < 0")
    ptc_synchronous_set_(ival)


def synchronousAfterPTC(ival):
    if ival >= 0:
        raise ValueError("synchronousAfterPTC requires ival < 0")
    ptc_synchronous_after_(ival)


def trackBunchThroughLatticePTC(lattice, bunch, phase_length):
    params_dict = {"bunch": bunch, "length": phase_length}
    for node in lattice.getNodes():
        node.track(params_dict)
