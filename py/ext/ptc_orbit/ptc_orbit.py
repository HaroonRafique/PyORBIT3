"""
Module. Includes classes for all PTC elements.
Built on the back of TEAPOT by J. Holmes.
"""

import sys
import os
import math

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# import some constants
from orbit.utils import consts

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import the teapot classes
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import BaseTEAPOT

# import the interface to PTC
from libptc_orbit import *


class PTC_Lattice(TEAPOT_Lattice):
    """
    PTC Subclass of the AccLattice class.
    Inherits from the TEAPOT lattice.
    """

    def __init__(self, name="no name"):
        super().__init__(name)

    def readPTC(self, PTC_File):
        """
        Reads the PTC file input and initializes all structures.
        Input PTC_File is the flat PTC file.
        """
        self.setName(PTC_File)
        length_of_name = len(PTC_File)
        ptc_init_(PTC_File, length_of_name - 1)
        
        (
            self.betax0, self.betay0, self.alphax0, self.alphay0,
            self.etax0, self.etapx0, self.etay0, self.etapy0,
            self.orbitx0, self.orbitpx0, self.orbity0, self.orbitpy0
        ) = ptc_get_twiss_init_()

        (self.nNodes, self.nHarm, self.lRing, self.gammaT) = ptc_get_ini_params_()

        for node_index in range(self.nNodes):
            (
                length, betax, betay, alphax, alphay, etax, etapx,
                etay, etapy, orbitx, orbitpx, orbity, orbitpy
            ) = ptc_get_twiss_for_node_(node_index)

            elem = PTC_Node("PTC_Node")
            elem.setparams(
                node_index, length, betax, betay, alphax, alphay,
                etax, etapx, etay, etapy, orbitx, orbitpx, orbity, orbitpy
            )
            self.addNode(elem)

        self.initialize()


class PTC_Node(BaseTEAPOT):
    """
    PTC element.
    """

    def __init__(self, name="ptc_node"):
        """
        Constructor. Creates a PTC element.
        """
        super().__init__(name)
        self.setType("ptc_node")

    def setparams(self, orbit_ptc_node_index, length,
                  betax, betay, alphax, alphay,
                  etax, etapx, etay, etapy, orbitx, orbitpx, orbity, orbitpy):
        """
        Sets element parameters.
        """
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
        """
        The PTC class implementation of the AccNodeBunchTracker class track(probe) method.
        """
        bunch = paramsDict["bunch"]
        PhaseLength = paramsDict["length"]
        orbit_ptc_node_index = self.getParam("node_index")
        action_type = -1
        # ptc_get_task_type_(orbit_ptc_node_index, action_type)

        if action_type == 1:
            print("===============================")
            print("PTC_Node.track.")
            print("Energy change actions have not been taken.")
            print("STOP.")
            raise SystemExit(1)

        ptc_trackBunch(bunch, PhaseLength, orbit_ptc_node_index)


def setBunchParamsPTC(bunch):
    """
    Sets the synchronous particle parameters of the bunch.
    """
    mass, charge, kin_energy = ptc_get_syncpart_()
    mass *= consts.mass_proton
    syncPart = bunch.getSyncParticle()
    syncPart.kinEnergy(kin_energy)
    bunch.charge(charge)
    bunch.mass(mass)


def readAccelTablePTC(Acc_File):
    """
    Gets the information for acceleration.
    """
    ptc_read_accel_table_(Acc_File, len(Acc_File) - 1)


def readScriptPTC(Script_File):
    """
    Reads a PTC Script file.
    """
    ptc_script_(Script_File, len(Script_File) - 1)


def updateParamsPTC(lattice, bunch):
    """
    Updates Twiss parameters of lattice.
    Updates element parameters.
    Updates synchronous particle parameters of the bunch.
    """
    (
        lattice.betax0, lattice.betay0, lattice.alphax0, lattice.alphay0,
        lattice.etax0, lattice.etapx0, lattice.etay0, lattice.etapy0,
        lattice.orbitx0, lattice.orbitpx0, lattice.orbity0, lattice.orbitpy0
    ) = ptc_get_twiss_init_()

    (lattice.nNodes, lattice.nHarm, lattice.lRing, lattice.gammaT) = ptc_get_ini_params_()

    for node in lattice.getNodes():
        node_index = node.getParam("node_index")
        (
            _, betax, betay, alphax, alphay, etax, etapx,
            etay, etapy, orbitx, orbitpx, orbity, orbitpy
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


def synchronousSetPTC(ival):
    """
    Calls ptc_synchronous_set_.
    """
    if ival >= 0:
        print("===============================")
        print("synchronousSetPTC requires ival < 0")
        print("STOP.")
        raise SystemExit(1)

    ptc_synchronous_set_(ival)


def trackBunchThroughLatticePTC(lattice, bunch, PhaseLength):
    """
    Tracks a bunch through the whole lattice.
    """
    paramsDict = {"bunch": bunch, "length": PhaseLength}
    for node in lattice.getNodes():
        node.track(paramsDict)


def trackBunchInRangePTC(lattice, bunch, PhaseLength, indexi, indexf):
    """
    Tracks a bunch from indexi through indexf, inclusive.
    """
    paramsDict = {"bunch": bunch, "length": PhaseLength}
    for node in lattice.getNodes():
        node_index = node.getParam("node_index")
        if indexi <= node_index <= indexf:
            node.track(paramsDict)
