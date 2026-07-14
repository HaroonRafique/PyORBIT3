"""
Module. Includes classes for the analytical space charge calculator nodes.
"""

import sys
import os
import math
import time

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

#import the base DirectForce AccNode class
from orbit.space_charge.scAccNodes import SC_Base_AccNode

class SCanalyticalAccNode(SC_Base_AccNode):
	"""
	The subclass of the AccNodeBunchTracker class. It uses SpaceChargeFrozen wrapper for the c++ space charge calculator.
	"""
	def __init__(self, sc_calculator, lattice_functions = None, name = "no name"):
		"""
		Constructor. Creates the SC accelerator node element.
		"""
		SC_Base_AccNode.__init__(self, sc_calculator, name)
		self.setType("FrozenSC")
		self.lattice_functions = lattice_functions or {}

	def set_lattice_functions(self, lattice_functions):
		self.lattice_functions = lattice_functions

	def track(self, paramsDict):
		"""
		It is tracking the bunch through the Space Charge calculator.
		"""
		bunch = paramsDict["bunch"]
		beta_x = self.lattice_functions["betax"]
		beta_y = self.lattice_functions["betay"]
		eta_x  = self.lattice_functions.get("etax", 0.0)
		eta_y  = self.lattice_functions.get("etay", 0.0)
		co_x   = self.lattice_functions.get("orbitx", 0.0)
		co_y   = self.lattice_functions.get("orbity", 0.0)
		self.sc_calculator.setLatticeParameters(beta_x, beta_y, eta_x, eta_y, co_x, co_y)
		self.sc_calculator.trackBunch(bunch, self.sc_length)
