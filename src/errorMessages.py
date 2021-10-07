#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana and Michael Johnston
#
#  Created: 2005-10-03 by Adrian L. Garcia-Lomana
#
#  This application is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public
#  License as published by the Free Software Foundation; either
#  version 2 of the License, or (at your option) any later version.
#
#  This application is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Library General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this library; if not, write to the Free
#  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 USA.
# 

# $Id: errorMessages.py,v 4.8 2008/12/16 21:35:06 alglomana Exp $

## \file
# This module deals with the error handling of the program.

import sys

class ClassByoDynException:

	'''
	The mother class of the exceptions of the program.
	'''

	def __init__(self, errorString):

		'''
		Constructor of the class.
		'''

		self.errorString = errorString

	def printExceptionInfo(self):

		'''
		This method prints the specific error string.
		'''

		print self.errorString


class ClassCentralException(ClassByoDynException):

	'''
	This class deals with the specific errors of the central module.
	'''

	def printExceptionInfo(self):

		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from Central module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassCheckerException(ClassByoDynException):

	'''
	This class deals with the specific errors of the checker module from benchmark directory.
	'''
	
	def printExceptionInfo(self):

		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from Checker module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassClusterException(ClassByoDynException):

	'''
	This class deals with the specific errors of the cluster module.
	'''

	def printExceptionInfo(self):

		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from Cluster module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassDynamicsReconstructerException(ClassByoDynException):

	'''
	This class deals with the specific errors of the dynamicsReconstructer module.
	'''

	def printExceptionInfo(self):
		
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''

		print 'ERROR from DynamicsReconstruncter module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassFormulasException(ClassByoDynException):

	'''
	This class deals with the specific errors of the formulas module of the library directory.
	'''

	def printExceptionInfo(self):
		
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''

		print 'ERROR from Formulas module of lib directory:',
		ClassByoDynException.printExceptionInfo(self)

		return None
	
class ClassIdentifiabilityAnalyzerException(ClassByoDynException):
	
	'''
	This class deals with the specific errors of the identifiabilityAnalyzer module.
	'''

	def printExceptionInfo(self):
	
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from IdentifiabilityAnalyzer module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassInitiatorException(ClassByoDynException):
	
	'''
	This class deals with the specific errors of the initiator module at the bin directory.
	'''

	def printExceptionInfo(self):
	
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from Initiator module:',
		ClassByoDynException.printExceptionInfo(self)

		return None
	
class ClassMatrixWorkerException(ClassByoDynException):

	'''
	This class deals with the specific errors of the formulas module of the library directory.
	'''

	def printExceptionInfo(self):
		
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''

		print 'ERROR from MatrixWorker module of lib directory:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassOptimalExperimentalDesignException(ClassByoDynException):
	
	'''
	This class deals with the specific errors of the optimalExperimentalDesign module.
	'''

	def printExceptionInfo(self):
		
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from OptimalExperimentalDesign module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassOptimiserException(ClassByoDynException):

	'''
	This class deals with the specific errors of the optimiser module.
	'''

	def printExceptionInfo(self):

		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from Optimiser module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassSamplerException(ClassByoDynException):
	
	'''
	This class deals with the specific errors of the sampler module.
	'''
	
	def printExceptionInfo(self):
		
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from sampler module:',
		ClassByoDynException.printExceptionInfo(self)

		return None
	
class ClassSBMLWorkerException(ClassByoDynException):

	'''
	This class deals with the specific errors of the sbmlWorker module.
	'''

	def printExceptionInfo(self):

		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from SBMLWorker module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassSimulatorEulerException(ClassByoDynException):
	
	'''
	This class deals with the specific errors of the simulatorEuler module.
	'''
	
	def printExceptionInfo(self):
		
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from simulatorEuler module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassSimulatorStochasticException(ClassByoDynException):

	'''
	This class deals with the specific errors of the simulatorStochastic module.
	'''

	def printExceptionInfo(self):

		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from simulatorStochastic module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassSensitivityAnalyzerException(ClassByoDynException):
	
	'''
	This class deals with the specific errors of the sensitivityAnalyzer module.
	'''

	def printExceptionInfo(self):
	
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from SensitivityAnalyzer module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassSimulatorException(ClassByoDynException):

	'''
	This class deals with the specific errors of the simulator module.
	'''

	def printExceptionInfo(self):
		
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from Simulator module:',
		ClassByoDynException.printExceptionInfo(self)

		return None
	
class ClassSimulatorXPPException(ClassByoDynException):

	'''
	This class deals with the specific errors of the simulator module.
	'''

	def printExceptionInfo(self):
		
		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''
		
		print 'ERROR from SimulatorXPP module:',
		ClassByoDynException.printExceptionInfo(self)

		return None
		

class ClassStarterException(ClassByoDynException):
	
	'''
	This class deals with the specific errors of the starter module.
	'''

	def printExceptionInfo(self):

		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''

		print 'ERROR from Starter module:',
		ClassByoDynException.printExceptionInfo(self)

		return None

class ClassSurfaceException(ClassByoDynException):

	'''
	This class deals with the specific errors of the surface module.
	'''

	def printExceptionInfo(self):

		'''
		This method prints the specific error string corresponding to the module and then the specific string of the error.
		'''

		print 'ERROR from Surface module:',
		ClassByoDynException.printExceptionInfo(self)		
		
		return None
		

def ErrorHandler(type, value, traceback):

	'''
	This function is responsible of providing a nice format of the error message when an exception happens.
	'''

	try:
		value.printExceptionInfo()	
	except AttributeError:	
		sys.exit()

	return None
