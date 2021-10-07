#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido, Miguel Hernandez-Sanchez, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Alex Gomez-Garrido
#
#  Created: 2007-07-10 by Alex Gomez-Garrido
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

# $Id: testers.py,v 4.33 2008/12/16 16:02:33 alglomana Exp $

## \file
# This module contais the different tests available for ByoDyn

import os, re, sys, libsbml, shutil, scipy
import central

class ClassTest:

    '''
    Class for the tests.
    '''

    def __init__(self):

	'''
	The constructor.
	'''

        self.expectedResult = None
	self.obtainedResult = None
	self.id = None
	self.runner = None

        return None

    def __execute(self):

	'''
	This method executes ByoDyn with the specific test.
	'''
	
        print 'RUNNING TEST %s ...'%(self.id)
	#/ 1.- removing first the output directory and make it again
	outputdir = os.environ.get('BYODYN_OUTPUT_DIR')
	shutil.rmtree(outputdir)
	os.mkdir(outputdir)
	os.mkdir(outputdir + '/scratch')
	#/ 2.- executing
	try:
            central.main(self.runner)
	except:
            return False 
        else:
            return True 

    def dataFilesChecker(self):

	'''
	This method compares if the numerical outputs of ByoDyn differ on 1 per cent.
	'''

	#/ 1.- reading
        r = open(self.obtainedResult, 'r')
        e = open(self.expectedResult, 'r')
        resultLines = r.readlines()
        r.close()
        expectedLines = e.readlines()
        e.close()
	#/ 2.- creating variables
	results = []
	expected = []
	for line in resultLines:
	    if line[0] != '#':
		resultsVector = line.split()
		for value in resultsVector:
		    results.append(float(value))
	for line in expectedLines:
	    if line[0] != '#':
		expectedVector = line.split()
		for value in expectedVector:
		    expected.append(float(value))
	#/ 3.- comparing the values
	for i in range(len(results)):
	    #/ 3.1.- in the case of being zero
	    if results[i] == 0.0:
		if expected[i] != 0.0:
		    return False
            elif expected[i] == 0.0:
                if results[i] != 0.0:
                    return False
	    #/ 3.2.- maximum difference, 1 per cent
	    else:
		ratio = results[i] / expected[i]
		difference = abs(1.0 - ratio)
		if difference >= 0.01:
		    return False
	return True

    def setOutputName(self, name):
    
	'''
	This method sets the path to the obtained output that is going to be checked.
	'''
	
	outputdir = os.environ.get('BYODYN_OUTPUT_DIR')

	return '%s/%s'%(outputdir, name)

    def run(self):

	'''
	This method evaluates if the test ran correctly'
	'''
	
        if self.__execute() == True:
            if self.analyseResults() == True:
                return 'OK' #/ message to be printed
            else:
                return 'WRONG RESULTS' #/ message to be printed
        else:
            return 'EXECUTION FAILED' #/ message to be printed   

    def setRunnerName(self, name):
	
	'''
	This method sets the option file for the test.
	'''
	
        byodynPath = os.environ.get('BYODYN_PATH')
        folderName = '%s/benchmark/runners' %byodynPath
        
	return folderName + '/' + name

    def setResultName(self, name):

	'''
	This method sets the output directory of the results of the tests.
	'''

        byodynPath = os.environ.get('BYODYN_PATH')
        folderName = '%s/benchmark/expectedResults' %byodynPath
        
	return folderName + '/' + name
            
#/
#/ From here on, we implement the tests derived from ClassTest.
#/  

class ClassCluster2DTest(ClassTest):

    '''
    Class for testing the 2D clustering.
    '''

    def analyseResults(self):

	'''
	This method analyses the results of the test for the 2D clustering.
	'''
	
	passedTest = True
	#/ 1.- checking that the results have the appropriated structure
	file = self.obtainedResult + '/clusterResults.txt'
	f = open(file, 'r')
	data = f.readlines()
	f.close()
	counter = 0
	for datum in data:
	    if datum[0] != '#':
		counter = counter + 1
		vector = datum.split()
		if len(vector) != 3:
		    passedTest = False
	if counter != 600:
	    passedTest = False
	#/ 2.- check for the number of centers
	counter = 0
	file = self.obtainedResult + '/clusterCenters.txt'
	f = open(file, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		counter = counter + 1
	if counter != 2:
	    passedTest = False
	#/ 3.- check for the existence of the figure
	file = self.obtainedResult + '/clusterResults.ps'
	size = os.stat(file)[6]
	if size == 0:
	    passedTest = False

	return passedTest

class ClassCluster3DTest(ClassTest):

    '''
    Class for testing the 3D clustering.
    '''

    def analyseResults(self):

	'''
	This method analyses the results of the test for the 3D clustering.
	'''
	
	passedTest = True
	#/ 1.- checking that the results have the appropriated structure
	file = self.obtainedResult + '/clusterResults.txt'
	f = open(file, 'r')
	data = f.readlines()
	f.close()
	counter = 0
	for datum in data:
	    if datum[0] != '#':
		counter = counter + 1
		vector = datum.split()
		if len(vector) != 4:
		    passedTest = False
	if counter != 4600:
	    passedTest = False
	#/ 2.- check for the number of centers
	counter = 0
	file = self.obtainedResult + '/clusterCenters.txt'
	f = open(file, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		counter = counter + 1
	if counter != 9:
	    passedTest = False
	#/ 3.- check for the existence of the figure
	file = self.obtainedResult + '/clusterResults.ps'
	size = os.stat(file)[6]
	if size == 0:
	    passedTest = False

	return passedTest

class ClassExportingTest(ClassTest):

    '''
    Class for testing exporting.
    '''

    def __init__(self):
	
	'''
	The constructor.
	'''
	
	self.id = 'exporting'
	self.runner = self.setRunnerName('exporting.rn')
	
	return None

    def analyseResults(self):

	'''
	This method analyses the results of the test for exporting.
	'''
	
	value = None
	#/ 1.- retrieving the parameters
	outputdir = os.environ.get('BYODYN_OUTPUT_DIR')
	file = outputdir + '/' + 'Kinetic_modelling_of_Amadori_degradation.xml'
	modelFile = libsbml.readSBML(file)
	sbmlModel = modelFile.getModel()
	reactions = sbmlModel.getListOfReactions()
	for reaction in reactions:
	    listOfLocalParameters = reaction.getKineticLaw().getListOfParameters()
	    for parameter in listOfLocalParameters:
		if parameter.getId() == 'k4':
		    value = parameter.getValue()
	#/ 2.- comparing with the proper result
	if value == 0.07:
	    return True
	else:
	    return False

class ClassFigureFormatTest(ClassTest):

    '''
    Class for testing the output format of graphics.
    '''
    
    def analyseResults(self):

	'''
	This method analyses the format of the output graphics.
	'''

	size = os.stat(self.obtainedResult)[6]
	if size == 0:
	    return False
	else:
	    return True

class ClassFitnessFunctionCalculationTest(ClassTest):
    
    '''
    Class for testing the fitness function calculation functionality.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses for the result of the fitness function calculation.
	'''
	
	f = open(self.obtainedResult, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		vector = datum.split()
	ratio = float(vector[len(vector)-1]) / .25
	difference = abs(1.0 - ratio)
	if difference >= 0.01:
	    return False
	else:
	    return True

class ClassFitnessFunctionSurfaceTest(ClassTest):
    
    '''
    Class for testing the fitness function surface functionality.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the fitness function surface functionality.
	'''
	
	passedTest = True
	#/ 1.- checking the existence of the figure
	file = self.obtainedResult + 'ps'
	if os.path.exists(file) == False:
	    passedTest = False
	#/ 2.- checking the values of the grid
	self.obtainedResult = self.obtainedResult + 'txt'
	result = self.dataFilesChecker()
	if result == False:
	    passedTest = False
	
	return passedTest

class ClassGeneticAlgorithmTest(ClassTest):
    
    '''
    Class for testing the results of the genetic algorithm.
    '''

    def analyseResults(self):
	
	'''
	This method analyses the results of the genetic algorithm.
	It checks that the number of iteration is 5, that the best value of the iterations are getting lower and lower and finally both directories named "parametersFromRandom" and "parametersFromGA" exist.
	'''
	
	#/ 1.- checking the iterations of the genetic algorithm.
	values = []
	passedTest = True
	f = open(self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.st', 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    vector = datum.split()
	    if vector[0] != '#':
		values.append(float(vector[3]))
	if len(values) != 5:
	    passedTest = False
	for i in range(len(values)):
	    if i != 0:
		if values[i-1] < values[i]:
		    passedTest = False	    
	#/ 2.- checking for the directory "parametersFromGA". Directory "parametersFromRandom" does not exist as described in Section 2.3.4.3.3
	if os.path.exists(self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation/parametersFromGA') == False:
	    passedTest = False
	
	return passedTest

class ClassHybridOnePhaseTest(ClassTest):
    
    '''
    Class for testing the results of the hybrid one phase optimisation.
    '''
    
    def analyseResults(self):

	'''
	This method analyses the results of the hybrid one phase optimisation.
	'''
	
	passedTest = True
	#/ 1.- checking that there are two generations
	file = self.obtainedResult + '/Kinetic_modelling_of_Amadori_degradation.st'
	f = open(file, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		generationNumber = int(datum.split()[1])
	if generationNumber != 2:
	    passedTest = False
	#/ 2.- checking that the parameter values of the xml are the global minimum
	file = self.obtainedResult + '/Kinetic_modelling_of_Amadori_degradation/Kinetic_modelling_of_Amadori_degradation.xml'
	modelFile = libsbml.readSBML(file)
	sbmlModel = modelFile.getModel()
	reactions = sbmlModel.getListOfReactions()
	for reaction in reactions:
	    listOfLocalParameters = reaction.getKineticLaw().getListOfParameters()
	    for parameter in listOfLocalParameters:
		if parameter.getId() == 'k10':
		    k10 = parameter.getValue()
		if parameter.getId() == 'k2':
		    k2 = parameter.getValue()
	ratio = k10 / 0.0707
	difference = abs(1.0 - ratio)
	if difference >= 0.01:
	    passedTest = False
	ratio = k2 / 0.0156
	difference = abs(1.0 - ratio)
	if difference >= 0.01:
	    passedTest = False	    

	return passedTest

class ClassHybridTwoPhasesTest(ClassTest):

    '''
    Class for testing the results of the hybrid two phases optimisation.
    '''

    def analyseResults(self):
	
	'''
	This method analyses the results of the hybrid two phases optimisation.
	'''

	def parameterValueObtainer(data):
	    
	    '''
	    This method retrieves the parameter value of the best fitness function result.
	    '''

	    fitnessValue = None
	    parameterValue = None
	    for datum in data:
		if datum[0] != '#':
		    if fitnessValue == None:
			fitnessValue = float(datum.split()[3])
			parameterValue = float(datum.split()[0])
		    else:
			if float(datum.split()[3]) < fitnessValue:
			    fitnessValue = float(datum.split()[3])
			    parameterValue = float(datum.split()[0])
	    
	    return parameterValue
	
	#/ 0.- starting the analyseResults method
	passedTest = True
	#/ 1.- checking that the value of the parameters are the global minimum
	k10Value = None
	k2Value = None
	file = self.obtainedResult + 'parametersFromHybridTwoPhases/k10v10'
	f = open(file, 'r')
	data = f.readlines()
	f.close()
	k10Value = parameterValueObtainer(data)
	file = self.obtainedResult + 'parametersFromHybridTwoPhases/k2v2'
	f = open(file, 'r')
	data = f.readlines()
	f.close()
	k2Value = parameterValueObtainer(data)
	ratio = k10Value / 0.0707
	difference = abs(1.0 - ratio)
	if difference >= 0.01:
	    passedTest = False
	ratio = k2Value / 0.0156
	difference = abs(1.0 - ratio)
	if difference >= 0.01:
	    passedTest = False
	#/ 2.- checking that the sbml contains the parameter values obtained
	file = self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.xml'
	modelFile = libsbml.readSBML(file)
	sbmlModel = modelFile.getModel()
	reactions = sbmlModel.getListOfReactions()
	for reaction in reactions:
	    listOfLocalParameters = reaction.getKineticLaw().getListOfParameters()
	    for parameter in listOfLocalParameters:
		if parameter.getId() == 'k10':
		    k10 = parameter.getValue()
		if parameter.getId() == 'k2':
		    k2 = parameter.getValue()
	ratio = k10 / 0.0707
	difference = abs(1.0 - ratio)
	if difference >= 0.01:
	    passedTest = False
	ratio = k2 / 0.0156
	difference = abs(1.0 - ratio)
	if difference >= 0.01:
	    passedTest = False

	return passedTest

class ClassIdentifiabilityTest(ClassTest):
    
    '''
    Class for testing the functionalities of identifiability.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the identifiability.
	'''
	
	passedTest = True
	#/ 1- checking for the ps figure
	file = self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.correlation.ps'
	if os.path.exists(file) == False:
	    passedTest = False
	else:
	    size = os.stat(file)[6]
	    if size == 0:
		passedTest = False
	#/ 2.- checking for the text files
	#/ 2.1.- FIM
	obtainedData = []
	expectedData = []
	f = open(self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.FIM.txt', 'r')
	obtainedResults = f.readlines()
	f.close()
	f = open(self.expectedResult + 'fim.txt', 'r')
	expectedResults = f.readlines()
	f.close()
	for datum in obtainedResults:
	    if datum[0] != '#':
		vector = datum.split()
		if len(vector) == 3:
		    obtainedData.append(float(vector[1]))
		    obtainedData.append(float(vector[2]))
	for datum in expectedResults:
	    if datum[0] != '#':
		vector = datum.split()
		if len(vector) == 3:
		    expectedData.append(float(vector[1]))
		    expectedData.append(float(vector[2]))
	for i in range(len(expectedData)):
	    ratio = expectedData[i] / obtainedData[i]
	    difference = abs(1.0 - ratio)
	    if difference >= 0.01:
		passedTest = False
	#/ 2.2.- COV
	obtainedData = []
	expectedData = []
	f = open(self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.COV.txt', 'r')
	obtainedResults = f.readlines()
	f.close()
	f = open(self.expectedResult + 'cov.txt', 'r')
	expectedResults = f.readlines()
	f.close()
	for datum in obtainedResults:
	    if datum[0] != '#':
		vector = datum.split()
		if len(vector) == 3:
		    obtainedData.append(float(vector[1]))
		    obtainedData.append(float(vector[2]))
	for datum in expectedResults:
	    if datum[0] != '#':
		vector = datum.split()
		if len(vector) == 3:
		    expectedData.append(float(vector[1]))
		    expectedData.append(float(vector[2]))
	for i in range(len(expectedData)):
	    ratio = expectedData[i] / obtainedData[i]
	    difference = abs(1.0 - ratio)
	    if difference >= 0.01:
		passedTest = False
	#/ 2.3.- correlation
	obtainedData = []
	expectedData = []
	f = open(self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.correlation.txt', 'r')
	obtainedResults = f.readlines()
	f.close()
	f = open(self.expectedResult + 'correlation.txt', 'r')
	expectedResults = f.readlines()
	f.close()
	for datum in obtainedResults:
	    if datum[0] != '#':
		vector = datum.split()
		if len(vector) == 3:
		    obtainedData.append(float(vector[1]))
		    obtainedData.append(float(vector[2]))
	for datum in expectedResults:
	    if datum[0] != '#':
		vector = datum.split()
		if len(vector) == 3:
		    expectedData.append(float(vector[1]))
		    expectedData.append(float(vector[2]))
	for i in range(len(expectedData)):
	    ratio = expectedData[i] / obtainedData[i]
	    difference = abs(1.0 - ratio)
	    if difference >= 0.01:
		passedTest = False
	#/ 2.4.- criteria
	obtainedData = []
	expectedData = []
	f = open(self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.criteria.txt', 'r')
	obtainedResults = f.readlines()
	f.close()
	f = open(self.expectedResult + 'criteria.txt', 'r')
	expectedResults = f.readlines()
	f.close()
	for datum in obtainedResults:
	    if datum[0] != '#':
		vector = datum.split()
		if vector[len(vector)-2] == '=':
		    obtainedData.append(float(vector[len(vector)-1]))
	for datum in expectedResults:
	    if datum[0] != '#':
		vector = datum.split()
		if vector[len(vector)-2] == '=':
		    expectedData.append(float(vector[len(vector)-1]))
	for i in range(len(expectedData)):
	    ratio = expectedData[i] / obtainedData[i]
	    difference = abs(1.0 - ratio)
	    if difference >= 0.01:
		passedTest = False

	return passedTest

class ClassLocalSearchOptimisationTest(ClassTest):

    '''
    Class for testing the results of local optimisations.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the local optimisations.
	'''

	f = open(self.obtainedResult, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		vector = datum.split()
		parameterValue = float(vector[0])
	ratio = parameterValue / 0.0707
	difference = abs(1.0 - ratio)
	if difference >= 0.01:
	    return False
	else:
	    return True

class ClassNumberOfSimulationsStopperTest(ClassTest):
    
    '''
    Class for testing the results of an optimisation constrained by numberOfSimulations variable.
    '''
    
    def analyseResults(self):

	'''
	This method analyses the result of an optimisation constrained by numberOfSimulations variable.
	'''
	
	f = open(self.obtainedResult, 'r')
	data = f.readlines()
	f.close()
	count = 0
	for datum in data:
	    if datum[0] != '#':
		count = count + 1
	if count == 10:
	    return True
	else:
	    return False

class ClassOEDTest(ClassTest):
    
    '''
    Class for testing the results of the optimal experimental design.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the optimal experimental design.
	'''
	
	passedTest = True
	criteria = ['MA', 'D', 'E', 'ME']
	#/ 1.- ps figures
	for criterium in criteria:
	    file = self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.oed.' + criterium + '.ps'
	    if os.path.exists(file) == False:
		passedTest = False
	    else:
		size = os.stat(file)[6]
		if size == 0:
		    passedTest = False
	#/ 2.- text files
	#/ 2.1.- summary file
	obtainedData = []
	expectedData = []
	obtainedResult = self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.oed.txt'
	expectedResult = self.expectedResult + 'oedSummary.txt'
	f = open(obtainedResult, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		obtainedData.append(float(datum.split()[2]))
	f = open(expectedResult, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		expectedData.append(float(datum.split()[2]))
	for i in range(len(obtainedData)):
	    ratio = obtainedData[i] / expectedData[i]
	    difference = abs(1.0 - ratio)
	    if difference >= 0.01:
		passedTest = False
	#/ 2.2.- trajectories files
	obtainedDir = self.obtainedResult
	expectedDir = self.expectedResult
	for criterium in criteria:
	    self.obtainedResult = obtainedDir + 'scratch/Kinetic_modelling_of_Amadori_degradation.' + criterium + '.txt'
	    self.expectedResult = expectedDir + criterium + '.txt'
	    result = self.dataFilesChecker()
	    if result == False:
		passedTest = False	
	
	return passedTest

class ClassOptionalOutputFormatTest(ClassTest):
    
    '''
    Class for testing the results of the optionalOutputFormat variable.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of a simulation with the output formatted in comma separated value.
	'''
	
	passedTest = True
	f = open(self.obtainedResult, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		vector = datum.split(',')
		if len(vector) != 15:
		    passedTest = False
	
	return passedTest

class ClassPlotKeysTest(ClassTest):

    '''
    Class for testing the results of the plotKeys variable.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of a simulation without plot keys on the gnuplot figure.
	'''
	
	passedTest = False
	f = open(self.obtainedResult, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		vector = datum.split()
		if vector[0] == 'set':
		    if vector[1] == 'nokey':
			passedTest = True
	
	return passedTest

class ClassRandomSearchOptimisationTest(ClassTest):
    
    '''
    Class for testing the results of the optimisations.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the optimisation.
	The best optimisation parameters should be stored on the xml file.
	'''
	
	def parameterValueObtainer(data):
	    
	    '''
	    This method retrieves the parameter value of the best fitness function result.
	    '''

	    fitnessValue = None
	    parameterValue = None
	    for datum in data:
		if datum[0] != '#':
		    if fitnessValue == None:
			fitnessValue = float(datum.split()[3])
			parameterValue = float(datum.split()[0])
		    else:
			if float(datum.split()[3]) < fitnessValue:
			    fitnessValue = float(datum.split()[3])
			    parameterValue = float(datum.split()[0])
	    
	    return parameterValue

	#/ 0.- starting the analyseResults method
	k10Value = None
	k2Value = None
	fitnessValue = None
	testResult = True
	#/ 1.- obtain parameter values
	#/ 1.1.- k10
	file = self.obtainedResult + 'parametersFromRandom/k10v10'
	f = open(file, 'r')
	data = f.readlines()
	f.close()
	k10Value = parameterValueObtainer(data)
	#/ 1.2.- k2
	file = self.obtainedResult + 'parametersFromRandom/k2v2'
	f = open(file, 'r')
	data = f.readlines()
	f.close()
	k2Value = parameterValueObtainer(data)
	#/ 2.- reading the xml file
	file = self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.xml'
	modelFile = libsbml.readSBML(file)
	sbmlModel = modelFile.getModel()
	reactions = sbmlModel.getListOfReactions()
	for reaction in reactions:
	    listOfLocalParameters = reaction.getKineticLaw().getListOfParameters()
	    for parameter in listOfLocalParameters:
		if parameter.getId() == 'k10':
		    k10 = parameter.getValue()
		if parameter.getId() == 'k2':
		    k2 = parameter.getValue()
	#/ 3.- comparing the values
	if str(k10Value) != str(k10):
	    print str(k10Value), str(k10)
	    testResult = False
	if str(k2Value) != str(k2):
	    print str(k2Value), str(k2)
	    testResult = False

	return testResult    

class ClassSciPyTest(ClassTest):

    '''
    Class for testing SciPy.
    '''
    
    def __init__(self):

	'''
	The constructor.
	'''

        self.id = 'scipyODE'
        self.runner = self.setRunnerName('scipyODE.rn')
        self.expectedResult = self.setResultName('scipyODE.out')
	self.obtainedResult = self.setOutputName('Kinetic_modelling_of_Amadori_degradation.0.out')

        return None
        
    def analyseResults(self):

	'''
	This method analyses the results of the SciPy test.
	'''

	#/ 1.- determining the output
	print 'Analising results obtained ...'
	#/ 2.- analysis of the results for the concentrations
	results = self.dataFilesChecker()
	#/ 3.- analysis of the results for the velocities
	self.obtainedResult = self.setOutputName('Kinetic_modelling_of_Amadori_degradation.veloc.0.out')
	self.expectedResult = self.setResultName('scipyVelocities.out')
	results = self.dataFilesChecker()

	return results

class ClassScoreStopperTest(ClassTest):
    
    '''
    Class for testing the results of an optimisation constrained by score stopper.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the optimisation constrained by score stopper.
	'''
	
	f = open(self.obtainedResult, 'r')
	data = f.readlines()
	f.close()
	for datum in data:
	    if datum[0] != '#':
		value = float(datum.split()[3])
	if value > 1e-4:
	    return False
	else:
	    return True

class ClassSensitivityTest(ClassTest):
    
    '''
    Class for testing the results of the sensitivity.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the sensitivity.
	'''
	
	passedTest = True
	#/ 1.- checking for the ps figures
	#/ 1.1.- time course
	file = self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.sens.global.timeCourse.ps'
	if os.path.exists(file) == False:
	    passedTest = False
	else:
	    size = os.stat(file)[6]
	    if size == 0:
		passedTest = False
	#/ 1.2.- relative
	file = self.obtainedResult + 'Kinetic_modelling_of_Amadori_degradation.sens.relative.ps'
	if os.path.exists(file) == False:
	    passedTest = False
	else:
	    size = os.stat(file)[6]
	    if size == 0:
		passedTest = False
	#/ 2.- checking for the text files
	obtainedDir = self.obtainedResult
	expectedDir = self.expectedResult
	#/ 2.1.- global
	self.obtainedResult = obtainedDir + 'Kinetic_modelling_of_Amadori_degradation.sens.global.txt'
	f = open(self.obtainedResult, 'r')
	data = f.readlines()
	f.close()
	values = []
	for datum in data:
	    if datum[0] != '#':
		values.append(float(datum.split()[1]))
	expectedValues = [0.30989602672599997, 0.089588046431099994]
	for i in range(len(values)):
	    ratio = values[i] / expectedValues[i]
	    difference = abs(1.0 - ratio)
	    if difference >= 0.01:
		passedTest = False
	#/ 2.2.- relative
	expectedFile = expectedDir + 'relative.txt'
	obtainedFile = obtainedDir + 'Kinetic_modelling_of_Amadori_degradation.sens.relative.txt'
	f = open(expectedFile, 'r')
	expectedData = f.readlines()
	f.close()
	f = open(obtainedFile, 'r')
	obtainedData = f.readlines()
	f.close()
	expectedValues = []
	obtainedValues = []
	for i in range(len(expectedData)):
	    if expectedData[i][0] != '#':
		values = expectedData[i].split()
		if values[0] == 'k2v2' or values[0] == 'k10v10':
		    for j in range(len(values)-1):
			expectedValues.append(float(values[j+1]))
	for i in range(len(obtainedData)):
	    if obtainedData[i][0] != '#':
		values = obtainedData[i].split()
		if values[0] == 'k2v2' or values[0] == 'k10v10':
		    for j in range(len(values)-1):
			obtainedValues.append(float(values[j+1]))
	for i in range(len(obtainedValues)):
	    if expectedValues[i] != 0.:
		ratio = obtainedValues[i] / expectedValues[i]
		difference = abs(1.0 - ratio)
		if difference >= 0.01:
		    passedTest = False
		
	return passedTest

class ClassSeparatedGraphsTest(ClassTest):
    
    '''
    Class for testing the results of a simulation with the option of a single file for each of the node trajectories.
    '''
    
    def analyseResults(self):

	'''
	This method analyses the results of a simulation asking for separated graphs.
	'''
	
	passedTest = True
	#/ 1.- cheking that the regular file does not exist
	file = self.obtainedResults + 'Kinetic_modelling_of_Amadori_degradation.ps'
	if os.path.exists(file) == True:
	    passedTest = False
	#/ 2.- checking for simulation graphics
	workingDir = self.obtainedResults + 'separatedGraphs/'
	nodes = ['AA', 'Cn', 'DFG', 'E1', 'E2', 'FA', 'Fru', 'Glu', 'Gly', 'MG', 'Man', 'Mel', '_1DG', '_3DG']
	for node in nodes:
	    file = workingDir + node + '.ps'
	    f = open(file, 'r')
	    data = f.readlines()
	    f.close()
	    if len(data) < 500: #/ totally arbitrary
		passedTest = False
	#/ 3.- checking for the velocity files
	workingDir = self.obtainedResults + 'separatedGraphs/velocity/'
	for node in nodes:
	    file = workingDir + node + '.ps'
	    f = open(file, 'r')
	    data = f.readlines()
	    f.close()
	    if len(data) < 500: #/ totally arbitrary
		passedTest = False

	return passedTest

class ClassSimulationMethodsTest(ClassTest):

    '''
    Class for testing the several methods of integration available from Scipy.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the simulations due to different integration methods.
	'''
	
	result = self.dataFilesChecker()
	
	return result

class ClassStochasticGeneralTest(ClassTest):
    
    '''
    Class for testing the general issues for the stochastic simulation.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the stochastic simulation in a general form.
	'''

	passedTest = False
	variables = []
	for i in range(300): #/ the number of stochastic runs
	    variables.append([])
	    #/ 1.- reading data
	    fileName = self.obtainedResult + 'isomerisation.0.' + 'stoch.' + str(i) + '.out'
	    file = open(fileName, 'r')
	    data = file.readlines()
	    file.close()
	    #/ 2.- selecting the variables
	    numberOfVariables = len(data[len(data)-1].split())
	    for j in range(numberOfVariables):
		variables[i].append([])
	    for datum in data:
		if datum[0] != '#':
		    vector = datum.split()
		    for j in range(numberOfVariables):
			variables[i][j].append(float(vector[j]))
	#/ 3.- calculating the standard deviation
	time = variables[0][0]
	means = []
	sds = []
	for j in range(len(time)):
	    values = []
	    for k in range(len(variables)):
		values.append(variables[k][1][j])
	    theMean = scipy.mean(values)
	    theStandardDeviation = scipy.std(values)
	    means.append(theMean)
	    sds.append(theStandardDeviation)   
	#/ 4.- final standard deviation
	finalSd = sds[40]
	if 4.3 <= finalSd <= 6.3:
	    passedTest = True
	
	return passedTest

class ClassStochasticLastStateTest(ClassTest):
    
    '''
    Class for testing the results of stochastic histograms.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the histograms resulting from stochastic simulations.
	'''
	
	passedTest = True
	separatedGraphsDirectory = self.obtainedResult + 'separatedGraphs/'
	if os.path.exists(separatedGraphsDirectory) == True:
	    separatedGraphs = self.obtainedResult + 'separatedGraphs/s'
	    for i in range(2):
		fileName = separatedGraphs + str(i+1) + '.ps'
		if os.path.exists(fileName) == False:
		    passedTest = False
		else:
		    size = os.stat(fileName)[6]
		    if size == 0:
			passedTest = False
	else:
	    singleFigure = self.obtainedResult + 'isomerisation.ps'
	    if os.path.exists(singleFigure) == False:
		passedTest = False
	    else:
		size = os.stat(singleFigure)[6]
		if size == 0:
		    passedTest = False	    

	return passedTest

class ClassStochasticSeparatedGraphsTest(ClassTest):

    '''
    Class for testing the results of stochastic simulations rendering separated graphs.
    '''
    
    def analyseResults(self):
	
	'''
	This method for testing the results of stochastic simulaitons rendering separated graphs.
	'''

	passedTest = True
	#/ 1.- checking that the trajectories are correct in size
	for i in range(3):
	    fileName = self.obtainedResult + 'isomerisation.0.' + 'stoch.' + str(i) + '.out'
	    file = open(fileName, 'r')
	    data = file.readlines()
	    file.close()
	    for datum in data:
		vector = datum.split()
	    if float(vector[0]) < 200: #/ last time should be 200
		passedTest = False
	#/ 2.- checking for the figures
	separatedGraphs = self.obtainedResult + 'separatedGraphs/s'
	for i in range(2):
	    fileName = separatedGraphs + str(i+1) + '.ps'
	    if os.path.exists(fileName) == False:
		passedTest = False
	    else:
		size = os.stat(fileName)[6]
		if size == 0:
		    passedTest = False
	
	return passedTest

class ClassStochasticSingleFigureTest(ClassTest):
    
    '''
    Class for testing the results of stochastic simulations creating a single figure.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the stochastic simulation creating a single figure.
	'''
	
	passedTest = True
	#/ 1.- checking that the trajectories are correct in size
	for i in range(3):
	    fileName = self.obtainedResult + 'isomerisation.0.' + 'stoch.' + str(i) + '.out'
	    file = open(fileName, 'r')
	    data = file.readlines()
	    file.close()
	    for datum in data:
		vector = datum.split()
	    if float(vector[0]) < 200: #/ last time should be 200
		passedTest = False
	#/ 2.- checking for the figure
	file = self.obtainedResult + 'isomerisation.ps'
	if os.path.exists(file) == False:
	    passedTest = False
	else:
	    size = os.stat(file)[6]
	    if size == 0:
		passedTest = False
	
	return passedTest

class ClassTagFormatTest(ClassTest):
    
    '''
    Class for testing the results of a multicellular simulation using an input model on tag format.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the simulation of a multicellular tag format model.
	'''
	
	passedTest = True
	#/ 1.- ps figures
	file = self.obtainedResult + 'otic9cell.ps'
	if os.path.exists(file) == False:
	    passedTest = False
	else:
	    size = os.stat(file)[6]
	    if size == 0:
		passedTest = False
	file = self.obtainedResult + 'otic9cell.grid.ps'
	if os.path.exists(file) == False:
	    passedTest = False
	else:
	    size = os.stat(file)[6]
	    if size == 0:
		passedTest = False
	#/ 2.- simulation text files
	self.obtainedResult = self.obtainedResult + 'otic9cell.0.out'
	result = self.dataFilesChecker()
	if result == False:
	    passedTest = False

	return passedTest
	

class ClassTrajectoriesReconstructionTest(ClassTest):
    
    '''
    Class for testing the trajectories reconstruction function.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of the trajectories reconstruction.
	'''
	
	nodes = ['AA', 'Cn', 'DFG', 'E1', 'E2', 'FA', 'Fru', 'Glu', 'Gly', 'MG', 'Man', 'Mel', '_1DG', '_3DG']
	passedTest = True
	#/ 1.- detecting pdf files
	directory = self.obtainedResult + 'reconstructedDynamics/'
	for node in nodes:
	    file = directory + node + '.0.pdf'
	    if os.path.exists(file) == False:
		passedTest = False
	#/ 2.- checking the trajectories
	expectedDir = self.expectedResult
	obtainedDir = self.obtainedResult + 'scratch/'
	for node in nodes:
	    self.expectedResult = expectedDir + node + '.0.sd'
	    self.obtainedResult = obtainedDir + node + '.0.sd'
	    result = self.dataFilesChecker()
	    if result == False:
		passedTest = False
		
	return passedTest

class ClassWithoutGraphicsTest(ClassTest):
    
    '''
    Class for testing the results of the withoutGraphics variable.
    '''
    
    def analyseResults(self):
	
	'''
	This method analyses the results of a simulation without graphics.
	'''
	
	if os.path.exists(self.obtainedResult) == False:
	    return True
	else:
	    return False
