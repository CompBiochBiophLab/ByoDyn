#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido and Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana, Alex Gomez-Garrido and David Sportouch
#
#  Created: 2006-02-10 by Adrian L. Garcia-Lomana
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

# $Id: sensitivityAnalyzer.py,v 4.7 2008/12/03 16:37:38 alglomana Exp $

## \file
# This module contains the algorithms necessary for the analysis of sensitivity of a given model.

import copy, math, sys, os, scipy, scipy.linalg
import errorMessages, optimiser, simulator, identifiabilityAnalyzer, matrixWorker, initiator

class ClassSensitivityAnalyzer:
    
    '''
    Class for the sensitivity analysis.
    '''

    def __init__(self, parameters):
        
	'''
	This is the constructor.
	'''

        self.parametersToStudy = parameters
        self.RSMatrix = matrixWorker.ClassMatrix()
        self.RSCoefficients = {}
        self.OSCoefficients = {}
        self.OS = {}
        self.identifiabilityCoefficients = {}
        self.restrictedTimePoints = []
        
	return None

    def __obtainSensitivity(self, model, metamodel, originalSimulation, parameterSimulation, h, parameterValue):

	'''
	This method calculates the sensitivity given the dynamics of both the system and an infinitesimally changed parameter system.
	It returns trajectories.
	'''
        
        originalNodes = copy.deepcopy(model.nodes)
        nodes = model.nodes
        for constant in model.constantNodes:
            if nodes.count(constant) != 0:
                nodes.remove(constant)
        for node in model.algebraicNodes:
            nodes.append(node)
	#/ 1.- initialising the variables
        relativeSensitivity = {}
        relativeSensitivityVector = []
	#/ 1.1.- initialising variables for the identifiability
        identifiabilitySensitivity = []
        identifiabilityVector = []
        #/ 2.- working for multicell
        if metamodel.target == {}:
            choosenNodes = range(len(nodes)*model.xlength*model.ywidth)
            targets = []
            for i in choosenNodes:
                if self.restrictedTimePoints == []:
                    targets.append(range(len(originalSimulation)))
                else:
                    targets.append([])
                    for t in self.restrictedTimePoints:
                        targets[i].append(int(float(t)/float(metamodel.simulationTimeStep)))
        else:
            choosenNodes = []
            targets = []
            for tp in metamodel.target.keys():
                for field in metamodel.target[tp]:
                    position = getPositionNode(field, model)
                    if choosenNodes.count(position) == 0:
                        choosenNodes.append(position)
                        targets.append([])
                        targets[choosenNodes.index(position)].append(int(float(tp)/float(metamodel.simulationTimeStep)))
                    else:
                        targets[choosenNodes.index(position)].append(int(float(tp)/float(metamodel.simulationTimeStep)))
		#/ Rescaling the target values when the simulation has an incorrect number of time steps.
        if len(originalSimulation) < int(float(metamodel.simulationTime)/float(metamodel.simulationTimeStep)):
            for i in range(len(targets)):
                for j in range(len(targets[i])):
                    targets[i][j] = targets[i][j]*len(originalSimulation)/int(float(metamodel.simulationTime)/float(metamodel.simulationTimeStep))
        #/ 3.- determining the formula
        identifiabilitySensitivity = []	
        for j in choosenNodes:
            relativeSensitivity[getNodeName(j, model)] = []
            relativeValue = 0
            for i in targets[choosenNodes.index(j)]:
		#/ 3.1.- calculating the derivative
                coefficient = math.fabs((float(originalSimulation[i][j]) - float(parameterSimulation[i][j]))/h)
		#/ 3.3.- acting in case of zero value of the original simulation
                if float(originalSimulation[i][j]) != 0:
		    normalisedSensitivity = coefficient * (parameterValue/float(originalSimulation[i][j]))
                    relativeSensitivity[getNodeName(j, model)].append(normalisedSensitivity)
                    relativeValue = relativeValue + normalisedSensitivity
		    identifiabilitySensitivity.append(normalisedSensitivity)
                else:
                    relativeSensitivity[getNodeName(j, model)].append(coefficient)
                    relativeValue = relativeValue + coefficient
		    identifiabilitySensitivity.append(coefficient)
            relativeSensitivityVector.append(relativeValue/len(originalSimulation)) # check that for the sensitivity, all time values are taken and not just hidden targets.
        model.nodes = copy.deepcopy(originalNodes)
	
	return relativeSensitivity, relativeSensitivityVector, identifiabilitySensitivity

    def calculateOverallSens(self):

	'''
	This method calculates single value sensitivities from the trajectories.
	'''
        
        for parameter in self.RSCoefficients.keys():
            self.OSCoefficients[parameter] = []
            self.OS[parameter] = 0
            for node in self.RSCoefficients[parameter].keys():
                for i in range(len(self.RSCoefficients[parameter][node])):
                    if self.RSCoefficients[parameter].keys().index(node) == 0: #/ the first one
                        self.OSCoefficients[parameter].append(math.sqrt(self.RSCoefficients[parameter][node][i]**2))
                    elif self.RSCoefficients[parameter].keys().index(node) == (len(self.RSCoefficients[parameter].keys())-1): #/ the last one
                        self.OSCoefficients[parameter][i] += math.sqrt(self.RSCoefficients[parameter][node][i]**2)
                        self.OSCoefficients[parameter][i] = self.OSCoefficients[parameter][i] / len(self.RSCoefficients[parameter].keys())
                        self.OS[parameter] += self.OSCoefficients[parameter][i]  
                    else: #/ all others, similar to the first one
                        self.OSCoefficients[parameter][i] += math.sqrt(self.RSCoefficients[parameter][node][i]**2)
            self.OS[parameter] = self.OS[parameter] / len(self.RSCoefficients[parameter][self.RSCoefficients[parameter].keys()[0]])
        
        return None
    
    def calculateSens(self, model, metamodel, outputfiles):

	'''
	This method runs the simulation for a given parameter values, 
	change the parameters infinitesimally and re-runs the simulation.
	With this data, the sensitivity is calculated.
	'''
        
        originalNodes = copy.deepcopy(model.nodes)
        #/ 1.- running the model simulation with original parameters and obtaining the values
        originalSimulationValues = simulator.obtainSimulationValues(metamodel, model, outputfiles, 'simulation')
        model.nodes = copy.deepcopy(originalNodes)
        #/ 2.1- varying the parameters
        epsilon = 1e-6
        originalParameters = copy.deepcopy(model.parameters)
        originalInitialConditions = copy.deepcopy(model.initialConditions)
        for parameter in self.parametersToStudy.keys():
            if model.parameters[parameter] == 0:
                h = epsilon
            else:
                h = epsilon * model.parameters[parameter]        
            model.parameters[parameter] = model.parameters[parameter] + h
            #/ If there are algebraic or assignment rules, we will need adjust the initial conditions
            for rule in model.rules:
                if rule.type == 'Assignment':
                    model = model.checkAssignmentRule(rule)
        #/ 2.2 running the simulation for each parameter and obtaining the values
            type = 'sensitivity.' + parameter
            parameterSimulationValues = simulator.obtainSimulationValues(metamodel, model, outputfiles, type)
		#/ 2.3 calculating the sensitivy for each parameter
            #/ Return the original parameters values
            model.parameters = copy.deepcopy(originalParameters)
            model.initialConditions = copy.deepcopy(originalInitialConditions)
            #/ calculate sensitivity
	    self.RSCoefficients[parameter], RSVector, self.identifiabilityCoefficients[parameter] = self.__obtainSensitivity(model, metamodel, originalSimulationValues, parameterSimulationValues, h, model.parameters[parameter])
            self.RSMatrix.m.append(RSVector)
        originalNodes = copy.deepcopy(model.nodes)
        nodes = model.nodes
        for constant in model.constantNodes:
            if nodes.count(constant) != 0:
                nodes.remove(constant)
        for node in model.algebraicNodes:
            nodes.append(node)        
        self.RSMatrix.colNames = nodes
        model.nodes = copy.deepcopy(originalNodes)
        self.RSMatrix.rowNames = self.parametersToStudy.keys()	
        
        return None

    def plotOSTimeCourse(self, metamodel, outputfiles):

        '''
        This method creates the sensitivity time course plot
        '''

        print 'The sensitivity has been calculated. Plotting ...'
        commandFile = outputfiles.sensitivitiesOSPlotCommands
        dataFile = outputfiles.sensitivitiesOSData
        plotFile = outputfiles.sensitivitiesOSPlotFile
	#/ 1.- creating the data file
        f = open(dataFile, 'w')
        #/ 1.1.- writing the headers
        f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
        f.write('# t')
        for parameter in self.OSCoefficients.keys():
            f.write('\t%s' %parameter)
        f.write('\n')
        time = 0.0
        for point in range(len(self.OSCoefficients[self.OSCoefficients.keys()[0]])):
            f.write('%s' %time)
            for parameter in self.OSCoefficients.keys():
                f.write('\t%s' %self.OSCoefficients[parameter][point])
            f.write('\n')
            time += metamodel.simulationTimeStep
        f.close()
        #/ 2. creating the file with the gnuplot commands.
        plot = open(commandFile, 'w')
        plot.write('#\n# generated by ByoDyn version %s\n#\n#\n# this is the input file for gnuplot\n#\n'%initiator.BYODYNVERSION)
        plot.write('set xlabel "time"\nset ylabel "Global sensitivity"\nset title "Globlal Sensitivity Time Course"\n')
        if metamodel.figureFormat == 'ps':
	    plot.write('set output "%s"\nset terminal postscript color\n' %(plotFile))
	else:
	    plot.write('set output "%s"\nset terminal png\n' %(plotFile))
        plot.write('plot \'%s\' using 1:2 title \'%s\' with lines' %(dataFile, self.OSCoefficients.keys()[0]))
        for i in range (3, len(self.OSCoefficients.keys()) + 2):
            plot.write(', \'%s\' using 1:%s title \'%s\' with lines' %(dataFile, i, self.OSCoefficients.keys()[i-2]))
        plot.close()
        #/ 3. invoking gnuplot
        os.system('gnuplot \"%s\"' %(commandFile))
        
        return None
        
    def writeOSTable(self, outputfiles):
        
        '''
        This method writes the global sensitivity values in output file
        '''

        f = open(outputfiles.sensitivitiesOS, 'w')
        f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
        f.write('# Parameter\tGlobal Sensitivity\n#\n')
        values = self.OS.values()
        values.sort(reverse=True)
        for value in values:
            for parameter in self.OS.keys():
                if self.OS[parameter] == value:
                    f.write('%s\t%s\n' %(parameter, value))
        f.close()

        return None

def central(model, metamodel, outputfiles):

    '''
    This is the central function of the module.
    It directs the flow of the program for the sensitivity analysis and the identifiability analysis.
    '''
	
    criterionValue = None  
    #/ 1.- getting parameters to study from metamodel
    parametersToStudy = componentsChecker(model, metamodel)
    #/ 2.- sensitivity analysis
    #/ 2.1.- calculation
    sensitivity = ClassSensitivityAnalyzer(parametersToStudy)
    sensitivity.calculateSens(model, metamodel, outputfiles)
    sensitivity.calculateOverallSens()   
    #/ 2.2.- plotting
    color = 'blackAndWhite'
    sensitivity.plotOSTimeCourse(metamodel, outputfiles)
    sensitivity.writeOSTable(outputfiles)
    message = 'Sensitivity values for each node with respect to a given parameter.'
    sensitivity.RSMatrix.fileWriter(message, outputfiles.sensitivitiesRSMatrix)
    sensitivity.RSMatrix.plotWriter(metamodel, outputfiles.sensitivitiesRSMatrixPlot, color)
    
    return criterionValue

def componentsChecker(model, metamodel):

    '''
    This function is an initial checking of the parameters we want to study:
    first that they exist on the model and second that they are constant during the simulation.
    Finally we select the parameters to study on a new working variable.
    At the very end we set the new parameters values set by "parameter" variable.
    '''

    #/ 1.- checking that parameters exist and are constant
    parametersToStudy = {}
    existing = False
    if metamodel.sensitivityParameters == []:
        metamodel.sensitivityParameters = model.parameters.keys()
    else:
        for sensitivityParameter in metamodel.sensitivityParameters:
            for modelParameter in model.parameters.keys():
                if sensitivityParameter == modelParameter:
                    existing = True
            for nonConstantParameter in model.nonConstantParameters.keys():
                if sensitivityParameter == nonConstantParameter:
                    existing = True
                    print 'WARNING: The parameter %s for which you want to study the sensitivity variates during the simulation. We can not study the sensitivity for this parameter.' %sensitivityParameter
                    metamodel.sensitivityParameters.remove(sensitivityParameter)
            if existing == False:
                raise errorMessages.ClassSensitivityAnalyzerException,  'the parameter %s for which you want to study the sensitivity does not exist at the model file.'%sensitivityParameter
            else:
                existing = False
    for parameter in metamodel.sensitivityParameters:
        parametersToStudy[parameter] = model.parameters[parameter]
    #/ 2.- determining the new value of the parameters
    if len(metamodel.parameters) != 0:
 	for parameter in metamodel.parameters.keys():
 	    if model.parameters.has_key(parameter) == True:
 	        model.parameters[parameter] = metamodel.parameters[parameter]
 	    else:
 		raise errorMessages.ClassSimulatorException, 'error at the runner file. Check the variable \"parameters\", because the parameter %s is unknown.' %parameter

    return parametersToStudy

def getNodeName(value, model):

    '''
    Given the position of a node in model.nodes list
    this function returns its name. 
    In the case of multicellular models the index positions
    are added in the format: "nodeName(yindex,xindex)".
    '''

    if model.xlength == 0 and model.ywidht == 0:
        return model.nodes[value]
    else:
        i = 0
        while (value >= len(model.nodes)):
            value -= len(model.nodes)
            i += 1
        xindex = 0
        yindex = 0
        for n in range(i):
            if xindex == model.xlength - 1:
                yindex += 1
                xindex = 0
            else:
                xindex += 1

        name = model.nodes[value] + '(' + str(yindex) + ',' + str(xindex) + ')'
        return name

def getPositionNode(field, model):

    '''
    Given a target field ("node/xindex,yindex/experimentalValue/variance"),
    this function returns the column position of this node 
    in the simulation output file.
    '''

    f = field.split('/')

    return model.nodes.index(f[0]) + (int(f[1].split(',')[1])*len(model.nodes)) + int(f[1].split(',')[0])*len(model.nodes)*model.xlength
