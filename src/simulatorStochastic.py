#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Pau Rue-Queralt
#
#  Created: 2008-06-05 by Pau Rue-Queralt
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

# $Id: simulatorStochastic.py,v 4.24 2008/12/16 21:35:06 alglomana Exp $

## \file 
# This module is contains the code for the stochastic simulation algorithms.

import os, re, sys, copy, scipy
import formulas, initiator, errorMessages, sbmlWorker
from affectors import *

stochasticdir = os.path.join(os.environ.get('BYODYN_PATH'), 'lib', 'stochastic')
sys.path.append('%s'%stochasticdir)

import gssa, tauleap

try:
    import matplotlib.pyplot
except ImportError:
    raise errorMessages.ClassSimulatorStochasticException, 'error while importing matplotlib.pyplot.'
     
class ClassSimulatorStochastic:

    '''
    Class for the Stochastic simulators.
    '''

    def __init__(self):

        '''
        The constructor.
        '''

        pass

    def __getVariables(self, model):

        '''
        This private method sets the initial conditions of the system.
        '''

        #/ 1.- dealing with initial conditions
        variables = []
        initialVariablesValues = []
        for n in model.nodes:
            if model.constantNodes.keys().count(n) == 0:
                variables.append(n)
                initialVariablesValues.append(model.initialConditions[model.nodes.index(n)])
        for n in model.algebraicNodes:
            variables.append(n)
            initialVariablesValues.append(model.algebraicNodes[n])
        for n in model.nonConstantParameters:
            variables.append(n)
            initialVariablesValues.append(model.nonConstantParameters[n])
        for n in model.nonConstantCompartments:
            variables.append(n) 
            initialVariablesValues.append(model.nonConstantCompartments[n])   
        #/ 2.- in the case of a multicellular system
        if model.xlength != 1 or model.ywidth != 1:
            initialVariablesValues = model.initialConditions
        
        return variables, initialVariablesValues

    def __loadIntegrationFunction(self, outputfiles):

        '''
        This private method loads the system of equations, that is the topology of the system.
        '''

        sys.path.append('%s' %outputfiles.scratchdir)
        import stochasticIntegrationFunction
        from stochasticIntegrationFunction import propensities

        return propensities


    def __plotSimulation(self, model, metamodel, outputfiles, run):

        '''
        This private method creates the graph of the trajectories.
        '''
	
	def exhaustivePlotter(variables, metamodel, outputfiles):
	    
	    '''
	    A method to make exhaustive plots for stochastic simulations.
	    This is the case no time step is defined.
	    '''
	    
	    numberOfPoints = 0
	    for i in range(len(variables)):
		numberOfPoints = numberOfPoints + len(variables[i])*len(variables[i][0])
	    print 'The second argument of the runner file "time&timestep" was set to 0. Therefore, all points will be plotted which involves %s points. If you want plots to be built faster, please specify an integer for the second argument of "time&timestep" and that will be the time step shown in the graphs. \nPlotting ...'%numberOfPoints
	    if metamodel.separatedGraphs == False:
		for i in range(len(variables)):
		    time = variables[i][0]
		    for j in range(len(variables[i])-1):
			matplotlib.pyplot.plot(time, variables[i][j+1], '.')
			matplotlib.pyplot.hold(True)
		matplotlib.pyplot.title(model.systemName) 
		matplotlib.pyplot.xlabel('Time') 
		matplotlib.pyplot.ylabel('Concentration') 
		startPlottingTime = -.05*metamodel.simulationTime
		endPlottingTime = metamodel.simulationTime + .05*metamodel.simulationTime
		matplotlib.pyplot.xlim(startPlottingTime, endPlottingTime)
		if metamodel.showingPlot == True:
		    matplotlib.pyplot.show()
		else:
		    matplotlib.pyplot.savefig(outputfiles.timePlot,dpi=1200)
		    matplotlib.pyplot.hold(False)
	    #/ separated graphs
	    else:
		separatedGraphsDir = outputfiles.outputdir + '/' + 'separatedGraphs'
		if os.path.exists(separatedGraphsDir) == False:
		    os.mkdir(separatedGraphsDir)
		#/ cheching the number of model.nodes is the same to the number of trajectories
		for i in range(len(variables)):
		    if len(model.nodes) != len(variables[i])-1:
			raise errorMessages.ClassSimulatorStochasticException, 'error from exhaustivePlotter. The dimensions of model nodes is not the same as the number of trajectories. Exiting at the plotting phase.'
		#/ determining the files
		for i in range(len(model.nodes)):
		    figureFile = separatedGraphsDir + '/' + model.nodes[i] + '.ps'
		    for j in range(len(variables)):
			time = variables[j][0]
			variable = variables[j][i+1]
			matplotlib.pyplot.plot(time, variable, '.')
			matplotlib.pyplot.hold(True)
		    matplotlib.pyplot.xlabel('Time') 
		    matplotlib.pyplot.ylabel('Concentration') 
		    matplotlib.pyplot.title(model.nodes[i]) 
		    startPlottingTime = -.05*metamodel.simulationTime
		    endPlottingTime = metamodel.simulationTime + .05*metamodel.simulationTime
		    matplotlib.pyplot.xlim(startPlottingTime, endPlottingTime)
		    #/ saving the figure
		    matplotlib.pyplot.savefig(figureFile,dpi=1200)
		    matplotlib.pyplot.hold(False)
	    
	    return None

	def lastStateGrapher(metamodel, outputfiles):
	    
	    '''
	    This function deals with the histograms for the last state of the variables
	    '''
	    
	    #/ 1.- reading variables
	    f = open(outputfiles.simulationResults, 'r')
	    data = f.readlines()
	    f.close()
	    variables = []
	    for datum in data:
		if datum[0] != '#':
		    numberOfVariables = len(datum.split())-1
		    break
	    if numberOfVariables != len(model.nodes):
		raise errorMessages.ClassSimulatorStochasticException, 'while building the histograms for the last state variables for the stochastic simulations, the number of nodes does not match the number of trajectories.'
	    for i in range(numberOfVariables):
		variables.append([])
	    for datum in data:
		if datum[0] != '#':
		    vector = datum.split()
		    for i in range(numberOfVariables):
			variables[i].append(float(vector[i+1]))
	    if metamodel.separatedGraphs == True:
		separatedGraphsDir = outputfiles.outputdir + '/' + 'separatedGraphs'
		if os.path.exists(separatedGraphsDir) == False:
		    os.mkdir(separatedGraphsDir)
		#/ building up the histograms
		for i in range(len(variables)):
		    figureFile = separatedGraphsDir + '/' + model.nodes[i] + '.ps'
		    if len(variables[i]) < 100:
			resolution = len(variables[i])
		    else:
			resolution = 100
		    matplotlib.pyplot.hist(scipy.array(variables[i]), resolution, normed=True)
		    matplotlib.pyplot.hold(True)
		    matplotlib.pyplot.title(model.nodes[i])
		    matplotlib.pyplot.xlabel('Number of particles')
		    matplotlib.pyplot.ylabel('p')
		    matplotlib.pyplot.savefig(figureFile,dpi=1200)
		    matplotlib.pyplot.hold(False)
	    else:
		for i in range(len(variables)):
		    if len(variables[i]) < 100:
			resolution = len(variables[i])
		    else:
			resolution = 100
		    matplotlib.pyplot.hist(scipy.array(variables[i]), resolution, normed=True)
		    matplotlib.pyplot.hold(True)
		    matplotlib.pyplot.title(model.systemName)
		    matplotlib.pyplot.xlabel('Number of particles')
		    matplotlib.pyplot.ylabel('p')
		matplotlib.pyplot.savefig(outputfiles.timePlot,dpi=1200)

	    return None
	
	def statisticalPlotter(variables, metamodel, outputfiles):
	    
	    '''
	    A method to make statistical plots for stochastic simulations.
	    This is the case when a time step is defined.
	    '''
	    
	    print 'Plotting ...'
	    if metamodel.separatedGraphs == False:
		time = variables[0][0]
		numberNodes = len(variables[0])-1
		maxValue = 0.0
		minValue = max(model.initialConditions)
		for i in range(numberNodes):
		    means = []
		    sds = []
		    for j in range(len(time)):
			values = []
			for k in range(len(variables)):
			    values.append(variables[k][i+1][j])
			theMean = scipy.mean(values)
			means.append(theMean)
			theStandardDeviation = scipy.std(values)
			sds.append(theStandardDeviation)
			if maxValue < theMean + theStandardDeviation + (theMean + theStandardDeviation)*.05:
			    maxValue = theMean + theStandardDeviation + (theMean + theStandardDeviation)*.05
			if minValue > theMean - theStandardDeviation - (theMean + theStandardDeviation)*.05:
			    minValue = theMean - theStandardDeviation - (theMean + theStandardDeviation)*.05
		    matplotlib.pyplot.errorbar(time, means, sds, fmt='.')
		matplotlib.pyplot.xlabel('Time') 
		matplotlib.pyplot.ylabel('Concentration')
		matplotlib.pyplot.title(model.systemName)
		startPlottingTime = -.05*metamodel.simulationTime
		endPlottingTime = metamodel.simulationTime + .05*metamodel.simulationTime
		matplotlib.pyplot.xlim(startPlottingTime, endPlottingTime)
		matplotlib.pyplot.ylim(minValue,maxValue)
		if metamodel.showingPlot == True:
		    matplotlib.pyplot.show()
		else:
		    matplotlib.pyplot.savefig(outputfiles.timePlot,dpi=1200)
	    #/ separated graphs
	    else:
		separatedGraphsDir = outputfiles.outputdir + '/' + 'separatedGraphs'
		if os.path.exists(separatedGraphsDir) == False:
		    os.mkdir(separatedGraphsDir)
		#/ cheching the number of model.nodes is the same to the number of trajectories
		for i in range(len(variables)):
		    if len(model.nodes) != len(variables[i])-1:
			raise errorMessages.ClassSimulatorStochasticException, 'error from statisticalPlotter. The dimensions of model nodes is not the same as the number of trajectories. Exiting at the plotting phase.'
		#/ determining the files
		time = variables[0][0]
		for i in range(len(model.nodes)):
		    figureFile = separatedGraphsDir + '/' + model.nodes[i] + '.ps'
		    means = []
		    sds = []
		    maxValue = 0.0
		    minValue = model.initialConditions[i]
		    for j in range(len(time)):
			values = []
			for k in range(len(variables)):
			    values.append(variables[k][i+1][j])
			theMean = scipy.mean(values)
			theStandardDeviation = scipy.std(values)
			means.append(theMean)
			sds.append(theStandardDeviation)
			if maxValue < theMean + theStandardDeviation + (theMean + theStandardDeviation)*.05:
			    maxValue = theMean + theStandardDeviation + (theMean + theStandardDeviation)*.05
			if minValue > theMean - theStandardDeviation - (theMean + theStandardDeviation)*.05:
			    minValue = theMean - theStandardDeviation - (theMean + theStandardDeviation)*.05
		    matplotlib.pyplot.figure()
		    matplotlib.pyplot.errorbar(time, means, sds, fmt='k.')
		    matplotlib.pyplot.xlabel('Time') 
		    matplotlib.pyplot.ylabel('Concentration') 
		    matplotlib.pyplot.title(model.nodes[i]) 
		    startPlottingTime = -.05*metamodel.simulationTime
		    endPlottingTime = metamodel.simulationTime + .05*metamodel.simulationTime
		    matplotlib.pyplot.xlim(startPlottingTime, endPlottingTime)
		    matplotlib.pyplot.ylim(minValue,maxValue)
		    #/ saving the figure
		    matplotlib.pyplot.savefig(figureFile,dpi=1200)
		    
	    return None
	    
	print 'Performing the statistical analysis for plotting ...'
	if metamodel.onlyLastState == True:
	    lastStateGrapher(metamodel, outputfiles)
	else:
	    variables = []
	    for i in range(metamodel.stochasticRuns):
		variables.append([])
		#/ 1.- reading data
		fileName = outputfiles.simulationResults[0:-3] + 'stoch.' + str(i) + '.out'
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
	    #/ 4.- actually plotting
	    if metamodel.simulationTimeStep == 0:
		exhaustivePlotter(variables, metamodel, outputfiles)
	    else:
		statisticalPlotter(variables, metamodel, outputfiles)

        return None


    def __storeResults(self, x, t, model, metamodel, outputfiles, run):

        '''
        This method stores the results of the simulation in a file in the output directory.
        '''
	
        if metamodel.onlyLastState == True:
            if run == 0:
                file = open(outputfiles.simulationResults, 'w')
                file.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
                file.write('# t\t%s\n'%("\t".join(model.nodes)))
		file.write('%s\t%s\t'%(str(t[-1]), "\t".join(map(str,x[-1,:]))))
		file.write('\n')
            else:
                file = open(outputfiles.simulationResults, 'a')
                file.write('%s\t%s\t'%(str(t[-1]), "\t".join(map(str,x[-1,:]))))
                file.write('\n')
            file.close()
        else:
	    fileName = outputfiles.simulationResults[0:-3] + 'stoch.' + str(run) + '.out'
            file = open(fileName, 'w')
            file.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
            var = []
            for j in range(len(x)):
                for i in range(len(x[j])):
                    if j == 0:
                        var.append([x[j][i]])
                    else:
                        var[i].append(x[j][i])
            for i in range(len(t)):
                file.write('%s\t'%t[i])
                for variable in var:
                    file.write('%s\t'%variable[i])
                file.write('\n')
            file.close()

        return None

    def __writePropensitiesInput(self, model, outputfiles, initialVariablesValues, variables):

        '''
        This private method writes a file with the topology of the system.
        '''

    #/ 1.- Writting some headers
	file = open(outputfiles.integrationInput, 'w')
        file.write('#/\n#/ generated by ByoDyn version %s\n#/\nimport numpy as np\n'%initiator.BYODYNVERSION)
    #/ 2.- Writting sbml defined functions
        for function in model.functions:
            file.write('def %s(' %(function.id))
            file.write('%s' %function.arguments[0])
            for argument in function.arguments[1:]:
                file.write(', %s' %argument)
            file.write('):\n')    
            squareRootDefinitions = re.findall('root\(2, [\w\[\]()/\+\-\*\s\d\.\^\,]*\)', function.output)
            if len(squareRootDefinitions) != 0:
                function.output = function.output.replace('root(2, ', 'sqrt(')
            file.write('\treturn %s\n' %function.output)
    #/ 3.- Starting the function
        file.write('def propensities(x):\n')
    #/ 4.- Defining required variables
        option = 'python'
        cellIndex = 0 
        equationNumber = 0
    #/ 5.- Parameters, compartments and constant nodes
        for parameter in model.parameters.keys():
            file.write('\t%s = %s\n'%(parameter, model.parameters[parameter]))
        for compartment in model.compartments.keys():
            file.write('\t%s = %s\n'%(compartment, model.compartments[compartment]['size']))
        for species in model.constantNodes.keys():
            file.write('\t%s = %s\n'%(species, model.constantNodes[specie]))
    #/ 6.- Defining the topology
        file.write("\t[" + ",".join(model.nodes) + "] = x\n")
        propensities = [reaction.propensity for reaction in model.reactions]
        file.write("\ta = np.array([" + ",".join(propensities) + "])\n")
        file.write('\treturn a\n')
        file.close()

        return self

    def run(self, metamodel, model, outputfiles):

        '''
        This method directs the simulation.
        It creates the initial conditions, 
        loads the system of equations, 
        simulates the system and
        plots the results.
        '''

        variables, initialVariablesValues = self.__getVariables(model)
	outputfiles.integrationInput = outputfiles.scratchdir + '/stochasticIntegrationFunction.py'
        self.__writePropensitiesInput(model, outputfiles, initialVariablesValues, variables)
        propensities = self.__loadIntegrationFunction(outputfiles)
        stoichiometry = model.stoichiometryMatrix
        time =  metamodel.simulationTime
        timestep = metamodel.simulationTimeStep
        for run in range(metamodel.stochasticRuns):
            x, t = self.simulate(propensities, stoichiometry, initialVariablesValues, time, timestep, metamodel.stochasticOption, metamodel.stochasticMethod)
	    self.__storeResults(x, t, model, metamodel, outputfiles, run)
	if metamodel.graphics == True:
	    self.__plotSimulation(model, metamodel, outputfiles, run)


    def simulate(self, propensities, stoichiometryMatrix, x_0, time, dt, option, method = None):

        '''
        This method simulate the system of equations.
        method     describes which method to use
        pyfile     python source file where the butcher tableau is found
        ''' 

        if option == 'ssa' or option == 'default':
            x, t = gssa.simulate (propensities, stoichiometryMatrix, x_0, time, dt)
        elif option == 'tau-leap':
            x, t = tauleap.simulate (propensities, stoichiometryMatrix, x_0, time, dt)
        else:
            #/ We should raise an error here instead
            print "Warning, unknown option: %s, defaulting to SSA"%(option)
            x, t = stoch.GSSA (propensities, stoichiometryMatrix, x_0, time, dt)

        return x, t
