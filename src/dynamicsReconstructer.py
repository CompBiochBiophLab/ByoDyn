#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana
#
#  Created: 2005-10-18 by Adrian L. Garcia-Lomana
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

# $Id: dynamicsReconstructer.py,v 4.14 2008/12/14 19:27:52 alglomana Exp $

## \file
# This module reconstructs the dynamics of a model given the parameter values.

import os,sys,math
import errorMessages, centralFunctions, initiator, simulatorXPP, simulatorOpenModelica, checker, simulator
try:
    import matplotlib.pyplot
except ImportError:
    raise errorMessages.ClassSimulatorStochasticException, 'error while importing matplotlib.pyplot.'


def central(metamodel, model, outputfiles):
    
    '''
    This is the central function of the module.
    It checks that the parameters introduced are correct, reads the solutions and runs the model.
    '''

    print 'Checking parameters ...'
    solutions = parametersDetector(metamodel, model)
    print 'Checking for initial conditions ...'
    initialConditions = initialConditionsDetector(metamodel, model)
    print 'Reconstructing trajectories ...'
    runner(metamodel, model, outputfiles, solutions, initialConditions)
    print 'Plotting ...'
    plotter(solutions,model,outputfiles,metamodel)
    
    return None

def initialConditionsDetector(metamodel, model):
    
    '''
    This function checks for the solutions of the initial conditions in case of
    '''
    
    if metamodel.initialConditionsSolutions == None:
	initialConditions = None
	print 'None detected.'
    else:
	print 'Detected.'
	if os.path.exists(metamodel.initialConditionsSolutions) == False:
	    raise errorMessages.ClassDynamicsReconstructerException, "the initial conditions file does not exist."
	#/ 2.2.- reading the file
	f=open(metamodel.initialConditionsSolutions,'r')
	lines = f.readlines()
	f.close()
        #/ 2.3.- determining the initial conditions
	initialConditions = []
	for line in lines:
	    if line[0] != '#':
		values = []
		vector = line.split()
		for element in vector:
		    values.append(float(element))
		#/ 2.4.- checking the size of the initial conditions
		if len(values) != len(model.nodes):
		    raise errorMessages.ClassDynamicsReconstructerException, 'the size of the initial conditions solutions does not match the model nodes.'
		initialConditions.append(values)
    
    return initialConditions

def modelDetermination(model, simulation, solutions, initialConditions):

    '''
    This function sets the model parameter values.
    '''

    #/ 1.- determining the model parameters
    for modelParameter in model.parameters:
        for simulationParameter in solutions.keys():
            if modelParameter == simulationParameter:
                model.parameters[modelParameter] = solutions[simulationParameter][simulation]
    #/ 2.- determining the model initial conditions, in case of
    if initialConditions != None:
	theInitialConditions = initialConditions[simulation]
	#/ 2.1.- detecting if there is a match on initial conditions
	if len(theInitialConditions) != len(model.initialConditions):
	    raise errorMessages.ClassDynamicsReconstructerException, 'the size of the initial conditions solutions does not match the model nodes.'
	model.initialConditions = theInitialConditions
                
    return model

def parametersDetector(metamodel, model):

    '''
    This function checks that the parameters' solutions' files located at the solutionsDirectory are part of the model.
    '''

    parameters={}
    #/ 1.- checking that the directory exists
    if os.path.exists(metamodel.solutionsDirectory) == False:
	raise errorMessages.ClassDynamicsReconstructerException, "the solutions's directory does not exist."
    for parameter in model.parameters:
	plausibleParameter = metamodel.solutionsDirectory + '/' + parameter
	if os.path.exists(plausibleParameter) == True:
	    print parameter,'solutions file detected.'
	    parameters[parameter]=[]
	    f=open(plausibleParameter, 'r')
	    lines = f.readlines()
	    f.close()
	    for line in lines:
		vector = line.split()
		if vector[0] != '#':
		    parameters[parameter].append(float(vector[0]))

    return parameters

def plotter(solutions, model, outputfiles, metamodel):

    '''
    This function plots the reconstructed trajectories.
    '''
    
    for i in range(model.xlength*model.ywidth):
	for node in model.nodes:
	    maxValue = 0.0
	    minValue = 0.0
	    startPlottingTime = 0.0
	    endPlottingTime = metamodel.simulationTime
	    solutionFile=outputfiles.scratchdir+'/'+node+'.%s.sd'%i
	    figureFile=outputfiles.statisticalDynamicsDirectory+'/'+node+'.%s.pdf'%i
	    f=open(solutionFile, 'r')
	    data=f.readlines()
	    f.close()
	    values=[]
	    time=data[1].split()
	    for j in range(len(data)-2):
		concentrations=[]
		stringConcentrations=data[j+2].split()
		for element in stringConcentrations:
		    concentrations.append(float(element))
		values.append(concentrations)
		if maxValue < max(concentrations):
		    maxValue = max(concentrations)
            matplotlib.pyplot.plot(time,values[0],'k .')
            matplotlib.pyplot.hold(True)
	    for j in range(len(values)-1):
		matplotlib.pyplot.plot(time,values[j+1])
		maxLocal=max(values[j+1])
		if maxLocal > maxValue:
		    maxValue=maxLocal
		matplotlib.pyplot.hold(True)
	    #/ setting the experimental points
	    for target in metamodel.target.keys():
		for point in metamodel.target[target]: 
		    vector=point.split('/')
		    gene=vector[0]
		    index=vector[1]
		    indexes=[int(index.split(',')[0]),int(index.split(',')[1])]
		    cell=model.xlength*indexes[0]+indexes[1]
		    value=float(vector[2])
		    standardDeviation=math.sqrt(float(vector[3]))
		    if cell==i and gene==node:
			matplotlib.pyplot.plot([float(target)],[value],'kD')
			matplotlib.pyplot.errorbar([float(target)],[value],[standardDeviation],ecolor='k')
			if value+standardDeviation > maxValue:
			    maxValue=value+standardDeviation+standardDeviation*.1
			if value-standardDeviation < minValue:
			    minValue=value-standardDeviation-standardDeviation*.1
			if float(target) < .05*metamodel.simulationTime:
			    startPlottingTime = -.05*metamodel.simulationTime
			if float(target) > .95*metamodel.simulationTime:
			    endPlottingTime = metamodel.simulationTime + .05*metamodel.simulationTime
	    matplotlib.pyplot.xlabel('Time') #/ set x-axis label
            matplotlib.pyplot.ylabel('Concentration') #/ set y-axis label
            matplotlib.pyplot.title(node+'.%s'%i) #/ set plot title
	    matplotlib.pyplot.xlim(startPlottingTime,endPlottingTime)
	    matplotlib.pyplot.ylim(minValue,maxValue)
	    matplotlib.pyplot.savefig(figureFile,dpi=1200)
	    matplotlib.pyplot.hold(False)

    return None

def run(model, metamodel, outputfiles):

    '''
    This function calls the model integrators to render the dynamics.
    '''
    if metamodel.integrationOption == 'automatic':
        #/ For the automatic integrationOption, we need to check if scipy, octave, xpp and openModelical are available on the current machine.
        tester = checker.ClassChecker()
        tester.detector()
        #/ Next we adjust the integrationOption depending of model characteristics.
        metamodel.integrationOption = simulator.setIntegrationOption(model, tester)
    if metamodel.integrationOption == 'python':
        centralFunctions.pythonIntegration(metamodel, model, outputfiles)
    elif metamodel.integrationOption == 'octave':
        centralFunctions.octaveIntegration(metamodel, model, outputfiles)
    elif metamodel.integrationOption == 'xpp':
	simuXPP = simulatorXPP.ClassSimulatorXPP()   
        simuXPP.createInput(metamodel, model, outputfiles)
        simuXPP.callSolver(outputfiles)
        simuXPP.createOutputs(model, outputfiles, metamodel)
    elif metamodel.integrationOption == 'openModelica':
        simulatorOM = simulatorOpenModelica.ClassSimulatorOpenModelica()      
	simulatorOM.createInput(metamodel, model, outputfiles)
	simulatorOM.callSolver(outputfiles)
	simulatorOM.createOutputs(model, outputfiles, metamodel)
    else:      
        raise errorMessages.ClassDynamicsReconstructerException, 'error at the runner file. Check the variable \"integrationMethod\".'

    return None

def runner(metamodel, model, outputfiles, solutions, initialConditions):

    '''
    This function, for each solution,  creates and simulates the corresponding model and stores the dynamics.
    '''
    
    for simulation in range(len(solutions[solutions.keys()[0]])):
        model = modelDetermination(model, simulation, solutions, initialConditions)
        run(model, metamodel, outputfiles)
        storeInfo(outputfiles, simulation, model)
 
    return None

def storeInfo(outputfiles, simulation, model):

    '''
    This function stores the dynamics of the solutions.
    It creates a file for each node of the system.
    '''


    #/ 1.- writting the headers and the time 
    if simulation == 0:
        times = []
        file = open(outputfiles.simulationResults, 'r')
        lines = file.readlines()
        #/ removing the commentaries
        commentaryIndexes = []
        for i in range(len(lines)):
            if lines[i][0] == '#':
                commentaryIndexes.append(i)
        for i in range(len(commentaryIndexes)):
            lines.pop(commentaryIndexes[len(commentaryIndexes) - i - 1]) # backward loop
        #/ writting the time
        for line in lines:
            times.append(line.split()[0])
        file.close()
        #/ writing a file for each of the nodes
        nodeFiles = []
        for node in model.nodes:
	    for i in range(model.xlength*model.ywidth):
		nodeFiles.append(outputfiles.scratchdir + '/' + node + '.%s.sd'%i)
        if os.path.exists(outputfiles.statisticalDynamicsDirectory) == False:
            os.mkdir(outputfiles.statisticalDynamicsDirectory)
        for nodeFile in nodeFiles:
            file = open(nodeFile, 'w')
            file.write('# generated by ByoDyn version %s\n'%initiator.BYODYNVERSION)
            for time in times:
                file.write('%s\t'%time)
            file.write('\n')
            file.close()
    # 2.- obtaining the data
    file = open(outputfiles.simulationResults, 'r')
    lines = file.readlines()
    file.close()
    # renmoving the commentaries
    commentaryIndexes = []
    for i in range(len(lines)):
        if lines[i][0] == '#':
            commentaryIndexes.append(i)
    for i in range(len(commentaryIndexes)):
        lines.pop(commentaryIndexes[len(commentaryIndexes) - i - 1]) #/ backward loop
    #/ obtaining the correct structure
    structure = {}
    for i in range(len(lines[0].split())-1):
	structure[i]=[]
    for i in range(len(lines)):
        values = lines[i].split('\t')
        values[len(values) - 1] = values[len(values) - 1].split('\n')[0] #/ removing the last \n 
        for j in range(len(values)-1):
            structure[j].append(float(values[j+1]))
    #/ 3.- storting the data
    for i in range(model.xlength*model.ywidth):
	for j in range(len(model.nodes)):
	    nodeFile=outputfiles.scratchdir + '/' + model.nodes[j] + '.%s.sd'%i
	    file=open(nodeFile,'a')
	    for value in structure[i*len(model.nodes)+j]:
		file.write('%s\t'%value)
	    file.write('\n')
	    file.close()
        
    return None
