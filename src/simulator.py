#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana and Alex Gomez-Garrido
#
#  Created: 2005-07-11 by Adrian L. Garcia-Lomana
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

# $Id: simulator.py,v 4.28 2008/12/16 13:33:33 alglomana Exp $

## \file 
# This module is responsible of the simulation of the model.

import sys, os, copy, string, libsbml
import checker, errorMessages, simulatorOpenModelica, simulatorXPP, simulatorRungeKutta, simulatorEuler, simulatorStochastic, initiator, parallel
from affectors import *

class ClassEpithelium:

    '''
    Class for the multicellular plots.
    '''

    def __init__(self):
	
	'''
	This is the constructer.
	'''

        self.color = []
        self.concentration = []
        self.length = 0
        self.nodes = []
        self.time = 0.0
        self.width = 0

        return None

    def dataObtainer(grid, model, outputfiles):
	
	'''
	This method initialises the object of the class ClassEpithelium.
	It takes the values from different sources, either the object "model" or from the simulationResults file.
	'''

        grid.nodes = model.nodes
        grid.length = model.xlength
        grid.width = model.ywidth
        for i in range(len(grid.nodes)):
            grid.color.append([])
        f = open(outputfiles.simulationResults, 'r')
        for line in f:
            fields = line.split()
        grid.time = float(fields[0])
        grid.concentration = fields[1:]
        for i in range(len(grid.concentration)):
            grid.concentration[i] = float(grid.concentration[i])
        f.close()

        return grid

def central(metamodel, model, outputfiles):

    '''
    This function controls the flow of the simulation depending on the simulation options.
    It checks the correctness of the input data, 
    it creates a pdf format file with the system of equations, 
    it integrates the system with the selected integrator and
    it finally builds the graphs of the simulation.
    '''
        
    model.summary(outputfiles) #/ 0.- Writing a short summary of the model description in a output file
    createPDFFormulae(model, metamodel, outputfiles)      
    if metamodel.integrationOption == 'automatic':
        #/ 1.- For the automatic integrationOption, we need to check if scipy, octave, xpp and openModelical are available on the current machine.
        tester = checker.ClassChecker()
        tester.detector()
        #/ 2.- Next we adjust the integrationOption depending of model characteristics.
        metamodel.integrationOption = setIntegrationOption(model, tester)   
    compatibilityChecker(metamodel, model) #/ checking some issues before the simulation
    if metamodel.integrationOption == 'python':
	if len(model.events) == 0:
	    centralFunctions.pythonIntegration(metamodel, model, outputfiles)
	else:
	    eventsDealer(metamodel, model, outputfiles)
    elif metamodel.integrationOption == 'octave':
        centralFunctions.octaveIntegration(metamodel, model, outputfiles)
    elif metamodel.integrationOption == 'openModelica':
	simulatorOM = simulatorOpenModelica.ClassSimulatorOpenModelica()	#/ Creating the simulator object
	simulatorOM.createInput(metamodel, model, outputfiles)
	simulatorOM.callSolver(outputfiles)
	simulatorOM.createOutputs(model, outputfiles, metamodel)
    elif metamodel.integrationOption == 'matlab':
	centralFunctions.matlabIntegrator(metamodel, model, outputfiles)
    elif metamodel.integrationOption == 'xpp':
        simuXPP = simulatorXPP.ClassSimulatorXPP()    #/ Creating the simulator object
        simuXPP.createInput(metamodel, model, outputfiles)
        simuXPP.callSolver(outputfiles)
        simuXPP.createOutputs(model, outputfiles, metamodel)
    elif metamodel.integrationOption == 'rungeKutta':
        simRK = simulatorRungeKutta.ClassSimulatorRungeKutta()
        simRK.run(metamodel, model, outputfiles)
        metamodel.graphics = False
    elif metamodel.integrationOption == 'euler':
        simEuler = simulatorEuler.ClassSimulatorEuler()
        simEuler.run(metamodel, model, outputfiles)
        metamodel.graphics = False
    elif metamodel.stochasticOption in ['ssa', 'default', 'tau-leap']:
        simStochastic = simulatorStochastic.ClassSimulatorStochastic()
        simStochastic.run(metamodel, model, outputfiles)
        metamodel.graphics = False
    else:
        raise errorMessages.ClassSimulatorException, 'error at the runner file. Check the variable "integrationMethod" and "stochasticMethod".'
    #/ 4.- Checking the results obtained
    if metamodel.stochasticOption == None:
	trajectoriesFile = outputfiles.simulationResults
	checkResults(metamodel, trajectoriesFile)
    else:
	if metamodel.onlyLastState == False:
	    for i in range(metamodel.stochasticRuns):
		trajectoriesFile = outputfiles.simulationResults[0:-3] + 'stoch.' + str(i) + '.out'
		checkResults(metamodel, trajectoriesFile)
    if metamodel.graphics == True:
	if metamodel.velocities == True:
		type = 'velocities'
		createGnuplot(model, metamodel, outputfiles, type)
	type = 'simulation' 
	createGnuplot(model, metamodel, outputfiles, type)
	if model.xlength * model.ywidth > 1:
	    createEpiplot(model, metamodel, outputfiles)

    return None

def checkResults(metamodel, trajectoriesFile):

    '''
    This function checks the integration output file. 
    It returns an error if the results are missing and a warning if the results are incomplete.
    '''

    dataFile = open(trajectoriesFile, 'r')
    rawInformation = []
    for line in dataFile:
        fields = line.split()
        if fields[0] != '#':
            rawInformation.append(fields)
    dataFile.close()
    if len(rawInformation) == 0:
        raise errorMessages.ClassSimulatorException, 'There was an error integrating the system. Please, change the integration options.'        
    if metamodel.simulationTimeStep != 0:
	if len(rawInformation) < (metamodel.simulationTime/metamodel.simulationTimeStep):
	    print 'WARNING: The integration seems to be truncated. Therefore, the results are incomplete.'

    return None

def compatibilityChecker(metamodel, model):

    '''
    This function checks the compatibility of some simulation options:
    if the model holds events has to be integrated by OpenModelica and
    that a certain parameter of the model exists for which its value is specified.
    '''

    if metamodel.integrationOption != 'xpp' and model.delayFunctions == True:
        raise errorMessages.ClassSimulatorException, 'because the model has delay functions, the integration option for models with delays is only xpp. Please substitute the integation option string by "xpp" at the runner file. Sorry for the inconvenience.'
    #/ 1.- checking that we use openModelica integration routines if the model has Events
    if metamodel.integrationOption != 'openModelica' and len(model.events) != 0:
	raise errorMessages.ClassSimulatorException, 'because the model holds "Events", the integration option for models with events is openModelica. Please substitute the integration option string by "openModelica" at the runner file. Sorry for the inconvenience.'
    #/ 2.- putting the parameters value specified in runner file inside the model
    if len(metamodel.parameters) != 0:
 	for parameter in metamodel.parameters.keys():
 	    if model.parameters.has_key(parameter) == True:
 	        model.parameters[parameter] = metamodel.parameters[parameter]
 	    else:
 		raise errorMessages.ClassSimulatorException, 'error at the runner file. Check the variable \"parameters\", because the parameter %s is unknown.' %parameter

    return None

def createEpiplot(model, metamodel, outputfiles):

    '''
    This function creates the multicellular plots.
    It calls the class ClassEpithelium, fills the object
    and creates a flat file data and the graph.
    '''

    grid = ClassEpithelium()
    grid = ClassEpithelium.dataObtainer(grid, model, outputfiles)
    flatFileWriter(outputfiles, grid)
    plotMaker(outputfiles, grid, model)
    #/ converting to png the grid figure if it is the case
    if metamodel.figureFormat == 'png':
	psFile = '%s/%s.grid.ps' %(outputfiles.outputdir, model.systemName)
	pngFile = '%s/%s.grid.png' %(outputfiles.outputdir, model.systemName)
	cmd = 'convert %s %s'%(psFile, pngFile)
	os.system(cmd)
	
    return None

def createGnuplot(model, metamodel, outputfiles, type):

    '''
    This function creates the postscript plots for the dynamic trajectories.
    It is responsible for the plotting of both the concentration versus time plots
    and the change of concentration versus time plots.
    '''

    #/ 0.- determining the files
    separatedGraphsDir = outputfiles.outputdir + '/' + 'separatedGraphs'
    if metamodel.separatedGraphs == True:
	if os.path.exists(separatedGraphsDir) == False:
	    os.mkdir(separatedGraphsDir)
    if type == 'velocities':
        commandFile = outputfiles.gnuplotVelocitiesInputFile
        dataFile = outputfiles.velocitiesFile
        plotFile = outputfiles.velocitiesPlot
	yLabel = 'Velocities'
	separatedGraphsDir = separatedGraphsDir + '/' + 'velocity'
	if metamodel.separatedGraphs == True:
	    if os.path.exists(separatedGraphsDir) == False:
		os.mkdir(separatedGraphsDir)
    elif type == 'simulation':
        commandFile = outputfiles.gnuplotInputFile
        dataFile = outputfiles.simulationResults
        plotFile = outputfiles.timePlot
	yLabel = 'Concentration'
    else:
        if len(type.split('.')) != 1:   
            dataFile = '%s/%s.sim.%s.out' %(outputfiles.outputdir, model.systemName, type.split('.')[1])
            plotFile = '%s/%s.sim.%s.ps' %(outputfiles.outputdir, model.systemName, type.split('.')[1])
            commandFile = '%s/scratch/%s.sim.%s.gnu' %(outputfiles.outputdir, model.systemName, type.split('.')[1])
        else:
            raise errorMessages.ClassSimulatorException, 'internal problem while creating the graphs.'
    variables = []
    for n in model.nodes:
        if model.constantNodes.keys().count(n) == 0:
            variables.append(n)
    for rule in model.rules:
        if rule.type == 'Rate' and variables.count(rule.variable) == 0:
            variables.append(rule.variable)           
    for n in model.algebraicNodes:
        if metamodel.integrationOption != 'xpp':
            variables.append(n)
    #/ 1.- creating the file with the gnuplot commands.
    plot = open(commandFile, 'w')
    plot.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
    plot.write('set xlabel "time"\nset ylabel "%s"\nset title "%s"\n'%(yLabel, model.systemName))
    if metamodel.figureFormat == 'ps':
	plot.write('set terminal postscript color\n')
    if metamodel.figureFormat == 'png':
	plot.write('set terminal png\n')
    if metamodel.separatedGraphs == False:
	plot.write('set output "%s"\n' %plotFile)
    plot.write('set pointsize 0.3\n')
    if metamodel.plottingKeys == False:
        plot.write('set nokey\n')
    numberOfTrajectories = model.xlength * model.ywidth * len(variables)
    if metamodel.separatedGraphs == True:
	plotFile = separatedGraphsDir + '/' + variables[0] + '.ps'
	plot.write('set output "%s"\n'%plotFile)
    plot.write('plot \'%s\' using 1:2 title \'%s\' with points' %(dataFile, variables[0]))
    for i in range (3, numberOfTrajectories + 2):
	if metamodel.separatedGraphs == False:
	    plot.write(', \'%s\' using 1:%s title \'%s\' with points' %(dataFile, i, variables[(i - 2) % len(variables)])) 
	else:
	    plotFile = separatedGraphsDir + '/' + variables[(i - 2) % len(variables)] + '.ps'
	    plot.write('\nset output "%s"'%plotFile)
	    plot.write('\nplot \'%s\' using 1:%s title \'%s\' with points' %(dataFile, i, variables[(i - 2) % len(variables)])) 
    plot.close()
    #/ 2.- invoking gnuplot
    os.system('gnuplot \"%s\"' %(commandFile))

    return None

def createPDFFormulae(model, metamodel, outputfiles):

    '''
    This function creates the latex format files with the system of equations. 
    It calls the different affector functions at the affectors module of the lib directory.
    '''

    option = 'latex'
    cellIndex = 0  
    #/ Creating latex file
    latex = open(outputfiles.latexFile, 'w')
    latex.write('%%\n%% generated by ByoDyn version %s \n%%\n%%\n%% This is a latex file with the equations of the model. Use "pdflatex" to build the .pdf file\n%%\n'%initiator.BYODYNVERSION) #/ the second % (%%) is for the python to distinguish between %% and %s. At the latex file, only one '%' is written
    latex.write('\\documentclass[leqno]{article}\n\\usepackage{lscape}\n\\usepackage{hyperref}\n')
    latex.write('\\oddsidemargin 0in\n\\evensidemargin 0in\n\\begin{document}\n\pagestyle{empty}\n')
    latex.write('\\tableofcontents\n')
    #/ Adding function definitions
    if len(model.functions) != 0:
        latex.write('\\section{Functions}\n')
        latex.write('\\tiny\n\\begin{eqnarray}\n')
        for f in model.functions:
            idLatex = f.id.replace('_', '\_')
            arguments = ''
            for arg in f.arguments:
                arguments += arg.replace('_', '\_')
                if arg != f.arguments[len(f.arguments) - 1]:
                    arguments += ', '            
            latex.write('%s(%s)' %(idLatex, arguments))
            ASTNode = libsbml.readMathMLFromString(f.mathAST)
            mathString = formulas.getMathExpression(ASTNode)[0]
            mathString = mathString.split('lambda')[1].replace('%s,' %arguments, '')    
            latex.write(' & = & %s\\nonumber \\\\\\nonumber & & \\\\\n' %mathString) 
        latex.write('\\nonumber\n\\end{eqnarray}\n')        
    #/ Adding Events
    if len(model.events) != 0:
        latex.write('\\section{Events}\n')
        latex.write('\\tiny\n\\begin{eqnarray}')
        for event in model.events:
            triggerLatex = event.trigger.replace('_', '\_')
            latex.write('when & %s & then \\nonumber \\\\\\nonumber & & \\\\' %(triggerLatex))
            if event.delay != None:
                latex.write('with\ delay & = & %s \\nonumber \\\\\\nonumber & & \\\\' %(event.delay))
            for a in event.assignment:
                aLatex = a.replace('_', '\_')
                ASTNode = libsbml.readMathMLFromString(event.assignmentAST[a])
                latex.write('\t%s & = & %s \\nonumber \\\\\\nonumber & & \\\\' %(aLatex, formulas.getMathExpression(ASTNode)[0]))
            latex.write('end\  when \\nonumber \\\\\\nonumber & & \\\\')
        latex.write('\\nonumber\n\\end{eqnarray}\n')
    #/ Adding Rules
    if len(model.rules) != 0:
        latex.write('\\section{Rules}\n')
        latex.write('\\tiny\n\\begin{eqnarray}')
        for rule in model.rules:
            ASTNode = libsbml.readMathMLFromString(rule.mathAST)
            if rule.id.startswith('delay') == False:
                if rule.type == 'Assignment':
                    variableLatex = rule.variable.replace('_', '\_')
                    latex.write('\n%s & = & %s \\nonumber \\\\\\nonumber & & \\\\' %(variableLatex, formulas.getMathExpression(ASTNode)[0]))
                if rule.type == 'Rate':
                    variableLatex = rule.variable.replace('_', '\_')
                    latex.write('\n\\frac{\mathrm{d}[%s]_{i}}{\mathrm{d}t} & = & %s \\nonumber \\\\\\nonumber & & \\\\' %(variableLatex, formulas.getMathExpression(ASTNode)[0]))                
                if rule.type == 'Algebraic':
                    latex.write('\n0 & = & %s \\nonumber \\\\\\nonumber & & \\\\' %(formulas.getMathExpression(ASTNode)[0]))
        latex.write('\\nonumber\n\\end{eqnarray}\n')
    #/ Adding kinetic laws
    latex.write('\\section{Kinetic laws}\n')
    latex.write('\\tiny\n\\begin{eqnarray}')
    reactions = {}
    for definition in model.topology:
        if metamodel.modelFormat == 'tags':
            reactions[definition.split('/')[1]] = definition
        else:
            if definition.split('/')[2] != 'RULE':
                reactions[definition.split('/')[3].split('.')[1]] = model.topologyAST[definition]
    for r in reactions:
        rLatex = r.replace('_', '\_')
        latex.write('\n%s & = & ' %str(rLatex))
        if metamodel.modelFormat == 'tags':
            fieldsDefinition = reactions[r].split("/")
            exec "%s(model, latex, option, cellIndex, reactions[r], fieldsDefinition)"%r
        else:
            reactions[r][0] = ''
            formulas.formulaLatex(reactions[r], latex, model)
        latex.write("\\\\")
    latex.write('\\nonumber\n\\end{eqnarray}\n')
    #/ Building ODEs
    latex.write('\\section{ODEs}\n')
    latex.write('\\tiny\n\\begin{eqnarray}')
    for node in model.nodes:
        nodeLatex = node.replace('_', '\_')
        latex.write('\n\\frac{\mathrm{d}[%s]_{i}}{\mathrm{d}t} & = & ' %str(nodeLatex))
        if node in model.constantNodes:
            latex.write('+ 0')
        else:
            for definition in model.topology.keys():
                fieldsDefinition = definition.split("/")
                if fieldsDefinition[0] == node:
                    if metamodel.modelFormat == 'tags':
                        latex.write('%s' %(definition.split('/')[1]))
                    else:
                        rLatex = definition.split('/')[3].split('.')[1].replace('_', '\_')
                        if model.topologyAST[definition][0] == '':
                            stoichiometry = model.topology[definition].split()[0]
                        else: 
                            stoichiometry = model.topologyAST[definition][0]
                        latex.write('%s %s' %(stoichiometry, rLatex))    
        latex.write(' \\nonumber \\\\\\nonumber & & \\\\')
    latex.write('\\nonumber\n\\end{eqnarray}\n')
    latex.write('\end{document}')
    latex.close()    
    
    return None

def eventsDealer(metamodel, model, outputfiles):

    '''
    This function simulates the model in the case it holds events.
    It runs the different simulations depending on the number of events and
    it joins them finally.
    '''

    #/ 1.- checking for the trigger
    triggerCheckPointValues = []
    for event in model.events:
	triggerVariable = event.trigger.split('(')[1].split(',')[0]
	if triggerVariable == 'time':
	    triggerCheckPoint = event.trigger.split('(')[1].split(',')[1].split(')')[0].split(' ')[1]
	    triggerCheckPointValue = model.parameters[triggerCheckPoint]
	    triggerCheckPointValues.append(triggerCheckPointValue)
	else:
	    raise errorMessages.SimulatorException,  'this version of ByoDyn only supports time triggered events. Sorry for the inconvenience.'
    #/ 2.- detecting that simulation time is larger than the time events
    for triggerCheckPointValue in triggerCheckPointValues:
	if triggerCheckPointValue >= metamodel.simulationTime:
	    raise errorMessages.ClassSimulatorException,  'simulation time (%s) is shorter that any of the Events (%s).'%(metamodel.simulationTime, triggerCheckPointValue)
	#/ 2.1.- detecting that the timestep is 100 times smaller than the first trigger
	if metamodel.simulationTimeStep >= triggerCheckPointValue/100:
	    print 'WARNING: simulation time step (%s) is equal or bigger than 1/100 times the time of one of the Events (%s).'%(metamodel.simulationTimeStep, triggerCheckPointValue)
    #/ 3.0.- remembering the total simulation time and the name of the SimulationResults file
    simulationEndingTime = metamodel.simulationTime
    simulationResultsFileName = outputfiles.simulationResults
    #/ 3.1.- sending the first simulation
    outputfiles.simulationResults = outputfiles.scratchdir + '/PartialSimulation' + '.0'
    metamodel.simulationTime = triggerCheckPointValues[0]
    centralFunctions.pythonIntegration(metamodel, model, outputfiles)
    #/ 3.2.- sending the simulations on the middle
    for i in range(len(triggerCheckPointValues)-1):
	model.initialConditions = lastStateRetriever(metamodel, model, outputfiles) # retrieving the concentration of the nodes after the simulation period, becoming the initial conditions of the next 
	outputfiles.simulationResults = outputfiles.scratchdir + '/PartialSimulation' + '.%s'%(i+1)
	metamodel.simulationTime = triggerCheckPointValues[i+1] - triggerCheckPointValues[i] # this is the time window to simulate
	centralFunctions.pythonIntegration(metamodel, model, outputfiles)
    #/ 3.3.- sending the last simulation
    model.initialConditions = lastStateRetriever(metamodel, model, outputfiles) 
    outputfiles.simulationResults = outputfiles.scratchdir + '/PartialSimulation' + '.%s'%((len(triggerCheckPointValues)))
    metamodel.simulationTime = simulationEndingTime - triggerCheckPointValues[len(triggerCheckPointValues)-1]
    centralFunctions.pythonIntegration(metamodel, model, outputfiles)
    #/ 3.4.- joining the partial simulations
    completeSimulations = []
    for i in range(len(triggerCheckPointValues)+1):
	file = outputfiles.scratchdir + '/PartialSimulation' + '.%s'%i
	f = open(file, 'r')
	partialSimulation = f.readlines()
	f.close()
	completeSimulations.append(partialSimulation)
    wholeSimulation = []
    for i in range(len(completeSimulations)):
	for j in range(len(completeSimulations[i])):
	    if i == 0:
		wholeSimulation.append(completeSimulations[i][j])
	    else:
		if completeSimulations[i][j].split()[0] != '#':
			simulationStep = completeSimulations[i][j].split()
			simulationStep[0] = str(float(simulationStep[0]) + triggerCheckPointValues[i-1])
			if simulationStep[0] != str(triggerCheckPointValues[i-1]):
			    wholeSimulation.append('\t'.join(simulationStep)+'\n')
    #/ 3.4.1.- writing the file and recovering the name
    outputfiles.simulationResults = simulationResultsFileName
    print simulationResultsFileName
    f = open(outputfiles.simulationResults, 'w', 1048576)
    for line in wholeSimulation:
	f.write(line)
    f.close()
    
    return None

def flatFileWriter(outputfiles, grid):

    '''
    This function writes a flat file with the information of the ClassEpithelium.
    '''

    f = open(outputfiles.epitheliumFlatFile, 'w')
    f.write('\ttime = %s\n\n\n' %grid.time)
    for i in range(len(grid.nodes)):
	f.write('%s\n\n' %(grid.nodes[i]))
	for j in range(grid.width):
	    for k in range(grid.length):
		value = grid.concentration[j * grid.length * len(grid.nodes) + k * len(grid.nodes) + i]
		f.write('\t%s' %value)
		grid.color[i].append(value)
	    f.write('\n')
    f.close()

    return None

def isAGene(node):

    '''
    This function returns a boolean determining if the node is or not a gene.
    The function's answer is based on whether the first character of the string is capital or not.
    This function is only used to determine the font of the node at the multicellular plot.
    '''
    if string.uppercase.find(node[0]) == -1:
	return True
    else:
	return False

def lastStateRetriever(metamodel, model, outputfiles):

    '''
    This function returns the values of the nodes at the last time of the simulation.
    '''

    lastState = []
    f = open(outputfiles.simulationResults, 'r')
    chunkSize = 100
    fd = f.fileno()
    size = os.fstat(fd)[6]
    f.seek(size)
    lines = []
    count = 1
    while len(lines) <= 1:
	bytes = -1*count*chunkSize #/ bytes has to be minus to indicate going backwards
	f.seek(bytes, 2)
	lines = f.readlines()
	count = count + 1
    lastLine = lines[len(lines) -1]
    elements = lastLine.split()
    for i in range(len(elements)-1):
	lastState.append(float(elements[i+1]))
    
    return lastState

def obtainSimulationValues(metamodel, model, outputfiles, type):
    
    '''
    This function obtains the simulation values of a given model.
    It is used by the optimiser and the sensitivityAnalyzer modules.
    '''
    simulationValues = []
    #/ 1.1.- running the simulation
    if len(type.split('.')) != 1:
        outputfiles.simulationResults = '%s/%s.%s.out' %(outputfiles.scratchdir, model.systemName, type.split('.')[1])
    else:      
        outputfiles.simulationResults = '%s/%s.%s.out' %(outputfiles.outputdir, model.systemName, parallel.currentProcessor())
    if metamodel.integrationOption == 'automatic':
        #/ For the automatic integrationOption, we need to check if scipy, octave, xpp and openModelical are available on the current machine.
        tester = checker.ClassChecker()
        tester.detector()
        #/ Next we adjust the integrationOption depending of model characteristics.
        metamodel.integrationOption = setIntegrationOption(model, tester)       
    if metamodel.integrationOption == 'python':
        centralFunctions.pythonIntegration(metamodel, model, outputfiles)
    elif metamodel.integrationOption == 'octave':
        centralFunctions.octaveIntegration(metamodel, model, outputfiles)
    elif metamodel.integrationOption == 'openModelica':
        simulator = simulatorOpenModelica.ClassSimulatorOpenModelica()	#/ Creating the simulator object
    elif metamodel.integrationOption == 'xpp':
        simulator = simulatorXPP.ClassSimulatorXPP()    #/ Creating the simulator object
    elif metamodel.integrationOption == 'rungeKutta':
        simulator = simulatorRungeKutta.ClassSimulatorRungeKutta()
    elif metamodel.integrationOption == 'euler':
        simulator = simulatorEuler.ClassSimulatorEuler()
    elif metamodel.integrationOption == 'stochastic':
        simulator = simulatorStochastic.ClassSimulatorStochastic()
    else:
        raise errorMessages.ClassSimulatorException,  'check the variable "integrationMethod".'
    if metamodel.integrationOption == 'openModelica' or metamodel.integrationOption == 'xpp':
        simulator.createInput(metamodel, model, outputfiles)
        simulator.callSolver(outputfiles)
        simulator.createOutputs(model, outputfiles, metamodel)
    #/ 1.2.- obtaining the simulation values from the file. 
    if type == 'simulation' or type == 'parameterEstimation':
        createGnuplot(model, metamodel, outputfiles, type)
    dataFile = open(outputfiles.simulationResults, 'r')
    originalNodes = copy.deepcopy(model.nodes)
    for line in dataFile:
        fields = line.split()
        if fields[0] != '#':
            fields.pop(0)
            simulationValues.append(fields)
    dataFile.close()

    return simulationValues

def plotMaker(outputfiles, grid, model):

    '''
    This function creates the postscript files for the multicellular plots.
    '''

    libdir = os.environ.get('BYODYN_PATH') + '/lib'
    #/ 1.- first obtain some information at the library necessary for the postscript file
    headerFile = open('%s/postscript/header.ps' %libdir, 'r')
    header =  headerFile.readlines()
    headerFile.close()
    #/
    boxFile = open('%s/postscript/box.ps' %libdir, 'r')
    box =  boxFile.readlines()    
    boxFile.close()
    #/
    frameFile = open('%s/postscript/frame.ps' %libdir, 'r')
    frame =  frameFile.readlines()
    frameFile.close()
    #/
    normalfontFile = open('%s/postscript/normalFont.ps' %libdir, 'r')
    normalfont = normalfontFile.readlines()
    normalfontFile.close()
    #/
    italicfontFile = open('%s/postscript/italicFont.ps' %libdir, 'r')
    italicfont = italicfontFile.readlines()
    italicfontFile.close()
    #/
    normalSmallFontFile = open('%s/postscript/normalSmallFont.ps' %libdir, 'r')
    normalSmallFont = normalSmallFontFile.readlines()
    normalSmallFontFile.close()
    #/ 2.- writing in the file
    f = open(outputfiles.epitheliumPlot, 'w')
    for line in header:
	f.write('%s' %line)
    for line in normalfont:
	f.write('%s' %line)
    for line in frame:
	f.write('%s' %line)
    for line in box:
	f.write('%s' %line)
    #/ 3.- calculating the displacement
    xDis = 0.0
    yDis = 0.0
    #/ 4.-calculating the dimentions of the boxes
    boxLength = 60.0 / grid.length
    boxWidth = 60.0 / grid.width
    xIni = 100.0
    yIni = 700.0 - boxWidth
    xFin = 100.0 + boxLength
    yFin = 700.0
    R = 0.0; G = 0.0; B = 0.0
    #/ 5.- calculating the concentrations
    localMax =[]
    for i in range(len(grid.nodes)):
	localMax.append(max(grid.color[i]))
    maxValue = max(localMax)
    nodeMatrix = 0 #/ just for counting. We'll show twelve matrices for page
    f.write('newpath 250 750  moveto (%s)show\n' %model.systemName) # setting the title
    for i in range(len(grid.nodes)):
	#/ painting the grid
	for j in range(grid.width):
	    for k in range(grid.length):
		R =  grid.color[i][(j * grid.length) + k] / maxValue
		G = B = 1.0 - R
		R = 1.0
		f.write('%s %s %s %s %s %s %s boxRGB\n' %(xIni + xDis, yIni + yDis, xFin + xDis, yFin + yDis, R, G, B))
		xIni = xIni + boxLength; xFin = xFin + boxLength
	    yIni = yIni - boxWidth; yFin = yFin - boxWidth
	    xIni = 100.0; xFin = 100.0 + boxLength
	yIni = 700.0 - boxWidth; yFin = 700.0
	f.write('%s %s %s %s %s %s 0.9 grid\n' %(100.0 + xDis, 640.0 + yDis, 160.0 + xDis, 700 + yDis, grid.length, grid.width)) #/ setting the grid
	if string.uppercase.find(grid.nodes[i][0]) != -1 :
	    for line in normalSmallFont:
		f.write('%s' %line)
	elif string.uppercase.find(grid.nodes[i][0]) == -1:
	    for line in italicfont:
		f.write('%s' %line)
	else:
	    raise errorMessages.ClassSimulatorException, 'there is an error with a name of a node while plotting. Check that all nodes start with a letter.'
	f.write('newpath %s %s  moveto (%s)show\n' %(100.0 + xDis, 620 + yDis, grid.nodes[i])) #/ setting the node
	nodeMatrix = nodeMatrix + 1
	xDis = xDis + 100.0
	if nodeMatrix % 4 == 0:
	    xDis = 0.0
	    yDis = yDis - 110.0
	if nodeMatrix == 20:
	    nodeMatrix = 0
	    f.write('showpage\n\n')
    f.write('newpath 250 100 moveto (t = %s)show\n' %(grid.time)) #/ setting the time
    f.close()

    return None

def setIntegrationOption(model, tester):
    '''
    This function set the integration option depending of model characteristics
    and the software installed on the current machine. 
    '''
    #/ Ranking the integration methods availables. Developers can change the order depending of robustness and velocity of each method.
    preferences = [] #/ a list sorted by our preferences about integration methods  
    if tester.ScipyFound == True:
        preferences.append('python')
    if tester.OpenModelicaFound == True:
        preferences.append('openModelica')
    if tester.OctaveFound == True:
        preferences.append('octave')
    if tester.XPPFound == True:
        preferences.append('xpp')
    if len(preferences) == 0:
        raise errorMessages.ClassSimulatorException, 'Error: The automatic integration method cannot detect any integration software in your system. Please, check the installation.'
    #/ Choosing integration method  
    if len(model.events) != 0:
        if 'openModelica' in preferences:
            return 'openModelica'
        else:
            raise errorMessages.ClassSimulatorException,  'Error: The SBML model have events and Byodyn only work with events using openModelica. But, you do not have openModelica installed in your system. Installation guide contains more detail information about the installation of optional requeriments.'
    elif model.delayFunctions == True:
        if 'xpp' in preferences:
            return 'xpp'
        else:
            raise errorMessages.ClassSimulatorException, 'Error: The SBML model have delay function and Byodyn only work with delay function using XPPAUT. But, you do not have openModelica installed in your system. Installation guide contains more detail information about the installation of optional requeriments.'
    elif len(model.rules) != 0:
        integrator = preferences[0]
        for rule in model.rules:
            if rule.type == 'Assignment':
                if 'python' in preferences:
                    preferences.remove('python')
                if len(preferences) == 0:
                    raise errorMessages.ClassSimulatorException, 'Error: The SBML model have assignment rules, you cannot use scipy but Byodyn work with assignment rules using openModelica, Octave or XPP. But, you only have installed scipy in your system. Installation guide contains more detail information about the installation of optional requeriments.'                    
            if rule.type == 'Algebraic':
                if 'python' in preferences:
                    preferences.remove('python')
                if 'xpp' in preferences:
                    preferences.remove('xpp')
                if len(preferences) == 0:
                    raise errorMessages.ClassSimulatorException, 'Error: The SBML model have algebraic rules and Byodyn only work with assignment rules using openModelica or Octave. But you do not have OpenModelica nor Octave installed scipy in your system. Installation guide contains more detail information about the installation of optional requeriments.'                    
        return preferences[0]
    else:
        return preferences[0]
    
