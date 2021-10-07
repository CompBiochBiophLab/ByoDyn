#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana
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

# $Id: optimiser.py,v 4.16 2008/12/03 16:37:38 alglomana Exp $

## \file
# This module is responsible of the optimisation of fitness function. 
# It directs the flow of the program to the genetic algorithm or the local search.
# It contains most of the parallel code.

from libsbml import *
import random, copy, os, math, sys, time
import central, centralFunctions, sbmlWorker, errorMessages, localOptimiser, simulatorOpenModelica, simulatorXPP, parallel, identifiabilityAnalyzer, sensitivityAnalyzer, simulator, checker, initiator

def calculateConfidenceIntervals(model, metamodel, outputfiles):

    '''
    This function calculates the covariance matrix and returns the object identifiability.
    This object is necessary to calculate the confidence intervals.
    The implementation of confidence intervals is still under development, not finished yet.
    '''
    
    sensitivity = sensitivityAnalyzer.ClassSensitivityAnalyzer(model.parameters)
    sensitivity.calculateSens(model, metamodel, outputfiles)   
    #/ 2.- Identifiability analysis
    identifiability = identifiabilityAnalyzer.ClassIdentifiabilityAnalyzer()
    identifiability.setSensitivity(sensitivity.identifiabilityCoefficients)
    identifiability.calculateFIM(model.parameters, metamodel.target)
    identifiability.setCOV(model.parameters)
    
    return identifiability

def central(metamodel, model, outputfiles, solutions):

    '''
    This function directs the optimisation to the different optimisation options.
    '''

    #/ PARALLEL: MAIN PROCESSOR AND IDENTATION
    if (parallel.mainProcessor()):
    # 1.- checking consistency of the parameters
        checkParametersToVary(model, metamodel)
        checkDataPoints(model, metamodel)
        checkTargetNodes(model, metamodel)
	checkInitialConditionsToVary(model, metamodel)
    # 2- starting the optimisation loop
    stopping = False
    iteration = 1
    #/ PARALLEL: Karyotype global 
    karyotype = None
    while stopping == False:
        if metamodel.optimisationMethod == 'geneticAlgorithm' or metamodel.optimisationMethod == 'hybridTwoPhases' or metamodel.optimisationMethod == 'hybridOnePhase':
            #/ PARALLEL: MAIN PROCESSOR AND IDENTATION
            if (parallel.mainProcessor()):
                if iteration == 1:
                    karyotype = gaInitialisePopulation(model, metamodel)
                #/ If the option runWithBackup is activated, we will write the GA information in a file
                if metamodel.runWithBackup == True:  
                    if iteration == 1:
                        #/ If the file exists, we will start the optimation from this point
                        if os.path.exists(outputfiles.runBackup) == True:
                            print 'WARNING: The file %s already exists in your output directory, we start the optimzation from this point..' %outputfiles.runBackup   
                            f = open(outputfiles.runBackup, 'r')
                            for line in f:
                                if line.split()[0] == 'GENERATION':                                                
                                    iteration = int(line.split()[1]) + 1
                                if line.split()[0] == 'CHROMOSOME':
                                    karyotype[int(line.split()[1])][line.split()[2]] = float(line.split()[3])
                            f.close()
    	    if metamodel.optimisationMethod == 'hybridOnePhase': #/ when we want to run a local for each of the genetic algorithm elements
                #/ PARALLEL: MAIN PROCESSOR AND IDENTATION
                if (parallel.mainProcessor()):
                    scores, karyotype = gaPopulationScoresObtainerLocal(iteration, karyotype, metamodel, model, outputfiles)
    	    else: #/ regular fitness evaluation
                scores = gaPopulationScoresObtainer(iteration, karyotype, metamodel, model, outputfiles)
            #/ PARALLEL: MAIN PROCESSOR AND IDENTATION
            if (parallel.mainProcessor()):
                bestChromosomes, sexChromosomes = gaNaturalSelection(karyotype, iteration, scores, metamodel, outputfiles)
                karyotype = gaSex(karyotype, bestChromosomes, sexChromosomes, model, metamodel)
                #/ If the option runWithBackup is activated, we will write the GA information in a file 
                if metamodel.runWithBackup == True:
                    f = open(outputfiles.runBackup, 'w')
                    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
                    f.write('GENERATION\t%s\n' %iteration)
                    for chromosome in karyotype:
                         for gene in karyotype[chromosome]:
                             f.write('CHROMOSOME\t%s\t%s\t%s\n' %(chromosome, gene, karyotype[chromosome][gene])) 
                    f.close()
        elif metamodel.optimisationMethod == 'randomSearch':
            #/ PARALLEL: MAIN PROCESSOR AND IDENTATION
            if (parallel.mainProcessor()):
                score = randomSearch(metamodel, model, outputfiles)
                evaluateSolution(iteration, score, metamodel, model, outputfiles)
                scores = []; scores.append(score) # needed for the "evaluateStopping" function
        elif metamodel.optimisationMethod == 'localSearch':
            #/ PARALLEL: MAIN PROCESSOR AND IDENTATION
            if (parallel.mainProcessor()):
                #/ 2.1.- determining the starting points
                #/ 2.1.1.- randomly if they are not defined
                for parameterToVary in metamodel.parametersToVary:
                    if parameterToVary.split('/')[3] == 'lin':
                        model.parameters[parameterToVary.split('/')[0]] = linearVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))
                    elif parameterToVary.split('/')[3] == 'log': #/ logarithmic scale
                        model.parameters[parameterToVary.split('/')[0]] = logarithmicVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))
                    else:
                        raise errorMessages.ClassOptimiserException, 'there is a problem with the way of exploring the parameters, either linear of logarithmic. Other value given.'
		#/ 2.1.2.- or from a given point from the runner
		for parameter in metamodel.parameters:
		    for modelParameter in model.parameters:
			if parameter == modelParameter:
			    model.parameters[parameter] = metamodel.parameters[parameter]
		#/ 2.1.3.- starting from a normally distributed initial condition, in case of
		if metamodel.initialConcentrationsToVary != []:
		    model = normalVariationForInitialCondition(metamodel, model)
		#/ 2.1.4.- storing the random search values
		score = scoreObtainer(metamodel, model, outputfiles)
		evaluateSolution(iteration, score, metamodel, model, outputfiles)
		#/ 2.1.5.- launching the calculation
		score = localOptimiser.central(model, metamodel, outputfiles) 
		#/ 2.1.6.- needed for the "evaluateStopping" function
		scores = [] #/ needed for the "evaluateStopping" function
	else:
	    raise errorMessages.ClassOptimiserException, 'second argument of variable "runningType" has not been recognized.'
        #/ PARALLEL: MAIN PROCESSOR AND IDENTATION
        if (parallel.mainProcessor()):
            stopping = evaluateStopping(metamodel, iteration, scores)
            iteration = iteration + 1
            parallel.sendAll(stopping, 0)
            parallel.sendAll(iteration, 1)
        stopping = parallel.receiveAny(0)
        iteration = parallel.receiveAny(1)
    
    return None

def checkDataPoints(model, metamodel):

    '''
    This function checks for the consistency of the target.
    It checks if the target constrain has a time value larger than the simulation time.
    '''

    for targetPoint in metamodel.target.keys():
         if float(targetPoint) > metamodel.simulationTime:
             raise errorMessages.ClassOptimiserException,  'the time you are targetting (%s) is longer than the simulation time (%s).'%(targetPoint, metamodel.simulationTime)
	 
    return None

def checkInitialConditionsToVary(model, metamodel):

    '''
    This function checks for the consistency of the initial conditions to vary.
    It checks that the values of the initial conditions are part of the nodes.
    Apart it checks for incompatibilities about model boundary conditions.
    '''

    #/ 1.- initial conditions to vary should be model nodes
    for element in metamodel.initialConcentrationsToVary:
	putativeNode = element.split('/')[0]
	if model.nodes.count(putativeNode) == 0 and model.algebraicNodes.keys().count(putativeNode)  == 0 and model.nonConstantParameters.keys().count(putativeNode) == 0 and model.nonConstantCompartments.keys().count(putativeNode) == 0:
	    raise errorMessages.ClassOptimiserException,  'you are trying to set an initial condition to vary node "%s" which is not part of the variables of the model: %s %s %s %s' %(evaluatingPoint.split('/')[0], model.nodes, model.algebraicNodes.keys(), model.nonConstantParameters.keys(), model.nonConstantCompartments.keys())
    #/ 2.- initial conditions to vary should not be boundary conditions on the model
    for constant in model.constantNodes.keys():
	for ic in metamodel.initialConcentrationsToVary:
	    if constant == ic.split('/')[0]:
		raise errorMessages.ClassOptimiserException, 'the species %s has a boundary condition and therefore the initial condition cannot be optimised. Please change the model accordingly.'%constant

    return None

def checkParametersToVary(model, metamodel):
    
    '''
    This function checks for the consistency of the introduced parameters to vary:
    the parameter to vary has to be one of the model's,
    the parameter to vary cannot be a constant parameter and
    the lower value of the exploring range has to be specified before the higher value.
    '''

    if metamodel.parametersToVary == [] and metamodel.target.keys().count(0.0) == 0:
	raise errorMessages.ClassOptimiserException,  'the field parametersToVary in runner file is empty.'
    else:
	for parameterToCheck in metamodel.parametersToVary:
	    consistency = False
	    for parameter in model.parameters.keys():
		if parameterToCheck.split('/')[0] == parameter:
		    consistency = True
            for parameter in model.nonConstantParameters.keys():
                if parameterToCheck.split('/')[0] == parameter:
                    consistency = True
	    if consistency == False:
		raise errorMessages.ClassOptimiserException,  'the parameter to vary %s is not found at the model\'s parameters: %s ' %(parameterToCheck, model.parameters)
	    if float(parameterToCheck.split('/')[1]) > float(parameterToCheck.split('/')[2]):
		raise errorMessages.ClassOptimiserException,  'the minimum value of the parameter to vary %s should be at the first position.' %(parameterToCheck.split('/')[0])
    
    return None
        
def checkTargetNodes(model, metamodel):

    '''
    This function checks that the nodes of the experimental data are part of the model.
    '''

    for targetTime in metamodel.target:
        for evaluatingPoint in metamodel.target[targetTime]:
            if model.nodes.count(evaluatingPoint.split('/')[0]) == 0 and model.algebraicNodes.keys().count(evaluatingPoint.split('/')[0])  == 0 and model.nonConstantParameters.keys().count(evaluatingPoint.split('/')[0]) == 0 and model.nonConstantCompartments.keys().count(evaluatingPoint.split('/')[0]) == 0:
                raise errorMessages.ClassOptimiserException,  'you are trying to target node "%s" which is not part of the variables of the model:%s %s %s %s' %(evaluatingPoint.split('/')[0], model.nodes, model.algebraicNodes.keys(), model.nonConstantParameters.keys(), model.nonConstantCompartments.keys())
    
    return None

def evaluateIteration(karyotype, model, metamodel, outputfiles, scores, iteration):

    '''
    This function helps the evaluation of the stopping of the program.
    It evaluates at the same time the stopper options and the stopper variables.
    It modifies the stopper options for the proper evaluation of the stopping by evaluateStopping.
    '''

    if metamodel.stopper[0] == 'iteration':
        if int(metamodel.stopper[1]) <= iteration: #/ that means that it is the last iteration and the stopper is not the score but the number of iterations. If so, we pick up the best element and evaluate if it is subject for a local search.
            metamodel.stopper[0] = 'iterationFinished'
            orderedScores = copy.deepcopy(scores)
            orderedScores.sort()
            score = orderedScores[0]
            model.parameters = karyotype[scores.index(orderedScores[0])]
            evaluateSolution(iteration, score, metamodel, model, outputfiles)
    elif metamodel.stopper[0] == 'score':
        for score in scores:
            if score <= float(metamodel.stopper[1]):
                 metamodel.stopper[0] = 'iterationFinished'
                 orderedScores = copy.deepcopy(scores)
                 orderedScores.sort()
                 score = orderedScores[0]
                 model.parameters = karyotype[scores.index(orderedScores[0])]
                 evaluateSolution(iteration, score, metamodel, model, outputfiles)
                 metamodel.stopper[0] = 'score'	
    elif metamodel.stopper[0] == 'numberOfSimulations':
	if int(metamodel.stopper[1]) <= metamodel.simulationNumber:
            metamodel.stopper[0] = 'iterationFinished'
            orderedScores = copy.deepcopy(scores)
            orderedScores.sort()
            score = orderedScores[0]
            model.parameters = karyotype[scores.index(orderedScores[0])]
            evaluateSolution(iteration, score, metamodel, model, outputfiles)
	    metamodel.stopper[0] = 'numberOfSimulations'	
    
    return None

def evaluateSolution(iteration, score, metamodel, model, outputfiles):
    
    '''
    This function evaluates if the solution is susceptible of being stored.
    It involves the creation of solution directories and files.
    '''

    #/ 1.- if the solution is comming from a random search, just store the new file    
    if metamodel.optimisationMethod == 'randomSearch' or metamodel.optimisationMethod == 'localSearch':
	#/ 0.- check if it is the first solution ever, and you have to create a solutions directory with the name of the system
        if os.path.exists(outputfiles.solutionDirectory) == False:
            os.mkdir(outputfiles.solutionDirectory)
	if iteration == 1:
	    metamodel.lastIterationScore = score
	    storeSolution(score, metamodel, model, outputfiles.solutionDirectory)
	else:
	    if score < metamodel.lastIterationScore:
		metamodel.lastIterationScore = score
		storeSolution(score, metamodel, model, outputfiles.solutionDirectory)
	#/ 1.- saving the initial conditions in case of
	if metamodel.initialConcentrationsToVary != []:
	    file = outputfiles.solutionDirectory + '/initialConditionsFromRandom'
	    #/ 1.1.1.- writing the header for the first time the file is created
	    if os.path.exists(file) == False:
		f = open(file, 'w')
		f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
		f.close()
            #/ 1.1.2.- appending information to the files
	    f = open(file, 'a')
	    for ic in model.initialConditions:
		f.write('%s\t'%ic)  
	    f.write('\n')
	    f.close()
	    #/ 1.2.- writing the score of the fitness function in the case of absence of parameters to vary
	    if metamodel.parametersToVary == [] and metamodel.optimisationMethod == 'randomSearch':
		file = outputfiles.solutionDirectory + '/scoresForInitialConcentrationsFromRandomSearch'
		#/ 1.2.1.- setting the headers
		if os.path.exists(file) == False:
		    f = open(file, 'w')
		    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
		    f.close()
		#/ 1.2.2.- appending the information
		f = open(file, 'a')
		f.write('%s\n'%score)
		f.close()
	    #/ 1.2.3.- for the case of the local search
	    if metamodel.parametersToVary == [] and metamodel.optimisationMethod == 'localSearch':
		file = outputfiles.solutionDirectory + '/scoresForInitialConcentrationsFromRandomSearch'
		#/ setting the headers
		if os.path.exists(file) == False:
		    f = open(file, 'w')
		    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
		    f.close()
		#/ appending the information
		f = open(file, 'a')
		f.write('%s\n'%score)
		f.close()
        #/ 1.3.- storing the parameters
        parametersDirectory = outputfiles.solutionDirectory + '/parametersFromRandom'
        if os.path.exists(parametersDirectory) == False:
            os.mkdir(parametersDirectory)
	#/ 1.4.- writing the header for the first time the file is created
        for parameter in metamodel.parametersToVary:
            file = parametersDirectory + '/%s' %parameter.split('/')[0]
	    if os.path.exists(file) == False:
		f = open(file, 'w')
		f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
		f.close()
	    #/ 1.5.- appending information to the files
	    f = open(file, 'a')
	    f.write('%s\t%s\t%s\t%s\n'%(model.parameters[parameter.split('/')[0]], parameter.split('/')[1], parameter.split('/')[2], score))
	    f.close()
    elif score <= metamodel.threshold or metamodel.stopper[0] == 'iterationFinished': 
        #/ 0.- check if it is the first solution ever, and you have to create a solutions directory with the name of the system
        if os.path.exists(outputfiles.solutionDirectory) == False:
            os.mkdir(outputfiles.solutionDirectory)
        #/ 2.- if the solution is comming from a genetic algorithm search ...
        if metamodel.optimisationMethod == 'geneticAlgorithm':
            #/ 2.1.- store the parameters
            parametersDirectory = outputfiles.solutionDirectory + '/parametersFromGA'
            if os.path.exists(parametersDirectory) == False:
                os.mkdir(parametersDirectory) 
	    #/ storing the last best value, check the user reference in case of doubt, section 2.3.4.3.3
            for parameter in metamodel.parametersToVary:
                file = parametersDirectory + '/%s' %parameter.split('/')[0]
                if os.path.exists(file) == False:
                    f = open(file, 'w')
                    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
                f = open(file, 'a')
		f.write('%s\t%s\t%s\t%s\n'%(model.parameters[parameter.split('/')[0]], parameter.split('/')[1], parameter.split('/')[2], score))	    
            #/ 2.2.- create a directory with the generation number
	    storeSolution(score, metamodel, model, outputfiles.solutionDirectory)
	    #/ 2.3.- optionally we can calculate the confidence intervals of the optimisation values.
        if metamodel.confidenceIntervals == True:
            identifiability = calculateConfidenceIntervals(model, metamodel, outputfiles)
            for parameter in metamodel.parametersToVary:
                file = parametersDirectory + '/%s' %parameter.split('/')[0]
                #/ writing the header for the first time the file is created
                if os.path.exists(file) == False:
                    f = open(file, 'w')
                    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
                f = open(file, 'a')
		if metamodel.confidenceIntervals == False:
		    f.write('%s\t%s\t%s\t%s\n'%(model.parameters[parameter.split('/')[0]], parameter.split('/')[1], parameter.split('/')[2], score))
		elif metamodel.confidenceIntervals == True:
		    ci = identifiability.confidenceIntervals(model.parameters[parameter.split('/')[0]], model.parameters.keys().index(parameter.split('/')[0]))
		    f.write('%s\t%s\t%s\t%s\t%s\n'%(model.parameters[parameter.split('/')[0]], ci, parameter.split('/')[1], parameter.split('/')[2], score))
                f.close()
        elif metamodel.optimisationMethod == 'hybridTwoPhases' or metamodel.optimisationMethod == 'hybridOnePhase':
    	    #/ defining the name of the directory
    	    if metamodel.optimisationMethod == 'hybridTwoPhases':
                parametersDirectory = outputfiles.solutionDirectory + '/parametersFromGA'
    	    elif metamodel.optimisationMethod == 'hybridOnePhase':
                parametersDirectory = outputfiles.solutionDirectory + '/parametersFromHybridOnePhase'
		#/ saving the best element in an SBML file
		sbmlWorker.sbmlWriter(SBMLWriter(), model, metamodel, outputfiles.solutionDirectory + '/%s.xml' %model.systemName)
	    if os.path.exists(parametersDirectory) == False:
		os.mkdir(parametersDirectory)
	    for parameter in metamodel.parametersToVary:
		file = parametersDirectory + '/%s' %parameter.split('/')[0]
	        #/ writing the header for the first time the file is created
		if os.path.exists(file) == False:
		    f = open(file, 'w')
		    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
		    f.write('%s\t%s\t%s\t%s\n'%(model.parameters[parameter.split('/')[0]], parameter.split('/')[1], parameter.split('/')[2], score)) #/ this line is necessary for the first time. The next is for appending in the next runs.
		else:
		    f = open(file, 'a')
		    f.write('%s\t%s\t%s\t%s\n'%(model.parameters[parameter.split('/')[0]], parameter.split('/')[1], parameter.split('/')[2], score))
		    f.close()	    
    	    if metamodel.stopper[0] == 'iterationFinished' and metamodel.optimisationMethod == 'hybridTwoPhases': #/ just in the case of the hybrid, in the case of the "hybridOnePhase" the optimisation is finnished already
                localOptimiser.central(model, metamodel, outputfiles) #/ launching the local search of the hybridTwoPhases
		sbmlWorker.sbmlWriter(SBMLWriter(), model, metamodel, outputfiles.solutionDirectory + '/%s.xml' %model.systemName) #/ storing the solution

	  
    return None
    
def evaluateStopping(metamodel, iteration, scores):
    
    '''
    This function evaluates if the program needs to stop due to the score value or the iteration number.
    '''
    
    if metamodel.optimisationMethod == 'localSearch' and metamodel.stopper[0] == 'score':
        stopping = True
    elif metamodel.optimisationMethod == 'parallelGA':
	stopping = True
    else:
        if metamodel.stopper[0] == 'score':
            stopping = False
            for score in scores:
                if float(metamodel.stopper[1]) >= score:
                    stopping = True
            return stopping
        elif metamodel.stopper[0] == 'iteration':
            if int(metamodel.stopper[1]) <= iteration:
                print 'Iteration %s completed from optimisation.' %(iteration)
		print 'Total number of model integrations has been %s.'%metamodel.simulationNumber
                return True
            else:
                return False
        elif metamodel.stopper[0] == 'iterationFinished':
            print 'Iteration %s completed from hybrid optimisation.'%(iteration)
	    print 'Total number of model integrations has been %s.'%metamodel.simulationNumber
            return True
	elif metamodel.stopper[0] == 'numberOfSimulations':
		if int(metamodel.stopper[1]) <= metamodel.simulationNumber:
			print 'The optimisation has arrived to the maximum number of simulations (%s).' %metamodel.simulationNumber
			return True
		else:
			return False
        else:
		raise errorMessages.ClassOptimiserException, 'Error in optimiser options.'

def gaInitialisePopulation(model, metamodel):

    '''
    This function create the initial random karyotype for the genetic algorithm.
    '''

    karyotype = {}
    for i in range(metamodel.gaPopulation):
        karyotype[i] = copy.deepcopy(model.parameters)
        for gene in karyotype[i]:
            for parameterToVary in metamodel.parametersToVary:
                if gene == parameterToVary.split('/')[0]:
		    if parameterToVary.split('/')[3] == 'lin':
			karyotype[i][gene] = linearVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))
		    elif parameterToVary.split('/')[3] == 'log': #/ logarithmic scale
			karyotype[i][gene] = logarithmicVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))    
		    else:
			raise errorMessages.ClassOptimiserException, 'there is a problem with the way of exploring the parameters, either linear of logarithmic. Other value given.'
    #/ 1.- setting back to the runner value in case the GA was initialised from an specific point
    for fixedParameter in metamodel.fixedParameters:
	for i in range(metamodel.gaPopulation):
	    for gene in karyotype[i].keys():
		if gene == fixedParameter[0]:
		    karyotype[i][gene] = fixedParameter[1]
    #/ 2.- appending the initial concentrations to calibrate
    if metamodel.initialConcentrationsToVary != []:
	for i in karyotype.keys():
	    for initialConcentrationToVary in metamodel.initialConcentrationsToVary:
		name = initialConcentrationToVary.split('/')[0]
		concentration = float(initialConcentrationToVary.split('/')[2])
		variance = float(initialConcentrationToVary.split('/')[3])
		value = -1
		while value < 0:
		    value = random.normalvariate(concentration,math.sqrt(variance))
		karyotype[i][name] = value
                                        
    return karyotype

def gaNaturalSelection(karyotype, iteration, scores, metamodel, outputfiles):

    '''
    This function selects the best fit elements from the last generation karyotype.
    '''
    
    #/ 0.- handling with the scores
    if metamodel.optimisationMethod != 'hybridOnePhase': #/ 0.1.- determining the global and the best score for the regular GA or hybrid
	globalScore = 0.0
	bestScore = 0.0
	for score in scores:
	    globalScore = globalScore + score
	globalScore = globalScore / len(scores)
	bestScore = min(scores)
        #/ 0.2.- writing the file
	if iteration == 1:
	    #/ writing the headers
	    stad = open(outputfiles.gaStatistics,"w")
	    stad.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
	    stad.write('Generation\t%s\t%s\t%s\n' %(iteration, globalScore, bestScore))
	    stad.close()
	elif iteration != 1 and metamodel.optimisationMethod != 'hybridOnePhase': #/ this line should not be written when a hybridOnePhase algorithm is running
	    stad = open(outputfiles.gaStatistics,"a")
	    stad.write('Generation\t%s\t%s\t%s\n' %(iteration, globalScore, bestScore))
	    stad.close()
    #/ 1.- choosing the best 10 % chromosomes and the 40 best for sex
    #/ 1.1.- selecting at least one for best
    bests = int(0.1 * metamodel.gaPopulation)
    if bests < 1:
	bests = 1
    sexuals = (metamodel.gaPopulation - 2 * bests) / 2
    #/ retrieving the correspondance of scores
    correspondance = {}
    for i in range(len(scores)):
            correspondance[i] = scores[i]
    scores.sort()
    scoresChoosen = scores[:sexuals]
    chromosomesChoosen = []
    for i in range(len(scoresChoosen)):
        for j in range(len(correspondance)):
            if scoresChoosen[i] == correspondance[j]:
                chromosomesChoosen.append(j)
                break #/ CRUCIAL
    #/ 3.- retrieving the chromosomes
    bestChromosomes = chromosomesChoosen[:bests]
    sexChromosomes = chromosomesChoosen[:sexuals]

    return bestChromosomes, sexChromosomes

def gaPopulationScoresObtainer(iteration, karyotype, metamodel, model, outputfiles): #/ PARALLEL: METAMODEL, MODEL & OUTPUTFILES ARE SHARED BY ALL PROCESSORS

    '''
    This function evaluates the score of the karyotype by calling scoreObtainer function.
    Special emphasis is given for the parallel code.
    '''

    #/ PARALLEL: ADDING CODE
    if (parallel.mainProcessor()):
        parallel.sendAll(iteration, 0)
        parallel.sendAll(karyotype, 1)
    iteration = parallel.receiveAny(0)
    karyotype = parallel.receiveAny(1)
    partitionSize = metamodel.gaPopulation / parallel.totalProcessors()
    partitionIndex = 0
    #/
    scores = []
    for chromosome in karyotype:
        #/ PARALLEL: PARTITION MUST BE IN THIS WAY BECAUSE EVALUATE_SOLUTION IS POSITION DEPENDING
        if (partitionSize * parallel.currentProcessor() <= partitionIndex and partitionIndex < partitionSize * (parallel.currentProcessor() + 1)):
            model.parameters = karyotype[chromosome]
            score = scoreObtainer(metamodel, model, outputfiles)
            evaluateSolution(iteration, score, metamodel, model, outputfiles)
            scores.append(score)
        partitionIndex = partitionIndex + 1
    parallel.send(scores, 0, 0)        
    #/ PARALLEL: ADDING CODE & IDENTATION
    if (parallel.mainProcessor()):
        fullScores = []
        for i in range (0, parallel.totalProcessors(), 1):
            tmpScores = parallel.receive(i, 0)
            fullScores.extend(tmpScores)
        evaluateIteration(karyotype, model, metamodel, outputfiles, fullScores, iteration)    
        
	return fullScores
    
    return None

def gaPopulationScoresObtainerLocal(iteration, karyotype, metamodel, model, outputfiles):

    '''
    This function evaluates the fitness function for the hybridOnePhase algorithm.
    For each element it launches a local search algorithm and takes as fitness function value the result of the local optimisation.
    '''
    
    #/ 1.- Running the calculations
    scores = [] 
    for chromosome in karyotype:
	model.parameters = karyotype[chromosome]
	score = float(localOptimiser.central(model, metamodel, outputfiles)) #/ be very careful with this function. localOptimiser writes the memory position of model.parameters. As karyotype elements point to each of the model.parameter, when launch the optmization, the karyotype is upgraded automatically
	scores.append(score)
    #/ 2.- Storing the data
    if iteration == 1:
	#/ writing the headers
	stad = open(outputfiles.gaStatistics,"w")
	stad.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
	if metamodel.optimisationMethod == 'hybridOnePhase':
	    stad.write('Generation\t%s' %iteration)
	    for score in scores:
		stad.write('\t%s'%score)
	    stad.write('\n')
	stad.close()
    else:
	stad = open(outputfiles.gaStatistics,"a")
	stad.write('Generation\t%s' %iteration)
	for score in scores:
	    stad.write('\t%s'%score)
	stad.write('\n')
	stad.close()
    evaluateIteration(karyotype, model, metamodel, outputfiles,  scores, iteration)

    return scores, karyotype

def gaSex(karyotype, bestChromosomes, sexChromosomes, model, metamodel):

    '''
    This function directs the crossing over and mutation from one to the next generation.
    '''

    #/ 1.- choosing the new generation
    bestKaryotype = {}
    sexKaryotype = {}
    randomKaryotype = {}
    #/ 1.1.- determining the best chromosomes
    for i in range(len(bestChromosomes)):
        bestKaryotype[i] = copy.deepcopy(karyotype[bestChromosomes[i]])
    #/ 1.2.- determining the sex chromosomes
    for i in range(len(sexChromosomes)):
        sexKaryotype[i] =  copy.deepcopy(karyotype[sexChromosomes[i]])
    #/ 1.3.- duplicating the sex chromosomes
    for i in range(len(sexChromosomes)):
        sexKaryotype[len(sexChromosomes) + i] =  copy.deepcopy(sexKaryotype[i])
    #/ 1.3.- determining 10 % of random chromosomes
    for i in range(len(bestChromosomes)):
        randomKaryotype[i] = copy.deepcopy(model.parameters)
        for gene in randomKaryotype[i]:
            for parameterToVary in metamodel.parametersToVary:
                if gene == parameterToVary.split('/')[0]:
                    if parameterToVary.split('/')[3] == 'lin':
                        randomKaryotype[i][gene] = linearVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))
                    elif parameterToVary.split('/')[3] == 'log': #/ logarithmic scale
                        randomKaryotype[i][gene] = logarithmicVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))
                    else:
                        raise errorMessages.ClassOptimiserException, 'there is a problem with the way of exploring the parameters, either linear of logarithmic. Other value given.'
    #/ 1.4.- killing one generation to start another one. That's life, dude !!!
    karyotype = {}
    #/ 1.5.- filing up the new generation
    for i in range(len(bestKaryotype)):
        karyotype[i] = bestKaryotype[i]
    for i in range(len(sexKaryotype)):
        karyotype[len(bestKaryotype) + i] = sexKaryotype[i]
    for i in range(len(randomKaryotype)):
        karyotype[len(bestKaryotype) + len(sexKaryotype) + i] = randomKaryotype[i]
    #/ 2.- mutation. Only sex chromosomes are susceptible of being mutated. And only genes that are parameters to vary.
    numberOfMutations = 0
    numberOfGenesForMutation = len(sexChromosomes) * 2 * len(metamodel.parametersToVary)
    for i in range(numberOfGenesForMutation):
        if random.random() <= metamodel.gaMutationRate:
            numberOfMutations = numberOfMutations + 1
    for i in range(numberOfMutations):
        #/ determining chromosome to mutate
        chromosomeIndex = int(len(bestKaryotype) + random.random() * len(sexKaryotype))
        #/ determining the gene to mutate
        mutatingGene = metamodel.parametersToVary[int(random.random() * len(metamodel.parametersToVary))]
        if mutatingGene.split('/')[3] == 'lin':
            karyotype[chromosomeIndex][mutatingGene.split('/')[0]] = linearVariation(float(mutatingGene.split('/')[1]), float(mutatingGene.split('/')[2]))
        elif mutatingGene.split('/')[3] == 'log':
            karyotype[chromosomeIndex][mutatingGene.split('/')[0]] = logarithmicVariation(float(mutatingGene.split('/')[1]), float(mutatingGene.split('/')[2]))
        else:
            raise errorMessages.ClassOptimiserException, 'there is a problem with the way of exploring the parameters, either linear of logarithmic. Other value given.'
    #/ 3.- translocation
    #/ 3.1.- number of translocations
    numberOfTranslocations = 0
    for i in range(len(sexChromosomes)):
        if random.random() <= metamodel.gaTranslocationRate:
            numberOfTranslocations = numberOfTranslocations + 1
    #/ 3.2.- translocating
    for i in range(numberOfTranslocations):
        maleChromosome = int(random.random() * (len(sexChromosomes) * 2)) + len(bestChromosomes)
        femaleChromosome = int(random.random() * (len(sexChromosomes) * 2)) + len(bestChromosomes)
        while maleChromosome == femaleChromosome:
            femaleChromosome = int(random.random() * len(sexChromosomes * 2)) + len(bestChromosomes)
        breakPoint =  1 +  int(random.random() * (len(karyotype[0]) - 1)) #/ get a pen and a paper and think about it. It's OK anyway ...
        #/ preparing the chromosomes: transforming them from dictionary to list
        male = []
        female = []
        for gene in karyotype[maleChromosome].keys():
            male.append([gene, karyotype[maleChromosome][gene]])
            female.append([gene, karyotype[femaleChromosome][gene]])
        #/ breaking the chromosomes
        maleReproducing = male[:breakPoint]
        maleComplementary = male[breakPoint:]
	femaleReproducing = female[breakPoint:]
        femaleComplementary = female[:breakPoint]
        #/ first joint ...
        maleReproducing.extend(femaleReproducing)
        #/ ... and second
        femaleComplementary.extend(maleComplementary)
        #/ new quimeras are in the world !
        male = maleReproducing
        female = femaleComplementary
        #/ transforming the lists back into dictionaries of the karyotype
        for i in range(len(male)):
           karyotype[maleChromosome][male[i][0]] = male[i][1]
           karyotype[femaleChromosome][female[i][0]] = female[i][1]
    
    return karyotype

def getSimulateValue(originalSimulationValues, node, target, model, metamodel):

    '''
    This function retrieves the simulation value of a single node at a single time point.
    '''

    column = model.nodes.index(node) - 1 # provisional
    row = int(float(target) / float(metamodel.simulationTimeStep))

    return float(originalSimulationValues[row][column])

def linearVariation(minimalValue, maximalValue):

    '''
    This function generates a random value for a parameter explored in a linear scale range.
    '''

    value = minimalValue + random.random() * (maximalValue - minimalValue)

    return value

def logarithmicVariation(minimalValue, maximalValue):

    '''
    This function generates a random value for a parameter explored in a logarithmic
    '''

    cocientOfInterval = maximalValue / minimalValue
    logValue = math.log(minimalValue, 10) + random.random() * math.log(cocientOfInterval, 10)
    value = 10**logValue

    return value

def normalVariationForInitialCondition(metamodel, model):

    '''
    This function variates the initial condition of the model based on a normal distribution.
    '''
    
    for newIC in metamodel.initialConcentrationsToVary:
	vector = newIC.split('/')
	node = vector[0]
	cell = vector[1]
	meanIC = float(vector[2])
	variance = float(vector[3])
	if len(model.initialConditions) != len(model.nodes):
	    raise errorMessages.ClassOptimiserException, 'detected a mismatch on the size of initial conditions and model nodes.'
	for i in range(len(model.nodes)):
	    if node == model.nodes[i]:
		newInitialCondition = random.normalvariate(meanIC,math.sqrt(variance))
		model.initialConditions[i] = newInitialCondition
    
    return model

def randomSearch(metamodel, model, outputfiles):

    '''
    This function evaluates the score of a new random point in the parameter space.
    '''

    #/ 1.- variate parameters
    for parameterToVary in metamodel.parametersToVary:
        for i in range(len(model.parameters.keys())):
            if parameterToVary.split('/')[0] == model.parameters.keys()[i]:
                if parameterToVary.split('/')[3] == 'lin':
                    model.parameters[parameterToVary.split('/')[0]] = linearVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))
                elif parameterToVary.split('/')[3] == 'log': #/ logarithmic scale
                    model.parameters[parameterToVary.split('/')[0]] = logarithmicVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))
                else:
                    raise errorMessages.ClassOptimiserException, 'there is a problem with the way of exploring the parameters, either linear of logarithmic. Other value given.'
    #/ 2.- recovering the fixed parameters
    for parameter in metamodel.parameters:
	    for modelParameter in model.parameters:
		if parameter == modelParameter:
		    model.parameters[parameter] = metamodel.parameters[parameter]
    #/ 3.- changing the initial conditions, in case of
    if model.xlength != 1:
	if model.ywidth != 1:
	    raise errorMessages.ClassOptimiserException,  'the possibility of estimate initial conditions is restricted to unicellular systems. Sorry for the inconvenience.' 
    if metamodel.initialConcentrationsToVary != []:
	model = normalVariationForInitialCondition(metamodel, model)
    #/ 4.- obtaining the score
    score = scoreObtainer(metamodel, model, outputfiles)
   
    return score

def scoreObtainer(metamodel, model, outputfiles):

    '''
    This function obtains the fitness function value given a model and the experimental target.
    It involves the simulation of the model.
    A set of target values are determined, with its corresponding importance, that is, the ponderation.
    Then a set of corresponding values at the simulation are searched. Because at the corresponding time there are several simulation values, depending on the time step, the simulation value of the earlier data point is taken. 
    The values are substracted from one list to the other.
    '''
    
    #/ 1.- Counting number of simulations
    metamodel.simulationNumber = metamodel.simulationNumber + 1
    #/ 2.- Checking assignment with the new values of the parameters
    for rule in model.rules:
        if rule.type == 'Assignment':
            model = model.checkAssignmentRule(rule)          
    #/ 3.- running the model
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
    elif metamodel.integrationOption == 'openModelica':
        simulatorObject = simulatorOpenModelica.ClassSimulatorOpenModelica()	#/ Creating the simulator object
    elif metamodel.integrationOption == 'xpp':
        simulatorObject = simulatorXPP.ClassSimulatorXPP()    #/ Creating the simulator object
    else:
        raise errorMessages.ClassOptimiserException,  'check the variable "integrationMethod".'
    if metamodel.integrationOption == 'openModelica' or metamodel.integrationOption == 'xpp':
        simulatorObject.createInput(metamodel, model, outputfiles)
        simulatorObject.callSolver(outputfiles)
        simulatorObject.createOutputs(model, outputfiles, metamodel)    
    #/ 4.- defining variables
    score = 0.
    simulationValues = []
    targetNodes = [] # necessary for the normalization
    targetValues = []
    targetValuesVariance = []
    #/ 5.- calculating the target values(per time, per cell) with their correspondant ponderation.
    targetTimePoints = []
    for i in range(len(metamodel.target.keys())):
        targetTimePoints.append(metamodel.target.keys()[i])
    targetTimePoints.sort()
    for i in range(len(targetTimePoints)):
        for singleEvaluation in metamodel.target[targetTimePoints[i]]:
            targetNodes.append(singleEvaluation.split('/')[0])
            targetValues.append(float(singleEvaluation.split('/')[2]))
            targetValuesVariance.append(float(singleEvaluation.split('/')[3]))
    #/ 5.- calculating the simulation values. THEY MUST BE CORRESPONDANT BECAUSE AFTERWARDS IS JUST A SUBSTRACTION BETWEEN VALUES, THERE'S NO FURTHER CHECKING
    #/ 5.1.- obtaining the values from the file. 
    dataFile = open(outputfiles.simulationResults, 'r')
    rawInformation = []
    for line in dataFile:
        fields = line.split()
        if fields[0] != '#':
            rawInformation.append(fields)
    dataFile.close()
    #/ 5.2.- appending all possible variables of the system of equations for the target: nodes, algebraic nodes, non constant parameters and non constant compartments
    variables = []
    for n in model.nodes:
        if model.constantNodes.keys().count(n) == 0:
            variables.append(n)
    for n in model.algebraicNodes:
        variables.append(n)
    for n in model.nonConstantParameters:
        variables.append(n)
    for n in model.nonConstantCompartments:
        variables.append(n)
    
    #/ 5.3.- obtaining the simulation values
    for i in range(len(targetTimePoints)):
        found = False
        for j in range(len(rawInformation)):
            if targetTimePoints[i] <= float(rawInformation[j][0]) and found == False:
                found = True
                for singleEvaluation in metamodel.target[targetTimePoints[i]]:
                    simulationValues.append(float(rawInformation[j][(int(singleEvaluation.split('/')[1].split(',')[0]) * model.xlength) + int(singleEvaluation.split('/')[1].split(',')[1]) * (len(variables)) + variables.index(singleEvaluation.split('/')[0]) + 1]))
    
    #/ 6.- subtracting the values
    scoreValues = []
    if len(targetValues) > len(simulationValues):
        print 'WARNING: Simulation seems to be truncated due to the combination of parameters. The objective function will be rescaled to a high value.'
        print 'Parameters values:'
        for parameter in metamodel.parametersToVary:
            if parameter in model.nonConstantParameters.keys():
                print '%s: %s' %(parameter.split('/')[0], model.nonConstantParameters[parameter.split('/')[0]])
            else:
		print '%s: %s' %(parameter.split('/')[0], model.parameters[parameter.split('/')[0]])
        for i in range(len(targetValues)):
            scoreValues.append(1e10*max(targetValues))
    else:
        for i in range(len(targetValues)):
            scoreValues.append((targetValues[i] - simulationValues[i])**2 / targetValuesVariance[i])
    #/ 7.- determining if the score is global or is for each data point
    if metamodel.hybridScore == False: #/ calculating the distance for the python modules
        for i in range(len(scoreValues)):
            score = score + scoreValues[i]
        score = score / len(scoreValues)
        return score 
    elif metamodel.hybridScore == True: #/ calculating the distance for the fortran modules. I am sending the SQUARED differences. Keep it in mind !!!
        for i in range(len(scoreValues)):
             scoreValues[i] = scoreValues[i] * (2./len(scoreValues))
        return scoreValues
    else:
        raise errorMessages.ClassOptimiserException, 'Error: it can be determined if the stopper is global, single or array.'
   
def storeSolution(score, metamodel, model, workingDirectory):
    
    '''
    This function saves a solution in the tags or SBML format.
    '''

    if metamodel.modelFormat == 'tags':
        print 'WARNING: tags solution storer is underconstruction.'
    elif metamodel.modelFormat == 'SBML':
        file = workingDirectory + '/%s.xml' %model.systemName
        #/ initializing the writer
        w = SBMLWriter()
        w.setProgramName('ByoDyn')
        w.setProgramVersion('with a score of %s' % score)
        sbmlWorker.sbmlWriter(w, model, metamodel, file)
    else:
        raise errorMessages.ClassOptimiserException, 'Error: it can be determinded if the model is in SBML format or in tags format'

    return None
