#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author:  Adrian L. Garcia-Lomana and Alex Gomez-Garrido
#
#  Created: 2004-08-30 by Adrian L. Garcia-Lomana
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

# $Id: central.py,v 4.17 2008/12/03 16:37:38 alglomana Exp $

## \file
# This module is the first module called from the initiator.
# This module contains the code to direct all the jobs of the execution.

import sys, os
import centralFunctions, starter, sbmlWorker, simulator, optimiser, tagParser, errorMessages, dynamicsReconstructer, surface, sensitivityAnalyzer, fitnessFunctionEvaluator, simulatorOpenModelica, identifiabilityAnalyzer, optimalExperimentalDesign, parallel, exporter, sampler, cluster

class ClassFile:

    ''' 
    This class specifies the paths for all required files that ByoDyn uses, both input and output files. 
    '''

    pass

def profile():

    ''' 
    This function determines the number of times ByoDyn is run for profiling. 
    '''

    for i in range(1): #/ number of runs of the program to perform profiling
        main()

##############################
#       Main Functions       #                                  
##############################

def modelReader(metamodel): 
    
    ''' 
    This function reads the model file (the name of the file is inside the metamodel object) and creates a model object with the model information and details. 
    It works for either formats, "tags" or "SBML". 
    '''

    if metamodel.runningType != 'clustering': #/ we do not need to build the model for the clustering
	if metamodel.modelFormat == 'SBML':
	    model = sbmlReader(metamodel)
	elif metamodel.modelFormat == 'tags':
	    model = tagsReader(metamodel)
	else:
	    raise errorMessages.ClassCentralException,  'check for "modelFormat" variable.'
	return model
    else:
	return sbmlWorker.ClassModelSBML() #/ just initialise the object to avoid errors
        
def optionReader(runnerFile): 

    ''' 
    This function reads the runner file of ByoDyn and creates the metamodel object that contains its information. 
    '''

    metamodel = starter.central(runnerFile)

    return metamodel

def runner(metamodel, model):

    ''' 
    This function discriminates the different main running options of ByoDyn and calls its appropiate functions. 
    The different possibilities are simulation, parameter estimation, sensitivity analysis, model format exporting, simulation trajectories reconstruction from a given parameter value set, creation of score surfaces of the fitness with respect of the vector of parameters and fitness function evaluation. 
    '''
    
    def fileInitializator(model, metamodel): 

	''' 
	This function defines the path to all files used by ByoDyn, both input files and output files.
	'''
	
        outputfiles = ClassFile()
        outputdir = os.environ.get('BYODYN_OUTPUT_DIR')
        scratchdir = os.environ.get('BYODYN_SCRATCH_DIR')
        outputfiles.outputdir = outputdir
        outputfiles.scratchdir = scratchdir
	outputfiles.clusterCentersTextFile = '%s/clusterCenters.txt'%outputdir
	outputfiles.clusterCode = '%s/clusteringCode.m'%scratchdir
	outputfiles.clusterResultsTextFile = '%s/clusterResults.txt'%outputdir
	outputfiles.clusterResultsFigure = '%s/clusterResults.ps'%outputdir
	outputfiles.clusterScaledInputData = '%s/clusterScaledInputData.txt'%scratchdir
        outputfiles.epitheliumFlatFile = '%s/%s.grid.data' %(outputdir, model.systemName)
        outputfiles.epitheliumPlot = '%s/%s.grid.ps' %(outputdir, model.systemName)
        outputfiles.exportFile = '%s/%s.xml' %(outputdir, model.systemName)
        outputfiles.identifiabilityCorrelationPlot = '%s/%s.correlation.ps' %(outputdir, model.systemName)
        outputfiles.integrationInput = '%s/integrationFunction.py' %(scratchdir)
        outputfiles.gaStatistics = '%s/%s.st' %(outputdir, model.systemName)
	outputfiles.gnuplot3DPlotter = '%s/%s.3D.gnu' %(scratchdir, model.systemName)
	outputfiles.gnuplot3DData = '%s/%s.3D.data' %(scratchdir, model.systemName)
	outputfiles.gnuplotErrorMessagesFile = '%s/%s.gnu.error.message.txt' %(scratchdir, model.systemName)
        outputfiles.gnuplotInputFile = '%s/%s.gnu' %(scratchdir, model.systemName)
        outputfiles.gnuplotVelocitiesInputFile = '%s/%s.veloc.gnu' %(scratchdir, model.systemName)
        outputfiles.identifiabilityCriteria = '%s/%s.criteria.txt' %(outputdir, model.systemName)
        outputfiles.identifiabilityFIM = '%s/%s.FIM.txt' %(outputdir, model.systemName)
        outputfiles.identifiabilityCOV = '%s/%s.COV.txt' %(outputdir, model.systemName)
        outputfiles.identifiabilityCorrelation = '%s/%s.correlation.txt' %(outputdir, model.systemName)
        outputfiles.latexFile = '%s/%s.tex' %(outputdir, model.systemName)
	outputfiles.logFile = '%s/%s.log'%(scratchdir, model.systemName)
        outputfiles.octaveInputFile = '%s/%s.%s.oc' %(scratchdir, model.systemName, parallel.currentProcessor())
        outputfiles.octaveIntegrationData = '%s/%s.%s.data' %(scratchdir, model.systemName, parallel.currentProcessor())
        outputfiles.octaveIntegrationTime = '%s/%s.%s.time' %(scratchdir, model.systemName, parallel.currentProcessor())
        outputfiles.oedTimePlot = '%s/%s.oed.' %(outputdir, model.systemName)
        outputfiles.oedResults = '%s/%s.oed.txt' %(outputdir, model.systemName)
        outputfiles.openModelicaModel = '%s/%s.mo' %(scratchdir, model.systemName)
        outputfiles.openModelicaInput = '%s/%s.mos' %(scratchdir, model.systemName)
        outputfiles.openModelicaOutput = '%s/%s_res.plt' %(scratchdir, model.systemName) #/ This name is mandatory
        outputfiles.PGAPackInput = '%s/pgapack.input.txt'%(scratchdir)
        outputfiles.pythonIntegrationFile = '%s/%s.%s.integ.py' %(scratchdir, model.systemName, parallel.currentProcessor())
        outputfiles.runBackup = '%s/running.bkp' %outputdir
        outputfiles.sampleResults = '%s/%s.sample.txt' %(outputdir, model.systemName)
        outputfiles.sampleNonAcceptedResults = '%s/%s.sample.nonAccepted.txt' %(scratchdir, model.systemName)
        outputfiles.sensitivitiesOSData = '%s/%s.sens.global.out' %(scratchdir, model.systemName)
        outputfiles.sensitivitiesOSPlotCommands = '%s/%s.sens.global.gnu' %(scratchdir, model.systemName)
        outputfiles.sensitivitiesOS = '%s/%s.sens.global.txt' %(outputdir, model.systemName)
        outputfiles.sensitivitiesRSMatrix = '%s/%s.sens.relative.txt' %(outputdir, model.systemName)
        outputfiles.sensitivitiesRSMatrixPlot = '%s/%s.sens.relative.ps' %(outputdir, model.systemName)
        outputfiles.simulationResults = '%s/%s.%s.out' %(outputdir, model.systemName, parallel.currentProcessor())
        outputfiles.simulationResultsCSV = '%s/%s.%s.csv' %(outputdir, model.systemName, parallel.currentProcessor())
        outputfiles.solutionDirectory = '%s/%s' %(outputdir, model.systemName)
        outputfiles.summaryModel = '%s/%s.description.txt' %(outputdir, model.systemName)
        outputfiles.statisticalDynamicsDirectory = '%s/reconstructedDynamics' %(outputdir)
        outputfiles.velocitiesFile = '%s/%s.veloc.%s.out' %(outputdir, model.systemName, parallel.currentProcessor())
        outputfiles.xppIntegrationFile = '%s/%s.%s.xpp' %(scratchdir, model.systemName, parallel.currentProcessor())
        outputfiles.xppOutputFile = '%s/%s.XPP.data' %(scratchdir, model.systemName) 
        outputfiles.xppOutputDetails = '%s/%s.XPP.log' %(scratchdir, model.systemName)
	#/ setting the output format of the files
	if metamodel.figureFormat == 'ps':
	    outputfiles.timePlot = '%s/%s.ps' %(outputdir, model.systemName)
	    outputfiles.velocitiesPlot = '%s/%s.veloc.ps' %(outputdir, model.systemName)
	    outputfiles.sensitivitiesOSPlotFile = '%s/%s.sens.global.timeCourse.ps' %(outputdir, model.systemName)
	    outputfiles.samplePlot = '%s/%s.sample.ps' %(outputdir, model.systemName)
	elif metamodel.figureFormat == 'png':
	    outputfiles.timePlot = '%s/%s.png' %(outputdir, model.systemName)
	    outputfiles.velocitiesPlot = '%s/%s.veloc.png' %(outputdir, model.systemName)
	    outputfiles.sensitivitiesOSPlotFile = '%s/%s.sens.global.timeCourse.png' %(outputdir, model.systemName)
	    outputfiles.samplePlot = '%s/%s.sample.png' %(outputdir, model.systemName)
	else:
	    raise errorMessages.ClassCentralException, 'output figure format not recognised.'

        return outputfiles
    
    outputfiles = fileInitializator(model, metamodel) #/ Initializing the output files
    #/ PARALLEL: INSERT
    parallelParameterEstimation = False
    if (parallel.mainProcessor()):
    #/ PARALLEL: IDENTATION
        if metamodel.runningType == 'simulation':
            print 'Simulation running ...'
            simulator.central(metamodel, model, outputfiles)
        elif metamodel.runningType == 'parameterEstimation':
            #/ PARALLEL: FLAG ADDED
            parallelParameterEstimation = True
        elif metamodel.runningType == 'exporting':
            print 'Exporting ...'
            exporter.central(metamodel, model, outputfiles)
        elif metamodel.runningType == 'dynamicsReconstruction':
            print 'Reconstructing the dynamics from the solutions ...'
            dynamicsReconstructer.central(metamodel, model, outputfiles)
        elif metamodel.runningType == 'scoreSurface':
            print 'Creating score surfaces ...'
            surface.central(metamodel, model, outputfiles)
        elif metamodel.runningType == 'sensitivityAnalysis':
            print 'Performing the sensitivity analysis ...'
            sensitivityAnalyzer.central(model, metamodel, outputfiles)
        elif metamodel.runningType == 'identifiabilityAnalysis':
            print 'Performing the identifiability analysis ...'
            identifiabilityAnalyzer.central(model, metamodel, outputfiles)
        elif metamodel.runningType == 'calculateFunction':
            print 'Calculating the value of the fitness function ...'
            fitnessFunctionEvaluator.central(model, metamodel, outputfiles)            
        elif metamodel.runningType == 'optimalExperimentalDesign':
            print 'Performing optimal experimental design analysis ...'
            optimalExperimentalDesign.central(model, metamodel, outputfiles)
        elif metamodel.runningType == 'sampleSurface':
            print 'Running Monte Carlo sampling of the fitness function surface ...'
            sampler.central(model, metamodel, outputfiles)
	elif metamodel.runningType == 'clustering':
	    print 'Clustering ...'
	    cluster.main(metamodel, outputfiles)
    #/ PARALLEL: CODE ADDED
    if parallel.mainProcessor():
        parallel.sendAll(parallelParameterEstimation, 1)
        
    parallelParameterEstimation = parallel.receiveAny(1)
    if parallelParameterEstimation == True:
        if parallel.mainProcessor():
            print 'Running optimisation ...'
        solutions = 0
        optimiser.central(metamodel, model, outputfiles, solutions)
        #/ Repeat the calculation in order to obtain several samples of the result
        while metamodel.statisticalOptimisation == True:
            solutions = solutions + 1
            optimiser.central(metamodel, model, outputfiles, solutions)
            if solutions >= metamodel.statisticalOptimisationSolutions:
                metamodel.statisticalOptimisation = False

    return None
    
def sbmlReader(metamodel):

    ''' 
    This function reads a model in SBML format and create the model object. 
    '''

    model = sbmlWorker.ClassModelSBML() #/ To create the object
    model.parser(metamodel) #/ Building the model object with the model infomation
    model.checkCompatibilities(metamodel) #/ Checking if the model can run with ByoDyn

    return model

def tagsReader(metamodel):

    '''
    This function reads the model in "tags" format and creates the model object. 
    '''

    model = tagParser.ClassModelTags()
    model = model.readInput(metamodel)

    return model

##############################
#       Main Program         #
##############################

def main(runnerFile):

    ''' 
    This is the main function of the entire program.
    It builds an object called "metamodel" with the options to run ByoDyn, another object called "model" with the structure of the system and finally the function calls the function "runner" to execute the instructions. 
    '''
    
    metamodel = None
    model = None
    if (parallel.mainProcessor()):
        metamodel = optionReader(runnerFile) #/ Reading the input runner with the details
        model = modelReader(metamodel) #/ Convert model in internal model object
        parallel.sendAll(metamodel, 0)
        parallel.sendAll(model, 1)
    metamodel = parallel.receiveAny(0)
    model = parallel.receiveAny(1)
    #/ do not run the program if we call PGAPack
    if metamodel.optimisationMethod == 'parallelGA':
	return model
    else:
	runner(metamodel, model) #/ Running the program
#/ necessary for profiling
if __name__ =='__main__':
    main()
