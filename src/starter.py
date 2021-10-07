#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author:  Adrian L. Garcia-Lomana and Alex Gomez-Garrido
#
#  Created: 2004-11-05 by Adrian L. Garcia-Lomana
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

# $Id: starter.py,v 4.33 2008/12/15 10:26:45 alexgomez Exp $

## \file
# This module is the parser for the running options of ByoDyn.

import getopt, os, shutil, time, sys
import errorMessages

class ClassMetaModel:

      '''
      Class for the running options of Byodyn.
      '''

      def __init__(self):

            '''
            The constructor.
            '''
            
            self.bandwidth = None
            self.clusterInputData = None
            self.confidenceIntervals = False
            self.figureFormat = 'ps' #/ by default
            self.fixedParameters = []
            self.gaPopulation = 10
            self.gaTranslocationRate = 0.8
            self.gaMutationRate = 0.3
            self.graphics = True
            self.onlyLastState = False
            self.hybridScore = False
            self.hybridPath = None
            self.identifiabilityCriteria = ['MA', 'E', 'D', 'ME']
            self.identifiabilityOutputs = ['FIM', 'COV', 'Correlation']
            self.initialConcentrationsToVary = []
            self.initialConditionsSolutions = None
            self.integrationMethod = None
            self.integrationOption = None
            self.integrationTolerance = []
            self.stochasticOption = None
            self.stochasticMethod = None
            self.stochasticRuns = None
            self.RKCustomButcherTableau = None
            self.RKDampedParameters = None
            self.lastIterationScore = 0.
            self.modelFile = None
            self.modelFormat = None
            self.optimisationMethod = None
            self.optionalOutputFormat = []
            self.oedTask = None
            self.oedTargetSpecies = None
            self.oedResolution = 0
            self.parameters = {}
            self.parametersToVary = []
            self.plottingKeys = None
            self.runWithBackup = False
            self.runningType = None
            self.sampleMethod = None
            self.sampleSize = 0
            self.sensitivityParameters = []
            self.separatedGraphs = False
            self.showingPlot = False #/ by default
            self.simulationNumber = 0
            self.simulationTime = 0.0
            self.simulationTimeStep = 0.0
            self.solutionsDirectory = None
            self.statisticalOptimisation = False
            self.statisticalOptimisationSolutions = 0
            self.stopper = []
            self.strictCheckSBML = False
            self.surfaceParameters = []
            self.surfaceResolution = 0
            self.target = {}
            self.threshold = -1
            self.velocities = False

            return None

      
      def metamodelParser(self, optionFile):

            '''
            This method gets the information of the running obtions 
            and it sets the object.
            '''
            
            def emptyLineChecker(fields):

                  '''
                  This function checks that there are no empty lines on the input tag file.            
                  '''
                  
                  if len(fields) == 0:
                        raise errorMessages.ClassStarterException, 'the runner file contains empty lines, please remove them.'

                  return None
            
            for line in optionFile:
                  fields = line.split()
                  emptyLineChecker(fields)
                  if fields[0] == 'integrationMethod':
                        self.integrationOption = fields[1]
                        if self.integrationOption != 'python' and self.integrationOption != 'octave' and self.integrationOption != 'openModelica' and self.integrationOption != 'matlab' and self.integrationOption != 'xpp' and self.integrationOption != 'rungeKutta' and self.integrationOption != 'euler' and self.integrationOption != 'automatic':
                              raise errorMessages.ClassStarterException, 'invalid setting of the integration method at the runner file. Make sure the variable "integrationMethod" takes as a first argument either "automatic", "python", "octave", "openModelica", "xpp", or "matlab".'
                        self.integrationMethod = fields[2]
                        #/ checking the correctness of the integration method
                        correctPossibilitiesForIntegrationMethod = ['default', 'adams', 'non-stiff', 'bdf', 'stiff']
                        if correctPossibilitiesForIntegrationMethod.count(self.integrationMethod) != 1:
                            raise errorMessages.ClassStarterException, 'wrong integration method at the runner file. Make sure the variable "integrationMethod" takes as a second argument one of these posibilities: default, adams, non-stiff, bdf, stiff'
                  elif fields[0] == 'integrationTolerance':
                      self.integrationTolerance = [float(fields[1]), float(fields[2])]
                  elif fields[0] == 'stochasticMethod':
                        self.stochasticOption = fields[1]
                        self.stochasticMethod = fields[2]
                        #/ checking self.stochasticOption and self.stochasticMethod
                        if self.stochasticOption == 'ssa':
                              if self.stochasticMethod != 'gillespie':
                                    raise errorMessages.ClassStarterException, 'invalid setting of the stochastic option at the runner file. Make sure that the second argument for the "stochasticMethod" variable takes as second argument "gillespie". No other option is currently available for exact stochastic simulations.'
                        elif self.stochasticOption == 'tau-leap':
                              if self.stochasticMethod != 'default':
                                    raise errorMessages.ClassStarterException, 'invalid setting of the stochastic option at the runner file. Make sure that the second argument for the "stochasticMethod" variable takes as second argument "default". No other option is currently available for the tau-leap algorithm.'
                        else:
                              raise errorMessages.ClassStarterException, 'invalid setting of the stochastic method at the runner file. Make sure that the variable "stochasticMethod" takes as a first argument either "ssa" or "tau-leap.'
                        self.stochasticRuns = int(fields[3])
                        if self.stochasticRuns <= 0:
                            raise errorMessages.ClassStarterException, 'wrong stochastic runs number. "stochasticRuns" should be a positive integer'
                  elif fields[0] == 'RKCustomButcherTableau':
                        self.RKCustomButcherTableau = map(float, fields[1:])
                        #/ this should be "s a11 a12 .. a1s a21 ..a2s ..as1 .. ass b1 .. bs" # so the number of fields is 1 + s + s**2 = s*(s+1) +1
                        #/ we check that the tableau sizes are consitent
                        if (len(self.RKCustomButcherTableau) != self.RKCustomButcherTableau[0]**2 + self.RKCustomButcherTableau[0] + 1 ):
                           raise errorMessages.ClassStarterException, 'invalid setting of the "RKCustomButcherTableau". Refer to the documetation for how to use this option.'                  
                  elif fields[0] == 'RKDampedParameters':
                        self.RKDampedParameters = map(float, fields[1:])
                        if (len(self.RKDampedParameters) != 2):
                           raise errorMessages.ClassStarterException, 'invalid setting of the "RKDampedParameters". This option takes 2 arguments, arg1 is the number of stages of the method and arg2 is the damping factor.'
                  elif fields[0] == 'modelFile':
                        self.modelFile = fields[1]
                  elif fields[0] == 'modelFormat':
                        self.modelFormat = fields[1]
                        if self.modelFormat != 'SBML' and self.modelFormat != 'tags':
                              raise errorMessages.ClassStarterException, 'invalid setting of the "modelFormat". Allowed values are "SBML" or "tags".'
                  elif fields[0] == 'plotKeys':
                        if fields[1] == 'YES':
                              self.plottingKeys = True
                        elif fields[1] == 'NO':
                              self.plottingKeys = False
                        else:
                              raise errorMessages.ClassStarterException, 'invalid setting of the "plottingKeys" variable at the runner file. The possible values are "YES" or "NO".'
                  elif fields[0] == 'runningType':
                        self.runningType = fields[1]
                        if self.runningType != 'simulation' and self.runningType != 'parameterEstimation' and self.runningType != 'exporting' and self.runningType != 'dynamicsReconstruction' and self.runningType != 'scoreSurface' and self.runningType != 'sensitivityAnalysis' and self.runningType != 'calculateFunction' and self.runningType != 'identifiabilityAnalysis' and self.runningType != 'optimalExperimentalDesign' and self.runningType != 'sampleSurface' and self.runningType != 'clustering':
                              raise errorMessages.ClassStarterException,  'check for "runningType" variable. Valid options are "simulation", "parameterEstimation", "exporting", "dynamicsReconstruction", "scoreSurface", "calculateFunction", "identifiabilityAnalysis", "optimalExperimentalDesign", "sampleSurface", "sensitivityAnalysis" or "clustering".'
                        if self.runningType == 'simulation':
                              if fields[2] == 'velocity':
                                    self.velocities = True
                              elif fields[2] == 'noVelocity':
                                    self.velocities = False
                              else:
                                    raise errorMessages.ClassStarterException, 'set as second argument for "runningType" either "velocity" or "noVelocity".'
                        if self.runningType == 'parameterEstimation':
                              if fields[2] == 'geneticAlgorithm':
                                    self.optimisationMethod = 'geneticAlgorithm'
                              elif fields[2] == 'randomSearch':
                                    self.optimisationMethod = 'randomSearch'
                              elif fields[2] == 'hybridTwoPhases':
                                    self.optimisationMethod = 'hybridTwoPhases'
                              elif fields[2] == 'localSearch':
                                    self.optimisationMethod = 'localSearch'
                              elif fields[2] == 'hybridOnePhase':
                                    self.optimisationMethod = 'hybridOnePhase'
                              elif fields[2] == 'parallelGA':
                                    self.optimisationMethod = 'parallelGA'
                              else:
                                    raise errorMessages.ClassStarterException, 'invalid setting of the second argument for the "runningType" variable. Acceptable values are "geneticAlgorithm", "parallelGA", "randomSearch", "localSearch", "hybridTwoPhases" or "hybridOnePhase".'
                        if self.runningType == 'scoreSurface':
                              self.surfaceResolution = int(fields[2])
                              self.surfaceParameters = fields[3:]
                        if self.runningType == 'optimalExperimentalDesign':
                            self.oedTask = fields[2]
                        if self.runningType == 'sampleSurface':
                            if len(fields) == 4:                               
                                if fields[2] == 'MonteCarlo':
                                    self.sampleMethod = fields[2]    
                                else:
                                    raise errorMessages.ClassStarterException, 'invalid setting of the argument for the "sampleSurface" method. Acceptable value is only "MonteCarlo".'                                
                                self.sampleSize = int(fields[3])
                            else:
                                raise errorMessages.ClassStarterException, 'invalid setting of the arguments for the "sampleSurface" method. You need indicate the method and the sample size. Ex: runningType    sampleSurface    [method]    [integer]'
                  elif fields[0] == 'time&timestep':
                        self.simulationTime = float(fields[1])
                        self.simulationTimeStep = float(fields[2])
                  elif fields[0] == 'statisticalOptimisation':
                        statisticalOptimisationOptions = fields[1:]
                        if statisticalOptimisationOptions[0] == 'YES':
                              self.statisticalOptimisation = True
                        elif statisticalOptimisationOptions[0] == 'NO':
                              self.statisticalOptimisation = False
                        else:
                              raise errorMessages.ClassStarterException, 'invalid setting of the optimisation values at the runner file. The variable "statisticalOptimisation" can only take as first argument "YES" or "NO".'
                        if self.statisticalOptimisation == True and len(statisticalOptimisationOptions) < 2:
                              raise errorMessages.ClassStarterException, 'if you want to run a statistical optimisation, you have to set the number of solutions you want.'
                        self.statisticalOptimisationSolutions = int(statisticalOptimisationOptions[1])
                  elif fields[0] == 'gaOptions':
                        gaOptions = fields[1:]
                        for option in gaOptions:
                              if option.split('/')[0] == 'populationSize':
                                    self.gaPopulation = int(option.split('/')[1])
                              elif option.split('/')[0] == 'translocationRate':
                                    self.gaTranslocationRate = float(option.split('/')[1])
                              elif option.split('/')[0] == 'mutationRate':
                                    self.gaMutationRate = float(option.split('/')[1])
                              else:
                                    raise errorMessages.ClassStarterException, 'invalid setting of the variable "gaOptions". The possible arguments are "populationSize", "translocationRate" or "muationRate" with the following strucature: "argument/value".'
                  elif fields[0] == 'stopper':
                        self.stopper = fields[1:]
                        if self.stopper[0] != 'iteration' and self.stopper[0] != 'score' and self.stopper[0] != 'numberOfSimulations':
                              raise errorMessages.ClassStarterException, 'invalid setting of the first argument of the variable "stopper". The possible values are "iteration", "score" or "numberOfSimulations".'
                        if self.stopper[0] == 'score' and len(self.stopper) < 2:
                              raise errorMessages.ClassStarterException, 'Fatal error: the number of arguments of the variable "stopper" is two for a score based stopping. %s given. Please provide the threshold of the stopper.' % len(self.stopper)
                  elif fields[0] == 'threshold':
                        self.threshold = float(fields[1])
                  elif fields[0] == 'parameter':
                        specifiedParameters = fields[1:]
                        for specifiedParameter in specifiedParameters:
                              unbound = specifiedParameter.split('=')
                              name = unbound[0] 
                              value = float(unbound[1])
                              self.fixedParameters.append([name,value]) #/ necessary for the GA at the gaInitialisePopulation function
                              if self.parameters.has_key(name) == True:
                                    print 'Warning: In the runner file you have the parameter %s specified more than one time. We only use the first value' %name
                              else:
                                    self.parameters[name] = value
                  elif fields[0] == 'parametersToVary':
                        self.parametersToVary = fields[1:]
                  elif fields[0] == 'sensitivityParameters':
                        self.sensitivityParameters = fields[1:]
                  elif fields[0] == 'target':
                        if self.target.has_key(float(fields[1])):
                              for element in fields[2:]:
                                    self.target[float(fields[1])].append(element)
                        else:
                              self.target[float(fields[1])] = fields[2:]
                  elif fields[0] == 'solutionsDirectory':
                        self.solutionsDirectory = fields[1]
                  elif fields[0] == 'runWithBackup':
                        self.runWithBackup = True
                  elif fields[0] == 'withoutGraphics':
                        self.graphics = False
                  elif fields[0] == 'onlyLastState':
                        self.onlyLastState = True
                  elif fields[0] == 'identifiabilityCriteria':
                        #/ checking elements
                        if len(fields[1:]) > 4:
                             raise errorMessages.ClassStarterException, 'no more than four different identifiability criteria are available: MA, D, E and ME. Remove the repeats.'
                        for element in fields[1:]:
                              if self.identifiabilityCriteria.count(element) == 0:
                                    raise errorMessages.ClassStarterException, 'no recognised identifiability criteria. Possible options are MA, D, E or ME.'
                        #/ incorporating the values
                        self.identifiabilityCriteria = fields[1:]
                  elif fields[0] == 'identifiabilityOutputs':
                        self.identifiabilityOutputs = fields[1:]
                  elif fields[0] == 'showingPlot':
                        self.showingPlot = True
                  elif fields[0] == 'targetSpecies':
                        self.oedTargetSpecies = fields[1:]
                  elif fields[0] == 'OEDResolution':
                        self.oedResolution = int(fields[1])
                  elif fields[0] == 'optionalOutputFormat':
                        self.optionalOutputFormat = fields[1:]
                  elif fields[0] == 'checkConsistencySBML':
                        self.strictCheckSBML = True
                  elif fields[0] == 'figureFormat':
                        self.figureFormat = fields[1]
                        if self.figureFormat != 'png' and self.figureFormat != 'ps':
                              raise errorMessages.ClassStarterException, '"figureFormat" variable can only take "ps" (as default) or "png" arguments.'
                  elif fields[0] == 'clusteringDataFile':
                        self.clusterInputData = fields[1]
                  elif fields[0] == 'bandwidth':
                        self.bandwidth = fields[1]
                  elif fields[0] == 'separatedGraphs':
                        self.separatedGraphs = True
                  elif fields[0] == 'initialConcentrationsToVary':
                        self.initialConcentrationsToVary = fields[1:]
                  elif fields[0] == 'initialConditionsSolutions':
                        self.initialConditionsSolutions = fields[1]
                  else:
                        if fields[0][0] != '#': #/ this is a commentary
                              print fields
                              raise errorMessages.ClassStarterException, 'tag "%s" is not recognised.'%fields[0]

            return self

def central(runnerFile):

      '''
      This function directs the reading of the running options of ByoDyn.
      Then it checks for incompatibilities.
      '''

      #/ 1.- Initializing objects
      metamodel = ClassMetaModel()
      #/ 2.- Obtainig the values for the running options and the metamodel
      workingDir = os.environ.get('WORKING_DIR')
      optionFile = open(runnerFile, 'r')
      metamodel = metamodel.metamodelParser(optionFile)
      optionFile.close()
      #/ 2.- looking for incompatibility
      incompatibilityChecker(metamodel)
      #/ 3.- compiling some fortran modules necesary for the hybrid optimisation
      if metamodel.optimisationMethod == 'hybridTwoPhases' or metamodel.optimisationMethod == 'localSearch' or metamodel.optimisationMethod == 'hybridOnePhase':
            metamodel = fortranModulesCreator(metamodel)
      if metamodel.runningType == 'parameterEstimation':
            if metamodel.stopper != []:
                  if metamodel.stopper[0] == 'score' and metamodel.threshold == -1:
                        metamodel.threshold = float(metamodel.stopper[1])

      return metamodel

def incompatibilityChecker(metamodel):

      '''
      This function checks that the options of the metamodel do not contain incompatibilities.
      '''

      if metamodel.integrationOption == 'octave' and metamodel.velocities == True:
            raise errorMessages.ClassStarterException, 'this version of ByoDyn does not support the analysis of the reaction velocities using Octave. Sorry for the inconvenience.'
      if metamodel.integrationOption == 'openModelica' and metamodel.velocities == True:
            raise errorMessages.ClassStarterException, 'this version of ByoDyn does not support the analysis of the reaction velocities using openModelica. Sorry for the inconvenience.'
      if metamodel.integrationOption == 'xpp' and metamodel.velocities == True:
            raise errorMessages.ClassStarterException, 'this version of ByoDyn does not support the analysis of the reaction velocities using xpp. Sorry for the inconvenience.'      
      if metamodel.runningType == 'parameterEstimation' and metamodel.stopper == []:
            if metamodel.optimisationMethod != 'localSearch':
                  raise errorMessages.ClassStarterException, 'you have to set the variable "stopper" for running a parameter estimationA.'
      if metamodel.runningType == 'parameterEstimation' and metamodel.parametersToVary == []:
            targetForTimeZero = False
            for element in metamodel.target.keys():
                  if element == 0.0:
                        targetForTimeZero = True
            if targetForTimeZero == False:
                  raise errorMessages.ClassStarterException, 'you have to set the variable "parametersToVary" for running a parameter estimation.'
      if metamodel.runningType == 'parameterEstimation' and metamodel.target == {}:
            raise errorMessages.ClassStarterException, 'you have to set the variable "target" for running a parameter estimation.'
      for parameterToVary in metamodel.parametersToVary:
            if len(parameterToVary.split('/')) < 4:
                  raise errorMessages.ClassStarterException, 'the number of arguments of the parameters to vary must be four: "parameter/lowerValue/upperValue/scale".'
            if parameterToVary.split('/')[3] != 'log' and parameterToVary.split('/')[3] != 'lin':
                  raise errorMessages.ClassStarterException, 'the fourth argument of the parameters to vary must be either "log" or "lin".'
      if metamodel.runningType == 'dynamicsReconstruction' and metamodel.solutionsDirectory == None:
            raise errorMessages.ClassStarterException, 'if you want to reconstruct the dynamics of the solutions, you have to set the variable "solutionsDirectory".'
      #/ 1.- checking that the starting points of the optimisation has their corresponding ranges
      for parameter in metamodel.parameters:
            for parameterToVary in metamodel.parametersToVary:
                  parameterToVaryName = parameterToVary.split('/')[0]
                  parameterToVaryRange = [float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2])]
                  if parameter == parameterToVaryName:
                        startingPointValue = float(metamodel.parameters[parameter])
                        if parameterToVaryRange[0] > startingPointValue or parameterToVaryRange[1] < startingPointValue:
                              raise errorMessages.ClassStarterException, 'the parameter value of %s you are suggesting is out of the parameter range.'%parameter
      #/ 2.- checking that the optionalOutputFormat value is csv
      if metamodel.optionalOutputFormat != []:
            if metamodel.optionalOutputFormat != ['csv']:
                  raise errorMessages.ClassStarterException, 'the only alternative output format is comma separated value. Please check that the value of "optionalOutputFormat" is "csv".'
      #/ 3.- detecting that using sensitivity, the integrators Euler and Runge-Kutta are not valid
      if metamodel.runningType == 'sensitivityAnalysis':
            if metamodel.integrationOption == 'euler' or metamodel.integrationOption == 'rungeKutta':
                  raise errorMessages.ClassStarterException, 'the sensitivity analysis cannot be handled using "euler" or "rungeKutta" integration methods.'
      #/ 4.- detecting that there are no target if we run sensitivity analysis
      if metamodel.runningType == 'sensitivityAnalysis':
            if metamodel.target != {}:
                  raise errorMessages.ClassStarterException, 'please remove the forbidden tag "target" when running the sensitivity analysis.'
      #/ 5.- imposing the target for identifiability analysis
      if metamodel.runningType == 'identifiabilityAnalysis' and metamodel.target == {}:
            raise errorMessages.ClassStarterException, 'you have to set the variable "target" for running an identifiability analysis.'
      #/ 6.- detecting repeated targets
      for time in metamodel.target.keys():
            constrainedNodes = []
            for constrain in metamodel.target[time]:
                constrainedNodes.append(constrain.split('/')[0]+constrain.split('/')[1]) #/ selecting the node and the cell position
            for node in constrainedNodes:
                  if constrainedNodes.count(node) != 1:
                        raise errorMessages.ClassStarterException, 'node %s is constrained several times at time %s. Please specify only one constrain.' %(node,time)
      #/ 7.- detecting convert command
      if metamodel.figureFormat == 'png':
            scratchdir = os.environ.get('BYODYN_SCRATCH_DIR')
            messagesFile = scratchdir + '/convertMessagesFile'
            os.system('convert > %s'%messagesFile)
            f = open(messagesFile, 'r')
            lines = f.readlines()
            f.close()
            if len(lines) == 0:
                  print 'WARNING:"convert" command was not found which is required to convert postscript figures on png files.\nOutput format will be postscript, as default.'
                  metamodel.figureFormat == 'ps'
      #/ 8.- sampling iterations
      if metamodel.runningType == 'sampleSurface':
          if metamodel.sampleSize <= 0:
              raise errorMessages.ClassStarterException, 'you have to indicate a size of sample higher to 0.'       
      #/ 9.- detecting the presence of octave for the clustering
      scratchdir = os.environ.get('BYODYN_SCRATCH_DIR')
      if metamodel.runningType == 'clustering':
            octaveOutputMessagesFile = scratchdir + '/' + 'octaveOutputMessages'
            octaveErrorMessagesFile = scratchdir + '/' + 'octaveErrorMessages'
            os.system('octave -version 1> %s 2> %s'%(octaveOutputMessagesFile, octaveErrorMessagesFile))
            f = open(octaveOutputMessagesFile, 'r')
            lines = f.readlines()
            f.close()
            if len(lines) > 5: 
                  vector = lines[0].split()
                  for i in range(len(vector)):
                        if vector[i] == 'version':
                              version = vector[i+1]
                  versionNumber = version.split('.')
                  if int(versionNumber[0]) < 3:
                        raise errorMessages.ClassStarterException, 'Octave version 3.0.1 or higher is required for clustering. Update Octave software, please.'                        
            else:
                  raise errorMessages.ClassStarterException, 'Octave is required for clustering and it seems that its executable cannot be found. Please make sure you have it accessible from the terminal.'
      #/ 10.- detecting the presence of initialConcentrationsToVary if some target refers to zero
      if metamodel.runningType != 'dynamicsReconstruction':
            targetedICs = []
            setICs = []
            if metamodel.target.keys().count(0.0) == 1:
                  for element in metamodel.target[0.0]:
                        targetedICs.append(element.split('/')[0])
            for element in metamodel.initialConcentrationsToVary:
                  setICs.append(element.split('/')[0])
            for element in targetedICs:
                  if setICs.count(element) != 1:
                        raise errorMessages.ClassStarterException, 'not all nodes targeted for time zero have been set by the "initialConcentrationsToVary" field. Please use it to define the starting initial conditions of the model.'
      #/ 11.- checking that time step is zero only if gillespie stochastic simulation is used
      if metamodel.simulationTimeStep == 0:
            if metamodel.stochasticOption != 'ssa' and metamodel.stochasticMethod != 'gillespie' and metamodel.runningType != 'exporting' and metamodel.runningType != 'clustering':
                  raise errorMessages.ClassStarterException, 'second argument for "time&timestep" should not be zero. Please change it.'
      #/ 12.- checking that integrationTolerance is only used with xpp integrator
      if len(metamodel.integrationTolerance) != 0 and metamodel.integrationOption != 'xpp':
          raise errorMessages.ClassStarterException, 'The tolerance setting is only available for XPP integrator. Remove the integrationTolerance option' 
 
                          
      return None

def fortranModulesCreator(metamodel):
      
      '''
      This function builds a module that can be called from Python from Fortran.
      It must be compiled now because there are some functions that are static and some arguments have to be declared before compilation.
      '''
       
      print 'Compiling Fortran code ...'
      #/ 1.- determining the value of the control variables
      rdim = 0
      for key in metamodel.target:
            rdim = rdim + len(metamodel.target[key])
      xdim = len(metamodel.parametersToVary) + len(metamodel.initialConcentrationsToVary)
      ddliv = 82 + 4*xdim
      ddlv = 105 + xdim*(rdim + 2*xdim + 21) + 2*rdim
      #/ 2.- writing the file of the control variables
      scratchdir = os.environ.get('BYODYN_SCRATCH_DIR')
      fortranModuleDir = os.environ.get('BYODYN_PATH') + '/lib/localSearch'
      metamodel.hybridPath = scratchdir + '/localSearch' + str(int(time.time()*100))
      os.mkdir(metamodel.hybridPath)
      currentDir = os.getcwd()
      controlFile = metamodel.hybridPath + '/auxVar.dc'
      f = open(controlFile, 'w')
      f.write('\tinteger RRDIM, XXDIM, ddliv, ddlv\n\n')
      f.write('\tparameter (RRDIM = %s)\n'%rdim)
      f.write('\tparameter (XXDIM = %s)\n'%xdim)
      f.write('\tparameter (ddliv = %s)\n'%ddliv)
      f.write('\tparameter (ddlv = %s)\n'%ddlv)
      f.close()
      #/ 3.- copying the files
      shutil.copy(fortranModuleDir + '/localSearch.f', metamodel.hybridPath  + '/localSearch.f')
      shutil.copy(fortranModuleDir + '/makefile', metamodel.hybridPath + '/makefile')
      #/ 4.- compiling the module
      os.chdir(metamodel.hybridPath)
      libportPath = fortranModuleDir + '/PORT'
      os.system('make libport=%s > %s/compilationMessages'%(libportPath, metamodel.hybridPath))
      os.chdir(currentDir)

      return metamodel
