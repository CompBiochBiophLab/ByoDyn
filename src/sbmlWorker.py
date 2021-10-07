#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana, Miguel A. Hernandez and Alex Gomez-Garrido
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

# $Id: sbmlWorker.py,v 4.25 2008/12/09 09:35:49 alexgomez Exp $

## \file
# This module is dedicated to the parsing of the models in SBML.
# All the information that ByoDyn requires is the object called "model".

import string, re, sys, scipy, os #/ python modules
import cPickle as pickle 
import errorMessages, formulas #/ ByoDyn modules
try:
    import libsbml 
except ImportError:
     raise errorMessages.ClassSBMLWorkerException, 'Error in libSBML library import.\n You needed to install libSBML with Python bindings correctly.\n Please visit "http://sbml.org/software/libsbml" for instructions.'

class ClassRules:
    
    '''
    Class for SBML rules.
    '''	
    
    def __init__(self):

	'''
	This is the constructer.
	'''
    
        self.id = None
        self.type = None
        self.variable = None
        self.math = None
        self.mathAST = None

        return None

class ClassEvent:
 
    '''
    Class for SBML events.
    '''
    
    def __init__(self):

	'''
	This is the constructor.
	'''

        self.id = None
        self.trigger = None
        self.delay = None
        self.assignment = {}
        self.assignmentAST = {}
        
        return None

class ClassReaction:

    '''
    Class for the reactions.
    Required for the stochastic simulations.
    '''

    def __init__(self):

	'''
	This is the constructor.
	'''

        self.id = None
        self.kineticLaw = None
        self.products = []
        self.reactants = []

        return None


class ClassFunction:
 
    '''
    Class for SBML function definitions.
    '''
    
    def __init__(self):

        '''
        This is the constructor.
        '''

        self.id = None
        self.arguments = []
        self.output = None
        self.mathAST = None
        
        return None

class ClassModelSBML:
 
    '''
    This is the classs with the model description and methods necessary to work with the model.
    '''
    
    def __init__(self):
	
	'''
	The constructor.
	'''
    
        self.systemName = None
        self.xlength = 1
        self.ywidth = 1
        self.nodes = [] #/ the nodes that are in topology. The algebraic are missing but not the constant nodes
        self.algebraicNodes = {} #/ specifying some of the nodes. The nodes that are modified for algebraic, assignment rules, or events
        self.constantNodes = {} #/ specifying some of the nodes. The nodes that are constant
        self.topology = {}
        self.topologyAST = {}
        self.initialConditions = []
        self.parameters = {} #/ the constant parameters
        self.nonConstantParameters = {} #/ the parameters that are non constant. These are not part of the variable parameters
        self.compartments = {} #/ the constant compartments
        self.nonConstantCompartments = {}
        self.rules = []
        self.functions = [] #/ a list of functions
        self.events = [] #/ a list of events.
        self.reactions = [] #/ a dictionary with the rate equations of the all reactions
        self.delayFunctions = False
        self.delays = [] #/ a list of delays
        self.stoichiometryMatrix = []
        
        return None

    def parser(self, metamodel): #/ method of ClassModelSBML	
        
        '''
	This method reads a SBML model and transfrom this into  the class ClassModelSBML.
	'''

        modelFile = libsbml.readSBML(metamodel.modelFile) #/ Reading the sbml file in sbml file object
        if metamodel.strictCheckSBML == True:
            modelFile.checkConsistency()
        if modelFile.getNumErrors() != 0:
            for i in range(modelFile.getNumErrors()):
                print 'Message of the fatal error %s is: %s' %(i, modelFile.getError(i).getMessage())
            raise errorMessages.ClassSBMLWorkerException,  'fatal errors on the SBML file.'
        if modelFile.getLevel() != 2:
            modelFile.setLevel(2)
        sbmlModel = modelFile.getModel() #/ Getting a sbml model from sbml file object
        #/ 1.- Obtaining the name of the model
        self.systemName = sbmlModel.getId()
    	if self.systemName == "":
            print 'WARNING: impossible to obtain model name from the SBML file, the model name will be "withoutId".'
            self.systemName = 'withoutId'
        #/ 2.- Obtaining model list of compartments
        self.__readCompartments(sbmlModel)
        #/ 3.- Obtaining model nodes, constant nodes and initial conditions
        self.__readNodes(sbmlModel)
        #/ 4.- Obtaining model parameters
        self.__readParameters(sbmlModel, metamodel)
        #/ 5.- Obtaining list of Function definitions
        self.__readFunctionDefinitions(sbmlModel)
        #/ 6.- Obtaining rules. Rules provide a way to create constraints on variables for cases in which the constraints cannot be expressed using reactions
        self.__readRules(sbmlModel)
        #/ Obtaining initial Assignments
        self.__readInitialAssignments(sbmlModel)
        #/ 7.- Obtaining events
        self.__readEvents(sbmlModel)
        #/ 8. Verifying that the algebraicNodes are modified by rules or events.
        self.__checkAlgebraicNodes()
        #/ 9.- Obtaining model topology
        self.__readTopology(sbmlModel)
        #/ 10.Searching delay function in topology
        self.__searchDelayFunction(sbmlModel)
        #/ 11. Obtaining model constraints
        self.__readConstraints(sbmlModel)
        #/ 12.- Putting the rate rules inside the model topology and checking the initialConditions of the assignment rules variables:
        if len(self.rules) != 0:
            self.__checkAndPutRaterules()
        if len(self.nonConstantParameters) != 0 or self.nonConstantCompartments != 0:
            self.__checkNonConstant()
        #/ 13. Verifying that the model can be integrate 
        if self.canBeIntegrate() == False:
            print 'WARNING: The SBML model does not have kinetic information, it means kinetic laws, rules or events. Therefore, the model cannot be integrate.\n'
        return self

    def canBeIntegrate(self):
        '''
        This method searches if the model contains kinetic information. It means kinetic laws, rules or events.
        It returns True or False.
        '''

        if self.topology == {} and self.rules == [] and self.events == []:       
            return False
        else:
            return True
        
    def checkCompatibilities(self, metamodel):
 
        '''
	This method searches for possible problems between the model and the ByoDyn version.
	'''
        
        #/ checking if the model has Events
        if len(self.events) != 0:
        	if metamodel.integrationOption != 'openModelica' and metamodel.integrationOption != 'automatic':
        		raise errorMessages.ClassSBMLWorkerException,  'Error: The SBML model have events and Byodyn only work with events using openModelica. You have to change the integrationOption.'
        for rule in self.rules:
            if rule.type == 'Assignment' and metamodel.integrationOption == 'python':
                raise errorMessages.ClassSBMLWorkerException,  'Error: The SBML model have assignment rules, you cannot use scipy but Byodyn work with assignment rules using openModelica, Octave or XPP. You have to change the integrationOption.'
            if rule.type == 'Algebraic':
                if metamodel.integrationOption == 'python' or metamodel.integrationOption == 'xpp':
                    raise errorMessages.ClassSBMLWorkerException,  'Error: The SBML model have algebraic rules and Byodyn only work with assignment rules using openModelica or Octave. You have to change the integrationOption.'   
        if len(self.compartments.keys()) == 0 and len(self.nodes) == 0:
        	raise errorMessages.ClassSBMLWorkerException,  'Error: If the SBML model declares any species (nodes), the model must also declare at least one compartment).'
        for compartment in self.compartments:
        	if self.compartments[compartment]['spatialDimensions'] != 3:
        		print 'WARNING: The SBML model have no-three-dimensional compartments.'	        
        if self.delayFunctions == True and metamodel.integrationOption != 'xpp' and metamodel.integrationOption != 'automatic':
            raise errorMessages.ClassSBMLWorkerException,  'Error: The SBML model have delay function and Byodyn only work with delay function using XPPAUT. Please change the integrationOption.'
        
        return None

    def checkAssignmentRule(self, rule):
        
        '''
	This method checks for the assignment rules.
	'''

        value = formulas.solveFormula(self, rule.math) 
        #/ In the case of a species
        if self.nodes.count(rule.variable) != 0: 
            for definition in self.topology.keys():
                if definition.split('/') == rule.variable:
                    if self.constantNodes.has_key(rule.variable) == False:
                        raise errorMessages.ClassSBMLWorkerException, 'ERROR in the SBML file: the specie %s have assignmentRule and is also modified in reactions, but the species cannot have both model definitions. ' %rule.variable
                    else:
                        del self.topology[definition]
            self.algebraicNodes[rule.variable] = self.initialConditions.pop(self.nodes.index(rule.variable)) 
            self.nodes.remove(rule.variable)
            if value != self.algebraicNodes[rule.variable]:
                self.algebraicNodes[rule.variable] = value
       #/ In the case of algebraic nodes
        elif self.algebraicNodes.keys().count(rule.variable) != 0: 
            if self.algebraicNodes[rule.variable] != value:
                self.algebraicNodes[rule.variable] = value
       #/ In the case of a non constant parameter
        elif self.nonConstantParameters.keys().count(rule.variable) != 0: 
            if self.nonConstantParameters[rule.variable] != value:
                self.nonConstantParameters[rule.variable] = value
       #/ In the case of a constant parameter
        elif self.parameters.keys().count(rule.variable) != 0: 
             print 'WARNING: Error in the SBML file, the parameter %s is not constant.' %rule.variable
             if  self.parameters[rule.variable] != value:
                 self.nonConstantParameters[rule.variable] = value
             else:
                 self.nonConstantParameters[rule.variable] = self.parameters[rule.variable]
             del self.parameters[rule.variable]
        #/ In the case of a non constant compartment
        elif self.nonConstantCompartments.keys().count(rule.variable) != 0: 
            if self.nonConstantCompartments[rule.variable]['size'] != value:
                self.nonConstantCompartments[rule.variable]['size'] = value
        #/ In the case of constant compartment
        elif self.compartments.keys().count(rule.variable) != 0:
             print 'WARNING: Error in the sbml file, the compartment %s is not constant.' %rule.variable
             if  self.compartments[rule.variable]['size'] != value:
                 self.nonConstantCompartments[rule.variable] = self.compartments[rule.variable]
                 self.nonConstantCompartments[rule.variable]['size'] = value
             else:
                 self.nonConstantCompartments[rule.variable] = self.compartments[rule.variable]
             del self.compartments[rule.variable]
        #/ Error the variable is unknown
        else:
			print 'WARNING: The variable %s of assignment rule is unknown, the rule will be ignored.' %rule.variable
        
        return self


    def __checkAlgebraicNodes(self):
        '''
        This method searches if all the species defined with boundary condition true and constant false
        are used in rules or events. Otherwise, it converts this species like constant species.
        '''
        all = False
        if len(self.rules) == 0 and len(self.events) == 0:
            all = True
        
        if all == True:
            for specie in self.algebraicNodes:
                self.constantNodes[specie] = self.algebraicNodes[specie]
            self.algebraicNodes = {}
            
        return self




    def __checkAndPutRaterules(self):
	
        '''
	This method checks the rules and puts rate rules into the model topology.
	'''

        for rule in self.rules:
            if rule.type == 'Rate': 
                self.__putRateruleInTopology(rule)
            elif rule.type == 'Assignment':
                self.checkAssignmentRule(rule)

        return self
    
    def __checkNonConstant(self):

        '''
        This internal method checks all the nonConstant parameters and compartments.
        If it finds that one is not used in the equations of the system,
        this one is moved to constant.
        A warning is showed because in this point we are modifing the information of SBML file,
        and it can be an error of SBML file.          
        '''

        remove = []
        for nc in self.nonConstantParameters:
            if self.isNonConstant(nc) == False:
                print 'WARNING: %s is a non constant parameter but it is not actually modified by any equation. ByoDyn will consider %s as a constant.' %(nc, nc)
                self.parameters[nc] = self.nonConstantParameters[nc]
                remove.append(nc)
        for r in remove:
            self.nonConstantParameters.pop(r)
        remove = []
        for nc in self.nonConstantCompartments:
            if self.isNonConstant(nc) == False:
                print 'WARNING: %s is a non constant compartment but it is not actually modified by any equation. ByoDyn will consider %s as a constant.' %(nc, nc)
                self.compartments[nc] = self.nonConstantCompartments[nc]
                remove.append(nc)
        for r in remove:
            self.nonConstantCompartments.pop(nc)            
        
        return self

    def isNonConstant(self, nc):
	
	'''
	This method checks if some nonConstant element (parameter or compartment)
	is modified by rules or events. It returns a boolean with the answer.
	'''

        answer = False
        for rule in self.rules:
            if rule.variable == nc:
                answer = True
            if rule.type == 'Algebraic':
                if rule.math.count(nc) != 0:
                    answer = True
        for event in self.events:
            for assignment in event.assignment:
                if assignment == nc:
                    answer = True
        
        return answer        

    def __putRateruleInTopology(self, rule):
        
        '''
	This method introduces rate rules in the model topology.
	'''
        
        #/ In the case of a species
        if self.nodes.count(rule.variable) != 0: 
            for definition in self.topology.keys():
                if definition.split('/') == rule.variable:
                    del self.topology[definition]
        #/ In the case of a parameter        
        elif self.parameters.keys().count(rule.variable) != 0: 
             print 'WARNING: Error in the sbml file, the parameter %s is not constant.' %rule.variable
             self.nonConstantParameter[rule.variable] = self.parameters[rule.variable]
             del model.parameters[rule.variable]
        #/ In the case of a compartment    
        elif self.compartments.keys().count(rule.variable) != 0:               
             print 'WARNING: Error in the sbml file, the compartment %s is not constant.' %rule.variable
             self.nonConstantCompartments[rule.variable] = self.compartments[rule.variable]
             del self.compartments[rule.variable]                   
        #/ Error the variable is unknown
        elif self.nonConstantParameters.has_key(rule.variable) == False and self.nonConstantCompartments.has_key(rule.variable) == False:
            print 'WARNING: The variable %s of rate rule is unknown, this rule will be ignored' %rule.variable
        #/ Adding in the topology
        definition = '%s/SBML/RULE/%s' %(rule.variable, rule.id)        
        self.topology[definition] = rule.math
        
        return self

    def summary(self, outputfiles): #/ method of ClassModelSBML        
        
	'''
	This method writes a short summary of the model description in an output file.
	'''
        
        f = open(outputfiles.summaryModel, 'w')
        f.write('Description of the current model:\n\n')
        f.write('\tSystem name: %s ;\n' %(self.systemName))
        f.write('\tSpecies (nodes): %s ;' %(len(self.nodes)))
	for node in self.nodes:
	    f.write('\t%s'%node)
	f.write('\n')
        f.write('\tConstant Species: %s ;' %(len(self.constantNodes.keys())))
	for constant in self.constantNodes.keys():
	    f.write('\t%s'%constant)
	f.write('\n')
        f.write('\tConstant Parameters: %s ;' %(len(self.parameters.keys())))
	for constant in self.parameters.keys():
	    f.write('\t%s=%s'%(constant, self.parameters[constant]))
	f.write('\n')
        f.write('\tNon-constant Parameters: %s ;' %(len(self.nonConstantParameters.keys())))
	for constant in self.nonConstantParameters.keys():
	    f.write('\t%s'%constant)
	f.write('\n')
        #/ Rules
        if self.rules:  #if len(self.rules) != 0:
            f.write('\tRules: %s\t Descriptions:\n\n' %(len(self.rules)))
            for rule in self.rules:
                f.write('\t\t\tType:%s\n' %(rule.type))                
                f.write('\t\t\tId:%s\n' %(rule.id))
                if rule.type == 'Assignment':
                    f.write('\t\t\t%s = %s\n' %(rule.variable, rule.math))
                elif rule.type == 'Rate':
                    f.write('\t\t\td%s/dt = %s\n' %(rule.variable, rule.math))             
                elif rule.type == 'Algebraic':
                    f.write('\t\t\t%s = 0\n' %(rule.math))                            
                f.write('\n')            
        #/ Events
        if self.events:    
            f.write('\tEvents: %s\t Description:\n' %(len(self.events)))
            for event in self.events:             
                f.write('\t\t\tId:%s\n' %(event.id))
                f.write('\t\t\tdelay = %s\n' %(event.delay))
                f.write('\t\t\ttrigger = %s\n' %(event.trigger))
                f.write('\t\t\tnumAssig = %s\n' %(len(event.assignment)))
                f.write('\n')
        f.close() 
	
        return None
	
    def __readCompartments(self, sbmlModel):
        
        '''
	This method reads the compartments from a SBML model.
	It discriminates between constant and non constant compartments.
	'''

        listOfCompartments = sbmlModel.getListOfCompartments()
        for compartment in listOfCompartments:
        	if compartment.getConstant() == False:
        		print 'WARNING: The compartment %s is not constant, this indicates the compartment\'s size can be changed by rules and the size is actually intended to be initial size of the compartment.' %compartment.getId()
        		self.nonConstantCompartments[compartment.getId()] = {} 
        		self.nonConstantCompartments[compartment.getId()]['size'] = compartment.getSize()
        		self.nonConstantCompartments[compartment.getId()]['spatialDimensions'] = compartment.getSpatialDimensions()
        	else:
        		self.compartments[compartment.getId()] = {} 
        		self.compartments[compartment.getId()]['size'] = compartment.getSize()
        		self.compartments[compartment.getId()]['spatialDimensions'] = compartment.getSpatialDimensions()
        
        return self

    def __readConstraints(self, sbmlModel):

        '''
        This method reads constraints from a SBML file
        '''

        listOfConstraints = sbmlModel.getListOfConstraints()
        if len(listOfConstraints) != 0:
            print 'WARNING: This model have SBML constraints and ByoDyn does not work using this feature. It will be ignored'

        return self

    def __readNodes(self, sbmlModel):
 
        '''
	This method reads of species from a SBML file.
	It sets the initial conditions of the model.
	'''
        
        boundaryConditionNodes = {}
        initialConditions = []
        listOfSpecies = sbmlModel.getListOfSpecies()
        #/ We are doing the initial conditions in two steps to sort the initialConditions in the same manner that list of nodes in self.nodes
        for species in listOfSpecies:
            self.nodes.append(species.getId())
            compartment = species.getCompartment()
            if compartment in self.compartments:
                amount = species.getInitialAmount()
                concentration = species.getInitialConcentration()
                if amount == concentration == 0:
                    initialConcentration = 0
                else:
                    if amount != 0:
                        initialConcentration = amount / float(self.compartments[compartment]['size'])
                    else:
                        initialConcentration = concentration
            elif compartment in self.nonConstantCompartments:
                amount = species.getInitialAmount()
                concentration = species.getInitialConcentration()
                if amount == concentration == 0:
                    initialConcentration = 0
                else:
                    if amount != 0:
                        initialConcentration = amount / float(self.nonConstantCompartments[compartment]['size'])
                    else:
                        initialConcentration = concentration
            else:
                raise errorMessages.ClassSBMLWorkerException,  'Error: %s is in unknown compartment %s' %(species.getId(), species.getCompartment())	
            self.initialConditions.append(initialConcentration)        	
            if species.getConstant() == True:
                if species.getBoundaryCondition() == True:
                    boundaryConditionNodes[species.getId()] = initialConcentration
                else:
                    boundaryConditionNodes[species.getId()] = initialConcentration
            else:
                if species.getBoundaryCondition() == True:
                    self.algebraicNodes[species.getId()] = initialConcentration #/ in this case it can be changed by rules
        self.constantNodes = boundaryConditionNodes

        return self

    def __readFunctionDefinitions(self, sbmlModel):
        
        '''
        This method reads the user function definitions from the SBML file.
        '''        

        listFunctionDefinitions = sbmlModel.getListOfFunctionDefinitions()
        for function in listFunctionDefinitions:
            currentFunction = ClassFunction()
            currentFunction.id = function.getId()
            currentFunction.mathAST = libsbml.writeMathMLToString(function.getMath())
            expr = re.findall('\([\w\[\]()/\+\-\*\s\d\.|,\^]*,[\w\[\]()/\+\-\*\s\d\.\^|,]*\)', libsbml.formulaToString(function.getMath()))[0]
            expr = expr[1:(len(expr)-1)]
            fields = expr.split(',')
            for i in range(function.getNumArguments()):
                currentFunction.arguments.append(fields[i].strip())
            currentFunction.output = ''
            for j in range(i+1, len(fields)):
                if currentFunction.output == '':
                    currentFunction.output += fields[j].strip()
                else:
                    currentFunction.output += ', ' + fields[j].strip()        
            self.functions.append(currentFunction)

        return self

    def __readEvents(self, sbmlModel):

        '''
	This method reads the events from a SBML file.
	'''
        
        listOfEvents = sbmlModel.getListOfEvents()
        for event in listOfEvents:
            currentEvent = ClassEvent() #/ initializing an instance of the class ClassEvent
            currentEvent.id = event.getId() #/ obtaining the "id" of the event
#            currentEvent.trigger = libsbml.formulaToString(event.getTrigger().getMath()) #/ obtaining the trigger of the event. We send directly the string
            currentEvent.trigger, op = getBooleanExpression(event.getTrigger().getMath())
            currentEvent.delay = event.getDelay() #/ checking if there is a delay
            if currentEvent.delay != None:
                currentEvent.delay = libsbml.formulaToString(currentEvent.delay.getMath())
#                raise errorMessages.ClassSBMLWorkerException,  'Error: This version of ByoDyn does not support Events with Delays. Sorry for the inconvenience.'
            assignments = event.getListOfEventAssignments()
            for assignment in assignments:
                currentEvent.assignment[assignment.getVariable()] = libsbml.formulaToString(assignment.getMath())
                currentEvent.assignmentAST[assignment.getVariable()] = libsbml.writeMathMLToString(assignment.getMath())                
            self.events.append(currentEvent)

        return self

    def __readInitialAssignments(self, sbmlModel):

        '''
        This method reads the initial Assignments from a SBML file.
        '''

        listOfAssignments = sbmlModel.getListOfInitialAssignments()
        if len(listOfAssignments) != 0:
            print 'Model with initial Assignments'
            sys.exit()

        return self

    def __readRules(self, sbmlModel):
        
	'''
	This method reads the rules from a SBML file.
	'''
        
        listOfRules = sbmlModel.getListOfRules()
        #/ There are three different possible functional forms: Algebraic rule, Assignment Rule, Rate Rule
        #/ Warning: The ordering of assignment rules is significant: they are always evaluated in the order given in SBML!
        for rule in listOfRules:
            objectRule = ClassRules()
            objectRule.id = rule.getMetaId()
            objectRule.math = rule.getFormula()
            objectRule.mathAST = libsbml.writeMathMLToString(rule.getMath())
            #/ Assignment rules
            if rule.isAssignment() == True: 
                objectRule.type = 'Assignment'
                objectRule.variable = rule.getVariable()
	    #/ Rate rules
            elif rule.isRate() == True: #/ 18 is the code for rate rules
                objectRule.type = 'Rate'
                objectRule.variable = rule.getVariable()            
	    #/ Algebraic rules
            elif rule.isAlgebraic() == True: 
                objectRule.type = 'Algebraic'
            self.rules.append(objectRule) 

        return self

    def __readParameters(self, sbmlModel, metamodel):
    	
        ''' 
	This method reads the parameters from a SBML file.
	'''

        listOfGlobalParameters = sbmlModel.getListOfParameters()
        for parameter in listOfGlobalParameters:
            if parameter.getConstant() == False:
                self.nonConstantParameters[parameter.getId()] = parameter.getValue() 
            else:
                self.parameters[parameter.getId()] = parameter.getValue() 
        reactions = sbmlModel.getListOfReactions()
        for reaction in reactions:
            if reaction.getKineticLaw() != None:
                listOfLocalParameters = reaction.getKineticLaw().getListOfParameters()
                for parameter in listOfLocalParameters:
                    if parameter.getConstant() == False:
                        self.nonConstantParameters['%s%s' %(parameter.getId(), reaction.getId())] = parameter.getValue() 
                    else:
                        self.parameters['%s%s' %(parameter.getId(), reaction.getId())] = parameter.getValue() 

        return self

    def __searchDelayFunction(self, sbmlModel):

	'''
	This internal method searchs delay functions in kinetic laws.
    If delay is found, it will create a non constant parameter 
    and an assignment rule in order to handle this delay in integrator solvers.
	'''

        delayParameters = []
        for definition in self.topology.keys():
            if self.topology[definition].rfind('delay') != -1:
                self.delayFunctions = True
                delays = re.findall('delay\([\w\[\]/\+\-\*\s\d\.\,\^]*,[\w\[\]/\+\-\*\s\d\.\^\,]*\)', self.topology[definition])  
                for delay in delays:
                    first = delay[6:len(delay)-1]
                    terms = re.findall('[\w\[\]/\+\-\*\s\d\.\^]*', first)
                    self.delays.append([terms[0].strip(), terms[2].strip()])
                    if len(terms[0].split()) > 1:   
                        if delayParameters.count(terms[0]) == 0:
                            delayParameters.append(terms[0])
                        self.topology[definition] = self.topology[definition].replace(delay, 'delay(delay%s,%s)' %(delayParameters.index(terms[0]), terms[2]))
        for i in range(len(delayParameters)):
            objectRule = ClassRules()
            objectRule.id = 'delay%s' %i
            objectRule.math = delayParameters[i]
            objectRule.type = 'Assignment'
            objectRule.variable = 'delay%s' %i
            self.rules.append(objectRule)
            self.nonConstantParameters['delay%s' %i] = formulas.solveFormula(self, objectRule.math)

        return self


    def __readTopology(self, sbmlModel):
 
        '''
	This method builds the network topology from SBML file.
	'''

        reactions = sbmlModel.getListOfReactions()
        reactionIndex = 0        

        self.stoichiometryMatrix = scipy.zeros([len(self.nodes), len(reactions)])
        for reaction in reactions:
            reactionIndex = reactionIndex + 1
            reactionObject = ClassReaction()
            reactionObject.id = reaction.getId()
            if reaction.getKineticLaw() != None:
                reactionObject.kineticLaw = reaction.getKineticLaw().getFormula()
                reactionObject.propensity = translateLocalParametersNamesInFormula(reaction.getKineticLaw().getFormula(), reaction)
            if reaction.getFast():
                print 'WARNING: The reaction %s has fast reactions but ByoDyn does not work using this features, so it will be ignored.' %reaction.getId()            
            stoichiometryReactants = {}
            stoichiometryProducts = {}  
            reactants = []
            products = []
            reactantList = reaction.getListOfReactants()
            productList = reaction.getListOfProducts()
            for reactant in reactantList:
                if reactant.isSetStoichiometryMath():
                    print 'WARNING: The stoichiometry in reactions cointains StoichiometryMath expressions. You need to check the ByoDyn behaviour.'
                    stoichiometryReactants[reactant.getSpecies()] = libsbml.formulaToString(reactant.getStoichiometryMath().getMath())
                    self.stoichiometryMatrix[self.nodes.index(reactant.getSpecies())][reactionIndex-1] -= reactant.getStoichiometry()
                else:
                    stoichiometryReactants[reactant.getSpecies()] = reactant.getStoichiometry()
                    self.stoichiometryMatrix[self.nodes.index(reactant.getSpecies())][reactionIndex-1] -= reactant.getStoichiometry()
                reactants.append(reactant.getSpecies())
                reactionObject.reactants.append(reactant.getSpecies())
            for product in productList:
                if product.isSetStoichiometryMath():
                    print 'WARNING: The stoichiometry in reactions cointains StoichiometryMath expressions. You need to check the ByoDyn behaviour.'
                    stoichiometryProducts[product.getSpecies()] = libsbml.formulaToString(product.getStoichiometryMath().getMath())
                    #/ stoichiometryMatrix cannot have math expressions, so stoichiometryMath expressions are not compatibles with stochastic simulations. 
                    self.stoichiometryMatrix[self.nodes.index(product.getSpecies())][reactionIndex-1] += product.getStoichiometry()
                else:
                    stoichiometryProducts[product.getSpecies()] = product.getStoichiometry()   
                    #/ stoichiometryMatrix cannot have math expressions, so stoichiometryMath expressions are not compatibles with stochastic simulations.
                    self.stoichiometryMatrix[self.nodes.index(product.getSpecies())][reactionIndex-1] += product.getStoichiometry()
                products.append(product.getSpecies())
                reactionObject.products.append(product.getSpecies())
            self.reactions.append(reactionObject)
            #/ checking that there is reactant and product or there is flux reaction
            if reactants == []:
                reactants.append('flux')
                stoichiometryReactants['flux'] = 1
            if products == []:
                products.append('flux')
                stoichiometryProducts['flux'] = 1
            for productNode in products:
                reactantNode = reactants[0]
                if reaction.getKineticLaw() != None:
                    ASTformula = ['', '', '']
                    formula =  reaction.getKineticLaw().getFormula()
                    formula = translateLocalParametersNamesInFormula(formula, reaction)
                    ASTformula[1] = libsbml.writeMathMLToString(reaction.getKineticLaw().getMath())
                    formula, ASTformula[2] = divCompartment(sbmlModel, formula, productNode) 
    	        #/ if the stoichiometry is different than 1, we need multiply the formula by the stoichiometry coefficient
                    if stoichiometryProducts[productNode] != 1:
                        self.topology['%s/SBML/%s/forward.%s' % (productNode, reactantNode, reaction.getId())] = '+ %s * (%s)' %(stoichiometryProducts[productNode], formula)
                        ASTformula[0] = '+ %s *' %stoichiometryProducts[productNode]
                        self.topologyAST['%s/SBML/%s/forward.%s' % (productNode, reactantNode, reaction.getId())] = ASTformula
                    else:
                        self.topology['%s/SBML/%s/forward.%s' % (productNode, reactantNode, reaction.getId())] = '+ %s' %formula
                        ASTformula[0] = '+'
                        self.topologyAST['%s/SBML/%s/forward.%s' % (productNode, reactantNode, reaction.getId())] = ASTformula
            for reactantNode in reactants: 
                productNode = products[0]
                if reaction.getKineticLaw() != None:
                    ASTformula = ['', '', '']
                    formula =  reaction.getKineticLaw().getFormula()
                    formula = translateLocalParametersNamesInFormula(formula, reaction)
                    ASTformula[1] = libsbml.writeMathMLToString(reaction.getKineticLaw().getMath())
                    formula, ASTformula[2] = divCompartment(sbmlModel, formula, reactantNode)
                    #/ if the stoichiometry is different than 1, we need multiply the formula by the stoichiometry coefficient
                    if stoichiometryReactants[reactantNode] != 1:
                        self.topology['%s/SBML/%s/backward.%s' % (reactantNode, productNode, reaction.getId())] = '- %s * (%s)' %(stoichiometryReactants[reactantNode], formula)
                        ASTformula[0] = '- %s *' %stoichiometryReactants[reactantNode]
                        self.topologyAST['%s/SBML/%s/backward.%s' % (reactantNode, productNode, reaction.getId())] = ASTformula
                    else:
                        self.topology['%s/SBML/%s/backward.%s' % (reactantNode, productNode, reaction.getId())] = '- (%s)' %formula  
                        ASTformula[0] = '-'
                        self.topologyAST['%s/SBML/%s/backward.%s' % (reactantNode, productNode, reaction.getId())] = ASTformula
        #/ Dealing with boundary conditions and constant nodes
        if self.constantNodes.keys() != []:
            formula = '+ 0'
            for constantNode in self.constantNodes.keys():
                for term in self.topology:
                    if term.split('/')[0] == constantNode and self.topology[term] != formula:
                        self.topology[term] = formula    
  
        return self
      
    def getDefaultValue(self, a):

	'''
	This method gets the values of the non constant parameters, 
	the non constant compartments and the algebraic nodes.
	This method is required from the simulatorOpenModelica for proper writing of the equations.
	'''

	if a in self.nonConstantParameters.keys():
            return self.nonConstantParameters[a]
	elif a in self.nonConstantCompartments.keys():
            return self.nonConstantCompartments[a]
	elif a in self.algebraicNodes.keys():
            return self.algebraicNodes[a]
	else:
            if a in self.parameters.keys():
                print 'WARNING: %s variable in event is a constant parameter\n' %a
                return self.parameters[a]
            elif a in self.nodes:
                return self.initialConditions[self.nodes.index(a)]
            else:
                raise errorMessages.ClassSBMLWorkerException,  'Error: Some variable in event is not found in the model.'


#/
#/ From this point functions are defined. There are no more methods of classModel.
#/
 
def divCompartment(sbmlModel, formula, node):

    '''
    This function divides the formula by compartment volume if it is needed.
    '''

    listOfSpecies = sbmlModel.getListOfSpecies()
    listOfCompartments = sbmlModel.getListOfCompartments()
    astDiv = ''
    for species in listOfSpecies:
        if node == species.getId():
            for compartment in listOfCompartments:
                if species.getCompartment() == compartment.getId():
                    if compartment.getSize() != 1 or compartment.getConstant() == False:
                        formula = '(%s)/%s' %(formula, compartment.getId())
                        astDiv = '/ %s' %compartment.getId()
                    break 
            break

    return formula, astDiv


def getBooleanExpression(ASTNode):
    
    '''
    This function reads an ASTNode (equation in Abstract Syntax Tree class)
    and builds the equation as a math string using recursive calls.
    It returns the math string and the operator.
    '''
    
    numChildren = ASTNode.getNumChildren()
    #/ 1.- this is a single value
    if (numChildren == 0):    
        if (ASTNode.isReal()):    return str(ASTNode.getReal()), None
        if (ASTNode.isInteger()):    return str(ASTNode.getInteger()), None
        return ASTNode.getName(), None;

    #/ 2.- the expression is ...
    eq = ""                                                       
    if ASTNode.isOperator():
        operator = ASTNode.getCharacter()
    else:
        operator = translateMathFactor(ASTNode.getName())
    for i in range (0, numChildren, 1):
        subeq, subop = getBooleanExpression(ASTNode.getChild(i))
        if (subop == None):
            if subeq == 't':
                subeq = 'time'
            eq += (subeq + operator)
        else:
            eq += ("("+subeq+")"+operator)
    eq = eq[0:len(eq)-len(operator)]
    
    return eq, operator

        
def translateMathFactor(mathFactor):

    '''
    This function converts some string logical affectors to Open Modelica format.
    '''

    if mathFactor == 'gt':    return '>'
    elif mathFactor == 'lt':    return '<'
    elif mathFactor == 'eq':    return '<=' #/ it is <= because openModelica don't can work with '=='. It is a problem for models where the trigger variable can have negative numbers
    elif mathFactor == 'neq':    return '<>'
    elif mathFactor == 'geq':    return '>='
    elif mathFactor == 'leq':    return '<='
    elif mathFactor == 'or':    return 'or'
    elif mathFactor == 'and':    return 'and'
    elif mathFactor == 'plus':    return '+'
    elif mathFactor == 'minus':    return '-'
    else:
	       return None

def translateLocalParametersNamesInFormula(formula, reaction):

    '''
    This function replaces in a formula all the names of parameters,
    because we need convert the local paramters to a new format ('parameterId''ReactionId')
    in order to work correctly with SBML definition of parameters
    where the local and global parameters can have the same ids.
    '''

    for parameter in reaction.getKineticLaw().getListOfParameters():
        newName = parameter.getId() + reaction.getId()
        prog = re.compile('\A%s\s|\s%s\Z' %(parameter.getId(), parameter.getId()))
        formula = prog.sub(newName, formula)
        prog = re.compile('\A%s\Z' %(parameter.getId()))
        formula = prog.sub(newName, formula)
        prog = re.compile('/%s\Z' %(parameter.getId()))
        formula = prog.sub(newName, formula)
        formula = formula.replace('-%s ' %parameter.getId(),'-%s ' %newName)    
        formula = formula.replace(' %s ' %parameter.getId(),' %s ' %newName)    
        formula = formula.replace('(%s ' %parameter.getId(),'(%s ' %newName)
        formula = formula.replace(' %s)' %parameter.getId(),' %s)' %newName)
        formula = formula.replace('-%s)' %parameter.getId(),'-%s)' %newName)
        formula = formula.replace(' %s/' %parameter.getId(),' %s/' %newName)
        formula = formula.replace('/%s ' %parameter.getId(),'/%s ' %newName)
        formula = formula.replace('(%s,' %parameter.getId(),'(%s,' %newName)
        formula = formula.replace(',%s)' %parameter.getId(),',%s)' %newName)
        formula = formula.replace('(%s/' %parameter.getId(),'(%s/' %newName)
        formula = formula.replace('/%s)' %parameter.getId(),'/%s)' %newName)
        formula = formula.replace(' %s,' %parameter.getId(),' %s,' %newName)

    return formula

#/
#/ writting SBML 
#/

def sbmlWriter(w, model, metamodel, file): 

    '''
    This function writes a SBML model given the object model.
    '''
    
    #/ 0.- Checking that the model is single cell
    if model.xlength != 1 or model.ywidth != 1:
        raise errorMessages.ClassSBMLWorkerException, 'sorry, this version of ByoDyn is only capable to export unicellular models into SBML.'
    #/ 1.- reading the file for the second time
    modelFile = libsbml.readSBML(metamodel.modelFile)
    sbmlModel = modelFile.getModel()
    #/ 2.- setting the new model
    #/ 2.1.- setting the new parameter values
    #/ 2.1.1.- for global parameters
    globalConstantParameters = []
    globalParameters = sbmlModel.getListOfParameters()
    for parameter in globalParameters:
	if parameter.getConstant() == True:
	    novelParameter = parameter.getId()
	    for modifiedParameter in model.parameters.keys():
		if modifiedParameter == novelParameter:
		    parameter.setValue(model.parameters[modifiedParameter])
    #/ 2.1.2.- for local parameters
    localConstantParameters = []
    for reaction in sbmlModel.getListOfReactions():
	for parameter in reaction.getKineticLaw().getListOfParameters():
	    if parameter.getConstant() == True:
		novelParameter = parameter.getId()+reaction.getId()
		for modifiedParameter in model.parameters.keys():
		    if modifiedParameter == novelParameter:
			parameter.setValue(model.parameters[modifiedParameter])
    #/ 2.2.- setting the new initial conditions, in case of
    listOfSpecies = sbmlModel.getListOfSpecies()
    for i in range(len(listOfSpecies)):
	listOfSpecies[i].setInitialConcentration(model.initialConditions[i])
    #/ 3.- converting the model into a new document
    doc = libsbml.SBMLDocument()
    doc.setModel(sbmlModel)
    #/ 4.- actual writting of the file    
    w.writeSBML(doc, file) #/ This line is not compatible for libSBML < 3.0.0

    return None

