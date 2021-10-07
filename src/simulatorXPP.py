#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Alex Gomez-Garrido
#
#  Created: 2007-10-01 by Alex Gomez-Garrido
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

# $Id: simulatorXPP.py,v 4.13 2008/12/16 21:35:06 alglomana Exp $

## \file
# This module simulates the model in the case the integration option xpp has been selected.

import os, re, sys, copy
import formulas, initiator, errorMessages, sbmlWorker
from affectors import *

class ClassSimulatorXPP:

    '''
    Class for the XPP-AUT simulator.
    '''    

    def __init__(self):

        '''
        The constructor.
        '''
        
        return None

    def __checkNames(self, model):

        '''
	XPPAUT converts all the strings into capittal letters, so we need to check for repeated nodes.
        '''

	#/ 1.- defining new variables
        names = []
        newNames = []
        conflictives = []
        modifyModel = copy.deepcopy(model)
        modifyModel.nodes = []
        modifyModel.algebraicNodes = {}
        modifyModel.compartments = {}
        modifyModel.nonConstantCompartments = {}
        modifyModel.nonConstantParameters = {}
        modifyModel.parameters = {}
        modifyModel.rules = []  
        modifyModel.constantNodes = {}
        modifyModel.topology = {}
        modifyModel.topologyAST = {}
        names.append('PI')
        newNames.append('PI')
        renamed = {}
	#/ 2.- model parameters
        for name in model.parameters:
            newName = convertName(name)
            if newName in newNames:
                newName += 'conf'
            modifyModel.parameters[newName] = model.parameters[name]
            newNames.append(newName)
            names.append(name)
	#/ 3.- model compartments
        for name in model.compartments:
            newName = convertName(name)
            if newName in newNames:
                newName += 'conf'
            modifyModel.compartments[newName] = model.compartments[name]
            newNames.append(newName)
            names.append(name)
	#/ 4.- model algebraic nodes
        for name in model.algebraicNodes:
            newName = convertName(name)
            if newName in newNames:
                newName += 'conf'
            modifyModel.algebraicNodes[newName] = model.algebraicNodes[name]
            newNames.append(newName)
            names.append(name)
	#/ 5.- model non constant parameters
        for name in model.nonConstantParameters:           
            newName = convertName(name)
            if newName in newNames:
                newName += 'conf'
            modifyModel.nonConstantParameters[newName] = model.nonConstantParameters[name]
            newNames.append(newName)
            names.append(name)
	#/ 6.- model non constant compartments
        for name in model.nonConstantCompartments:
            newName = convertName(name)
            if newName in newNames:
                newName += 'conf'
            modifyModel.nonConstantCompartments[newName] = model.nonConstantCompartments[name]
            newNames.append(newName)
            names.append(name)
	#/ 7.- model nodes
        for name in model.nodes:
            newName = convertName(name)	    
            if newName in newNames:
                newName += 'conf'
            modifyModel.nodes.append(newName)
            if name in model.constantNodes:
                modifyModel.constantNodes[newName] = model.constantNodes[name]
            newNames.append(newName)
            names.append(name)
	#/ 8.- model rules
        for rule in model.rules:
            modifyRule = sbmlWorker.ClassRules()
            ruleVariable = convertName(rule.variable)
            modifyRule.id = rule.id
            modifyRule.type = rule.type
            modifyRule.math = replaceNames(rule.math, names, newNames, conflictives)
            modifyRule.variable = ruleVariable
            modifyModel.rules.append(modifyRule)
	#/ 9.- model topology
        for t in model.topology:
            topology = model.topology[t]
            fields = t.split('/')
            for i in range(len(names)):
                if fields[0] == names[i]:
                    fields[0] = newNames[i]
                if fields[2] == names[i]:
                    fields[2] = newNames[i]
            newT = fields[0] + '/' + fields[1]+ '/' + fields[2]+ '/' + fields[3]
            topology = replaceNames(topology, names, newNames, conflictives)
            modifyModel.topology[newT] = topology
            modifyModel.topologyAST[newT] = None

        return modifyModel

    def createInput(self, metamodel, model, outputfiles):

        '''
        This method creates the input file for the XPPAUT simulator.
        '''

        #/ 0.- XPP converts all the strings in capittal letters, so we need to check repetitions
        modifyModel = self.__checkNames(model)
        #/ 1.- Some initial variables necessary for the affectors module.
        option = 'xpp'
        cellIndex = 0
        #/ 2.- writing the input file        
        file = open(outputfiles.xppIntegrationFile, 'w')        
        file.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
        file.write('pow(x,y)=x^y\n')
        #/ 2.1.- writting sbml defined functions
        for function in model.functions:
            file.write('%s(' %(function.id))
            file.write('%s' %function.arguments[0])
            for argument in function.arguments[1:]:
                file.write(',%s' %argument)
            file.write(')=')    
            squareRootDefinitions = re.findall('root\(2, [\w\[\]()/\+\-\*\s\d\.\^\,]*\)', function.output)
            if len(squareRootDefinitions) != 0:
                function.output = function.output.replace('root(2, ', 'sqrt(')
            file.write('%s\n' %function.output)
        #/ 2.2.- writting the nodes and its initial conditions
        for i in range(len(modifyModel.initialConditions)):
            if modifyModel.nodes[i] in modifyModel.constantNodes.keys():
                file.write('par %s=%s\n'%(modifyModel.nodes[i], modifyModel.initialConditions[i]))
            else:
                file.write('init %s=%s\n' %(modifyModel.nodes[i], modifyModel.initialConditions[i]))
	#/ 2.3.- checking rules
        for rule in modifyModel.rules:
            if rule.type == 'Rate':
                if modifyModel.nodes.count(rule.variable) == 0:
                    file.write('init %s=%s\n' %(rule.variable, modifyModel.getDefaultValue(rule.variable)))
	#/ 2.4.- writing the parameters
        for parameter in modifyModel.parameters.keys():
            file.write('par %s=%s\n'%(parameter, modifyModel.parameters[parameter]))
        for compartment in modifyModel.compartments.keys():
            file.write('par %s=%s\n'%(compartment, modifyModel.compartments[compartment]['size']))
        #/ 2.5.- writting the equations
        for i in range (modifyModel.xlength):
            for j in range (modifyModel.ywidth):
               for node in modifyModel.nodes:
                    if modifyModel.constantNodes.keys().count(node) == 0:
                        file.write('d%s/dt = 0 ' %(node))
                        for definition in modifyModel.topology.keys():
                            fieldsDefinition = definition.split('/')
                            if fieldsDefinition[0] == node:
				#print modifyModel.topology[definition]
                                exec "%s(modifyModel, file, option, cellIndex, definition, fieldsDefinition)"%fieldsDefinition[1]
                        file.write('\n')
               for rule in modifyModel.rules:
                    if rule.type == 'Rate':
                        file.write('d%s/dt = ' %(rule.variable))
                        formulas.writeXPPFormula(rule.math, file)
                    if rule.type == 'Algebraic':       
                        raise errorMessages.ClassSimulatorXPPException, 'XPP cannot handle algebraic rules. Please use Octave or OpenModelica solvers.'
                        sys.exit()
                    if rule.type == 'Assignment':
                        file.write('%s = ' %rule.variable)
                        formulas.writeXPPFormula(rule.math, file)
                    file.write('\n')                   
        #/ 2.6- defining some simulation options
	#/ 2.6.1.- selecting the biggest delay
	if model.delayFunctions == True:
	    uniqueDelays = []
	    for delay in model.delays:
		if uniqueDelays.count(delay[1]) == 0:
		    uniqueDelays.append(delay[1])
	    delayValues = []
	    for delay in uniqueDelays:
		for parameter in model.parameters.keys():
		    if parameter == delay:
			delayValues.append(model.parameters[delay])
	    maxDelay = max(delayValues)	    
	    file.write('@ delay=%s\n'%maxDelay)
	#/ 2.6.2.- other simulation parameters
        if len(metamodel.integrationTolerance) == 0:
            file.write('@ meth=cvode, tol=1e-6, atol=1e-8\n')
        else:
            tol = float(metamodel.integrationTolerance[0])
            atol = float(metamodel.integrationTolerance[1])
            if tol < 6.62554e-16:
                raise errorMessages.ClassSimulatorXPPException, 'XPP cannot use tolerance (first field of tolerance) of integrator lower than 6.62554e-16.'
            if atol < 6.62554e-18:
                raise errorMessages.ClassSimulatorXPPException, 'XPP cannot use absolute tolerance (second field of tolerance) of integrator lower than 6.62554e-18.'
            file.write('@ meth=cvode, tol=%s, atol=%s\n' %(tol, atol))
        file.write('@ total=%s, dt=%s\n' %(metamodel.simulationTime, metamodel.simulationTimeStep))
        file.write('@ maxstor=%s, bound=1000200000\n' %(int(metamodel.simulationTime / metamodel.simulationTimeStep  + 10.0)))
	file.write('@ output=%s.XPP.data\n'%model.systemName)
        file.write('done')
        file.close()

        return None
        
    def callSolver(self, outputfiles):
    
        '''
        This method calls xppaut to simulate the model
        '''
	
	currentDirectory = os.getcwd()
        os.chdir(outputfiles.scratchdir)
        os.chdir(outputfiles.scratchdir)
        os.system('xppaut -silent %s > %s' %(outputfiles.xppIntegrationFile.split('/')[len(outputfiles.xppIntegrationFile.split('/'))-1], outputfiles.xppOutputDetails))
	os.chdir(currentDirectory)

        return None
    
    def createOutputs(self, model, outputfiles, metamodel):

        '''
        This method converts the format of the XPP output file into a more adequate for ByoDyn.
        '''    

        #/ 1.- setting the headers
        stringNodes = '# t'
        for node in model.nodes:
            if model.constantNodes.keys().count(node) == 0:
                stringNodes = stringNodes + ' ' + node 
        for rule in model.rules:
            if rule.type == 'Rate':
                if model.nodes.count(rule.variable) == 0:
                    stringNodes = stringNodes + ' ' + rule.variable 
        stringNodes = stringNodes + '\n'
        f = open(outputfiles.xppOutputFile, 'r')
        o = open(outputfiles.simulationResults, 'w')
	o.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
        o.write(stringNodes)
        for line in f:
	    vector = line.split()
	    o.write('%s'%vector[0])
	    for i in range(len(vector)-1):
		o.write('\t%s'%vector[i+1])
            o.write('\n')
        f.close()
        o.close()
	#/ 2.- for the case of comma separated value format
        if 'csv' in metamodel.optionalOutputFormat:
            removeColumn = []
            f = open(outputfiles.xppOutputFile, 'r')
            o = open(outputfiles.simulationResultsCSV, 'w')            
            o.write('time')
            for s in stringNodes.split()[2:]:
                if s.startswith('delay') == False:
                    o.write(',' + s)
                else:
                    removeColumn.append(stringNodes.split()[2:].index(s))                    
            o.write('\n')
            for line in f:
                o.write(line.split()[0])
                for i in range(len(line.split())):
                    if i in removeColumn:
                        print 'Warnning: The column %s is not included in csv output file' %i
                    else:
                        o.write(',' + line.split()[i])
                o.write('\n')
            o.close()
            f.close()        
        
        return None
    
def convertName(name):

    '''
    This function prepares the format of the names suitable for XPPAUT format.
    Mainly string names of more than 7 letters are converted into a new one compoused of the first 5 and the last 2.
    Also underscores are removed.
    '''
    
    #/ 1.- string nodes of more than 7 strings are converted into a new string compoused of the first 5 and the last 2.
    if len(name) > 7:
        newName = name.upper()[:5] + name.upper()[len(name)-2:]
    else:
        newName = name.upper()
    #/ 2.- removing underscores
    if len(newName.split('_')) > 1:
        nameWithoutUnderscore = ''
        for n in newName.split('_'):
            nameWithoutUnderscore += n
        newName = nameWithoutUnderscore

    return newName

def replaceNames(expr, names, newNames, conflictives):
    
    '''
    This function replace from an expression the original names by the XPPAUT-adequate format.
    '''

    for i in range(len(names)):
        prog = re.compile('\A%s\s|\s%s\Z' %(names[i], names[i]))
        expr = prog.sub(' %s ' %newNames[i], expr)
        prog = re.compile('\A%s\Z' %(names[i]))
        expr = prog.sub(newNames[i], expr)
        prog = re.compile('/%s\Z' %(names[i]))
        expr = prog.sub('/ ' + newNames[i], expr)
        prog = re.compile('\*%s\Z' %(names[i]))
        expr = prog.sub('* ' + newNames[i], expr)
	#/ plausible replacing positions
        expr = expr.replace(' %s ' %names[i],' %s ' %newNames[i])    
        expr = expr.replace('(%s ' %names[i],'(%s ' %newNames[i])
        expr = expr.replace(' %s)' %names[i],' %s)' %newNames[i])
        expr = expr.replace(' %s/' %names[i],' %s/' %newNames[i])
        expr = expr.replace('/%s ' %names[i],'/%s ' %newNames[i])
        expr = expr.replace('*%s ' %names[i],'*%s ' %newNames[i])
        expr = expr.replace('(%s*' %names[i],'(%s*' %newNames[i])        
        expr = expr.replace('(%s,' %names[i],'(%s,' %newNames[i])
        expr = expr.replace(',%s)' %names[i],',%s)' %newNames[i])
        expr = expr.replace('(%s/' %names[i],'(%s/' %newNames[i])
        expr = expr.replace('/%s)' %names[i],'/%s)' %newNames[i])
	expr = expr.replace(' %s,' %names[i],' %s,' %newNames[i])

    return expr
