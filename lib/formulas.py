#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez Garrido and Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Alex Gomez Garrido and Adrian L. Garcia-Lomana
#
#  Created: 2006-01-30 by Alex Gomez Garrido
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

# $Id: formulas.py,v 4.13 2008/12/03 16:37:38 alglomana Exp $

## \file
# This module contains different functions necessary to interconvert the formula string formats.

import sys, string, re, copy
import libsbml
from math import *
from scipy import piecewise as pcw #/ piecewise function is used in MathML definition
#/ Importing comparison functions and changing the name to the MathML name
from scipy import less as lt
from scipy import less_equal as leq
from scipy import greater as gt
from scipy import greater_equal as geq
from scipy import equal as eq
from scipy import not_equal as neq 
#/ ByoDyn modules
import errorMessages

def formatPowers(formula, option):

    '''
    This function checks for the powers and it substitute 'pow(k, n)' to 'k**n'.
    '''

    powList = re.findall('pow\([\w\[\]()/\+\-\*\s\d\.\^]*,[\w\[\]()/\+\-\*\s\d\.\^]*\)', formula)
    for pow in powList:
        thePower = re.split(',', pow)
        b = re.split('pow', thePower[0])
        if len(b) > 2:
            insidePow = b[1] + 'pow' + b[2] + ',' + thePower[1]
            newInsidePow = formatPowers(insidePow, option)
            formula = formula.replace(insidePow, newInsidePow)
        powList2 = re.findall('pow\([\w\[\]()/\+\-\*\s\d\.\^]*,[[\w\[\]()/\+\-\*\s\d\.\^]*\)', formula)
        if len(powList2) != 0:
            pow = powList2[0]
            thePower = re.split(',', pow)
            b = re.split('pow', thePower[0])
            basis = b[1] + ')'
            exponent = thePower[1][1:]
            if option == 'octave':
                newFormat = '%s**(%s' %(basis, exponent)
            elif option == 'openModelica':
                newFormat = '%s^(%s' %(basis, exponent)
            formula = formula.replace(pow, newFormat)

    return formula

       
def translateMathFactor(mathFactor):

    '''
    This function converts operators to python format.
    '''

    if mathFactor == 'gt':    return '>'
    elif mathFactor == 'lt':    return '<'
    elif mathFactor == 'eq':    return '=='
    elif mathFactor == 'neq':    return '<>'
    elif mathFactor == 'geq':    return '>='
    elif mathFactor == 'leq':    return '<='
    elif mathFactor == 'or':    return 'or'
    elif mathFactor == 'and':    return 'and'
    elif mathFactor == 'plus':    return '+'
    elif mathFactor == 'minus':    return '-'
    elif mathFactor == 'power':    return '^'
    else:
        return mathFactor.replace('_', '\_')

def getMathExpression(ASTNode):

    '''
    This function return a string with the math expression in latex syntax and the operator of a ASTNode.
    ASTNode is a object codifying a formula using a abstract syntax tree object.
    It does recursive calls to itself in order to build the math expression.
    '''

    numChildren = ASTNode.getNumChildren()
    if (numChildren == 0):    
        if (ASTNode.isReal()):    return str(ASTNode.getReal()), None
        if (ASTNode.isInteger()):    return str(ASTNode.getInteger()), None
        name = ASTNode.getName()
        if name.find('_') != -1:
            name = name.replace('_', '\_')
        return name, None;
    eq = ""                                                       
    if ASTNode.isOperator():
        operator = ASTNode.getCharacter()
    else:
        operator = translateMathFactor(ASTNode.getName())
    if operator=="/":
        eq = "\\frac"
    for i in range (0, numChildren, 1):
        subeq, subop = getMathExpression(ASTNode.getChild(i))
        if (subop == None):
            if (operator=="/"):    eq += "{" + subeq + "}"
            elif operator == '^':    eq += ("{" + subeq + "}" + operator)
            else:
                if len(operator) < 2:
                    eq += (subeq + operator)
                else:
                    if eq != '':
                        eq += (', ' + subeq)
                    else:
                        eq += subeq
        else:
            if (ASTNode.getPrecedence() >= ASTNode.getChild(i).getPrecedence()):
                if (operator=="/"):    eq += "{" + subeq + "}"
                elif (operator=="^"):    eq += ("{"+subeq+"}"+operator)
                elif operator == 'exp': eq += subeq
                else:
                    if len(operator) < 2:
                        eq += (subeq + operator)
                    else:
                        if eq != '':
                            eq += (', ' + subeq)
                        else:
                            eq += subeq
            else:
                if (operator=="/"):    eq += "{\\left(" + subeq + "\\right)}"
                elif (operator=="^"):    eq += ("{\\left("+subeq+"\\right)}"+operator)
                else: 
                    eq += ("\\left("+subeq+"\\right)"+operator)
    if (operator != "/" and len(operator) < 2):
        eq = eq[0:len(eq)-len(operator)]
    elif len(operator) > 2:
         eq = '%s\\left(%s\\right)' %(operator, eq)
    return eq, operator


def formulaLatex(xmlFormula, file, model):

    '''
    This function converts the MathML of SBML into latex format.
    '''

    ASTNode = libsbml.readMathMLFromString(xmlFormula[1])
    formula, op = getMathExpression(ASTNode)
    if xmlFormula[2] != '':
        compartmentLatex = xmlFormula[2].replace('_', '\_')
        formula = xmlFormula[0] + ' \\frac{' + formula + '}{'+ compartmentLatex +'}'
        formula = formula.replace('}{/ ','}{')
    else:
        formula = xmlFormula[0]+ ' ' + formula
    formula = formula + ' \\nonumber \\\\\\nonumber & &'    
    file.write('%s '%formula)

    return None

def piecewise(a,b,c):
    '''
    This function converts the format of piecewise function in MathML to the piecewise function in Scipy. 
    '''
    return pcw([a],b,[c])[0]
    
def readWriteFormula(model, file, cellIndex, option, formula):

    '''
    This funtion inputs the SBML string of the formula and writes out on Python or Octave formats.
    '''

    formula = checkBrackets(formula)
    terms = re.split('[\W]*', formula)
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
    for node in variables:
        for term in terms:
            if term == node:
                if option == 'python':
                    nodeEquation = 'x[%s]' % str(variables.index(node) + (len(variables) * cellIndex))
                elif option == 'octave':                   
                    nodeEquation = 'x(%s)' % str(variables.index(node) + (len(variables) * cellIndex) + 1)
                detailFormula = re.split('[\s]', formula)      
                i = 0
                for detail in detailFormula:
                    detail2 = re.split('[\-\+\*\\(),/]', detail)
                    if detail2.count(node) != 0:
                        if detail2.count(node) != detailFormula[i].count(node):
                            detailFormula[i] = detailFormula[i].replace(node, nodeEquation, detail2.count(node))             
                        else:
                            detailFormula[i] = detailFormula[i].replace(node, nodeEquation)
                    i = i + 1
                formula = string.join(detailFormula)
                break
    #/ replacing other variables known: Time = t
    if formula.count('Time') != 0:
        formula = formula.replace('Time', 't')
    file.write('%s '%formula)

    return None
    
def writeOpenModelicaFormula(formula, file):
    
    '''
    This function inputs the SBML string of the formula and writes out on OpenModelica format.
    '''

    formula = checkBrackets(formula)
    if formula.count('Time') != 0:
	formula = formula.replace('Time', 'time')
    if formula.count(' t ') != 0:
        formula = formula.replace(' t ', ' time ')
    if formula.count('(t ') != 0:
        formula = formula.replace('(t ', '(time ')
    if formula.count(' t)') != 0:
        formula = formula.replace(' t)', ' time)')
    if formula.count('(t,') != 0:
        formula = formula.replace('(t,', '(time,')
    if formula.count('flow') !=0:
        formula = formula.replace('flow', 'flo') #/ because 'flow' is a restricted character in Modelica
    if formula.count('+ -') != 0:
        formula = formula.replace('+ -', '- ')
    if formula.count('* -') != 0:
        formula = formula.replace('* -', '* (-1) * ')        
    file.write('%s '%formula)

    return None

def writeXPPFormula(formula, file):

    '''
    This function inputs the SBML string of the formula and writes out on XPP format.
    '''

    if formula.count('+ -') != 0:
        formula = formula.replace('+ -', '- ')
    file.write('%s '%formula)

    return None

def replaceConstants(model, constant):

    '''
    This functions replaces a parameter for its value or a variable for its initial condition value.
    '''

    value = ''
    if model.nodes.count(constant) != 0:
        value = str(model.initialConditions[model.nodes.index(constant)])
    elif model.algebraicNodes.has_key(constant) == True:
        value = str(model.algebraicNodes[constant])
    elif model.parameters.has_key(constant) == True:
        value = str(model.parameters[constant])
    elif model.nonConstantParameters.has_key(constant) == True:
        value = str(model.nonConstantParameters[constant])
    elif model.compartments.has_key(constant) == True:
        value = str(model.compartments[constant]['size'])
    elif model.nonConstantCompartments.has_key(constant) == True:
        value = str(model.nonConstantCompartments[constant]['size'])
    if value == '':
        value = constant      

    return value
        
def solveFormula(model, formula):

    '''
    This function evaluates a algebraic string, replacing all the parameters and variables for its values in time 0 and returns the value.
    '''

    terms = re.split('\s', formula)
    newFormula = ''
    for term in terms:  
        details = re.split('[(,)\-\+]', term)
        for detail in details:
            constant = replaceConstants(model, detail)
            term = term.replace(detail, constant)    
        newFormula = newFormula + term + ' '    
    newFormula = includeFunctions(model, newFormula)
    if newFormula.count('Time') != 0:        
        newFormula = newFormula.replace('Time', '0')
    if newFormula.count('time') != 0:        
        newFormula = newFormula.replace('time', '0')
    if newFormula.count(' t ') != 0:        
        newFormula = newFormula.replace(' t ', ' 0 ')
    if newFormula.count('(t ') != 0:
        newFormula = newFormula.replace('(t ', '(0 ')
    if newFormula.count(' t)') != 0:
        newFormula = newFormula.replace(' t)', ' 0)')
    if newFormula.count('(t,') != 0:
        newFormula = newFormula.replace('(t,', '(0,')
    newFormula = checkBrackets(newFormula) 
    try:
        value = eval(newFormula)
    except ZeroDivisionError:
        value = 0.
    return value

def includeFunctions(model, formula):

    '''
    This function get a formula and return a new formula without external functions used.
    A external function is a non-python function like the user-defined functions in SBML. 
    It replaces the function definitions with the function output in order to create a string without external functions.
    '''
    
    for function in model.functions:
        #/ Searching definitions of this function in formula.
        definitions = re.findall('%s\([\w\[\]/\+\-\*\s\d\.\^\,]*\)' %function.id, formula)
        #/ For each definition translate to python syntax
        if len(definitions) != 0:
            #/ Getting arguments and modifying function output
            argumentsString = definitions[0][len(function.id)+1:len(definitions[0])-1]
            arguments = []
            for argument in argumentsString.split(','):
                arguments.append(argument.strip())
            output = function.output
            for i in range(len(arguments)):
                output = output.replace(function.arguments[i], arguments[i])
            #/ Searching square root definitions like root(2, X) and changing to sqrt(X)
            squareRootDefinitions = re.findall('root\(2, [\w\[\]()/\+\-\*\s\d\.\^\,]*\)', output)
            if len(squareRootDefinitions) != 0:
                output = output.replace('root(2, ', 'sqrt(')
            #/ Recursive calls because sometimetimes inside the function output there are used other functions.
            output = includeFunctions(model, output)
            #/ Replacing the function definitions with the function output in order to create a string without external functions.
            formula = formula.replace(definitions[0], output)
        
    return formula

def checkBrackets(formula):

    '''
    This function checks if the use of brackets is correct in a math expression.
    It returns the formula modified if the number of brackets is not correct and prints a warning in the standard output.
    '''

    while formula.count('(') != formula.count(')'):
        print 'WARNING: Problems with the brackets in rule math expression. We are adding a close bracket.'
        print 'original Formula: ', formula
        formula = formula + ')'
        print 'new Formula: ', formula
        
    return formula
