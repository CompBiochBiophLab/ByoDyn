#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana
#
#  Created: 2005-07-15 by Adrian L. Garcia-Lomana
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

# $Id: affectors.py,v 4.4 2008/12/03 16:37:38 alglomana Exp $

## \file
# This is the affectors module. 
# Here you will find all the possible biochemical reactions and what they do.
# This library is prepared to handle python, octave and latex formats. 

import centralFunctions, formulas #/ "formulas" module is at the lib directory

def complexExtraBack(model, file, option, cellIndex, definition, fieldsDefinition): 

    '''
    This function writes one of the mathematical terms (the other is done by complexExtraFwd) of a biological extracellular complex formation.
    This term affects the nodes defined as receptor and ligand.
    d[RECEPTOR]i / dt = - 1 / 64 k_binding_RECEPTOR.LIGAND [RECEPTOR]i <Sum>(j = 1;j = k) [LIGAND]j
    d[LIGAND]i / dt = - 1 / 64  k_binding_RECEPTOR.LIGAND [LIGAND]i <Sum>(j = 1;j = k) [RECEPTOR]j
    '''

    #/ 1. withdrawing information
    constant = model.topology[definition][0].split('/')[0]
    receptor = fieldsDefinition[2]
    receptorEquationIndex = model.nodes.index(receptor) + (len(model.nodes) * cellIndex)
    ligand = fieldsDefinition[3]
    ligandEquationIndex = model.nodes.index(ligand) + (len(model.nodes) * cellIndex)
    #/ 1.1. withdrawing more information: determining neighbours
    neighbourCellsPositions = [] 
    neighbourReceptorEquationIndexes = []
    neighbourLigandEquationIndexes = []
    neighbourCellsPositions = centralFunctions.neighboursFinder(model, cellIndex)
    #/ 1.2. determining the indexes of the neighbour nodes
    for i in range(len(neighbourCellsPositions)):
        neighbourReceptorEquationIndexes.append(model.nodes.index(receptor) + ((neighbourCellsPositions[i][0] * model.xlength + neighbourCellsPositions[i][1]) * len(model.nodes)))
        neighbourLigandEquationIndexes.append(model.nodes.index(ligand) + ((neighbourCellsPositions[i][0] * model.xlength + neighbourCellsPositions[i][1]) * len(model.nodes)))
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex':#/ the formula is different depending on being the receptor or the ligand
        if fieldsDefinition[0] == fieldsDefinition[2]: #/ d[RECEPTOR]i / dt = - 1 / 64 k_binding_RECEPTOR.LIGAND [RECEPTOR]i <Sum>(j = 1;j = k) [LIGAND]j
            file.write(' - \\frac{1}{64}\;%s_{\mathrm{%s}}^{\mathrm{%s}}\; [%s]_{i}\;\\sum_{j = 1}^{j = k} [%s]_{j}\\\\\\nonumber & &'%(constant.split('_')[0], constant.split('_')[1], constant.split('_')[2], receptor, ligand))
        elif fieldsDefinition[0] == fieldsDefinition[3]: #/ d[LIGAND]i / dt = - 1 / 64  K_binding_RECEPTOR.LIGAND [LIGAND]i <Sum>(j = 1;j = k) [RECEPTOR]j
            file.write(' - \\frac{1}{64}\;%s_{\mathrm{%s}}^{\mathrm{%s}}\; [%s]_{i}\;\\sum_{j = 1}^{j = k} [%s]_{j}\\\\\\nonumber & &'%(constant.split('_')[0], constant.split('_')[1], constant.split('_')[2], ligand, receptor))
    #/ 2.2.- OCTAVE
    elif option == 'octave':#/ the formula is different depending on being the receptor or the ligand
        if fieldsDefinition[0] == fieldsDefinition[2]: #/ d[RECEPTOR]i / dt = - 1 / 64 k_binding_RECEPTOR.LIGAND [RECEPTOR]i <Sum>(j = 1;j = k) [LIGAND]j
            file.write(' - %s * (1/64.0 * x(%s) * x(%s)' %(constant, receptorEquationIndex + 1, neighbourLigandEquationIndexes[0] + 1))
            for i in range(len(neighbourLigandEquationIndexes) - 1):
                file.write(' + 1/64.0 * x(%s) * x(%s)'%(receptorEquationIndex + 1, neighbourLigandEquationIndexes[i + 1] + 1))
            file.write(')')
        elif fieldsDefinition[0] == fieldsDefinition[3]:#/ d[LIGAND]i / dt = - 1 / 64  k_binding_RECEPTOR.LIGAND [LIGAND]i <Sum>(j = 1;j = k) [RECEPTOR]j
            file.write(' - %s * (1/64.0 * x(%s) * x(%s)' %(constant, ligandEquationIndex + 1, neighbourReceptorEquationIndexes[0] + 1))
            for i in range(len(neighbourReceptorEquationIndexes) - 1):
                file.write(' + 1/64.0 * x(%s) * x(%s)'%(ligandEquationIndex + 1, neighbourReceptorEquationIndexes[i + 1] + 1))
            file.write(')')
    #/ 2.3.- PYTHON
    elif option == 'python': #/ the formula is different depending on being the receptor or the ligand
        if fieldsDefinition[0] == fieldsDefinition[2]:#/ d[RECEPTOR]i / dt = - 1 / 64 k_binding_RECEPTOR.LIGAND [RECEPTOR]i <Sum>(j = 1;j = k) [LIGAND]j
            file.write(' - %s * (1/64.0 * x[%s] * x[%s]' %(constant, receptorEquationIndex, neighbourLigandEquationIndexes[0]))
            for i in range(len(neighbourLigandEquationIndexes) - 1):
                file.write(' + 1/64.0 * x[%s] * x[%s]'%(receptorEquationIndex, neighbourLigandEquationIndexes[i + 1]))
            file.write(')')
        elif fieldsDefinition[0] == fieldsDefinition[3]:#/ d[LIGAND]i / dt = - 1 / 64  k_binding_RECEPTOR.LIGAND [LIGAND]i <Sum>(j = 1;j = k) [RECEPTOR]j
            file.write(' - %s * (1/64.0 * x[%s] * x[%s]' %(constant, ligandEquationIndex, neighbourReceptorEquationIndexes[0]))
            for i in range(len(neighbourReceptorEquationIndexes) - 1):
                file.write(' + 1/64.0 * x[%s] * x[%s]'%(ligandEquationIndex, neighbourReceptorEquationIndexes[i + 1]))
            file.write(')')
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option
    
    return None

def complexExtraFwd(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This function writes one of the mathematical terms (the other is done by complexExtraBack) of biological extracellular complex formation.
    This term affects the node defined as the complex
    d[COMPLEX]i / dt = +  1 / 64 * k_binding_RECEPTOR.LIGAND [RECEPTOR]i  <Sum>(j = 1;j = k) [LIGAND]j
    '''

    #/ 1. withdrawing information
    constant = model.topology[definition][0].split('/')[0]
    receptor = fieldsDefinition[2]
    receptorEquationIndex = model.nodes.index(receptor) + (len(model.nodes) * cellIndex)
    ligand = fieldsDefinition[3]
    ligandEquationIndex = model.nodes.index(ligand) + (len(model.nodes) * cellIndex)
    #/ 1.1. withdrawing more information: determining neighbours
    neighbourCellsPositions = [] 
    neighbourLigandEquationIndexes = []
    neighbourLigandValues = []
    neighbourCellsPositions = centralFunctions.neighboursFinder(model, cellIndex)
    #/ 1.2. determining the indexes of the neighbour nodes
    for i in range(len(neighbourCellsPositions)):
        neighbourLigandEquationIndexes.append(model.nodes.index(ligand) + ((neighbourCellsPositions[i][0] * model.xlength + neighbourCellsPositions[i][1]) * len(model.nodes)))
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex': 
            file.write(' + \\frac{1}{64}\;%s_{\mathrm{%s}}^{\mathrm{%s}}\; [%s]_{i}\;\\sum_{j=1}^{j=k} [%s]_{j}\\\\\\nonumber & &'%(constant.split('_')[0], constant.split('_')[1], constant.split('_')[2], receptor, ligand))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
        file.write(' + %s * (1/64.0 * x(%s) * x(%s)' %(constant, receptorEquationIndex + 1, neighbourLigandEquationIndexes[0] + 1))
        for i in range(len(neighbourLigandEquationIndexes) - 1):
            file.write(' + 1/64.0 * x(%s) * x(%s)'%(receptorEquationIndex + 1, neighbourLigandEquationIndexes[i + 1] + 1))
        file.write(')')
    #/ 2.3.- PYTHON
    elif option == 'python':
        file.write(' + %s * (1/64.0 * x[%s] * x[%s]' %(constant, receptorEquationIndex, neighbourLigandEquationIndexes[0]))
        for i in range(len(neighbourLigandEquationIndexes) - 1):
            file.write(' + 1/64.0 * x[%s] * x[%s]'%(receptorEquationIndex, neighbourLigandEquationIndexes[i + 1]))
        file.write(')')
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option
    
    return None

def constant(model, file, option, cellIndex, definition, fieldsDefinition): 

    '''
    This function writes the mathematical term of a non changing biological node.
    d [node] / dt =  0
    '''

    #/ 1. withdrawing information
    rate = 0
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex': 
        file.write(' + %s\\\\\\nonumber & &'%(rate))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
        file.write(' + %s'%(rate))
    #/ 2.3.- PYTHON
    elif option == 'python':
        file.write(' + %s'%(rate))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option

    return None

def constitutive(model, file, option, cellIndex, definition, fieldsDefinition): 

    '''
    This function writes the mathematical term of a biological constitutive expression.
    d [node] / dt = constant > 0
    '''

    #/ 1. withdrawing information
    rate = model.topology[definition][0].split('/')[0]
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex': 
        file.write(' + %s^{\mathrm{%s}}\\\\\\nonumber & &'%(rate.split('_')[0], rate.split('_')[1]))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
        file.write(' + %s'%(rate))
    #/ 2.3.- PYTHON
    elif option == 'python':
        file.write(' + %s'%(rate))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option

    return None

def degradation(model, file, option, cellIndex, definition, fieldsDefinition): 

    '''
    This function writes the mathematical term of a biological degradation.
    d [node] / dt = -k * [node]
    '''

    #/ 1. withdrawing information
    constant =  model.topology[definition][0].split('/')[0]
    node = fieldsDefinition[2]
    nodeEquationIndex = model.nodes.index(node) + (len(model.nodes) * cellIndex)
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
        file.write(' - %s_{\mathrm{%s}}^{\mathrm{%s}}\;[%s]_{i}\\\\\\nonumber & &' %(constant.split('_')[0], constant.split('_')[1], constant.split('_')[2], node))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
        file.write(' - %s * x(%s)' %(constant, nodeEquationIndex + 1))
    #/ 2.3.- PYTHON
    elif option == 'python':
        file.write(' - %s * x[%s]' %(constant, nodeEquationIndex))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option

    return None

def dissociationExtraBack(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This function writes one of the mathematical terms (the other is done by dissociationExtraFwd) of a biological extracellular complex dissociation.
    This term affects the node defined as the complex.
    d[COMPLEX]i / dt = -  1 / 8 k_dissociation_COMPLEX <Sum>(j = 1;j = k) [COMPLEX]i
    '''

    #/ 1. withdrawing information
    constant = model.topology[definition][0].split('/')[0]
    complex = fieldsDefinition[0]
    complexEquationIndex = model.nodes.index(complex) + (len(model.nodes) * cellIndex)
    #/ 1.1. withdrawing more information: determining neighbours
    neighbourCellsPositions = [] 
    neighbourComplexEquationIndexes = []
    neighbourCellsPositions = centralFunctions.neighboursFinder(model, cellIndex)
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
        file.write(' - %s_{\mathrm{%s}}^{\mathrm{%s}}\;\\sum_{j = 1}^{j = k} \\frac{1}{8}\;[%s]_{i}\\\\\\nonumber & &'%(constant.split('_')[0], constant.split('_')[1], constant.split('_')[2], complex))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
        file.write(' - %s * (1/8.0 * x(%s)' %(constant, complexEquationIndex + 1))
        for i in range(len(neighbourComplexEquationIndexes) - 1):
            file.write(' - 1/8.0 * x(%s)' %(complexEquationIndex + 1)) 
        file.write(')')
    #/ 2.3.- PYTHON
    elif option == 'python':
        file.write(' - %s * (1/8.0 * x[%s]' %(constant, complexEquationIndex))
        for i in range(len(neighbourComplexEquationIndexes) - 1):
            file.write(' - 1/8.0 * x[%s]' %(complexEquationIndex)) 
        file.write(')')
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option

    return None

def dissociationExtraFwd(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This function writes one of the mathemathical terms (the other is done by dissociationExtraFwd) of a biological extracellular complex dissociation.
    This term affects the nodes defined as receptor and ligand.
    The formula is different depending on being the receptor or the ligand
    d[RECEPTOR]i / dt = + 1 / 8 k_dissociation_COMPLEX <Sum>(j = 1;j = k) [COMPLEX]i
    d[LIGAND]i / dt = + 1/8 k_dissociation_COMPLEX <Sum>(j = 1;j = k) [COMPLEX]j
    '''

    #/ 1. withdrawing information
    constant = model.topology[definition][0].split('/')[0]
    receptorActiv = fieldsDefinition[2]
    receptorActivEquationIndex = model.nodes.index(receptorActiv) + (len(model.nodes) * cellIndex)
    ligandActiv = fieldsDefinition[3]
    ligandActivEquationIndex = model.nodes.index(ligandActiv) + (len(model.nodes) * cellIndex)
    complex = fieldsDefinition[4]
    complexEquationIndex = model.nodes.index(complex) + (len(model.nodes) * cellIndex)
    #/ 1.1. withdrawing more information: determining neighbours
    neighbourCellsPositions = [] 
    neighbourComplexEquationIndexes = []
    neighbourComplexValues = []
    neighbourCellsPositions = centralFunctions.neighboursFinder(model, cellIndex)
    #/ 1.2. determining the indexes of the neighbour nodes
    for i in range(len(neighbourCellsPositions)):
        neighbourComplexEquationIndexes.append( model.nodes.index(complex) + ((neighbourCellsPositions[i][0] * model.xlength + neighbourCellsPositions[i][1]) * len(model.nodes)))
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex':#/ the formula is different depending on being the receptor or the ligand
        if fieldsDefinition[0] == fieldsDefinition[2]: #/ d[RECEPTOR]i / dt = + 1 / 8 k_dissociation_COMPLEX <Sum>(j = 1;j = k) [COMPLEX]i
            file.write(' + %s_{\mathrm{%s}}^{\mathrm{%s}}\;\\sum_{j = 1}^{j = k} \\frac{1}{8}\;[%s]_{i}\\\\\\nonumber & &'%(constant.split('_')[0], constant.split('_')[1], constant.split('_')[2], complex))
        elif fieldsDefinition[0] == fieldsDefinition[3]: #/ d[LIGAND]i / dt = + 1/8 k_dissociation_COMPLEX <Sum>(j = 1;j = k) [COMPLEX]j
            file.write(' + %s_{\mathrm{%s}}^{\mathrm{%s}}\;\\sum_{j = 1}^{j = k} \\frac{1}{8}\;[%s]_{j}\\\\\\nonumber & &'%(constant.split('_')[0], constant.split('_')[1], constant.split('_')[2], complex))
    #/ 2.2.- OCTAVE
    elif option == 'octave':#/ the formula is different depending on being the receptor or the ligand
        if fieldsDefinition[0] == fieldsDefinition[2]: #/ d[RECEPTOR]i / dt = + 1 / 8 k_dissociation_COMPLEX <Sum>(j = 1;j = k) [COMPLEX]i
            file.write(' + %s * (1/8.0 * x(%s)' %(constant, complexEquationIndex + 1))
            for i in range(len(neighbourComplexEquationIndexes) - 1):
                file.write(' + 1/8.0 * x(%s)' %(complexEquationIndex + 1)) 
            file.write(')')
        elif fieldsDefinition[0] == fieldsDefinition[3]: #/ d[LIGAND]i / dt = + 1/8 k_dissociation_COMPLEX <Sum>(j = 1;j = k) [COMPLEX]j
            file.write(' + %s * (1/8.0 * x(%s)' %(constant, neighbourComplexEquationIndexes[0] + 1))
            for i in range(len(neighbourComplexEquationIndexes) - 1):
                file.write(' + 1/8.0 * x(%s)' %(neighbourComplexEquationIndexes[i + 1] + 1))
            file.write(')')
    #/ 2.3.- PYTHON
    elif option == 'python': #/ the formula is different depending on being the receptor or the ligand
        if fieldsDefinition[0] == fieldsDefinition[2]:#/ d[RECEPTOR]i / dt = + 1 / 8 k_dissociation_COMPLEX <Sum>(j = 1;j = k) [COMPLEX]i
            file.write(' + %s * (1/8.0 * x[%s]' %(constant, complexEquationIndex))
            for i in range(len(neighbourComplexEquationIndexes) - 1):
                file.write(' + 1/8.0 * x[%s]' %(complexEquationIndex)) 
            file.write(')')
        elif fieldsDefinition[0] == fieldsDefinition[3]:#/ d[LIGAND]i / dt = + 1/8 k_dissociation_COMPLEX <Sum>(j = 1;j = k) [COMPLEX]j
            file.write(' + %s * (1/8.0 * x[%s]' %(constant, neighbourComplexEquationIndexes[0]))
            for i in range(len(neighbourComplexEquationIndexes) - 1):
                file.write(' + 1/8.0 * x[%s]' %(neighbourComplexEquationIndexes[i + 1]))
            file.write(')')
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option

    return None

def inhibition(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This function writes the mathematical term of a biological inhibition.
    d [PROTEIN] / dt = + k_basal / 1 + ([INHIBITOR]/k_inhibition)^s
    '''
    
    #/ 1. withdrawing information
    basalRate = model.topology[definition][0].split('/')[0]
    inhibitionConstant = model.topology[definition][1].split('/')[0]
    cooperativeCoefficient =  model.topology[definition][2].split('/')[0]
    inhibitor = fieldsDefinition[2]
    inhibitorEquationIndex = model.nodes.index(inhibitor) + (len(model.nodes) * cellIndex)
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
        file.write(' + \\frac{%s^{\mathrm{%s}}}{1 + \\left(\\frac{[%s]_{i}}{%s_{\mathrm{%s}}^{\mathrm{%s}}}\\right)^{%s}}\\\\\\nonumber & &'%(basalRate.split('_')[0], basalRate.split('_')[1], inhibitor, inhibitionConstant.split('_')[0], inhibitionConstant.split('_')[1], inhibitionConstant.split('_')[2], cooperativeCoefficient.split('_')[0]))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
        file.write(' + %s / (1 + (x(%s)/%s)^%s)'%(basalRate, inhibitorEquationIndex + 1, inhibitorConstant, cooperativeCoefficient))
    #/ 2.3.- PYTHON
    elif option == 'python':
        file.write(' + %s / (1 + (x[%s]/%s)**%s)'%(basalRate, inhibitorEquationIndex, inhibitionConstant, cooperativeCoefficient))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option

    return None

def SBML(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This function converts the SBML formula format into the integrators' format.
    It calls different functions from the formulas module.
    '''
    
    #/ 1. withdrawing information
    productNode = fieldsDefinition[0]
    reactantNode = fieldsDefinition[2]
    reactionIndex = fieldsDefinition[3]
    formula = model.topology['%s/SBML/%s/%s' %(productNode, reactantNode, reactionIndex)]
    #/ 2. writting different files
    #/ 2.1.- LATEX
    if option == 'latex':#/ the formula is different depending on being the reactant of the product
        formulaAST = model.topologyAST['%s/SBML/%s/%s' %(productNode, reactantNode, reactionIndex)]
        formulas.formulaLatex(formulaAST, file, model);
    #/ 2.2.- OCTAVE
    elif option == 'octave':
        formulas.readWriteFormula(model, file, cellIndex, option, formula)
    #/ 2.3.- PYTHON
    elif option == 'python':
        formulas.readWriteFormula(model, file, cellIndex, option, formula)
    elif option == 'openModelica':
        formulas.writeOpenModelicaFormula(formula, file)
    elif option == 'xpp':
        formulas.writeXPPFormula(formula, file)
    #/ 2.4 OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option

    return None

def transcription(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This function writes the mathematical term of a biological transcription.
    d[mRNA] / dt = v * [ACTIVATOR]^m / k^m + [ACTIVATOR]^m
    '''

    #/ 1. withdrawing information
    constant1 = model.topology[definition][0].split('/')[0]
    constant2 = model.topology[definition][1].split('/')[0]
    cooperativityCoefficient = model.topology[definition][2].split('/')[0]
    activator = fieldsDefinition[2]
    activatorEquationIndex = model.nodes.index(activator) + (len(model.nodes) * cellIndex)
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
        file.write(' + \\frac{%s^{\mathrm{%s}}\;[%s]_{i}^%s}{{%s_{\mathrm{%s}}^{%s} + [%s]_{i}^{%s}}}\\\\\\nonumber & &'%(constant1.split('_')[0], constant1.split('_')[1], activator, cooperativityCoefficient.split('_')[0],constant2.split('_')[0], constant2.split('_')[1], cooperativityCoefficient.split('_')[0], activator, cooperativityCoefficient.split('_')[0]))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
        file.write(' + (%s * x(%s)^%s) / (x(%s)^%s + %s^%s)' %(constant1, activatorEquationIndex + 1, cooperativityCoefficient, activatorEquationIndex + 1, cooperativityCoefficient, constant2, cooperativityCoefficient))
    #/ 2.3.- PYTHON
    elif option == 'python':
        file.write(' + (%s * x[%s]**%s) / (x[%s]**%s + %s**%s)' %(constant1, activatorEquationIndex, cooperativityCoefficient, activatorEquationIndex, cooperativityCoefficient, constant2, cooperativityCoefficient))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option

    return None

def translation(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This function writes the mathematical term of a biological translation.
    d[PROTEIN] / dt = + k_translation_gene * [gene]
    '''

    #/ 1. withdrawing information
    constant = model.topology[definition][0].split('/')[0]
    gene = fieldsDefinition[2]
    geneEquationIndex = model.nodes.index(gene) + (len(model.nodes) * cellIndex)
    #/ 2. writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
        file.write(' + %s_{\mathrm{%s}}^{\mathrm{%s}}\;[%s]_{i}\\\\\\nonumber & &' %(constant.split('_')[0], constant.split('_')[1], constant.split('_')[2], gene))  
    #/ 2.2.- OCTAVE
    elif option == 'octave':
        file.write(' + %s * x(%s)' %(constant, geneEquationIndex + 1))
    #/ 2.3.- PYTHON
    elif option == 'python':
        file.write(' + %s * x[%s]' %(constant, geneEquationIndex))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option

    return None

##########################################
# From here below, non dimension affectors
##########################################

def NonDimConstitutiveDegradation(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This affector compiles the constitutive gene expression and degradation.
    dnode/dtau = k_deg*(rate - node)
    '''

    #/ 1.- withdrawing information
    degradationConstant = model.topology[definition][0].split('/')[0]
    constitutiveConstant = model.topology[definition][1].split('/')[0]
    node = fieldsDefinition[2]
    nodeEquationIndex = model.nodes.index(node) + (len(model.nodes) * cellIndex)
    #/ 2.- writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
	file.write(' T_0\;%s_{\mathrm{%s}}^{\mathrm{%s}}\;\left(%s_{\mathrm{%s}}-%s_{i}(\\tau)\\right)\\\\\\nonumber & &'%(degradationConstant.split('_')[0],degradationConstant.split('_')[1],degradationConstant.split('_')[2],constitutiveConstant.split('_')[0],constitutiveConstant.split('_')[1],node))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
	file.write('%s * (%s - x(%s))' %(degradationConstant, constitutiveConstant, nodeEquationIndex+1))
    #/ 2.3.- PYTHON
    elif option == 'python':
        file.write('%s * (%s - x[%s])' %(degradationConstant, constitutiveConstant, nodeEquationIndex))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
        raise 'The integration format file %s is not supported.'%option
    
    return None

def NonDimTranslationDegradationBinding(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This affector compiles the translation, degradation and ligand binding.
    dNODE/dtau = k_deg(node - NODE) - (k/n^2) * k_bind [LIGAND]_0 RECEPTOR_i \sum LIGAND_j
    '''
    
    #/ 1.- withdrawing information
    degradationConstant = model.topology[definition][0].split('/')[0]
    bindingConstant = model.topology[definition][1].split('/')[0]
    characteristicConcentration = model.topology[definition][2].split('/')[0]
    receptor = fieldsDefinition[0]
    gene = fieldsDefinition[2]
    ligand = fieldsDefinition[3]
    geneEquationIndex = model.nodes.index(gene) + (len(model.nodes) * cellIndex)
    receptorEquationIndex = model.nodes.index(receptor) + (len(model.nodes) * cellIndex)
    #/ 1.1. withdrawing more information: determining neighbours
    neighbourCellsPositions = [] 
    neighbourLigandEquationIndexes = []
    neighbourLigandValues = []
    neighbourCellsPositions = centralFunctions.neighboursFinder(model, cellIndex)
    #/ 1.2. determining the indexes of the neighbour nodes
    for i in range(len(neighbourCellsPositions)):
        neighbourLigandEquationIndexes.append(model.nodes.index(ligand) + ((neighbourCellsPositions[i][0] * model.xlength + neighbourCellsPositions[i][1]) * len(model.nodes)))
    #/ 2.- writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
	file.write('T_0\;%s_{\mathrm{%s}}^{\mathrm{%s}}\;\left(%s_{i}(\\tau)-%s_{i}(\\tau)\\right)\\\\\\nonumber & &-\\frac{k}{n^2}\;%s_{\mathrm{%s}}^{\mathrm{%s}}\;[%s]_%s\;%s_{i}(\\tau)\;\\sum_{j=1}^{k}%s_{j}(\\tau)\\\\\\nonumber & &'%(degradationConstant.split('_')[0],degradationConstant.split('_')[1],degradationConstant.split('_')[2],gene,receptor,bindingConstant.split('_')[0],bindingConstant.split('_')[1],bindingConstant.split('_')[2],characteristicConcentration.split('_')[0],characteristicConcentration.split('_')[1],receptor,ligand))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
	file.write('%s * (x(%s) - x(%s)) - (%s/64.)*%s*%s*x(%s)*(x(%s)'%(degradationConstant,geneEquationIndex+1,receptorEquationIndex+1,len(neighbourCellsPositions),bindingConstant,characteristicConcentration,receptorEquationIndex+1,neighbourLigandEquationIndexes[0]+1))
	for i in range(len(neighbourLigandEquationIndexes)-1):
	    file.write('+x(%s)'%(neighbourLigandEquationIndexes[i+1]+1))
	file.write(')')
    #/ 2.3.- PYTHON
    elif option == 'python':
	file.write('%s * (x[%s] - x[%s]) - (%s/64.)*%s*%s*x[%s]*(x[%s]'%(degradationConstant,geneEquationIndex,receptorEquationIndex,len(neighbourCellsPositions),bindingConstant,characteristicConcentration,receptorEquationIndex,neighbourLigandEquationIndexes[0]))
	for i in range(len(neighbourLigandEquationIndexes)-1):
	    file.write('+x[%s]'%(neighbourLigandEquationIndexes[i+1]))
	file.write(')')
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
	raise 'The integration format file %s is not supported.'%option

    return None

def NonDimInhibitionDegradation(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This affector compiles an inhibition and a degradation.
    dnode/dtau = k_deg * (1 - (INHIBITOR^m/kappa^m+INHIBITION^m) - node)
    '''

    #/ 1.- withdrawing information
    degradationConstant = model.topology[definition][0].split('/')[0]
    kappa = model.topology[definition][1].split('/')[0]
    cooperativityCoefficient = model.topology[definition][2].split('/')[0]
    gene = fieldsDefinition[0]
    inhibitor = fieldsDefinition[2]
    geneEquationIndex = model.nodes.index(gene) + (len(model.nodes) * cellIndex)
    inhibitorEquationIndex = model.nodes.index(inhibitor) + (len(model.nodes) * cellIndex)
    #/ 2.- writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
	file.write('T_0\;%s_{\mathrm{%s}}^{\mathrm{%s}}\;\left(1-\\frac{%s_{i}^{%s}(\\tau)}{\\kappa_{\mathrm{%s}}+%s_{i}^{%s}(\\tau)}-%s_{i}(\\tau)\\right)\\\\\\nonumber & &'%(degradationConstant.split('_')[0],degradationConstant.split('_')[1],degradationConstant.split('_')[2],inhibitor,cooperativityCoefficient,kappa.split('_')[1],inhibitor,cooperativityCoefficient,gene))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
	file.write('%s * (1 - x(%s)^%s/(%s^%s+x(%s)^%s) - x(%s))'%(degradationConstant,inhibitorEquationIndex+1,cooperativityCoefficient,kappa,cooperativityCoefficient,inhibitorEquationIndex+1,cooperativityCoefficient,geneEquationIndex+1))
    #/ 2.3.- PYTHON
    elif option == 'python':
	file.write('%s * (1 - x[%s]**%s/(%s**%s+x[%s]**%s) - x[%s])'%(degradationConstant,inhibitorEquationIndex,cooperativityCoefficient,kappa,cooperativityCoefficient,inhibitorEquationIndex,cooperativityCoefficient,geneEquationIndex))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
	raise 'The integration format file %s is not supported.'%option

    return None

def NonDimBindingDegradation(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This affector compiles an complex binding and a degradation.
    dcomplex/dtau = (k/n^2) * k_bind [LIGAND]_0 RECEPTOR_i \sum LIGAND_j - k_deg COMPLEX_i
    '''
    
    #/ 1.- withdrawing information
    degradationConstant = model.topology[definition][0].split('/')[0]
    bindingConstant = model.topology[definition][1].split('/')[0]
    characteristicConcentration = model.topology[definition][2].split('/')[0]
    complex = fieldsDefinition[0]
    receptor = fieldsDefinition[2]
    ligand = fieldsDefinition[3]
    complexEquationIndex = model.nodes.index(complex) + (len(model.nodes) * cellIndex)
    receptorEquationIndex = model.nodes.index(receptor) + (len(model.nodes) * cellIndex)
    #/ 1.1.- withdrawing more information: determining neighbours
    neighbourCellsPositions = [] 
    neighbourLigandEquationIndexes = []
    neighbourLigandValues = []
    neighbourCellsPositions = centralFunctions.neighboursFinder(model, cellIndex)
    #/ 1.2.- determining the indexes of the neighbour nodes
    for i in range(len(neighbourCellsPositions)):
        neighbourLigandEquationIndexes.append(model.nodes.index(ligand) + ((neighbourCellsPositions[i][0] * model.xlength + neighbourCellsPositions[i][1]) * len(model.nodes)))
    #/ 2.- writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
	file.write('T_0\;\left(\\frac{k}{n^2}\;%s_{\mathrm{%s}}^{\mathrm{%s}}\;[%s]_%s\;%s_{i}(\\tau)\;\sum_{j=1}^{k}%s_{j}(\\tau)-%s_{\mathrm{%s}}^{\mathrm{%s}}\;%s_{i}(\\tau)\\right)'%(bindingConstant.split('_')[0],bindingConstant.split('_')[1],bindingConstant.split('_')[2],characteristicConcentration.split('_')[0],characteristicConcentration.split('_')[1],receptor,ligand,degradationConstant.split('_')[0],degradationConstant.split('_')[1],degradationConstant.split('_')[2],complex))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
	file.write('(%s/64.)*%s*%s*x(%s)*(x(%s)'%(len(neighbourCellsPositions),bindingConstant,characteristicConcentration,receptorEquationIndex+1,neighbourLigandEquationIndexes[0]+1))
	for i in range(len(neighbourLigandEquationIndexes)-1):
	    file.write('+x(%s)'%(neighbourLigandEquationIndexes[i+1]+1))
	file.write(')')
	file.write('-%s*x(%s)'%(degradationConstant,complexEquationIndex+1))
    #/ 2.3.- PYTHON
    elif option == 'python':
	file.write('(%s/64.)*%s*%s*x[%s]*(x[%s]'%(len(neighbourCellsPositions),bindingConstant,characteristicConcentration,receptorEquationIndex,neighbourLigandEquationIndexes[0]))
	for i in range(len(neighbourLigandEquationIndexes)-1):
	    file.write('+x[%s]'%neighbourLigandEquationIndexes[i+1])
	file.write(')')
	file.write('-%s*x[%s]'%(degradationConstant,complexEquationIndex))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
	raise 'The integration format file %s is not supported.'%option
    
    return None

def NonDimTranscriptionDegradation(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This affector compiles transcription and degradation.
    dnode/dtau = k_deg(ACTIVATOR^m/(kappa^m+ACTIVATOR^m) - node)
    '''

    #/ 1.- withdrawing information
    degradationConstant = model.topology[definition][0].split('/')[0]
    kappa = model.topology[definition][1].split('/')[0]
    cooperativityCoefficient = model.topology[definition][2].split('/')[0]
    gene = fieldsDefinition[0]
    activator = fieldsDefinition[2]
    geneEquationIndex = model.nodes.index(gene) + (len(model.nodes) * cellIndex)
    activatorEquationIndex = model.nodes.index(activator) + (len(model.nodes) * cellIndex)
    #/ 2.- writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
	file.write('T_0\;%s_{\mathrm{%s}}^{\mathrm{%s}}\;\left(\\frac{%s_{i}^{%s}(\\tau)}{\\kappa_{\mathrm{%s}}+%s_{i}^{%s}(\\tau)}-%s_{i}(\\tau)\\right)'%(degradationConstant.split('_')[0],degradationConstant.split('_')[1],degradationConstant.split('_')[2],activator,cooperativityCoefficient,kappa.split('_')[1],activator,cooperativityCoefficient,gene))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
	file.write('%s * (x(%s)^%s/(%s^%s+x(%s)^%s) - x(%s))'%(degradationConstant,activatorEquationIndex+1,cooperativityCoefficient,kappa,cooperativityCoefficient,activatorEquationIndex+1,cooperativityCoefficient,geneEquationIndex+1))
    #/ 2.3.- PYTHON
    elif option == 'python':
	file.write('%s * (x[%s]**%s/(%s**%s+x[%s]**%s) - x[%s])'%(degradationConstant,activatorEquationIndex,cooperativityCoefficient,kappa,cooperativityCoefficient,activatorEquationIndex,cooperativityCoefficient,geneEquationIndex))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
	raise 'The integration format file %s is not supported.'%option

    return None

def NonDimTranslationDegradation(model, file, option, cellIndex, definition, fieldsDefinition):

    '''
    This affector compiles translation and degradation.
    dnode/dtau = k_deg(gene -PROTEIN)
    '''
    
    #/ 1.- withdrawing information
    degradationConstant = model.topology[definition][0].split('/')[0]
    gene = fieldsDefinition[2]
    protein = fieldsDefinition[0]
    geneEquationIndex = model.nodes.index(gene) + (len(model.nodes) * cellIndex)
    proteinEquationIndex = model.nodes.index(protein) + (len(model.nodes) * cellIndex)
    #/ 2.- writing different files
    #/ 2.1.- LATEX
    if option == 'latex':
	file.write('T_0\;%s_{\mathrm{%s}}^{\mathrm{%s}}\;\left(%s_{i}(\\tau)-%s_{i}(\\tau)\\right)'%(degradationConstant.split('_')[0],degradationConstant.split('_')[1],degradationConstant.split('_')[2],gene,protein))
    #/ 2.2.- OCTAVE
    elif option == 'octave':
	file.write('%s * (x(%s) - x(%s))'%(degradationConstant,geneEquationIndex+1,proteinEquationIndex+1))
    #/ 2.3.- PYTHON
    elif option == 'python':
	file.write('%s * (x[%s] - x[%s])'%(degradationConstant,geneEquationIndex,proteinEquationIndex))
    #/ 2.4.- OR THE FORMAT IS NOT CORRECT
    else:
	raise 'The integration format file %s is not supported.'%option

    return None
