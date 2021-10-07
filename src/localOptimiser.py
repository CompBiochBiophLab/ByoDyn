#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana
#
#  Created: 2005-11-15 by Adrian L. Garcia-Lomana
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

## \file
#  This module minimises the fitness function based on the Fortran program dn2fb from the PORT Mathematical Subroutine Library

# $Id: localOptimiser.py,v 4.6 2008/12/03 16:37:38 alglomana Exp $

import os, scipy, math, copy, sys
import libsbml
import centralFunctions, optimiser, errorMessages, initiator, sbmlWorker

def central(model, metamodel, outputfiles):

    '''
    This is the main function of the module.
    First it appends the compiling directory to the system path.
    Then it makes global the three main objects of the program (model, metamodel and outputfiles) and the module localSearch.
    It calls the optimisation routine, obtains the new value of the parameters and the fitness function value.
    It updates the model and sends the very last fitness function.
    '''

    #/ 1.- defining variables as global, necessary for the python function called by fortran
    threshold = globalDefiner(model, metamodel, outputfiles)
    #/ 2.- running the actual optimisation
    value, popt = optimisation(threshold)
    #/ 3.- recovering the calibrated model
    #/ 3.1.- in the case that all free parameters are model parameters
    if len(popt) == len(metamodel.parametersToVary):
	for i in range(len(metamodel.parametersToVary)):
	    model.parameters[metamodel.parametersToVary[i].split('/')[0]] = popt[i]
    #/ 3.2.- in the case of free initial conditions
    else:
	for i in range(len(metamodel.parametersToVary)):
	    model.parameters[metamodel.parametersToVary[i].split('/')[0]] = popt[i]
	for i in range(len(metamodel.initialConcentrationsToVary)):
	    node = metamodel.initialConcentrationsToVary[i].split('/')[0]
	    for j in range(len(model.nodes)):
		if node == model.nodes[j]:
		    model.initialConditions[j] = popt[len(metamodel.parametersToVary)+i]
    #/ 4.- in the case of a pure local search, store the results
    if metamodel.optimisationMethod != 'hybridOnePhase':
        savingResults(value, model, metamodel, outputfiles)

    return value

def globalDefiner(model, metamodel, outputfiles):

    '''
    This function adds the compilation directory to the system path.
    Then it makes global the model, metamodel and outputfiles objects and the localSearch fortran function.
    It finally defines the metamodel.hybridScore variable as True for the correct calculation of the score: the Fortran function returns a list and Python functions work with a float.
    '''

    locatingModuleDir(metamodel)
    global MODEL, METAMODEL, OUTPUTFILES, LOCAL_SEARCH #/ global variables are expressed with capital characters
    import localSearch as LOCAL_SEARCH
    MODEL = model
    METAMODEL = metamodel
    OUTPUTFILES = outputfiles
    METAMODEL.hybridScore = True
    #/ 1.- dealing with the score thresholds
    if metamodel.stopper[0] == 'score':
	threshold = float(metamodel.stopper[1])
    else:
	threshold = 1e-20

    return threshold

def locatingModuleDir(metamodel):

    '''
    This function appends the compiling directory to the system path.
    '''

    sys.path.append('%s'%metamodel.hybridPath)

    return None

def optimisation(threshold): #/ default stopping fitness function value is 1e-20.

    '''
    This function interacts with the Fortran code directly.
    It calls the main routine of the Fortran program.
    It sends the variables: x for the parameter values and b for the parameter ranges. And optionally it can send the value for the stopper threshold (v31).
    Finally it collects the optimised parameter vector and the fitness function value.
    '''
    print 'Local optimisation running ...'
    #/ 1.- preparing some variables necessary for the fortran routine ...
    x = []; b = [[], []] 
    s = threshold
    for parameter in METAMODEL.parametersToVary:
        x.append(MODEL.parameters[parameter.split('/')[0]])
        b[0].append(float(parameter.split('/')[1]))
        b[1].append(float(parameter.split('/')[2]))  
    #/ 1.1.- in case of initial conditions optimisation
    for ic in METAMODEL.initialConcentrationsToVary:
	node = ic.split('/')[0]
	for i in range(len(MODEL.nodes)):
	    if node == MODEL.nodes[i]:
		x.append(MODEL.initialConditions[i])
		#/ 1.1.1.- setting arbitrarily large numbers
		b[0].append(1e-200)
		b[1].append(1e200)
    #/ 2.- calling the fortran program
    LOCAL_SEARCH.principal(x, b, s, pcalcr) 
    #/ 3.- retrieving the values
    popt = LOCAL_SEARCH.xxshared.xx
    value = LOCAL_SEARCH.vvshared.vv
    
    return value, popt

def pcalcr(x, p):

    '''
    This function is called by the Fortran program localSearch.
    This function obtains the fitness function value for the new parameter position given by the Fortran program.
    '''
    
    #/ 1.- in the case all free parameters are model parameters
    if len(x) == METAMODEL.parametersToVary:
	for i in range(len(x)):
	    MODEL.parameters[METAMODEL.parametersToVary[i].split('/')[0]] = x[i]
    #/ 2.- in the case some free parameters are model initial conditions
    else:
	#/ 2.1.- filling up the model parameters
	for i in range(len(METAMODEL.parametersToVary)):
	    MODEL.parameters[METAMODEL.parametersToVary[i].split('/')[0]] = x[i]
	#/ 2.2.- setting the new initial condition defined by the fortran routine
	for i in range(len(METAMODEL.initialConcentrationsToVary)):
	    node = METAMODEL.initialConcentrationsToVary[i].split('/')[0]
	    for j in range(len(MODEL.nodes)):
		if node == MODEL.nodes[j]:
		    MODEL.initialConditions[j] = x[len(METAMODEL.parametersToVary)+i]
    #/ 3.- obtaining the score and placing it on a common place for fortran
    LOCAL_SEARCH.rrshared.rr = optimiser.scoreObtainer(METAMODEL, MODEL, OUTPUTFILES)
    
    return None

def savingResults(value, model, metamodel, outputfiles):

    '''
    This function saves the results of the local optimisation.
    '''
    
    #/ 1.- saving the parameters
    if metamodel.optimisationMethod == 'localSearch':
        if os.path.exists(outputfiles.solutionDirectory) == False:
            os.mkdir(outputfiles.solutionDirectory)
        parametersDirectory = outputfiles.solutionDirectory + '/parametersFromLocalSearch'
    else:
        parametersDirectory = outputfiles.solutionDirectory + '/parametersFromHybridTwoPhases'
    if os.path.exists(parametersDirectory) == False:
        os.mkdir(parametersDirectory)
    for parameter in metamodel.parametersToVary:
	file = parametersDirectory + '/%s' %parameter.split('/')[0]
        #/ 1.1.- writing the header for the first time the file is created
	if os.path.exists(file) == False:
	    f = open(file, 'w')
	    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
	    f.close()
	#/ 1.2.- appending information to the files
	f = open(file, 'a')
	f.write('%s\t%s\t%s\t%s\n'%(model.parameters[parameter.split('/')[0]], parameter.split('/')[1], parameter.split('/')[2], value))     
	f.close()
    #/ 2.- saving the initial conditions, in case of
    if metamodel.initialConcentrationsToVary != []:
	if metamodel.optimisationMethod == 'localSearch':
	    file = outputfiles.solutionDirectory + '/initialConditionsFromLocalSearch' 
	else:
	    file = outputfiles.solutionDirectory + '/initialConditionsFromHybridTwoPhases'
	#/ 2.1.- writing the header for the first time the file is created
	if os.path.exists(file) == False:
	    f = open(file, 'w')
	    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
	    f.close()
	#/ 2.2.- appending information to the files
	f = open(file, 'a')
	for ic in model.initialConditions:
	    f.write('%s\t'%ic)  
	f.write('\n')
	f.close()
	#/ 2.3.- writing the score of the fitness function in the case of absence of parameters to vary
	if metamodel.parametersToVary == []:
	    file = outputfiles.solutionDirectory + '/scoresForInitialConcentrationsFromLocalSearch'
	    #/ 2.3.1.- setting the headers
	    if os.path.exists(file) == False:
		f = open(file, 'w')
		f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
		f.close()
	    #/ 1.2.2.- appending the information
	    f = open(file, 'a')
	    f.write('%s\n'%value)
	    f.close()
    #/ 3.- saving the SBML file
    sbmlWorker.sbmlWriter(libsbml.SBMLWriter(), model, metamodel, outputfiles.solutionDirectory + '/%s.xml' %model.systemName)
    #/ 4.- ending up the options of the hybrid algorithm
    METAMODEL.hybridScore = False 

    return None
