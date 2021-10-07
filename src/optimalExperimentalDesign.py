#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido and Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Alex Gomez-Garrido, David Sportouch and Adrian L. Garcia-Lomana
#
#  Created: 2008-01-28 by Alex Gomez-Garrido
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

# $Id: optimalExperimentalDesign.py,v 4.9 2008/12/14 19:27:52 alglomana Exp $

## \file
# This module contains the algorithms necessary for the optimal experimental design protocols of a given model.

import sensitivityAnalyzer, identifiabilityAnalyzer, errorMessages, initiator
import copy, sys, scipy
try:
    import matplotlib.pyplot
except ImportError:
    raise errorMessages.ClassSimulatorStochasticException, 'error while importing matplotlib.pyplot.'


def central(model, metamodel, outputfiles):

    '''
    This is the central function of the module.
    It directs the flow of the program for the optimal experimental design.
    '''
    
    #/ 1.- Setting the new parameter values of the of the model, if it is that was asked in the runner
    for parameter in metamodel.parameters:
        for modelParameter in model.parameters:
            if parameter == modelParameter:
        	model.parameters[parameter] = metamodel.parameters[parameter]	    
    #/ 2.- Getting parameters to study from metamodel
    parametersToStudy = sensitivityAnalyzer.componentsChecker(model, metamodel)
    sensitivity = sensitivityAnalyzer.ClassSensitivityAnalyzer(parametersToStudy)
    sensitivity.calculateSens(model, metamodel, outputfiles)   
    #/ 3.- Identifiability analysis
    identifiability = identifiabilityAnalyzer.ClassIdentifiabilityAnalyzer()
    identifiability.setSensitivity(sensitivity.identifiabilityCoefficients)
    identifiability.calculateFIM(parametersToStudy, metamodel.target)
    identifiability.setCOV(parametersToStudy)
    #/ 4.- choosing task
    if metamodel.oedTask == 'addNewPoint':
        print 'Protocol in use: "%s". Determining the most informative position ...'%metamodel.oedTask
        addNewPoint(identifiability, sensitivity, metamodel, model, outputfiles, parametersToStudy)
    elif metamodel.oedTask == 'rankTargets':
        print 'Protocol in use: "%s". Determining the most informative target ...'%metamodel.oedTask
        rankTargets(identifiability, sensitivity, metamodel, model, outputfiles, parametersToStudy)
    return None
  
def addNewPoint(identifiability, sensitivity, metamodel, model, outputfiles, parametersToStudy):

    '''
    This function runs the addNewPoint optimal experimental design protocol.
    In this task we evaluate the criteria value at different timePoints after adding a putative target of the given species at that specific time point. 
    We can also run the protocol for a list of species.  
    '''

    file = open(outputfiles.oedResults, 'w')
    file.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
    file.write('# This file contains, for each optimal experimental design criteria, for each species, the time point at which the criteria is optimum.\n#\n')
    
    for criteria in metamodel.identifiabilityCriteria:
	totalModelValues = {}
	dataFile = outputfiles.scratchdir + '/' + model.systemName + '.' + criteria + '.txt'
        initialCriteria = getCriteria(criteria, identifiability)
        results = {}
	startPlottingTime = 0.0
	endPlottingTime = metamodel.simulationTime
	timePoints = timePointsDeterminer(metamodel)
        for specie in metamodel.oedTargetSpecies:
	    totalModelValues[specie] = []
	    OEDvalues = []
            for tp in timePoints:
                original = copy.deepcopy(metamodel.target)
                identifiability = identifiabilityAnalyzer.ClassIdentifiabilityAnalyzer()
                for yindex in range(model.ywidth):
                    for xindex in range(model.xlength):
			#/ this is very interesting: the concentration value of the target does not matter the resut, it is just a sensitivity value at a given point
                        if metamodel.target.keys().count(tp) == 0:
                            metamodel.target[tp] = ['%s/%s,%s/1.0/1.0' %(specie, yindex, xindex)]
                        else:
                            metamodel.target[tp].append('%s/%s,%s/1.0/1.0' %(specie, yindex, xindex))
                sensitivity.calculateSens(model, metamodel, outputfiles)
                identifiability.setSensitivity(sensitivity.identifiabilityCoefficients)
		identifiability.calculateFIM(parametersToStudy, metamodel.target)
                identifiability.setCOV(parametersToStudy)
                newCriteria = getCriteria(criteria, identifiability)
		OEDvalues.append(newCriteria)
                metamodel.target = original
	    if criteria == 'MA' or criteria == 'D' or criteria == 'E':
		results[specie] = [timePoints[OEDvalues.index(chooseValue(criteria, OEDvalues))], chooseValue(criteria, OEDvalues)]
	    elif criteria == 'ME':
		results[specie] = [timePoints[OEDvalues.index(chooseValue(criteria, OEDvalues))], chooseValue(criteria, OEDvalues)]
	    else:
		raise errorMessages.ClassOptimalExperimentalDesignException, 'non recognised optimal criteria.'
	    totalModelValues[specie] = OEDvalues
	    matplotlib.pyplot.plot(timePoints, OEDvalues, 'o')
	    matplotlib.pyplot.hold(True)
	matplotlib.pyplot.legend(metamodel.oedTargetSpecies, loc=0)
	matplotlib.pyplot.xlabel('Time') #/ set x-axis label
	matplotlib.pyplot.ylabel('%s criteria' %criteria) #/ set y-axis label
	matplotlib.pyplot.title(model.systemName) #/ set plot title
	#/ setting x range for plotting
	for time in timePoints:
	    if time < .05*metamodel.simulationTime:
		startPlottingTime = -.05*metamodel.simulationTime
	    if time > .95*metamodel.simulationTime:
		endPlottingTime = metamodel.simulationTime + .05*metamodel.simulationTime
	matplotlib.pyplot.xlim(startPlottingTime, endPlottingTime)
	#/ choosing figure format
	if metamodel.figureFormat == 'ps':
	    timePlotFileName = outputfiles.oedTimePlot + '%s'%criteria + '.ps'
	elif metamodel.figureFormat == 'png':
	    timePlotFileName = outputfiles.oedTimePlot + '%s'%criteria + '.png'
	matplotlib.pyplot.savefig(timePlotFileName,dpi=1200)
	matplotlib.pyplot.hold(False)
	#/ writing data file
	f = open(dataFile, 'w')
	f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
	f.write('# Time')
	for specie in metamodel.oedTargetSpecies:
	    f.write('\t%s'%specie)
	f.write('\n')
	for i in range(len(timePoints)):
	    f.write('%s\t'%timePoints[i])
	    for species in metamodel.oedTargetSpecies:
		f.write('%s\t'%totalModelValues[species][i])
	    f.write('\n')
	f.close()
	#/ writing the total value for each criteria
        listResults = []
        file.write('# Species\tTime\tCriteria (%s)\n' %criteria)
        for r in results.keys():
            file.write('%s\t%s\t%s\n' %(r, results[r][0], results[r][1]))
            listResults.append(results[r][1])
	optimumSpecies=results.keys()[listResults.index(chooseValue(criteria, listResults))]
	optimumTime=results[results.keys()[listResults.index(chooseValue(criteria, listResults))]][0]
	optimumValue=chooseValue(criteria, listResults)
	print 'Species: %s; Time: %s; Criteria value(%s): %s'%(optimumSpecies,optimumTime,criteria,optimumValue)
    file.close()
    
    return None

def rankTargets(identifiability, sensitivity, metamodel, model, outputfiles, parametersToStudy):

    '''
    This function is still in development status. It sorts of targets indicating the more and the less important. 
    Less important means with low effect to improve the identifiability criteria.
    '''

    file = open(outputfiles.oedResults, 'w')
    for criteria in metamodel.identifiabilityCriteria:
        initialCriteria = getCriteria(criteria, identifiability)
        results = {}
        originalTargets = copy.deepcopy(metamodel.target)
        for specie in originalTargets[metamodel.target.keys()[0]]:
            for target in metamodel.target:
                for field in metamodel.target[target]:
                    if field.split('/')[0] == specie.split('/')[0]:
                        metamodel.target[target].remove(field)
                        break
            print metamodel.target
            sensitivity.calculateSens(model, metamodel, outputfiles)
            identifiability.setSensitivity(sensitivity.identifiabilityCoefficients)
            identifiability.calculateFIM(parametersToStudy, metamodel.target)
            identifiability.setCOV(parametersToStudy)
            newCriteria = getCriteria(criteria, identifiability)
            print specie.split('/')[0], newCriteria
            metamodel.target = copy.deepcopy(originalTargets)
#            metamodel.target.pop(target)
#            sensitivity.calculateSens(model, metamodel, outputfiles)
#            identifiability.setSensitivity(sensitivity.identifiabilityCoefficients)
#            identifiability.calculateFIM(parametersToStudy, metamodel.target)
#            identifiability.setCOV(parametersToStudy)
#            newCriteria = getCriteria(criteria, identifiability)
#
#            print newCriteria
#            metamodel.target = originalTargets


#            results[specie] = [timePoints[OEDvalues.index(min(OEDvalues))], chooseValue(criteria, OEDvalues)]
#        listResults = []
#        file.write('Species\tTime\tCriteria %s\n' %criteria)
#        for r in results.keys():
#            file.write('%s\t%s\t%s\n' %(r, results[r][0], results[r][1]))
#            listResults.append(results[r][1])
#        print 'Species: %s; Time: %s; Criteria value: %s' %(results.keys()[listResults.index(chooseValue(criteria, listResults))], results[results.keys()[listResults.index(chooseValue(criteria, listResults))]][0], chooseValue(criteria, listResults))
    file.close()
    
    return None

def getCriteria(criteria, identifiability):
    
    '''
    This function returns the criteria value specified by the user in runner. 
    '''

    if criteria == 'ME':
        return identifiability.getMEcriteria()
    elif criteria == 'MA':
        return identifiability.getMAcriteria()
    elif criteria == 'D':
        return identifiability.getDcriteria()
    elif criteria == 'E':
        return identifiability.getEcriteria()

def chooseValue(criteria, values):

    '''
    This function returns the maximmum or minimmun value of values lists depending on the criteria type. 
    '''
    
    if criteria == 'MA' or criteria == 'D' or criteria == 'E':
        return max(values)
    elif criteria == 'ME':
	return min(values)
    else:
        raise errorMessages.ClassOptimalExperimentalDesignException, 'non recognised optimal criteria.'

def timePointsDeterminer(metamodel):

    '''
    This function returns the time points of study of the optimal experimental design.
    It is based on the metamodel arguments of time&timestep and OEDResolution.
    '''

    
    if metamodel.oedResolution == 0:
	raise errorMessages.ClassOptimalExperimentalDesignException, 'OEDResolution should be a non zero integer.'
    else:
	maxResolution = int(metamodel.simulationTime / metamodel.simulationTimeStep)
	if metamodel.oedResolution > maxResolution:
	    raise errorMessages.ClassOptimalExperimentalDesignException, 'maximal OEDResolution accepted is %s. %s provided, please change it accordingly.'%(maxResolution, metamodel.oedResolution)
	else:
	    #/ calculating the timePoints based on time step
	    resolution = metamodel.simulationTimeStep * metamodel.oedResolution
	    timePoints = scipy.arange(0, metamodel.simulationTime + metamodel.simulationTimeStep, resolution)

    return timePoints
