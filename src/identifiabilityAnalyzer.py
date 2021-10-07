#
#  Projec: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido and Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Alex Gomez-Garrido and David Sportouch
#
#  Created: 2007-07-18 by Alex Gomez-Garrido
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

# $Id: identifiabilityAnalyzer.py,v 4.5 2008/12/03 16:37:38 alglomana Exp $

## \file
# This module contains the algorithms necessary for the analysis of identifiability of a given model.

import math, sys
import scipy.linalg
import matrixWorker, simulator, sensitivityAnalyzer, initiator, errorMessages

class ClassIdentifiabilityAnalyzer:
    
    '''
    Class for the identifiability analysis.
    '''
    
    def __init__(self):
        
        '''
        This is the constructor.
        '''
        
        self.sensitivityCoefficients = {}
        self.FIM = matrixWorker.ClassMatrix()
        self.COV = matrixWorker.ClassMatrix()
        self.CorrelationMatrix = matrixWorker.ClassMatrix()
       
	return None
    
    def setSensitivity(self, coefficients):
        
        '''
        This method set the sensitivity coefficients. 
        '''
        
        self.sensitivityCoefficients = coefficients

        return None
    
    def calculateFIM(self, parametersToStudy, metamodelTarget):
        
        '''
        This method calculates of Fisher Information Matrix (FIM) from sensitivity coefficients.
        '''
        
	#/ 1.- starting parameters
        k = 0
        sortParameters = parametersToStudy.keys()
        sortParameters.sort()
        #/ 2.- setting names
        self.FIM.colNames = sortParameters
        self.FIM.rowNames = sortParameters
	#/ 3.- detecting the weights. It should be the same way as in __obtainSensitivity method of sensitivityAnalyser module
	weights = []
	targetedNodes = []
	for target in metamodelTarget:
	    for point in metamodelTarget[target]:
		pointFields = point.split('/')
		plausibleNode = pointFields[0]
		if targetedNodes.count(plausibleNode) == 0:
		    targetedNodes.append(plausibleNode)  
	for targetedNode in targetedNodes:
	    for target in metamodelTarget:
		for point in metamodelTarget[target]:
		    pointFields = point.split('/')
		    plausibleNode = pointFields[0]
		    plausibleWeight = float(pointFields[3])
		    if targetedNode == plausibleNode:
			weights.append(plausibleWeight)
        #/ 3.- filling up the matrix 
        for i in range(len(sortParameters)):
            self.FIM.m.append([])
	    for j in range(len(sortParameters)):
		value = 0.
		for k in range(len(weights)):
		    value = value + (1/weights[k])*self.sensitivityCoefficients[sortParameters[i]][k]*self.sensitivityCoefficients[sortParameters[j]][k]
		self.FIM.m[i].append(value)

        return None

    def setCOV(self, parametersToStudy):
        
        '''
        This mehtod set the covariance Matrix and calculates its EigenValues, it is the inverse of FIM.
        '''

        #/ 1.- setting names
	sortParameters = parametersToStudy.keys()
	sortParameters.sort()
	self.COV.colNames = sortParameters
	self.COV.rowNames = sortParameters
	#/ 2.- converting the matrix
        if self.FIM.m != []:            
            if self.FIM.isSingular() == False:
                self.COV.m = scipy.linalg.inv(self.FIM.m)
                self.COV.setEigenValues()
            else:        
                raise errorMessages.ClassIdentifiabilityAnalyzerException,  'Error when we are calculating the FIM. It is singular.'
        else:
            print 'You need to calculate before the FIM matrix'

        return None

    def setCorrelationMatrix(self, parametersToStudy):
        
        '''
        This method calculate the Correlation Marix from the COV.
        '''
        
        i = 0
        if self.COV.m == []:
            self.setCOV()
        for row in self.COV.m:
            j = 0
            self.CorrelationMatrix.m.append([])
            for element in row:
                element = self.COV.m[i][j] / math.sqrt(abs(self.COV.m[i][i]*self.COV.m[j][j])) 
                if element < -1.0  or element > 1.0:
                    raise errorMessages.ClassIdentifiabilityAnalyzerException,  'Error when we are calculating the correlation matrix. It has values out of ranges (1, -1). It could happen if the covariance matrix has values with a big difference of magnitud order.'
                self.CorrelationMatrix.m[i].append(element)
                j+=1
            i += 1
        sortParameters = parametersToStudy.keys()
        sortParameters.sort()
        self.CorrelationMatrix.colNames = sortParameters
        self.CorrelationMatrix.rowNames = sortParameters

        return None
    
    def getMAcriteria(self):
        
        '''
        This method returns the modified A-optimal design criteria, which it is to maximise the trace of FIM.
        '''
        
        return self.FIM.getTrace()

    def getDcriteria(self):
        
        '''
        This method returns the D-optimal design criteria, which it is to maximise the determinant of FIM.
        '''

        return scipy.linalg.det(self.FIM.m)
    
    def getEcriteria(self):
        
        '''
        This method returns the E-optimal design criteria, which it is to maximise the lowest FIM eigenvalue.
        '''
        
        return self.FIM.getMinEigenValue()
    
    def getMEcriteria(self):

        '''
        This method returns the modified E-optimal design criteria, which it is to minimise the ratio among the largest to lowest eigenvalue.
        '''
        
        return self.FIM.getMaxEigenValue() / self.FIM.getMinEigenValue()

    def confidenceIntervals(self, i, threshold):
        
        '''
        This method calculates the confidence interval for one parameter.
        '''
        
        if threshold == 0.05:
            u = 1.96
        elif threshold == 0.5:
            u = 0.6745
        else:
            print 'value not handled yet, threshold = 0.05 by default'
            u = 1.96
        I = u*math.sqrt(self.COV.m[i][i])

        return I
    
    def residualsComputation(expValues, model, metamodel, outputfiles):

        '''
        This method calculates the residual for one experimental value.
        '''

	#/ This function has to compute the residuals values at each time point
	#/ for each specie and then to draw a plot like the one Fig.5 (Julio Banga)
	#/ see function sensitivityObtainer

        originalSimulationValues = simulator.obtainSimulationValues(metamodel, model, outputfiles, 'simulation')
        timePoints = []
        residuals = {}
        for t in expValues.keys():    
            timePoint = int(float(t)/float(metamodel.simulationTimeStep))
            for node in expValues[t].keys():
                if expValues.keys().index(t) == 0:
                    residuals[node] = []            
                residuals[node].append(float(expValues[t][node]) - float(originalSimulationValues[timePoint][model.nodes.index(node)]))
        
        return None
    
    def criteriaWriter(self, outputfiles, MA, D, E, ME):
    
        '''
        This method write the criteria information in its output file.
        '''
        
        f = open(outputfiles.identifiabilityCriteria, 'w')
        f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
	f.write('# Values of the different identifiability criteria.\n#\n')
        if MA != None:
            f.write('Modified A-optimal design criteria: max trace FIM.\n\t trace FIM = %s\n' %MA)
        if D != None:
            f.write('D-optimal design criteria: max det FIM.\n\t det FIM = %s\n' %D)
        if E != None:
            f.write('E-optimal design criteria: max lowest FIM eigenvalue.\n\t lowest FIM eigenvalue = %s\n' %E)
        if ME != None:
            f.write('Modified E-optimal design criteria: min ratio(largest FIM eigenvalue / lowest FIM eigenvalue).\n\t ratio(largest FIM eigenvalue / lowest FIM eigenvalue) = %s\n' %ME)
        f.close()

        return None

def central(model, metamodel, outputfiles):

    '''
    This is the central function of the module.
    It directs the flow of the program for the identifiability analysis.
    '''

    #/ 1.- Setting the parameters of the model, if it is necessary
    for parameter in metamodel.parameters:
        for modelParameter in model.parameters:
            if parameter == modelParameter:
                model.parameters[parameter] = metamodel.parameters[parameter] 
    #/ 1.- Getting parameters to study from metamodel
    parametersToStudy = sensitivityAnalyzer.componentsChecker(model, metamodel)
    sensitivity = sensitivityAnalyzer.ClassSensitivityAnalyzer(parametersToStudy)
    sensitivity.calculateSens(model, metamodel, outputfiles)   
    #/ 2.- Identifiability analysis
    identifiability = ClassIdentifiabilityAnalyzer()
    identifiability.setSensitivity(sensitivity.identifiabilityCoefficients)
    identifiability.calculateFIM(parametersToStudy, metamodel.target)
    identifiability.setCOV(parametersToStudy)
    identifiability.setCorrelationMatrix(parametersToStudy)
    #/ 3.- Identifiability criteria
    if metamodel.identifiabilityCriteria != []:
        MA = None
        D = None
        E = None
        ME = None
        for criterion in metamodel.identifiabilityCriteria:
            if criterion == 'MA':
                MA = identifiability.getMAcriteria()
            elif criterion == 'D':
                D = identifiability.getDcriteria()
            elif criterion == 'E':
                E = identifiability.getEcriteria()
            elif criterion == 'ME':
                ME = identifiability.getMEcriteria()
        identifiability.criteriaWriter(outputfiles, MA, D, E, ME)
    #/ 4.- outputting the results of identifiability
    if metamodel.identifiabilityOutputs != []:
        for output in metamodel.identifiabilityOutputs:
            if output == 'FIM':
		message = 'Fisher information matrix values.'
                identifiability.FIM.fileWriter(message, outputfiles.identifiabilityFIM)
            elif output == 'COV':
		message = 'Values of the covariance matrix from the Fisher information matrix.'
                identifiability.COV.fileWriter(message, outputfiles.identifiabilityCOV)
            elif output == 'Correlation':
		message = 'Correlation values of the covariance matrix from the Fisher information matrix.'
                identifiability.CorrelationMatrix.fileWriter(message, outputfiles.identifiabilityCorrelation)
                identifiability.CorrelationMatrix.plotWriter(metamodel, outputfiles.identifiabilityCorrelationPlot, 'blueRed')
    
    return None
