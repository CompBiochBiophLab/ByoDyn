#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana
#
#  Created: 2005-11-09 by Adrian L. Garcia-Lomana
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

# $Id: surface.py,v 4.5 2008/12/03 16:37:38 alglomana Exp $

## \file
# This module builds the fitness function surface for combination of 2 parameters.

import math, copy, os, sys, fpformat
import errorMessages, centralFunctions, optimiser, initiator

def central(metamodel, model, outputfiles):

    '''
    This function directs the creation of the surface of fitness function values.
    It checks that the selected parameters are in model.
    It variates the parameters and calls the scoreObtainer function from the 
    Optimiser module to calculate the fitness function value.
    '''
    
    #/ 0.-  modifying model parameters in case it was defined on the model
    for parameter in metamodel.parameters:
        for modelParameter in model.parameters:
            if parameter == modelParameter:
                model.parameters[parameter] = metamodel.parameters[parameter] 
    #/ 1.- getting a backup of a variable
    modelValues = copy.deepcopy(model.parameters)
    #/ 2.- checking introduced data
    checkSurfaceParameters(metamodel, model)
    #/ 3.- score obtaining
    for plot in metamodel.surfaceParameters:
	print 'Calculating %s surface ...'%plot
	surface = []
	surfaceComponent = []
        #/ 3.1.- obtaining the grid of values
        xScale = 'lin'; yScale = 'lin'
        xParameter = plot.split('/')[0]; yParameter = plot.split('/')[1]
        for parameterToVary in metamodel.parametersToVary:
            if xParameter == parameterToVary.split('/')[0]:
                xMax = float(parameterToVary.split('/')[2])
                xMin = float(parameterToVary.split('/')[1])
                #/ checking the lin/log scale
                if parameterToVary.split('/')[3] == 'log':
                    xMax = math.log(xMax, 10)
                    xMin = math.log(xMin, 10)
                    xScale = 'log'
            if yParameter == parameterToVary.split('/')[0]:
                yMax = float(parameterToVary.split('/')[2])
                yMin = float(parameterToVary.split('/')[1])
                #/ checking the lin/log scale
                if parameterToVary.split('/')[3] == 'log':
                    yMax = math.log(yMax, 10)
                    yMin = math.log(yMin, 10)
                    yScale = 'log'
        #/ 3.2.- checking that parameterMax is bigger than parameterMin
        if xMax < xMin or yMax < yMin:
            raise errorMessages.ClassSurfaceException, 'the range of the parameters is not at the correct order. First must be the minimum value to explore and second the maximum value to explore.'
        #/ 3.3.- obtaining the grid components
        xGrid = []; yGrid = []
        xInterval = xMax - xMin; yInterval = yMax - yMin
        xGridStep = xInterval / float(metamodel.surfaceResolution); yGridStep = yInterval / float(metamodel.surfaceResolution)
        for i in range(metamodel.surfaceResolution):
            xValue = xMin + (i * xGridStep) + xGridStep/2.; yValue = yMin + (i * yGridStep) + yGridStep/2.
            xGrid.append(xValue); yGrid.append(yValue)
        #/ 3.4.- converting again to proper values for the simulation
        if xScale == 'log':
            newXGrid = []
            for i in range(len(xGrid)):
                xValue = 10**xGrid[i]
                newXGrid.append(xValue)
            xGrid = newXGrid
        if yScale == 'log':
            newYGrid = []
            for i in range(len(yGrid)):
                yValue = 10**yGrid[i]
                newYGrid.append(yValue)
            yGrid = newYGrid
        for xPoint in xGrid:
	    surfaceComponent = []
            for yPoint in yGrid:
                for modelParameter in model.parameters:
                    if modelParameter == xParameter:
                        model.parameters[xParameter] = xPoint
                    if modelParameter == yParameter:
                        model.parameters[yParameter] = yPoint
                #/ obtaining the score
                score = optimiser.scoreObtainer(metamodel, model, outputfiles)
		surfaceComponent.append(math.log10(score))
	    surface.append(surfaceComponent)
        model.parameters = copy.deepcopy(modelValues)
	#/ 3.5.- creating the text file and the plot
	surfaceTextSaver(outputfiles, metamodel, plot, surface)
	surfacePlotter(metamodel, outputfiles, plot, surface, xGrid, yGrid)

    return None

def checkSurfaceParameters(metamodel, model):

    '''
    This function checks that the selected parameters for the surface plot do exist in the model.
    '''

    existing = False
    introducedParameters = []
    modelParameters = []
    for introducedParameter in metamodel.surfaceParameters:
        introducedParameters.append(introducedParameter.split('/')[0])
        introducedParameters.append(introducedParameter.split('/')[1])
    for modelParameter in model.parameters:
        modelParameters.append(modelParameter)
    for i in range(len(introducedParameters)):
        for j in range(len(modelParameters)):
            if introducedParameters[i] == modelParameters[j]:
                existing = True
        if existing == False:
            raise errorMessages.ClassSurfaceException, 'there is a parameter for the score surface that is not in the model.'
        existing == False

    return None

def rainbowColorCalculation(value, minValue, maxValue):

    '''
    This function calculates the color based on a rainbow scale
    '''
    
    #/ 1.- normalising from floats to 0. to 4.
    absoluteValue = value - minValue
    range = maxValue - minValue
    normalisedValue = (absoluteValue/range)*4. 
    value = normalisedValue
    #/ 2.- calculating the RGB colour
    if 0 <= value <= 1:
	R = 0.
	G = value
	B = 1.
    elif 1 < value <= 2:
	R = 0.
	G = 1.
	B = 2. - value
    elif 2 < value <= 3:
	R = value - 2.
	G = 1.
	B = 0.
    elif 3 < value <= 4:
	R = 1.
	G = 4. - value
	B = 0.
    else:
	print 'Value of %s can not be computed\nError during colour calculation.\nExiting ...'%value
	sys.exit()

    return R, G, B

def surfacePlotter(metamodel, outputfiles, plot, surface, xGrid, yGrid):
    
    '''
    This function creates the postscript surface plot.
    '''
    
    #/ 1.- opening the file
    fileName = outputfiles.outputdir + '/' + plot.split('/')[0] + '.' + plot.split('/')[1] + '.ps'
    f = open(fileName, 'w')
    #/ 2.- byodyn header
    f.write('%%\n%% generated by ByoDyn version %s\n%%\n'%initiator.BYODYNVERSION)
    #/ 2.- introducing the postscript functions
    libdir = os.environ.get('BYODYN_PATH') + '/lib'
    headerFile = open('%s/postscript/header.ps' %libdir, 'r')
    header =  headerFile.readlines()
    headerFile.close()
    for line in header:
	f.write('%s' %line)
    boxFile = open('%s/postscript/box.ps' %libdir, 'r')
    box =  boxFile.readlines()
    boxFile.close()
    for line in box:
	f.write('%s' %line)
    normalfontFile = open('%s/postscript/normalFont.ps' %libdir, 'r')
    normalfont = normalfontFile.readlines()
    normalfontFile.close()
    for line in normalfont:
	f.write('%s' %line)
    #/ 3.- calculating the positioning
    xMargin = 60.
    yMargin = 200.
    boxDimention = 500.
    resolution = boxDimention / metamodel.surfaceResolution
    maxValue = 0.
    for vector in surface:
	for element in vector:
	    if element > maxValue:
		maxValue = element
    minValue = maxValue
    for vector in surface:
	for element in vector:
	    if element < minValue:
		minValue = element
    #/ 4.- positioning the scale referece, labels and title
    #/ 4.1.- title
    f.write('newpath 200 750 moveto (Fitness function surface) show\n')
    #/ 4.2.- labels
    parameterX = plot.split('/')[0]
    parameterY = plot.split('/')[1]
    for parameter in metamodel.parametersToVary:
	vector = parameter.split('/')
	if vector[0] == parameterX:
	    if vector[3] == 'lin':
		xRange = [float(vector[1]), float(vector[2])]
		labelX = parameterX
	    elif vector[3] == 'log':
		xRange = [math.log10(float(vector[1])), math.log10(float(vector[2]))]
		labelX = 'log ' + parameterX 
	if vector[0] == parameterY:
	    if vector[3] == 'lin':
		yRange = [float(vector[1]), float(vector[2])]
		labelY = parameterY
	    elif vector[3] == 'log':
		yRange = [math.log10(float(vector[1])), math.log10(float(vector[2]))]
		labelY = 'log ' + parameterY
    f.write('newpath 250 160 moveto (%s) show\n'%labelY)
    f.write('newpath 50 710 moveto (%s) show\n'%labelX)
    #/ 4.3.- ranges
    f.write('newpath 30 180 moveto (%s) show\n'%fpformat.fix(yRange[0], 2))
    f.write('newpath 510 180 moveto (%s) show\n'%fpformat.fix(yRange[1], 2))
    f.write('newpath 10 200 moveto (%s) show\n'%fpformat.fix(xRange[1], 2))
    f.write('newpath 10 680 moveto (%s) show\n'%fpformat.fix(xRange[0], 2))
    #/ 4.4.- scale
    xIniScale = 60.
    xFinScale = 65. 
    yIniScale = 50.
    yFinScale = 100.
    for i in range(100):
	R, G, B = rainbowColorCalculation(i+.5, 0., 100.)
	f.write('%s %s %s %s %s %s %s boxRGB\n'%(xIniScale, yIniScale, xFinScale, yFinScale, R, G, B))
	xIniScale = xIniScale + 5.
	xFinScale = xFinScale + 5.
    f.write('newpath 40 30 moveto (%s) show\n'%fpformat.fix(minValue, 2))
    f.write('newpath 545 30 moveto (%s) show\n'%fpformat.fix(maxValue, 2))
    #/ 5.- writing the boxes
    xIni = xMargin
    xFin = xMargin + resolution
    yIni = yMargin + boxDimention - resolution
    yFin = yMargin + boxDimention
    for i in range(metamodel.surfaceResolution):
	for j in range(metamodel.surfaceResolution):
	    R, G, B = rainbowColorCalculation(surface[i][j], minValue, maxValue)
	    f.write('%s %s %s %s %s %s %s boxRGB\n'%(xIni, yIni, xFin, yFin, R, G, B))
	    xIni = xIni + resolution
	    xFin = xFin + resolution
	xIni = xMargin
	xFin = xMargin + resolution
	yIni = yIni - resolution
	yFin = yFin - resolution
    f.close()
    #/ 6.- converting to png if it is the case
    if metamodel.figureFormat == 'png':
	pngFile = fileName.replace('ps', 'png')
	cmd = 'convert %s %s'%(fileName, pngFile)
	os.system(cmd)

    return None

def surfaceTextSaver(outputfiles, metamodel, plot, surface):

    '''
    This function writes in a text file the results of the fitness function evaluation along the grid.
    '''

    #/ 1.- opening the file
    fileName = outputfiles.outputdir + '/' + plot.split('/')[0] + '.' + plot.split('/')[1] + '.txt'
    f = open(fileName, 'w')
    #/ 2.- byodyn header and further information
    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
    f.write('# fitness function values in the grid\n#\n')
    #/ 3.- writing the values
    for i in range(metamodel.surfaceResolution):
	for j in range(metamodel.surfaceResolution):
	    f.write('%s\t'%10**surface[i][j]) #/ we received the log values, we write the original fitness function values
	f.write('\n')
    f.close()
    

    return None
