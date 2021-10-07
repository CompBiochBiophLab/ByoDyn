#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Alex Gomez-Garrido
#
#  Created: 2008-07-09 by Alex Gomez-Garrido
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

# $Id: sampler.py,v 4.14 2008/12/14 19:27:52 alglomana Exp $

## \file 
# This module is responsible of sample methods of fitness function surface.

import math, random, scipy, numpy, sys, os
import optimiser, simulator, initiator, errorMessages
try:
    import matplotlib.pyplot
except ImportError:
    raise errorMessages.ClassSimulatorStochasticException, 'error while importing matplotlib.pyplot.'

class ClassSamplerMonteCarlo:

    '''
    Class for the Monte Carlo sampler. The method implemented is based in Metropolis-Hastings algorithm.
    '''    

    def __init__(self, metamodel):

        '''
        The constructor.
        '''

        self.acceptedSample = []
        self.nonAcceptedSample = []
        self.sampleSize = metamodel.sampleSize
        self.T = 100000
        
        return None

    def sampling(self, metamodel, model, outputfiles):

        '''
        This method sampling the fitness function surface using a Monte Carlo Metropolis-Hastings algorithm.
        It creates a list of accepted sampled elements and not accepted sampled elements.
        '''

	print 'Sampling the surface ...'
        sample = ClassSample() 
        sample.createRandomSetOfParameters(metamodel)
        sample.calculateFitnessFunction(metamodel, model, outputfiles)
        self.acceptedSample.append(sample)
	#/ 1.- Running the algorithm iterations
	sampledPoints = 0
        while len(self.acceptedSample) < self.sampleSize:
	    sampledPoints = sampledPoints + 1
            sample = ClassSample() 
            sample.createRandomSetOfParameters(metamodel)
            sample.calculateFitnessFunction(metamodel, model, outputfiles)            
            deltaFF = sample.fitnessFunction - self.__lastFF()
            if deltaFF/self.T > 100: #/ In order to solve numerical problems
                probability = 0
            elif deltaFF/self.T < 0:
                probability = 1 
            else:
                probability = math.exp(-deltaFF/self.T)      
            if probability >= random.random():
                self.acceptedSample.append(sample)
            else:
                self.nonAcceptedSample.append(sample)
	ratio = float(len(self.acceptedSample))/float(sampledPoints)
	percentage = ratio*100
	print 'Rate of accepted points: '+ '%.2f'%percentage + ' %'

        return self

    def setT(self, metamodel, model, outputfiles):

        '''
        This method sets the class variable T from a random exploration of Fitness Function surface.
        '''

	print 'Heuristic determination of sampling temperature ...'
        #numOfSample = 100
	numOfSample = 5
        FFvalues = []       
        for i in range(numOfSample):
            sample = ClassSample() 
            sample.createRandomSetOfParameters(metamodel)
            sample.calculateFitnessFunction(metamodel, model, outputfiles)            
            FFvalues.append(sample.fitnessFunction)
        sd = scipy.std(FFvalues)
        self.T = sd/10

        return self
        
    def __lastFF(self):

        '''
        This private method returns the fitness function of the last sampled element.
        '''

        return self.acceptedSample[len(self.acceptedSample)-1].fitnessFunction

class ClassSample:

    '''
    Class for the sample object.
    '''    

    def __init__(self):

        '''
        The constructor.
        '''

        self.parameters = {}
        self.fitnessFunction = None
        
        return None

    def createRandomSetOfParameters(self, metamodel):

        '''
        This method creates a random set of parameters.
        The parameters to explore and its range is indicated by the parametersToVary variable in the runner file.
        '''
        
        #/ 1.- variate parameters
        for parameterToVary in metamodel.parametersToVary:
            if parameterToVary.split('/')[3] == 'lin':
                self.parameters[parameterToVary.split('/')[0]] = optimiser.linearVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))
            elif parameterToVary.split('/')[3] == 'log': #/ logarithmic scale
                self.parameters[parameterToVary.split('/')[0]] = optimiser.logarithmicVariation(float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2]))
            else:
                raise errorMessages.ClassSamplerException, 'there is a problem with the way of exploring the parameters, either linear of logarithmic. Other value given.'            

        return self        
        
    def calculateFitnessFunction(self, metamodel, model, outputfiles):

        '''
        This method calculate the fitness function for a set of parameters.
        '''

        for parameter in self.parameters:
            model.parameters[parameter] = self.parameters[parameter]        
        self.fitnessFunction = optimiser.scoreObtainer(metamodel, model, outputfiles)        

        return self
    
def central(model, metamodel, outputfiles):

    '''
    This is the central function of the module.
    It directs the flow of the Monte Carlo sampling depending on the options.
    '''

    simulator.compatibilityChecker(metamodel, model)
    optimiser.checkDataPoints(model, metamodel)
    optimiser.checkTargetNodes(model, metamodel)
    optimiser.checkParametersToVary(model, metamodel)
    if metamodel.sampleMethod == 'MonteCarlo':
        samplerObject = ClassSamplerMonteCarlo(metamodel)
    samplerObject.setT(metamodel, model, outputfiles)
    samplerObject.sampling(metamodel, model, outputfiles)
    print 'Saving the sampled points ...'
    storeResults(samplerObject.acceptedSample, outputfiles.sampleResults, metamodel)
    storeResults(samplerObject.nonAcceptedSample, outputfiles.sampleNonAcceptedResults, metamodel)
    plotResults(metamodel, outputfiles, samplerObject.acceptedSample)

    return None

def rainbowColourCalculator(value): 

    '''
    Given a normalised value on the 0-4 range, this function returns the corresponding RGB colour code based on the rainbow scale. 
    '''

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
	raise errorMessages.ClassSamplerException, 'Value of %s can not be computed\nError during colour calculation.'%value

    return R, G, B

def storeResults(listOfSample, file, metamodel):

    '''
    This method stores the sampling solution on its specific file.
    '''

    #/ 1.- determining the scales
    scales = {}
    for parameterToVary in metamodel.parametersToVary:
	vector = parameterToVary.split('/')
	scales[vector[0]] = vector[3]
    #/ 2.- writing the file
    if len(listOfSample) != 0:    
        parameterNames = listOfSample[0].parameters.keys()
	f = open(file, 'w')
        f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
	f.write('#\t')
        for parameter in parameterNames:	    
            f.write('%s(%s)\t' %(parameter,scales[parameter]))
        f.write('fitnessFunction\n')
        for sample in listOfSample:
            for parameter in parameterNames:
                f.write('%s\t' %sample.parameters[parameter])
            f.write('%s\n' %sample.fitnessFunction)
    
    return None

def plot1D(metamodel, outputfiles, sample):
    
    '''
    This function creates a histogram of the distribution of the accepted sampled points.
    '''
    
    #/ 1.- initialising the values
    x = []
    dimensions = sample[0].parameters.keys()
    xlabel = None
    #/ 2.- determining the scale and appending the values
    if metamodel.parametersToVary[0].split('/')[3] == 'log':
	xLim = [math.log10(float(metamodel.parametersToVary[0].split('/')[1])),math.log10(float(metamodel.parametersToVary[0].split('/')[2]))]
	xlabel = 'log ' + dimensions[0]
	for point in sample:
	    value = point.parameters[dimensions[0]]
	    x.append(math.log10(value))
    else:
	xLim = [float(metamodel.parametersToVary[0].split('/')[1]),float(metamodel.parametersToVary[0].split('/')[2])]
	xlabel = dimensions[0]
	for point in sample:
	    value = point.parameters[dimensions[0]]
	    x.append(value)
    #/ 3.- plotting
    matplotlib.pyplot.hist(x, 100, normed=1)
    matplotlib.pyplot.xlabel(xlabel)
    matplotlib.pyplot.ylabel('P')
    matplotlib.pyplot.title('Sample Distribution')
    matplotlib.pyplot.xlim(xLim[0],xLim[1])
    matplotlib.pyplot.savefig(outputfiles.samplePlot,dpi=1200)
    
    return None

def plot2D(metamodel, outputfiles, sample):
    
    '''
    This function plots a 2D graphs with the accepted sampled points.
    '''

    #/ 1.- initialising the values
    x = []; y = []
    dimensions = sample[0].parameters.keys()
    xlabel = None; ylabel = None
    xLim = []; yLim = []
    #/ 2.- determining the scale 
    scales = {}
    ranges = {}
    for parameterToVary in metamodel.parametersToVary:
	scales[parameterToVary.split('/')[0]] = parameterToVary.split('/')[3]
	ranges[parameterToVary.split('/')[0]] = [float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2])]
    #/ 3.- appending the values, defining labels and limits
    #/ 3.1.- for the x
    if scales[dimensions[0]] == 'log':
	xlabel = 'log ' + dimensions[0]
	xLim = [math.log10(ranges[dimensions[0]][0]), math.log10(ranges[dimensions[0]][1])]
	for point in sample:
	    xValue = point.parameters[dimensions[0]]
	    x.append(math.log10(xValue))
    else:
	xlabel = dimensions[0]
	xLim = [ranges[dimensions[0]][0], ranges[dimensions[0]][1]]
	for point in sample:
	    xValue = point.parameters[dimensions[0]]
	    x.append(xValue)
    #/ 3.2.- for the y
    if scales[dimensions[1]] == 'log':
	ylabel = 'log ' + dimensions[1]
	yLim = [math.log10(ranges[dimensions[1]][0]), math.log10(ranges[dimensions[1]][1])]
	for point in sample:
	    yValue = point.parameters[dimensions[1]]
	    y.append(math.log10(yValue))
    else:
	ylabel = dimensions[1]
	yLim = [ranges[dimensions[1]][0], ranges[dimensions[1]][1]]
	for point in sample:
	    yValue = point.parameters[dimensions[1]]
	    y.append(yValue)
    #/ 3.- detecting the colour
    z = []
    colours = []
    sizes = []
    for point in sample:
	z.append(math.log10(point.fitnessFunction))
    #/ 3.1.- normalising the values
    minValue = 0
    #/ 3.1.1.- detecting negative values
    for value in z:
	if value < 0:
	    minValue = min(z)
	    break
    #/ 3.1.2.- reconverting the values
    for i in range(len(z)):
	z[i] = z[i] + abs(minValue)
    maxValue = max(z)
    for i in range(len(z)):
	z[i] = z[i]/maxValue*4.
    #/ 3.1.3.- calculating the RGB colours
    for i in range(len(z)):
	R, G, B = rainbowColourCalculator(z[i])
	theColour = (R,G,B,1.)
	colours.append(theColour)
	sizes.append(.1)
    #/ 4.- plotting
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(1,1.125,1) 
    matplotlib.pyplot.scatter(x, y, s=sizes, c=colours, linewidths=0)
    matplotlib.pyplot.title('Sample Distribution')
    matplotlib.pyplot.xlabel(xlabel)
    matplotlib.pyplot.ylabel(ylabel)
    #/ 4.1.- setting the ranges
    matplotlib.pyplot.xlim(xLim[0],xLim[1])
    matplotlib.pyplot.ylim(yLim[0],yLim[1])
    #/ 4.2.- color bar
    a = []
    resolution = 1000
    ax = fig.add_subplot(1,100,100)
    for i in range(resolution):
	a.append([i])
    a = numpy.array(a)
    #/ 4.2.1.- setting the ticks
    matplotlib.pyplot.xticks([])
    labels = []
    theRange = maxValue - minValue
    step = theRange/5.
    value = minValue
    labels.append(str('%.2f'%value))
    for i in range(5):
	value = value + step
	labels.append(str('%.2f'%value))
    matplotlib.pyplot.yticks([0,200,400,600,800,1000],labels)
    matplotlib.pyplot.text(3,500,'log F')
    #/ 4.2.2.- creating the palette
    im = ax.pcolor(a) 
    matplotlib.pyplot.savefig(outputfiles.samplePlot,dpi=1200)

    return None

def plot3D(metamodel, outputfiles, sample):
    
    '''
    This function creates a gnuplot script to plot in 3D the Monte Carlo results.
    '''
    
    #/ 1.- data converter
    f = open(outputfiles.sampleResults, 'r')
    data = f.readlines()
    f.close()
    f = open(outputfiles.gnuplot3DData, 'w')
    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
    scales = {}
    ranges = {}
    for parameterToVary in metamodel.parametersToVary:
	scales[parameterToVary.split('/')[0]] = parameterToVary.split('/')[3]
	ranges[parameterToVary.split('/')[0]] = [float(parameterToVary.split('/')[1]), float(parameterToVary.split('/')[2])]
    dimensions = sample[0].parameters.keys()
    localScales = []
    for dimension in dimensions:
	if scales[dimension] == 'log':
	    localScales.append('log')
	else:
	    localScales.append('lin')
    localScales.append('log') #/ for the fitness function
    for datum in data:
	if datum[0] != '#':
	    vector = datum.split()
	    for i in range(len(localScales)):
		if localScales[i] == 'log':
		    vector[i] = str(math.log10(float(vector[i])))
	    for element in vector:
		f.write('%s\t'%element)
	    f.write('\n')
    f.close()
    #/ 2.- detecting the ranges, labels 
    xLabel = None; yLabel = None; zLabel = None
    xLim = []; yLim = []; zLim = []
    #/ 2.1.1.- for the x
    if scales[dimensions[0]] == 'log':
	xLabel = 'log ' + dimensions[0]
	xLim = [math.log10(ranges[dimensions[0]][0]), math.log10(ranges[dimensions[0]][1])]
    else:
	xLabel = dimensions[0]
	xLim = [ranges[dimensions[0]][0], ranges[dimensions[0]][1]]
    #/ 2.1.2.- for the y
    if scales[dimensions[1]] == 'log':
	yLabel = 'log ' + dimensions[1]
	yLim = [math.log10(ranges[dimensions[1]][0]), math.log10(ranges[dimensions[1]][1])]
    else:
	yLabel = dimensions[1]
	yLim = [ranges[dimensions[1]][0], ranges[dimensions[1]][1]]
    #/ 2.1.3.- for the z
    if scales[dimensions[2]] == 'log':
	zLabel = 'log ' + dimensions[2]
	zLim = [math.log10(ranges[dimensions[2]][0]), math.log10(ranges[dimensions[2]][1])]
    else:
	zLabel = dimensions[2]
	zLim = [ranges[dimensions[2]][0], ranges[dimensions[2]][1]]
    #/ 3.- creating the code for the gnuplot figure
    f = open(outputfiles.gnuplot3DPlotter, 'w')
    f.write('#\n# generated by ByoDyn version %s\n#\n'%initiator.BYODYNVERSION)
    f.write('set term postscript color\n')
    f.write('set output "%s"\n'%outputfiles.samplePlot)
    f.write('set nokey\n')
    f.write('set pm3d\n')
    f.write('set palette rgbformulae 30,-13,-23\n')
    f.write('set pointsize 1\n')
    f.write('set xtics border out scale 2,2 mirror norotate offset character 0, -.5, 0\n')
    f.write('set ytics border out scale 2,2 mirror norotate offset character 0, -.5, 0\n')
    f.write('set label "log F" at screen .9, screen .4\n')
    f.write('set xrange [%s:%s]\n'%(xLim[0],xLim[1]))
    f.write('set yrange [%s:%s]\n'%(yLim[0],yLim[1]))
    f.write('set zrange [%s:%s]\n'%(zLim[0],zLim[1]))
    f.write('set title "Sample Distribution"\n')
    f.write('set xlabel "%s"\n'%xLabel)
    f.write('set ylabel "%s"\n'%yLabel)
    f.write('set zlabel "%s"\n'%zLabel)
    f.write('splot "%s" u 1:2:3:4  w points palette pt 7 \n'%outputfiles.gnuplot3DData)
    f.close()
    #/ 4.- calling the gnuplot
    os.system('gnuplot 2> %s %s'%(outputfiles.gnuplotErrorMessagesFile,outputfiles.gnuplot3DPlotter))
    
    return None

def plotResults(metamodel, outputfiles, sample):

    '''
    This function creates a 1D, 2D or 3D graphs if the parameter space is 1, 2 or 3 dimensional.
    Otherwise it gives a message expressing the impossibility of the graph.
    '''
    
    #/ 1.- detecting the number of parameters
    if len(metamodel.parametersToVary) == 0:
        raise errorMessages.ClassSamplerException, 'the number of parameters to vary cannot be zero.'
    elif len(metamodel.parametersToVary) == 1:
        plot1D(metamodel, outputfiles, sample)
    elif len(metamodel.parametersToVary) == 2:
        plot2D(metamodel, outputfiles, sample)
    elif len(metamodel.parametersToVary) == 3:
        plot3D(metamodel, outputfiles, sample)
    else:
        print 'The number of dimensions of the explored space is %s and therefore it cannot be plotted.'%len(metamodel.parametersToVary)
    
    
    return None
