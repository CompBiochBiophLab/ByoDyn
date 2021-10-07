#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana, Alex Gomez-Garrido and Jordi Villa-Freixa
#
#  Authors: Adrian L. Garcia-Lomana and Alex Gomez-Garrido
#
#  Created: 2004-11-10 by Adrian L. Garcia-Lomana
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

# $Id: initiator.py,v 4.9 2008/12/17 00:55:44 alglomana Exp $

## \file initiator.py
# ByoDyn is an open source computational package aimed at studying the dynamical behaviour of small to massive biochemical networks. 
# Models are input in the standard format of systems biology markup language (SBML).
# The model can be simulated, the sensitivity and the identifiability of the system with respect to the parameters can be analysed and finally kinetic parameters can be estimated using experimental time course data.
# Several state of the art optimisation algorithms have been implemented for this purpose. 
# ByoDyn can run in parallel for the some of them using MPI. 


#/ Importing python modules 
import os, sys, shutil
from getopt import *
#/ Adding byodyn folders to the PATH
srcdir = os.environ.get('BYODYN_PATH') + '/src'
benchmarkdir = os.environ.get('BYODYN_PATH') + '/benchmark'
libdir = os.environ.get('BYODYN_PATH') + '/lib'
sys.path.append('%s'%benchmarkdir)
sys.path.append('%s'%libdir)
sys.path.append('%s'%srcdir)
#/ Importing byodyn modules
import central, checker, parallel, errorMessages
#/
#/ Functions 
#/
def createExamples():

    ''' 
    The function creates the examples of ByoDyn at the examples directory. 
    '''

    #/ hard coded variables
    dirExamples = 'examples'
    examples = []
    #/
    if os.path.exists(dirExamples) == True:
        print 'The directory "examples" already exists. Remove it first if you want to start from the begining.'
    else:
        print 'Creating examples ...', 
        examplesDir = os.environ.get('BYODYN_PATH') + '/examples/quickStart'
        localExamplesDir = 'examples'
        shutil.copytree(examplesDir, localExamplesDir)
	shutil.rmtree(localExamplesDir+'/CVS') #/ removing the 'CVS' directory
        print 'done.'

    return None

def printingHelp(): 

    ''' 
    This function prints the command line help. 
    '''

    print 'Create the documentation at the docs directory for further information.'

    return None

def printingVersion(version):

    ''' 
    This function prints the current ByoDyn version. 
    '''

    print 'You are running ByoDyn version', version

    return None

def versionDefinitor():
    
    '''
    This funciton defines the ByoDyn version globally.
    '''

    cvsVersion = "$Id: initiator.py,v 4.9 2008/12/17 00:55:44 alglomana Exp $"
    global BYODYNVERSION
    BYODYNVERSION = cvsVersion.split()[2]
    
    return None

#/
#/ end functions
#/

#/
#/ starting main fuction of the file
#/

def initial(runnerFile = None):

    ''' 
    This function reads the command line options, it configures the ByoDyn options and checks its consistency and calls the modules of the source directory of ByoDyn. 
    '''
    
    #/ 1.- initializing the error handling. Required for stable releases.
    sys.excepthook = errorMessages.ErrorHandler
    #/ 2.- defining the ByoDyn version
    versionDefinitor()
    #/ 3.- first message to the world 
    if (parallel.mainProcessor()):
        print 'ByoDyn version', BYODYNVERSION, 'is running ...'
    #/ 4.- initializing getoptions
    pair = getopt(sys.argv[1:], 'o:tevhl', ('output=', 'testing', 'examples', 'version', 'help', 'limit'))
    options = pair[0]
    #/ 5.- variables and its default values
    output = ''
    testing = False
    examples = False
    version = False
    help = False
    limited = False
    in_parameters = sys.argv[1:] #/ In this variable we add all input arguments, later we remove the getopts, and the finish only should have the runner file
    #/ PARALLEL: ONE PROCESSOR
    parallelFlag = False
    if (parallel.mainProcessor()):
    #/ PARALLEL: IDENT ALL
        if runnerFile != None:
    	    in_parameters.append(runnerFile)
        #/ 6.- Reading getoptions
        if len(in_parameters) == 0:		# If we don't have a runner file, print help
            print 'Run byodyn -h for help.'
        else:
            for opt in options:
                if opt[0] == '-o' or opt[0] == '--output':
                    output = opt[1]
                    in_parameters.pop(in_parameters.index(opt[0]))
                    in_parameters.pop(in_parameters.index(opt[1]))
                if opt[0] == '-t' or opt[0] == '--testing': 
                    testing = True
                    in_parameters.pop(in_parameters.index(opt[0]))
                if opt[0] == '-e' or opt[0] == '--examples': 
                    examples = True
                    in_parameters.pop(in_parameters.index(opt[0]))
                if opt[0] == '-v' or opt[0] == '--version':
                    version = True
                if opt[0] == '-h' or opt[0] == '--help':
                    help = True	
                if opt[0] == '-l' or opt[0] == '--limit':
                    limited = True    
                    in_parameters.pop(in_parameters.index(opt[0]))
            #/        
            if version == True:
                printingVersion(BYODYNVERSION)
            elif help == True:
                printingHelp()
            else:
		#/ 7.- setting the output and scratch directories
		if output == '':
		    outputDir = os.getcwd() + '/output'
		    scratchDir = outputDir + '/scratch'
		else:
		    if output[0] == '/':
			outputDir = output + '/output'
			scratchDir = outputDir + '/scratch'
		    else:
			outputDir = os.getcwd() + '/' + output
			scratchDir = outputDir + '/scratch'
		os.environ['BYODYN_OUTPUT_DIR'] = outputDir
		os.environ['BYODYN_SCRATCH_DIR'] = scratchDir
                #/ 8.- running the examples without creating output nor scratch directories
                if examples == True:
                    createExamples()
                else:
		    #/ 9.1- removing the output directories in the case of testing
		    if testing == True:
			if os.path.exists(outputDir) == True:
			    shutil.rmtree(outputDir)
		    #/ 9.2.- creating the output directories
                    if os.path.exists(outputDir) == False:
                        os.mkdir(outputDir)
                    if os.path.exists(scratchDir) == False:
                        os.mkdir(scratchDir)
                    #/ 9.3 - sending to the checking modules
                    if testing == True:
                        checker.main()
                    else:
                        #/ PARALLEL: INCLUDE FLAG; COMMENT FUNCTION
                        parallelFlag = True                  
    #/ PARALLEL: SPLIT AND CALL FUNCTION
    if (parallel.mainProcessor()):
        parallel.sendAll(parallelFlag, 0)
    if (parallel.mainProcessor() and parallelFlag == True):
        parallel.sendAll(in_parameters, 1)
        parallel.sendAll(outputDir, 3)
        parallel.sendAll(scratchDir, 4)
	#/
    parallelFlag = parallel.receiveAny(0)
    if (parallelFlag == True):
        outputDir = parallel.receiveAny(3)
        scratchDir = parallel.receiveAny(4)
        os.environ['BYODYN_OUTPUT_DIR'] = outputDir
        os.environ['BYODYN_SCRATCH_DIR'] = scratchDir
	#/
        in_parameters = parallel.receiveAny(1)
        central.main(in_parameters[0])
    if (parallel.mainProcessor()):
	#/ 10.- last message to the world
        print '... exiting from ByoDyn.' 

    return None
