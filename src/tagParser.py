#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author: Adrian L. Garcia-Lomana
#
#  Created: 2005-08-11 by Adrian L. Garcia-Lomana
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

# $Id: tagParser.py,v 4.3 2008/12/03 16:37:38 alglomana Exp $

## \file 
# This module is the paser for the tag format files.
# A system model is obtained.

class ClassModelTags:

    '''
    Class for the model.
    '''

    def __init__(self):

	'''
	The constructor.
	'''

        self.systemName = None
        self.xlength = 0
        self.ywidth = 0
        self.nodes = []
        self.topology = {}
        self.initialConditions = {}
        self.parameters = {}
        #/ Next fields are not defined at the tag format, but they are necessary for ClassModelSBML compatibility
        self.delayFunctions = False
        self.algebraicNodes = {}
        self.constantNodes = {}
        self.nonConstantParameters = {}
        self.compartments = {}
        self.nonConstantCompartments = {}
        self.functions = {}
        self.rules = []
        self.events = []
               
        return None     
        
    def readInput(self, metamodel):

	'''
	This method gets the information of the input file and it sets the object.
	'''
    
        #/ 1.- reading file
        modelFile = open(metamodel.modelFile, 'r')
        #/ 2.- obtaining information
        for line in modelFile:
            fields = line.split()
            if fields[0] == "systemName":
                self.systemName = fields[1]
            elif fields[0] == "xlength":
                self.xlength = int(fields[1])
            elif fields[0] == "ywidth":
                self.ywidth = int(fields[1])
            elif fields[0] == "nodes":
                self.nodes = fields[1:]
            elif fields[0] == "topology":
                self.topology[str(fields[1])] = fields[2:]
            elif fields[0] == "initialCondition":
                self.initialConditions[str(fields[1])] = fields[2:] 
        #/ 3.- closing file
        modelFile.close()
        #/ 4.- avoiding redundancy on the parameters
        rawParameters = []
        for definition in self.topology.keys():
            for parameter in self.topology[definition]:
                fieldsParameter = parameter.split("/")
                rawParameters.append(fieldsParameter)
        repeated = 0
        parameters = []
        for i in range(0, int(len(rawParameters))):
            for j in range(int(i+1), int(len(rawParameters))):
                if (rawParameters[i] == rawParameters[j]):
                    repeated = repeated + 1
            if (repeated == 0):
                parameters = rawParameters[i]
                self.parameters['%s'%(parameters[0])] = float(parameters[1])
                repeated = 0
            else:
                repeated = 0
        #/ 5.- obtaining a list of the initial conditions ordered by node
        valuesOrderedByNode = []
        valuesOrderedByCell = []
        for node in self.nodes:
            for definition in self.initialConditions.keys():
                if definition == node:
                    for parameter in self.initialConditions[definition]:
                        valuesOrderedByNode.append(float(parameter.split("/")[1]))
        for i  in range (0, self.xlength * self.ywidth): # inverting the matrix
            for j in range (0, len(self.nodes)):
                valuesOrderedByCell.append(valuesOrderedByNode[(j * self.xlength * self.ywidth) + i])
        self.initialConditions = valuesOrderedByCell
        
        return self
    
    def summary(self, outputfiles):
	
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
        f.write('\tConstant Parameters: %s ;' %(len(self.parameters.keys())))
	for constant in self.parameters.keys():
	    f.write('\t%s=%s'%(constant, self.parameters[constant]))
	f.write('\n')
	f.close()
	
	return None
