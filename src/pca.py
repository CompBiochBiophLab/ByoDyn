#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Alex Gomez-Garrido and Adrian L. Garcia-Lomana
#
#  Created: 2008-12-03 by Alex Gomez-Garrido
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

# $Id: pca.py,v 1.3 2008/12/03 16:37:38 alglomana Exp $

# \file 
# This module holds the principal component analysis algorithms.

import matrixWorker
import scipy, math
from scipy import linalg, array

class ClassPCA:
    '''
    Class for the Principal Component Analysis
    '''  
    
    def __init__(self):
        
        '''
        This is the constructor.
        '''
        self.X = matrixWorker.ClassMatrix()
        self.B = matrixWorker.ClassMatrix()
        self.COV = matrixWorker.ClassMatrix()
        self.eigenValues  = []
        self.eigenVectors = matrixWorker.ClassMatrix()
        self.eigenList = []
        
    def setDataSet(self, dataMatrix):
        self.X.m = dataMatrix
        return None

    def setDeviationsFromMean(self):
        x_mean = scipy.array(self.X.getEmpiricalMean()).reshape([len(self.X.m), 1])
        print x_mean
        ones = scipy.ones((len(self.X.m[0]),1), float).reshape([1,len(self.X.m[0])])
        self.B.m = scipy.array(self.X.m) - scipy.dot(x_mean, ones)
        print self.B.m
        return None

    def setCOV(self):
        self.COV.m = self.B.getCovarianceMatrix()
        #self.COV.setEigenValues()
        return None
    
    def setEigen(self):
        self.COV.setEigenValues()
        self.eigenValues = self.COV.eigenValues
        self.eigenVectors.m = self.COV.eigenVectors
        self.eigenList = zip(self.eigenValues, self.eigenVectors.m.transpose())
        self.eigenList.sort( reverse = True )
        return None
        
def central():
    #data = [[1.,5., 2.],[2.,4.,1.], [3.,6.,2.], [1.,5., 9.]]
    data = [[2.5,0.5,2.2, 1.9, 3.1, 2.3, 2., 1., 1.5, 1.1],
            [2.4, 0.7, 2.9, 2.2, 3.0, 2.7, 1.6, 1.1, 1.6, 0.9]]
    
    pcaCalculator = ClassPCA()
    pcaCalculator.setDataSet(data)
    pcaCalculator.setDeviationsFromMean()
    pcaCalculator.setCOV()
    pcaCalculator.setEigen()
    #print pcaCalculator.eigenVectors.m
    #print pcaCalculator.eigenVectors.m[:,0]
    for l in pcaCalculator.eigenList:
        print l
    
central()
