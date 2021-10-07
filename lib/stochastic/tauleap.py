#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Pau Rue-Queralt
#
#  Created: 2008-06-05 by Pau Rue-Queralt
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

# $Id: tauleap.py,v 1.4 2008/12/11 15:00:39 paurue Exp $

## \file 
# This module is contains the the stochastic simulation algorithms for the tau-leap method.

import numpy as np

def simulate (evalPropensities, stoichiometry, x_0, time, tau, seed = None):

    '''
    Gillespie's Poisson tau-leap stochastic simulation algorithm
        evalPropensities  list of propensity functions
        stoichiometry     stoichiometrty matrix, each column corresponds to a reaction
        x_0               initial conditions
        time              simulation time
        tau               the time step for which the system state is written (if set 
                          to None all reactions are written in the output)
        seed              for the pseudo-random number generator
    '''

    def poisson(x):

        '''
        This function takes a list of (nonnegative) real numbers and samples a 
        Poisson random number for each of them with mean and variance given by 
        this. If any number in the list is negative, it is treated as zero.
        '''

        xt = np.maximum(x,0)

        return np.array(map(np.random.poisson, xt)).T

    if seed != None:
        np.random.seed(seed)

    nsteps = np.ceil(time/tau) + 1
    x = np.array(x_0, dtype = np.int64)
    tHist = tau * np.arange(nsteps)
    xHist = np.zeros([x.size, nsteps], dtype = np.int64)
    xHist[:,0] = x
    i = 1
    while(i < nsteps and sum(abs(x)) > 0):
        k = poisson(tau * evalPropensities(x))
        update = np.dot(stoichiometry, k)
        x += update
        xHist[:,i] = x
        i = i + 1

    return xHist.T, tHist

