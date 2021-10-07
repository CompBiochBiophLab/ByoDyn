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

# $Id: gssa.py,v 1.11 2008/12/12 17:34:59 alglomana Exp $

## \file 
# This module is contains the the stochastic simulation algorithms for the Guillespie Stochastic Simulation Algorithm (SSA).

import numpy as np

def simulate (evalPropensities, stoichiometry, x_0, time, hurdle, seed = None):

    '''
    Gillespie's Stochastic Simulation Algorithm (SSA) 
        evalPropensities  list of propensity functions
        stoichiometry     stoichiometrty matrix, each column corresponds to a reaction
        x_0               initial conditions
        time              simulation time
        hurdle            the time step for which the system state is written (if set 
                          to None all reactions are written in the output)
        seed              for the pseudo-random number generator
    '''
    if seed != None:
        np.random.seed(seed)
    if not hurdle:
        nreactions = 0
        nreactionsAllocated = 10000
        #/ each time we reach the allocated space we add nallocationBlock reactions
        nallocationBlock = 5000
        x = np.array(x_0, dtype = np.int64)
        tHist = np.zeros(nreactionsAllocated)
        xHist = np.zeros([x.size, nreactionsAllocated], dtype = np.int64)
        t = 0
        i = 1
        while(t <= time and np.sum(x) >= 0):
            propensities = evalPropensities(x)
            totalPropensity = np.sum(propensities)
            if totalPropensity == 0:
                print "No more reactions are likely to occur."
                break
            relativePropensities = np.array(propensities) / totalPropensity
            cumrelPropensities = np.cumsum(propensities) / totalPropensity
            tau = np.random.exponential(1. / totalPropensity)
            r = np.random.random()
            j = np.where(cumrelPropensities > r)[0][0]
            t = t + tau
            x += stoichiometry[:,j]
            xHist[:,i] = x
            tHist[i] = t
            i = i + 1
            if i == nreactionsAllocated:
                nreactionsAllocated = nreactionsAllocated + nallocationBlock
                tHist.resize(nreactionAllocated)
                xHist.resize([x.size, nreactionAllocated])
        tHist.resize(i)
        xHist = xHist[:,0:i]
    else:
        nsteps = np.ceil(time/hurdle) + 1
        x = np.array(x_0, dtype = np.int64)
        tHist = hurdle * np.arange(nsteps)
        xHist = np.zeros([x.size, nsteps], dtype = np.int64)
        xHist[:,0] = x
        nextHurdle = hurdle
        i = 1
        t = 0
        while(t <= time and np.sum(x) >= 0):
            propensities = evalPropensities(x)
            totalPropensity = np.sum(propensities)
            if totalPropensity == 0:
                print "No more reactions are likely to occur."
                break
            relativePropensities = np.array(propensities) / totalPropensity
            cumrelPropensities = np.cumsum(propensities) / totalPropensity
            tau = np.random.exponential(1. / totalPropensity)
            r = np.random.random()
            j = np.where(cumrelPropensities > r)[0][0]
            t = t + tau
            while t >= nextHurdle and i < nsteps:
                #/ if t reaches nextHurdle, the system state at t=nextHurdle is
                #/ the system state before updating. That applies for all the 
				#/ following hurdles t reaches.
                xHist[:,i] = x
                i = i + 1
                nextHurdle = nextHurdle + hurdle
            x += stoichiometry[:,j]
        tHist.resize(i)
        xHist = xHist[:,0:i]

    return xHist.T, tHist

