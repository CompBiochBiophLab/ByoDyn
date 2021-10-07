#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Alex Gomez-Garrido and Adrian L. Garcia-Lomana
#
#  Created: 2007-12-11 by Alex Gomez-Garrido
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

# $Id: simulatorRungeKutta.py,v 4.2 2008/12/03 16:37:38 alglomana Exp $

## \file 
# This module is responsible of the Runge-Kutta method for the numerical integrations.

import os, re, sys, copy
import formulas, initiator, errorMessages, sbmlWorker, simulatorEuler
from affectors import *

class ClassSimulatorRungeKutta(simulatorEuler.ClassSimulatorEuler):

    '''
    Class for the Runge-Kutta simulator.
    It heritages from ClassSimulatorEuler.
    '''    

    def integrate(self, f, x_0, time, dt):

	'''
	This method directs the integration
	'''
	
	#/ 1.- Definition of the initial variables
        t = [0.0]
        x = []
        x.append(x_0)
	#/ 2.- Runge-Kutta method up to 4 orders
        for n in range(int(time/dt)+1):
            t.append(t[n] + dt)
            xdot = f(x[n])
            k1 = []
            k1x = []
            for i in range(len(xdot)):
                k1.append(dt*xdot[i])
                k1x.append(x[n][i] + k1[i]/2.0) 

            k2xdot = f(k1x)
            k2 = []
            k2x = []
            for i in range(len(k2xdot)):
                k2.append(dt/2*k2xdot[i])
                k2x.append(x[n][i] + k2[i]/2.0) 
            
            k3xdot = f(k2x)
            k3 = []
            k3x = []
            for i in range(len(k3xdot)):
                k3.append(dt/2*k3xdot[i])
                k3x.append(x[n][i] + k3[i])
            
            k4xdot = f(k3x)
            k4 = [] 
            for i in range(len(k4xdot)):
                k4.append(dt*k4xdot[i])                        
            
            x.append([])
            for i in range(len(x[n])):
                x[n+1].append(x[n][i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6)

        return x, t
