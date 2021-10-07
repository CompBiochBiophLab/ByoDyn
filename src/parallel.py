#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Alex Gomez-Garrido, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Miguel A. Hernandez
#
#  Created: 2007-06-13 by Miguel A. Hernandez
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

# $Id: parallel.py,v 4.2 2008/12/03 16:37:38 alglomana Exp $

# \file 
# This module holds the parallel functions called by the program.

import pickle
try:
    from Scientific import MPI
except ImportError:    
    communicator = None
else:
    try:
        communicator = MPI.world.duplicate() 
        if (communicator.rank == 0):
            communicator.send("AABBCCDDDWWXXXYYYZZ", 0, 0)
            communicator.receiveString(0, 0)
    except Exception:
        communicator = None 
    else:
        communicator = MPI.world.duplicate() 
        if (communicator.size == 1):
            communicator = None
internalBuffer = {}

def runsOnParallelMachine():
    
    '''
    This function determines if the program runs in parallel or serial environment.
    It returns a boolean variable.
    '''

    if (communicator == None):
        return False
    else:
        return True

def mainProcessor():

    '''
    This function discriminates if the current processor is the main processor, that is, processor 0
    It returns a boolean variable.
    '''
    
    if (communicator == None):
        return True
    if (communicator.rank == 0):
        return True
    return False

def currentProcessor():
    
    '''
    This function provides which is the current processor.
    It returns an integer.
    '''
    
    if (communicator == None):
        return 0
    return communicator.rank

def totalProcessors():
    
    '''
    This function provides the amount of processors.
    It returns an integer.
    '''

    if (communicator == None):
        return 1

    return communicator.size

def sendAll(Object, Tag):

    '''
    This function send a given object to all processors.
    It takes as arguments:
    Object, the object to be sent, and
    Tag, an aditional identifier of the object.
    It returns nothing.
    '''
    
    strObject = pickle.dumps(Object)
    if (communicator != None):
        for i in range (0, communicator.size, 1):
            communicator.send(strObject, i, Tag)
    else:
        internalBuffer[Tag] = strObject    

    return None

def send(Object, Destination, Tag):
    
    '''
    This function sends an object to a defined processor.
    It takes as arguments:
    Object, the object, 
    Destination, the destination processor, and
    Tag, an aditional identifier of the object.
    It returns nothing.
    '''

    strObject = pickle.dumps(Object)
    if (communicator != None):
        communicator.send(strObject, Destination, Tag)
    else:
        internalBuffer[Tag] = strObject    
    
    return None

def waitProcessors():
    
    '''
    This function makes all the processors to wait up to that point.
    '''

    if (communicator != None):
        communicator.barrier()
    
    return None

def receiveAny(Tag):

    '''
    This function receives an object from any processor.
    It takes as argument, Tag, an aditional identifier of the object.
    It returns the object.
    '''

    if (communicator != None): 
        strObject, source, tag = communicator.receiveString(None, Tag)
    else:
        strObject = internalBuffer[Tag]
    Object = pickle.loads(strObject)
    
    return Object

def receive(Source, Tag):
    
    '''
    This function receives an object from a determined processor.
    It takes as arguments:
    Source, the source processor, and
    Tag, an aditional identifier of the object.
    It returns the object.
    '''

    if (communicator != None): 
        strObject, source, tag = communicator.receiveString(Source, Tag)
    else:
        strObject = internalBuffer[Tag]
    Object = pickle.loads(strObject)
    
    return Object
