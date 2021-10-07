#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
#
#  Author:  Adrian L. Garcia-Lomana and Alex Gomez-Garrido
#
#  Created: 2004-08-30 by Adrian L. Garcia-Lomana
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

# $Id: exporter.py,v 4.2 2008/12/03 16:37:38 alglomana Exp $

## \file
# This module contains the code for exporting SBML files.

import libsbml
import sbmlWorker

def central(metamodel, model, outputfiles):

    '''
    This central function directs the flow of the exporting module.
    '''

    #/ 1.- Determining new parameter values
    for parameter in metamodel.parameters:
        for modelParameter in model.parameters:
            if parameter == modelParameter:
        	model.parameters[parameter] = metamodel.parameters[parameter]
    #/ 2.- writing the model
    w = libsbml.SBMLWriter()
    file = outputfiles.outputdir + '/%s.xml' %model.systemName
    sbmlWorker.sbmlWriter(w, model, metamodel, file)
    

    return None
