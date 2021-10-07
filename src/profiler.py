#! /usr/bin/env python

#
#  Project: ByoDyn
#
#  Copyright (C) 2008 Michael Johston, Adrian L. Garcia-Lomana and Jordi Villa-Freixa
#
#  Author: Michael Johston
#
#  Created: 2004-12-15 by Michael Johston
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

# $Id: profiler.py,v 4.3 2008/12/03 16:37:38 alglomana Exp $

## \file
# This module profiles the program.
# In order to run use Profile.py [sciptname] [scriptargs]

import sys, profile, pstats
import central

scratchDir = os.environ.get('BYODYN_SCRATCH_DIR')
    
SCRIPTNAME = central
ERROR_FILE = '%s/profilingErrorFile' % str(scratchDir)
PROFILE_FILE = '%s/profilingOutputFile' % str(scratchDir)
STAT_FILE = '%s/profilingStatFile' % str(scratchDir)

def main():

	'''
	This is the main function of the profiling. 
	''' 

	theString = '%s.Profile()' % SCRIPTNAME.__name__
	code = compile(theString, ERROR_FILE, 'exec')
	print 'Beginning execution of %s for profiling...' % SCRIPTNAME
	profile.run(code, PROFILE_FILE)
	extraRuns = 0 # number of groups - 1 for the profiling
	fileList = []
	for i in range(1, extraRuns+1):
		filename = "%s%d" % (PROFILE_FILE, i)
		fileList.append(filename)
		profile.run(code, filename)
	StatMaker = pstats.Stats(PROFILE_FILE)
	for filename in fileList:
		StatMaker.add(filename)
	StatMaker.sort_stats('time')	
	StatMaker.print_stats()

if __name__ == '__main__':
	main()
