# coding=UTF-8
"""
  MD_Voronoi_Headers.py
  
  Authors:  Ardiani, Franco
            Manelli, Andrés
  
  Copyright 2015 Andres <andres@magic>
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
  MA 02110-1301, USA.
  
  ------------------------------------------------------------------

  Python Module MD_Voronoi

"""  

import time
import os

from MD_Voronoi_Defs import *

def voro_dump_header(f):
	ret = 	'# \n' + \
			'# Basic Voronoi analysis dump file for: \n' + \
			'# ' + f + '\n' + \
			'# \n' + \
			'# ' + time.strftime("%d/%m/%Y") + '\n' + \
			'# ' + time.strftime("%H:%M:%S") + '\n' + \
			'# Created with module MD_Voronoi \n' + \
			'# https://github.com/andresmanelli/MD_Scripts.git \n' + \
			'# \n' + \
			'# Reference: \n'
			
	for i in range(len(arrVoroIndices)):
		ret += ('# Type %i: %s\n' % (i+1, tuple(arrVoroIndices[i])))
	
	ret += '# Type 0: OTHER\n'
	ret += '# \n\n'

	return ret

def voro_hist_header():
	ret = 	'# \n' + \
			'# Voronoi type histogram \n' + \
			'# \n' + \
			'# ' + time.strftime("%d/%m/%Y") + '\n' + \
			'# ' + time.strftime("%H:%M:%S") + '\n' + \
			'# Created with module MD_Voronoi \n' + \
			'# https://github.com/andresmanelli/MD_Scripts.git \n' + \
			'# \n' + \
			'# Reference: \n'
			
	for i in range(len(arrVoroIndices)):
		ret += ('# Type %i: %s\n' % (i+1, tuple(arrVoroIndices[i])))
	
	ret += '# Type 0: OTHER\n'
	ret += '# \n# frame '
	for i in range(len(arrVoroIndices)+1):
		ret += ('n_type_%i ' % (i))

	ret += '\n\n'

	return ret

def voro_changes_header():
	ret = 	'# \n' + \
			'# Voronoi type-change histogram \n' + \
			'# \n' + \
			'# ' + time.strftime("%d/%m/%Y") + '\n' + \
			'# ' + time.strftime("%H:%M:%S") + '\n' + \
			'# Created with module MD_Voronoi \n' + \
			'# https://github.com/andresmanelli/MD_Scripts.git \n' + \
			'# \n' + \
			'# Reference: AB means from type A to type B\n' + \
			'# e.g.: 34 --> from type 3 to type 4\n' + \
			'# frame '
	for i in range(100):
		ret += ('%i ' % (i))

	ret += '\n\n'

	return ret
