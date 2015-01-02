#  MD_Voronoi_Headers.py
#  coding: UTF-8
#  
#  Authors:  Ardiani, Franco
#            Manelli, Andrés
#  
#  Copyright 2015 Andres <andres@magic>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  ------------------------------------------------------------------
#
#  Python Module MD_Voronoi
#  
#  

import time

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
			'# \n\n' + \
			'# Reference: \n'
			
	for i in range(len(arrVoroIndices)):
		ret += ('# Type %i: %s\n' % (i, tuple(arrVoroIndices[i])))
		
	ret += '# \n\n'

	return ret