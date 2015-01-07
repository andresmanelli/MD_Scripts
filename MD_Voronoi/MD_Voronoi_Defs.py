# coding=UTF-8
"""
  MD_Voronoi_Defs.py
  
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

import numpy

# Defines voronoi cells
# [0, 0, 0, 0, 12, 0] => 12 faces with 5 edges each.
# [0, 0, 0, 3, 6, 3]  => 3 faces with 4 edges, 6 faces with 5 edges, 3 faces with 6 edges
# Note: The first two indices are always zero because there are no faces with less than three edges.

arrVoroIndices_A = numpy.array([[0,0,5,2,6,0], \
								[0,0,4,4,6,0]])

arrVoroIndices_B = numpy.array([[0,0,0,0,12,0], \
								[0,0,0,2,8,2],	\
								[0,0,0,2,8,1],	\
								[0,0,0,3,6,3],	\
								[0,0,0,3,6,4], 	\
								[0,0,0,1,10,2],	\
								[0,0,0,2,8,3],	\
								[0,0,0,4,4,3],	\
								[0,0,1,0,9,3]])

arrVoroIndices = arrVoroIndices_A
								
# Creates a 'view' of array <<arr>> that points to the same data, but represent it as int32
# TODO: Verify, is this architecture dependant?
VoroIndices_Int32 = numpy.ascontiguousarray(arrVoroIndices).view([('', numpy.int32)] * 6)

# TODO: in init set these from info about lammpstrj file
CU = 1
ZR = 2
EL_NAMES = ['_DUMMY','_CU','_ZR']
