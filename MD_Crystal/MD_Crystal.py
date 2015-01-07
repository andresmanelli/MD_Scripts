# coding=UTF-8
"""
  MD_Crystal.py
  
  Authors:  Ardiani, Franco
            Manelli, Andrés
  
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

  Python Module MD_Crystal
  
"""

from ovito.io import *			# Import OVITO modules.
from ovito.modifiers import *

from MD_Crystal_Defs import *

"""Global variables"""
lammpstrj_file = '' # Keep track of lammpstrj file
node = 0 # Used for importing lammpstrj files
numOfParticles = 0 # Used for storing the number of particles
workDir = '~/' # Used for writing files
ids = [] # Particle ID
types = [] # Particle type
positions = [-1,[]] # Particle position
struct_type = [-1,[]] # Particle structure type (FCC...)

import sys
import numpy
import os.path

def insc_sphere(par, tol, struct, frame):
	"""Finds the biggest sphere with given crystal structure with center
	in [cx = par[0], cy = par[1], cz = par[2]], with the specified 
	tolerance. Initial radius is set to r = par[3]
	
	Keyword arguments:
	par -- Sphere parameters. [cx, cy, cz, radius]
	tol -- The particles in sphere will be (tol)% type (struct)
	struct -- The type of crystal, e.g. MD_Cystal.FCC
	frame -- Frame to analyze
	"""
	global struct_type
	
	if len(par) is not 4:
		print 'Parameters of sphere are not ok..'
		return 0
		
	cx = par[0]
	cy = par[1]
	cz = par[2]
	r = par[3]
	
	if frame is not struct_type[0]:
		baaMod = BondAngleAnalysisModifier()
		node.modifiers.append(baaMod)
		ovito.dataset.anim.current_frame = frame
		
		node.compute()
		struct_type = [frame,node.output['Structure Type'].array]
		
	tolAct = 0
	local_positions = node.output['Position'].array
	filtered = struct_type[1]
	while tolAct < tol:		
			
		f2 = []
		lp2 = []
		for i, p in enumerate(local_positions):
			if ((p[0] - cx)**2 + (p[1] - cy)**2 + (p[2] - cz)**2) <= r**2:
				f2.insert(len(f2),filtered[i])
				lp2.insert(len(lp2),p)
				
		local_positions = lp2
		filtered = f2
		
		N = len(filtered)
		
		if N == 0:
			print 'insc_sphere(',par,',',tol,',',struct,',',frame,'):'
			print 'Reached limit of N == 0. No particles left. Aborting..'
			break
		
		n_struct = numpy.bincount(filtered)
		# If last types are not present, length of n_struct is smaller than 
		# length of STRUCT_TYPES
		for i in range(len(STRUCT_TYPES)-len(n_struct)+1):
			n_struct = numpy.insert(n_struct,len(n_struct),[0])
		
		n_struct = n_struct[struct]
			
		tolAct = 100*n_struct/float(N)
		#print 'N:',N
		#print 'Radius:',r
		#print 'Number of FCC atoms:', n_struct
		#print 'Number of OTHER atoms:', N-n_struct
		#print 'Porcentage of',struct,'structure:',tolAct,'%'
		r -= 0.5
		
	return r	
	
def init(f):
	"""Sets the file for analysis"""
	# TODO: Init also the voronoi reference indices
	global lammpstrj_file, workDir, node, numOfParticles, positions
	
	lammpstrj_file = os.path.realpath(f)	
	if os.path.isfile(lammpstrj_file) is not True:
		print "The file",lammpstrj_file,"does not exists..."
		sys.exit()
	
	# Working directory. Not the one where we executed MD_Voronoi,
	# the path where we read lammpstrj_file
	workDir = os.path.dirname(os.path.realpath(lammpstrj_file))
		
	# Load a simulation snapshot of a Cu-Zr metallic glass.
	node = import_file(lammpstrj_file)
				
	# Number of particles
	numOfParticles = node.source.data['Position'].size
