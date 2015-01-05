# coding=UTF-8
"""
  MD_Voronoy.py
  
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

  Python Module MD_Voronoi
  
"""  

# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import *

import sys
import numpy
import os.path

# Import Voronoi files headers
import MD_Voronoi_Headers

# Import Voronoi analysis definitions
from MD_Voronoi_Defs import *

"""Global variables"""
# Keep track of lammpstrj file
lammpstrj_file = ''
# Used for importing lammpstrj files
node = 0
# Used for storing the number of particles
numOfParticles = 0
# Used for writing files
workDir = '~/'
# Voronoi types per particle [FRAME, [TYPES]]
voro_types = [-1,[]]
# Voronoi indices per particle (Reference to see changes) [FRAME, [TYPES]]
voro_types_ref = [-1,[]]
# Particle ID
ids = []
# Particle type
types = []
# Particle position
positions = []
# Variables for histogram analysis
unique = []
counts = []
# Voronoi histogram
voroHistogram = []
#
arrays = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]

def set_file(f):
	"""Sets the file for analysis"""
	global lammpstrj_file
	lammpstrj_file = f

def type_deduction(cell):
	"""Returns the type of Voronoi cell.
	
	Keyword arguments:
	cell -- Voronoi cell to analyse
	"""
	for i in range(len(arrVoroIndices)):
		if numpy.array_equal(arrVoroIndices[i],cell):
			return i+1			
	return 0

# For each particle, checks what voronoi cell type it belongs to.
# Returns an array of cell type for all particles, and an array of cell type for Cu-centered cells
# 
# a: Voronoi indices for each particle in frame
# index: Particles IDs
# ar: Voronoi cell reference (VoroIndices_Int32)
# typ: Particles types
#
# TODO: Delete argument and use global variable << VoroIndices_Int32 >>
def id_deduction(a,index,ar,typ):
	global numOfParticles
	# This way, we get uniqueness by row, not by each element
	ca = numpy.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
	# Get references of ordering with IDs. (This doesn'y modify the array, but returns an array with the ordering positions)
	indexes = numpy.argsort(index)[::1]
	# Cu particles types
	temp = []
	# All particles types
	temp2 = []
	for elemento in range(numOfParticles):
		# Type of each particle, respecting order with increasing ID
		currentType = type_deduction(ca[indexes[elemento]],ar)
		temp2.append(currentType)
		# 1 stands for a Cu particle
		# TODO: DEFINE COPPER 1
		# type == 0 <=> not Cu
		if typ[indexes[elemento]]==1:
			temp.append(currentType)
		else:
			temp.append(0)
	return (temp2,temp)

def first_line(lines):
	for i, line in enumerate(lines):
		if line[0] != '#':
			return i
			
	return len(lines)

def voro_histogram(frame, force=False, writeFile=True):	
	"""Generates a histogram of Voronoi types for the specified frame.
	
	Keyword arguments:
	frame -- Frame number to analyse, starting from 0
	force -- If last frame was not <<frame>>, re-voro_dump() frame.
	writeFile -- Write histogram to disk
	
	TODO: Not for frame, but for given array of voro_types,
	they will contain frame number now!
	"""
	global voroHistogram
	global unique
	global voro_types
	global counts
	
	if (voro_types[0] != frame) or (voro_types[0] == -1):
		print 'You have to first voro_dump() the requested frame'
		print 'Last analysed frame: ',voro_types[0]
		print 'Requested frame: ',frame
		return
	
	# Returns the sorted unique elements of an array.
	unique = numpy.unique(voro_types[1])
	# TODO: Delete this and add return_counts option in numpy.unique (numpy 1.9.0)
	counts = numpy.bincount(voro_types[1])
	# Convention: Initial frame = 0	
	if(len(voroHistogram) <= frame):
		for i in range(len(voroHistogram),frame+1):
			voroHistogram.insert(i,[])

	for i in range(len(arrVoroIndices)):
		if (unique[i] == i):
			voroHistogram[frame].insert(i,[i,counts[i],(counts[i]/float(numOfParticles))*100])
		else:
			voroHistogram[frame].insert(i,[i,0,0])
	
	if(writeFile):	
		histFile = os.path.realpath(workDir+'/hist_MD_Voronoi')	
		if os.path.isfile(histFile) is not True:
			hFile = open(histFile, 'w')
			hFile.write(MD_Voronoi_Headers.voro_hist_header())
			hFile.close()
		
		hFile = open(histFile, 'r')
		lines = hFile.readlines()
		hFile.close()
		first = first_line(lines)
		if(len(lines) < first+frame+1):
			for i in range(len(lines),first+frame+1):
				lines.insert(i,' \n')
		lines[first+frame] = '%u ' % (frame)
		for i in range(len(arrVoroIndices)):
			lines[first+frame] += '%u ' % (voroHistogram[frame][i][1])
		lines[first+frame] += '\n'
		
		hFile = open(histFile,'w')
		for line in lines:
			hFile.write(line)
		hFile.close()
		
	return

def calculo_ovito(frame):
	ovito.dataset.anim.current_frame = frame

	global node
	node.compute()

	#voro_indices = node.output["Voronoi Index"].array
	#indices = node.output["Particle Identifier"].array
	#type_at = node.output["Particle Type"].array
	#unique_indices, counts = row_histogram(voro_indices,VoroIndices_Int32)
	#tempI,tempCUI = id_deduction(voro_indices, indices, VoroIndices_Int32, type_at)

	return #(tempI,tempCUI,unique_indices,counts,len(voro_indices))

def escribe_archivo(fdR, fdW, arrs):
	"""Cuenta cuantos cambios ocurrieron del tipo i al tipo j respecto de otro frame"""
	j=0.005
	for lines in fdR:
		fdW.write('%.3f\t' % (j))
		for i in arrs:
			fdW.write('%i\t' % (lines.count(" "+str(i)+" ")))
		fdW.write('\n')
		j = j + 0.005
		
def voro_dump(frame, writeFile=False):
	"""Analyse per-particle Voronoi type and write dump file
	
	Keyword arguments:
	frame -- Frame to analyse, starting from 0
	writeFile -- Write to disk dump info
	"""
	global lammpstrj_file, workDir, node, numOfParticles
	global ids, types, positions
	
	lammpstrj_file = os.path.realpath(lammpstrj_file)	
	if os.path.isfile(lammpstrj_file) is not True:
		print "The file",lammpstrj_file,"does not exists..."
		sys.exit()
	
	# Working directory. Not the one where we executed MD_Voronoi,
	# the path where we read lammpstrj_file
	# TODO: workDir initliaized in module, not in this function
	workDir = os.path.dirname(os.path.realpath(lammpstrj_file))
		
	# Load a simulation snapshot of a Cu-Zr metallic glass.
	node = import_file(lammpstrj_file)
	
	# Check if requested frame exists
	if (frame >= node.source.num_frames):
		print 'Requested frame for analysis doesn\'t exists in node. (Check your lammpstrj files)'
		return [-1,[]]
		
	# Set up the Voronoi analysis modifier.
	node.modifiers.append(VoronoiAnalysisModifier(
		compute_indices = True,
		use_radii = False,
		edge_count = 6,
		edge_threshold = 0,
		use_cutoff = True,
		cutoff = 6.0
	))
	
	# Number of particles
	numOfParticles = node.source.data['Position'].size
	
	ovito.dataset.anim.current_frame = frame
	node.compute()

	voro_indices = node.output["Voronoi Index"].array
	ids = node.output["Particle Identifier"].array
	types = node.output["Particle Type"].array
	positions = node.output["Position"].array
	
	voro_types_local = [frame,[]]
	
	for i in range(numOfParticles):
		voro_types_local[1].insert(i,type_deduction(voro_indices[i]))
	
	if(writeFile):
		dump_file = open(workDir+'/dump_MD_Voronoi_Frame_'+frame, 'w')
		dump_file.write(MD_Voronoi_Headers.voro_dump_header(lammpstrj_file))
		dump_file.write('id particle_type x y z voro_type\n')	
		for i in range(numOfParticles):
			dump_file.write('%u %i %f %f %f %u\n' % (	ids[i], 			\
														types[i],			\
														positions[i][0], 	\
														positions[i][1],  	\
														positions[i][2], 	\
														voro_types_local[1][i]))
												
		dump_file.close()
		
	return voro_types_local

def voro_change(frameRef, frame, justCu=False, writeFile=True):
	"""Analyses changes in the voronoi cell type for each particle and
	saves an histogram of these changes.
	
	Keyword arguments:
	frameRef -- The reference frame
	frame -- Frame to analyse changes
	justCu -- Only for Cu-centered cells
	writeFile -- Write results to file.
	"""
	global voro_types, voro_types_ref, lammpstrj_file, types, numOfParticles
	
	if(frameRef == frame):
		# Same frame, no change!
		return
	# Prepare data
	elif(frameRef == voro_types_ref[0]):
		# Ref already there
		if(frame == voro_types[0]):
			# This frame already there
			pass
	elif(frameRef == voro_types[0]):
		# Ref already there ...
		voro_types_ref = voro_types
		# ... but we have to voro_dump this frame
		voro_types = voro_dump(frame)
	else:
		# Ref not there..
		voro_types_ref = voro_dump(frameRef)
		# Check this frame
		if(frame == voro_types[0]):
			# This frame already there
			pass
		else:
			voro_types = voro_dump(frame)
	
	# Go
	changes = []
	n_changes = []

	for i in range(numOfParticles):
		if((justCu is True and types[i] == 1) or (justCu is not True)):
			changes.insert(i,voro_types_ref[1][i]*10+voro_types[1][i])
	
	unique = numpy.unique(changes)
	# TODO: Delete this and add return_counts option in numpy.unique (numpy 1.9.0)
	counts = numpy.bincount(changes)

	for i in range(100):
		if len(unique) > i and unique[i] == i:
			n_changes.insert(i,[i,counts[i],(counts[i]/float(numOfParticles))*100])
		else:
			n_changes.insert(i,[i,0,0])
	
	if(writeFile):
		print 	'Warning: Voronoi type-change histogram file is only useful \
				when frame of reference remains constant. Keep in mind that running \
				full analysis will write down this file with respect of frame 0'
		which = 'all'
		if justCu is True:
			which = 'Cu'
		changesFile = os.path.realpath(workDir+'/changesHist_'+which+'_MD_Voronoi')	
		if os.path.isfile(changesFile) is not True:
			cFile = open(changesFile, 'w')
			cFile.write(MD_Voronoi_Headers.voro_changes_header())
			cFile.close()
		
		cFile = open(changesFile, 'r')
		lines = cFile.readlines()
		cFile.close()
		first = first_line(lines)
		if(len(lines) < first+frame+1):
			for i in range(len(lines),first+frame+1):
				lines.insert(i,' \n')
		lines[first+frame] = '%u ' % (frame)
		for i in range(len(n_changes)):
			lines[first+frame] += '%u ' % (n_changes[i][1])
		lines[first+frame] += '\n'
		
		cFile = open(changesFile,'w')
		for line in lines:
			cFile.write(line)
		cFile.close()
	
	return n_changes

# TODO: Just dump voronoi cell info per particle with position and ID
#       for later analysis.
#       The rest of processing like changes with respect of last frame,
#       should be in separate functions for modularity
def dump(lammps_file):

	global node

	# Load a simulation snapshot of a Cu-Zr metallic glass.
	node = import_file(lammpstrj_file)

	# Set up the Voronoi analysis modifier.
	node.modifiers.append(VoronoiAnalysisModifier(
		compute_indices = True,
		use_radii = False,
		edge_count = 6,
		edge_threshold = 0
	))

	# Number of particles
	global numOfParticles
	numOfParticles = node.source.data['Position'].size

	# Working directory. Not the one where we executed MD_Voronoi,
	# the path where we read lammpstrj_file
	global workDir
	workDir = os.path.dirname(os.path.realpath(lammpstrj_file))

	f_Ini = open(workDir+'/Cambio_inicial', 'w')
	f_Ant = open(workDir+'/Cambio_anterior', 'w')
	f_cuI = open(workDir+'/Cambio_CUI','w')
	f_cuA = open(workDir+'/Cambio_CUA','w')
	f = []
	f.insert(0,open(workDir+'/Tipo1', 'w'))
	f.insert(1,open(workDir+'/Tipo2', 'w'))
	f.insert(2,open(workDir+'/Tipo3', 'w'))
	f.insert(3,open(workDir+'/Tipo4', 'w'))
	f.insert(4,open(workDir+'/Tipo5', 'w'))
	f.insert(5,open(workDir+'/Tipo6', 'w'))
	f.insert(6,open(workDir+'/Tipo7', 'w'))
	f.insert(7,open(workDir+'/Tipo8', 'w'))
	f.insert(8,open(workDir+'/Tipo9', 'w'))

	tempInicial,tempCUInicial,unique_indices,counts,size_voro = calculo_ovito(0)

	temp2 = tempInicial
	tempCU2 = tempCUInicial

	for i in range(len(unique_indices)):
		f[i].write('%.3f\t%s\tTipo%i\t%i\t(%.1f %%)\n' % (0.000, tuple(unique_indices[i]), i+1, counts[i], 100.0*float(counts[i])/size_voro))

	for frame in range(1, ovito.dataset.anim.last_frame + 1):

		temp,tempCU,unique_indices,counts,size_voro = calculo_ovito(frame)

		for i in range(len(temp)):
			f_Ini.write('%i ' % (tempInicial[i]*10+temp[i]))
			f_Ant.write('%i ' % (temp2[i]*10+temp[i]))
			f_cuI.write('%i ' % (tempCUInicial[i]*10+tempCU[i]))
			f_cuA.write('%i ' % (tempCU2[i]*10+tempCU[i]))
		f_Ini.write('\n')
		f_Ant.write('\n')
		f_cuI.write('\n')
		f_cuA.write('\n')

		temp2 = temp
		tempCU2 = tempCU

		for i in range(len(unique_indices)):
			f[i].write('%.3f\t%s\tTipo%i\t%i\t(%.1f %%)\n' % (frame*0.005, tuple(unique_indices[i]), i+1, counts[i], 100.0*float(counts[i])/size_voro))


	for i in range (len(f)):
		f[i].close()
		
	f_Ini.close()
	f_Ant.close()
	f_cuI.close()
	f_cuA.close()

	f_Ini = open(workDir+'/Cambio_inicial', 'r')
	f_Ant = open(workDir+'/Cambio_anterior', 'r')
	f_cuI = open(workDir+'/Cambio_CUI','r')
	f_cuA = open(workDir+'/Cambio_CUA','r')

	fIni = open(workDir+'/Cambios_FrameInicial','w')
	fAnt = open(workDir+'/Cambios_FrameAnterior', 'w')
	fcuI = open(workDir+'/Cambios_CUFrameInicial','w')
	fcuA = open(workDir+'/Cambios_CUFrameAnterior','w')
		
	escribe_archivo(f_Ini,fIni,arrays)
	escribe_archivo(f_Ant,fAnt,arrays)
	escribe_archivo(f_cuI,fcuI,arrays)
	escribe_archivo(f_cuA,fcuA,arrays)

	f_Ini.close()
	f_Ant.close()
	f_cuI.close()
	f_cuA.close()

	fIni.close()
	fAnt.close()
	fcuI.close()
	fcuA.close()

if __name__ == "__main__":
	import sys
	lammpstrj_file = os.path.realpath(sys.argv[1])

	if os.path.isfile(lammpstrj_file	) is not True:
		print "The file",lammpstrj_file,"does not exists..."
		sys.exit()
    
	print "Running module, just dump functionality. For other outputs, run the different functions from your own script"
	voro_dump(0)
