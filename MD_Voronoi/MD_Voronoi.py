#  MD_Voronoy.py
#  coding: UTF-8
#  
#  Authors:  Ardiani, Franco
#            Manelli, Andrés
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

# Global variables
# Used for importing lammpstrj files
node = 0
# Keep track of last analysed frame
last_frame = -1
# Used for storing the number of particles
numOfParticles = 0
# Used for writing files
workDir = '~/'
# Voronoi indices per particle
voro_indices = []
# Particle ID
ids = []
# Particle type
types = []
# Particle position
positions = []
# Variables for histogram analysis
unique = []
counts = []
voro_types = []
# Voronoi histogram
voroHistogram = []
#
arrays = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]

# Returns the type of Voronoi cell.
# cell:  Voronoi cell to analyse
def type_deduction(cell):
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

# 
def row_histogram(frame):	
	global voroHistogram
	global unique
	global voro_types
	global counts
	
	if (last_frame != frame) or (last_frame == -1):
		print 'You have to first voro_dump() the requested frame'
		print 'Last analysed frame: ',last_frame
		print 'Requested frame: ',frame
		return
	
	# Returns the sorted unique elements of an array.
	unique = numpy.unique(voro_types)
	# TODO: Delete this and add return_counts option in numpy.unique (numpy 1.9.0)
	counts = numpy.bincount(voro_types)
	#TODO: Insert empty rows in histogram to keep the order. This will
	# allow to analyse any specified frame
	print 'Frame: ',frame,' len hist: ',len(voroHistogram)
	# Convention: Initial frame = 0
	if(len(voroHistogram) <= frame):
		for i in range(len(voroHistogram),frame-len(voroHistogram)+1):
			print i
			voroHistogram.insert(i,[])
	
	print 'After, Frame: ',frame,' len hist: ',len(voroHistogram)
	
	for i in range(len(arrVoroIndices)):
		if (unique[i] == i):
			voroHistogram[frame].insert(i,[i,counts[i],(counts[i]/float(numOfParticles))*100])
		else:
			voroHistogram[frame].insert(i,[i,0,0])
						
	return

def calculo_ovito(frame):
	ovito.dataset.anim.current_frame = frame

	global node
	node.compute()

	#voro_indices = node.output["Voronoi Index"].array
	#indices = node.output["Particle Identifier"].array
	#type_at = node.output["Particle Type"].array
	particlePositions = node.output["Position"].array
	#unique_indices, counts = row_histogram(voro_indices,VoroIndices_Int32)
	#tempI,tempCUI = id_deduction(voro_indices, indices, VoroIndices_Int32, type_at)

	return #(tempI,tempCUI,unique_indices,counts,len(voro_indices))

def escribe_archivo(fdR, fdW, arrs):
	j=0.005
	for lines in fdR:
		fdW.write('%.3f\t' % (j))
		for i in arrs:
			fdW.write('%i\t' % (lines.count(" "+str(i)+" ")))
		fdW.write('\n')
		j = j + 0.005

def voro_dump(lammpstrj_file, frame):
	
	lammpstrj_file = os.path.realpath(lammpstrj_file)	
	if os.path.isfile(lammpstrj_file) is not True:
		print "The file",lammpstrj_file,"does not exists..."
		sys.exit()
	
	global node		
	# Load a simulation snapshot of a Cu-Zr metallic glass.
	node = import_file(lammpstrj_file)
	
	# Check if requested frame exists
	if (frame >= node.source.num_frames):
		print 'Requested frame for analysis doesn\'t exists in node. (Check your lammpstrj files)'
		return
		
	global last_frame
	last_frame = frame
		
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
	global numOfParticles
	numOfParticles = node.source.data['Position'].size
	
	ovito.dataset.anim.current_frame = frame
	node.compute()

	voro_indices = node.output["Voronoi Index"].array
	ids = node.output["Particle Identifier"].array
	types = node.output["Particle Type"].array
	positions = node.output["Position"].array
	
	# Working directory. Not the one where we executed MD_Voronoi,
	# the path where we read lammpstrj_file
	# TODO: workDir initliaized in module, not in this function
	global workDir
	workDir = os.path.dirname(os.path.realpath(lammpstrj_file))
	
	global voro_types
	voro_types = []
	
	dump_file = open(workDir+'/dump_MD_Voronoi', 'w')
	dump_file.write(MD_Voronoi_Headers.voro_dump_header(lammpstrj_file))
	dump_file.write('id particle_type x y z voro_type\n')	
	for i in range(numOfParticles):
		voro_types.insert(i,type_deduction(voro_indices[i]))
		dump_file.write('%u %i %f %f %f %u\n' % (	ids[i], 			\
													types[i],			\
													positions[i][0], 	\
													positions[i][1],  	\
													positions[i][2], 	\
													voro_types[i]))
											
	dump_file.close()
	return

# TODO: Just dump voronoi cell info per particle with position and ID
#       for later analysis.
#       The rest of processing like changes with respect of last frame,
#       should be in separate functions for modularity
def dump(lammpstrj_file):

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
	voro_dump(lammpstrj_file)
