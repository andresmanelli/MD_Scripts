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

# Import NumPy module.
import numpy

# Used globally later for importing lammpstrj files
node = 0
# Used globally later for storing the number of particles
numofParticles = 0

# Defines voronoi cells
# [0, 0, 0, 0, 12, 0] => 12 faces with 5 edges each.
# [0, 0, 0, 3, 6, 3]  => 3 faces with 4 edges, 6 faces with 5 edges, 3 faces with 6 edges
# Note: The first two indices are always zero because there are no faces with less than three edges.
arrVoroIndices = numpy.array([[0,0,0,0,12,0],[0,0,0,2,8,2],[0,0,0,2,8,1],[0,0,0,3,6,3],[0,0,0,3,6,4],[0,0,0,1,10,2],[0,0,0,2,8,3],[0,0,0,4,4,3],[0,0,1,0,9,3]])
# Creates a 'view' of array <<arr>> that points to the same data, but represent it as int32
# TODO: Verify, is this architecture dependant?
VoroIndices_Int32 = numpy.ascontiguousarray(arrVoroIndices).view([('', numpy.int32)] * 6)
#
arrays = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]

# Returns the type of Voronoi cell.
# a:  This cell
# ar: Voronoi cell reference (VoroIndices_Int32)
#
# TODO: Delete argument and use global variable << VoroIndices_Int32 >>
def type_deduction(a,ar):
	for i in range(len(ar)):
		temp = ar[i][0]
		for j in range(7):
			if j==6:
				return i+1
			if temp[j]!=a[0][j]:
				break
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
	global numofParticles
	# This way, we get uniqueness by row, not by each element
	ca = numpy.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
	# Get references of ordering with IDs. (This doesn'y modify the array, but returns an array with the ordering positions)
	indexes = numpy.argsort(index)[::1]
	# Cu particles types
	temp = []
	# All particles types
	temp2 = []
	for elemento in range (numofParticles):
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
# a: Voronoi indices for each particle in frame
# ar: Voronoi cell reference (VoroIndices_Int32)
#
# TODO: Delete argument and use global variable << VoroIndices_Int32 >>
def row_histogram(a,ar):
	# This way, we get uniqueness by row, not by each element
	ca = numpy.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
	# Returns the sorted unique elements of an array.
	unique, indices, inverse = numpy.unique(ca, return_index=True, return_inverse=True)
	# TODO: Delete this and add return_counts option in numpy.unique (numpy 1.9.0)
	counts = numpy.bincount(inverse)
	# vector1 stores the voronoi indices of those cells found in this frame
	# and that matches one of those defined in arrVoroIndices
	vector1 = []
	# vector2 stores the cell count for the corresponding cell type (same position
	# in vector1)
	vector2 = []
	# For each voronoi reference type
	for i in range(len(ar)):
		# TODO: Is it necessary to convert arrVoroIndices into VoroIndices_Int32?
		#		Yes, because ovito output is int32
		# This convertion adds a needed depth in the selection ([i][j][k])
		# ar[0][0] = (0, 0, 0, 0, 12, 0)
		temp = ar[i][0]
		# For each unique voronoi type in frame
		for w in range (len(unique)):
			if numpy.array_equal(temp,unique[w]):
				vector1.insert(i,a[indices[w]])
				vector2.insert(i,counts[w])
				break
	return (vector1, vector2)

def calculo_ovito(frame):
	ovito.dataset.anim.current_frame = frame

	global node
	node.compute()

	voro_indices = node.output["Voronoi Index"].array
	indices = node.output["Particle Identifier"].array
	type_at = node.output["Particle Type"].array
	# TODO: Verify
	particlePositions = node.source.data["Position"].array
	unique_indices, counts = row_histogram(voro_indices,VoroIndices_Int32)
	tempI,tempCUI = id_deduction(voro_indices, indices, VoroIndices_Int32, type_at)

	return (tempI,tempCUI,unique_indices,counts,len(voro_indices))

def escribe_archivo(fdR, fdW, arrs):
	j=0.005
	for lines in fdR:
		fdW.write('%.3f\t' % (j))
		for i in arrs:
			fdW.write('%i\t' % (lines.count(" "+str(i)+" ")))
		fdW.write('\n')
		j = j + 0.005

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
	global numofParticles
	numofParticles = node.source.data['Position'].size

	f_Ini = open('Cambio_inicial', 'w')
	f_Ant = open('Cambio_anterior', 'w')
	f_cuI = open('Cambio_CUI','w')
	f_cuA = open('Cambio_CUA','w')
	f = []
	f.insert(0,open('Tipo1', 'w'))
	f.insert(1,open('Tipo2', 'w'))
	f.insert(2,open('Tipo3', 'w'))
	f.insert(3,open('Tipo4', 'w'))
	f.insert(4,open('Tipo5', 'w'))
	f.insert(5,open('Tipo6', 'w'))
	f.insert(6,open('Tipo7', 'w'))
	f.insert(7,open('Tipo8', 'w'))
	f.insert(8,open('Tipo9', 'w'))

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

	f_Ini = open('Cambio_inicial', 'r')
	f_Ant = open('Cambio_anterior', 'r')
	f_cuI = open('Cambio_CUI','r')
	f_cuA = open('Cambio_CUA','r')

	fIni = open('Cambios_FrameInicial','w')
	fAnt = open('Cambios_FrameAnterior', 'w')
	fcuI = open('Cambios_CUFrameInicial','w')
	fcuA = open('Cambios_CUFrameAnterior','w')
		
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
	import os.path
	lammpstrj_file = sys.argv[1]

	if os.path.isfile(lammpstrj_file) is not True:
		print "The file",lammpstrj_file,"does not exists..."
		sys.exit()
    
	print "Running module, just dump functionality. For other outputs, run the different functions from your own script"
	dump(lammpstrj_file)
