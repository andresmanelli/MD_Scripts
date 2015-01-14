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

from ovito.io import *   # Import OVITO modules.
from ovito.modifiers import *

import sys
import numpy
import os.path

# Import Voronoi files headers
import MD_Voronoi_Headers

# Import Voronoi analysis definitions
from MD_Voronoi_Defs import *

"""Global variables"""
group = 'all' # 'all' or array of particle ids
lammpstrj_file = '' # Keep track of lammpstrj file
node = 0 # Used for importing lammpstrj files
numOfParticles = 0 # Used for storing the number of particles
workDir = '~/' # Used for writing files
voro_indices = [-1,[]] # Voronoi indices per particle of last analysed frame
voro_types = [-1,[]] # Voronoi types per particle [FRAME, [TYPES]]
voro_types_ref = [-1,[]] # Voronoi indices per particle (Reference to see changes) [FRAME, [TYPES]]
ids = [] # Particle ID
types = [] # Particle type
positions = [-1,[]] # Particle position
counts = [] # Variable for histogram analysis
voroHistogram = [] # Voronoi histogram

def group(shape,par,frame=0):
	"""Returns the indexes of the particles contained in this region for the
	current frame, for accesing the arrays (ids, positions, etc)
	
	Keyword arguments:
	shape -- 'sphere' or [...]
	par --  [cx, cy, cz, radius] for 'sphere'
			[xlo, xhi, ylo, yhi, zlo, zhi] for 'cube'
			[x0, y0, z0, a, b, c] for 'ellipsoid'
			[...] for [...]
	frame -- frame
	"""
	global positions, ids, voro_types
	
	if(frame == positions[0]):	
		g = []
		if shape == 'sphere':
			if len(par) is not 4:
				print 'Parameters of sphere are not ok..'
				return []
			cx = par[0]
			cy = par[1]
			cz = par[2]
			r = par[3]
			for i,p in enumerate(positions[1]):
				if ((p[0] - cx)**2 + (p[1] - cy)**2 + (p[2] - cz)**2) <= r**2:
					g.insert(len(g),i)
					
			print 'Selected spherical shape of',len(g),'atoms'
			return g
		elif shape == 'cube':
			if len(par) is not 6: # xlo xhi, ylo yhi, zlo zhi
				print 'Parameters of cube are not ok..'
				return []
			xlo, xhi, ylo, yhi, zlo, zhi = par[0],par[1],par[2],par[3],par[4],par[5]
			for i,p in enumerate(positions[1]):
				if 	(p[0] >= xlo and p[0] <= xhi) and \
					(p[1] >= ylo and p[1] <= yhi) and \
					(p[2] >= zlo and p[2] <= zhi):
					g.insert(len(g),i)				
			
			print 'Selected cubic shape of',len(g),'atoms'
			
			return g
		elif shape == 'ellipsoid':
			if len(par) is not 6:
				print 'Parameters of ellipsoid are not ok..'
				return []
			x0, y0, z0, a, b, c = par[0],par[1],par[2],par[3],par[4],par[5]
			for i,p in enumerate(positions[1]):
				if 	((p[0] - x0)**2)/float(a)**2 + \
					((p[1] - y0)**2)/float(b)**2 + \
					((p[2] - z0)**2)/float(c)**2 <= 1:
					g.insert(len(g),i)	
			
			print 'Selected ellipsoidal shape of',len(g),'atoms'
			
			return g
		else:
			print 'Shape not recognized'
			return []
	else:
		voro_types = voro_dump(frame)
		return group(shape,par,frame)
		
def init(f):
	"""Sets the file for analysis"""
	# TODO: Init also the voronoi reference indices
	global lammpstrj_file, workDir, node, numOfParticles
	
	lammpstrj_file = os.path.realpath(f)	
	if os.path.isfile(lammpstrj_file) is not True:
		print "The file",lammpstrj_file,"does not exists..."
		sys.exit()
	
	# Working directory. Not the one where we executed MD_Voronoi,
	# the path where we read lammpstrj_file
	workDir = os.path.dirname(os.path.realpath(lammpstrj_file))
		
	# Load a simulation snapshot of a Cu-Zr metallic glass.
	node = import_file(lammpstrj_file)
			
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

def type_deduction(cell):
	"""Returns the type of Voronoi cell.
	
	Keyword arguments:
	cell -- Voronoi cell to analyse
	"""
	for i in range(len(arrVoroIndices)):
		if numpy.array_equal(arrVoroIndices[i],cell):
			return i+1			
	return 0

def first_line(lines):
	"""Returns the first writable line in file after header"""
	for i, line in enumerate(lines):
		if line[0] != '#':
			return i
			
	return len(lines)

def filter_group(group,prop):
	"""Returns a filtered array with the indexes given by group.
	
	Keyword arguments:
	group -- indexes of particles in group
	prop -- property to return for this group
	"""
	ret = []
	for i in group:
		ret.insert(len(ret),prop[i])
		
	return ret

def voro_freq(frame, just=False, group='all'):
	"""Returns an array with 9 more frequent voronoi cells in this frame
	
	Keyword arguments:
	frame -- Frame number to analyse, starting from 0
	"""
	global voro_indices, voro_types, types
	
	voro_filtered = []
	voro_pre_filtered = []
	if(frame == voro_indices[0]):
		if group is not 'all':
			voro_pre_filtered = filter_group(group, voro_indices[1])
			
		if just is not False:
			for i, index in enumerate(group):
				if types[index] == just:
					voro_filtered.insert(len(voro_filtered), voro_pre_filtered[i])
		else:
			voro_filtered = voro_pre_filtered
		
		voro_filtered = numpy.array(voro_filtered)
		
		# In order to find unique rows
		filtTemp = numpy.ascontiguousarray(voro_filtered).view(numpy.dtype((numpy.void, voro_filtered.dtype.itemsize * voro_filtered.shape[1])))
		_, idx, inverse = numpy.unique(filtTemp, return_index=True, return_inverse=True)
		
		unique = voro_filtered[idx]
		counts = numpy.bincount(inverse)
		indices = numpy.argsort(counts)[::-1]
		
		for i in range(10):
			print unique[indices[i]],' : ',counts[indices[i]]
			
		return [unique[indices],counts[indices]]
	else:
		voro_types = voro_dump(frame)
		return voro_freq(frame, group=group, just=just)

def voro_histogram(frame, force=False, writeFile=True, group='all', just=False):	
	"""Generates a histogram of Voronoi types for the specified frame.
	
	Keyword arguments:
	frame -- Frame number to analyse, starting from 0
	force -- If last frame was not <<frame>>, re-voro_dump() frame.
	writeFile -- Write histogram to disk
	group -- Only this particles. See functon group()
	justCu -- Only for Cu-centered cells
	
	TODO: Not for frame, but for given array of voro_types,
	they will contain frame number now!
	"""
	global voroHistogram, voro_types, counts
	
	print 'Histogram analysis for frame: ',frame
	
	if (voro_types[0] != frame) or (voro_types[0] == -1):
		print 'You have to first voro_dump() the requested frame'
		print 'Last analysed frame: ',voro_types[0]
		print 'Requested frame: ',frame
		return
	
	voro_filtered = []
	voro_pre_filtered = voro_types[1]
	suffix = 'all'
	if group is not 'all':
		suffix = 'group'
		voro_pre_filtered = filter_group(group, voro_types[1])
		
	if just is not False:
		suffix += EL_NAMES[just]
		for i, t in enumerate(types):
			if t == just:
				voro_filtered.insert(len(voro_filtered), voro_pre_filtered[i])
	else:
		voro_filtered = voro_pre_filtered
	
	counts = numpy.bincount(voro_filtered)
	# Convention: Initial frame = 0	
	if(len(voroHistogram) <= frame):
		for i in range(len(voroHistogram),frame+1):
			voroHistogram.insert(i,[])
	
	# If last types are not present, length of counts is smaller than 
	# length of arrVoroIndices
	for i in range(len(arrVoroIndices)-len(counts)+1):
		counts = numpy.insert(counts,len(counts),[0])
	
	for i in range(len(counts)):
		voroHistogram[frame].insert(i,[i,counts[i],(counts[i]/float(numOfParticles))*100])
	
	if(writeFile):
		histFile = os.path.realpath(workDir+'/hist_MD_Voronoi_'+suffix)	
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
		for i in range(len(counts)):
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
		
def voro_dump(frame, writeFile=False):
	"""Analyse per-particle Voronoi type and write dump file
	
	Keyword arguments:
	frame -- Frame to analyse, starting from 0
	writeFile -- Write to disk dump info
	"""
	global lammpstrj_file, workDir, node, numOfParticles
	global ids, types, positions, voro_indices, voro_types
	
	# Check if requested frame exists
	if (frame >= node.source.num_frames):
		print 'Requested frame for analysis doesn\'t exists in node. (Check your lammpstrj files)'
		return [-1,[]]
	
	# Check if already dump-ed
	if voro_types[0] == frame:
		return voro_types
	
	print 'Dumping frame ',frame
	
	ovito.dataset.anim.current_frame = frame
	node.compute()

	voro_indices[0] = frame
	voro_indices[1] = node.output["Voronoi Index"].array
	ids = node.output["Particle Identifier"].array
	types = node.output["Particle Type"].array
	positions[0] = frame
	positions[1] = node.output["Position"].array
	
	voro_types_local = [frame,[]]
	
	for i in range(numOfParticles):
		voro_types_local[1].insert(i,type_deduction(voro_indices[1][i]))
	
	if(writeFile):
		dump_file = open(workDir+'/dump_MD_Voronoi_Frame_'+`frame`, 'w')
		dump_file.write(MD_Voronoi_Headers.voro_dump_header(lammpstrj_file))
		dump_file.write('id particle_type x y z voro_type\n')	
		for i in range(numOfParticles):
			dump_file.write('%u %i %f %f %f %u\n' % (	ids[i], 			\
														types[i],			\
														positions[1][i][0], 	\
														positions[1][i][1],  	\
														positions[1][i][2], 	\
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
	
	print 'Change analysis for ref: ',frameRef,'and frame: ',frame
	
	if(frameRef == frame):
		# Same frame, no change!
		voro_types_ref = voro_types
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
	
	counts = numpy.bincount(changes)

	# Length of counts is different than length of arrVoroIndices^2 
	# if there aren't some kinds of changes
	# +1 because type 0 is not defined in MD_Voronoi_Defs
	for i in range((len(arrVoroIndices)+1)**2-len(counts)):
		counts = numpy.insert(counts,len(counts),[0])

	for i in range(len(counts)):
		n_changes.insert(i,[i,counts[i],(counts[i]/float(numOfParticles))*100])
	
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

def voro_combo_A():
	""" - Dump for all particles
		- Histogram for all particles
		- Changes with respect to frame 0
		- Write all files
	"""
	global node, voro_types
	for frame in range(node.source.num_frames):
		voro_types = voro_dump(frame,writeFile=True)
		voro_histogram(frame)
		voro_change(0,frame)

def voro_combo_B(group=False):
	""" - Dump for all particles (Don't write file)
		- Histogram for all particles (Write file)
		- Histogram for given group (Write file)
	"""
	if group is False:
		print 'This combo requies a group to be specified'
		return
		
	global node, voro_types
	frames = range(node.source.num_frames)
	print 'frames:',frames
	# Avoid re-dumping
	if voro_types[0] is not -1:
		frames = numpy.roll(frames,(len(frames)-voro_types[0])%4)
		print 'frames:',frames
		print 'Rolling',(len(frames)-voro_types[0])%4,'frames to frame',frames[0]
	for frame in frames:
		voro_types = voro_dump(frame)
		voro_histogram(frame) # all
		voro_histogram(frame,group=group) # group

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
