#!/usr/bin/python

import numpy as np 
from progress.bar import FillingCirclesBar 

class WXYZ(object):
	
	def __init__(self, traj, atom_types):
		
		self.types = atom_types
		self.atomT = traj.atomT
		self.atomC = traj.atomC
		self.nFrames = traj.nFrames
		self.nAtoms = traj.nAtoms 
		self.boxSize = traj.boxSize

	def count_atoms(self):
		"""count the number of atoms for given types"""
		type_dic = {}
		n_atoms = 0
		for each_type in self.types:
			n = 0
			for j in range(self.nAtoms):
				if self.atomT[0][j] == each_type:
					n += 1
			n_atoms += n
			type_dic[str(each_type)] = n
			
		return (type_dic, n_atoms)
		
	def get_coords_and_v(self):
		"""get the coordinates and velocites for given types"""
		n = self.count_atoms()[1]
		coords_dic = {}
#        v = np.empty([self.nFrames, n, 3], dtype=np.float)
		for each_type in self.types:
			coords = np.empty([self.nFrames, n, 3], dtype=np.float)
			for i in range(self.nFrames):
				j = 0
				for k in range(self.nAtoms):
					if self.atomT[i][k] == each_type:
						coords[i][j] = self.atomC[i][k]
#                   v[i][j] = self.atomV[i][k]
						j += 1
			coords_dic[str(each_type)] = coords
			
#        return (coords, v)
		return coords_dic
		
	def out_put(self):
		"""output xyz fild, get from lammpstraj"""
		coords_dic = self.get_coords_and_v()
		raw_data = "Type" + ' '.join(str(e) for e in self.types)+ ".xyz"
		n = self.count_atoms()[1]
		with open(raw_data, 'w') as f:
			bar = FillingCirclesBar('Processing', max=self.nFrames)
			for i in range(self.nFrames):
				f.write(str(n) + "\n")
				f.write('{0} {1} {2}\n'.format(self.boxSize[i][0], self.boxSize[i][1], self.boxSize[i][2]))
				for each_type in self.types:
					for j in range(self.count_atoms()[0][str(each_type)]):
						f.write('{0} {1} {2} {3}\n'.format(each_type, coords_dic[str(each_type)][i][j][0],
								coords_dic[str(each_type)][i][j][1], coords_dic[str(each_type)][i][j][2]))
				
				bar.next()
			bar.finish()
								
				
					
