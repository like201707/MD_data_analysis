#!/usr/bin/python
import numpy as np
from progress.bar import FillingCirclesBar

class Distance(object):
	"""
	Get the distance between atom1 and atom2
	atom 1 is the central atom
	"""
	def __init__(self,traj1, traj2):
		
		self.Frames = traj1.nFrames 
		self.coords1 = traj1.atomC
		self.coords2 = traj2.atomC
		self.n_atoms1 = traj1.nAtoms
		self.n_atoms2 = traj2.nAtoms
		self.bzx, self.bzy, self.bzz = traj1.boxSize
		
	def distance_one(self, coords1, coords2):
		"""
		Get the distance between atom1 and atom2 in one frame
		"""
		dist = np.zeros([self.n_atoms1, self.n_atoms2], dtype=np.float)
		
		# periodic boundary conditions
		for i in range(self.n_atoms1):
			atom1 = coords1[i]
			for j in range(self.n_atoms2):
				atom2 = coords2[j]
				dr = atom1 - atom2
				if dr[0] > self.bzx * 0.5:
					dr[0] -= self.bzx
				elif dr[0] <= -self.bzx * 0.5:
					dr[0] += self.bzx
				if dr[1] > self.bzy * 0.5:
					dr[1] -= self.bzy
				elif dr[1] <= -self.bzy * 0.5:
					dr[1] += self.bzy
				if dr[2] > self.bzz * 0.5:
					dr[2] -= self.bzz
				elif dr[2] <= -self.bzz * 0.5:
					dr[2] += self.bzz
					
				dist[i][j] = np.sqrt(dr[0]**2 + dr[1]**2 + dr[2]**2)
		
		return dist
		
	def distance_all(self):
		"""
		Get the distance between atom1 and atom1 in all frames
		"""
		bar = FillingCirclesBar('Processing', max=self.Frames)
		dist_all = np.empty([self.Frames, self.n_atoms1, self.n_atoms2], dtype=np.float)
		
		for i in range(self.Frames):
			dist_all[i] = self.distance_one(self.coords1[i], self.coords2[i])
		
			bar.next()
		bar.finish()
		
		return dist_all
			

				
				
				