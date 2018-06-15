#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as pl
import os
from progress.bar import FillingCirclesBar
	
class Diffusion_coefficent(object):

	def __init__(self, traj, dt_max):
		
		self.c = traj.atomC
		self.v = traj.atomV
		self.nFrames, self.nAtoms, random = self.c.shape
		base = os.path.basename(traj.filename)
		self.fprefix = os.path.splitext(base)[0]
		self.dt_max = dt_max 

	def v_auto_correlation(self):
		"""
		 velocity veloccity autocorrelation function
		 one way to calculate the diffusion coefficent
		"""
		ctt = np.empty([self.dt_max], dtype=float)
		for t in range(self.dt_max):
			ct = 0.0
			for i in range(self.nFrames):
				if i + t < self.nFrames:
					for j in range(self.nAtoms):
						ct += np.dot(self.v[i][j], self.v[i+t][j])
			ct = ct / (self.nAtoms * (self.nFrames - t))
			ctt[t] = ct
		fig_vv = pl.figure()
		X = np.arange(0, self.dt_max, 1)
		Y = ctt
		pl.plot(X, Y, 'g-', linewidth=2)
		pl.xlabel("t(fs)", size=12)
		pl.ylabel(r"$<V(0)\cdot V(t)>$", size=12)
		pl.xticks(size=12)
		pl.yticks(size=12)
		fig_vv.savefig('vv.png', dpi=fig_vv.dpi)

	def vv_diffusion(self):
		d = 1.0 / 3.0 * np.trapz(self.ctt)
		return d

	def root_mean_square_displacement(self):
		"""
		mean square displacement
		another way to calculate the diffusionc coefficent
		"""
		bar = FillingCirclesBar('Processing', max=self.dt_max)
		rmsd = np.empty([self.dt_max], dtype=float)
		for t in range(self.dt_max):
			msd = 0.0
			mid = 0.0
			for i in range(self.nFrames):
				if i + t < self.nFrames:
					for j in range(self.nAtoms):
						mid += np.square(np.subtract(self.c[i][j],
										 self.c[i+t][j]))
						msd = np.sum(mid)
			rmsd[t] = msd / (self.nAtoms * (self.nFrames - t))
			bar.next()
		bar.finish()
		return rmsd
	
	def out_put(self, taskname = 'diffusion'):
		"""output raw data"""
		rmsd = self.root_mean_square_displacement()
		raw_data = self.fprefix + '_' + taskname + '.dat'
		with open(raw_data, 'w') as f:
			for i,j in enumerate(rmsd):
				f.write('{}\t{}\n'.format(i, j))
		fig_msd = pl.figure()
		X = np.arange(0, self.dt_max, 1)
		Y = rmsd
		pl.plot(X, Y, 'g-', linewidth=2)
		pl.xlabel("t(ps)", size=12)
		pl.ylabel("MSD")
		pl.xticks(size=12)
		pl.yticks(size=12)
		fig_msd.savefig('msd.png', dpi=fig_msd.dpi)
			



