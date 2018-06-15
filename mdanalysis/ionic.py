#!/usr/bin/python

import numpy as np
import pylab as pl
import yaml
import os 
from scipy.optimize import curve_fit
from progress.bar import FillingCirclesBar

class Ionic_Conductivity(object):
	
	def __init__(self, traj, inputfile):
		
		with open(inputfile, 'r') as f:
			self.input = yaml.load(f)	
			
		self.E = float(self.input['efield'])
		self.ti = self.input['t_start']
		self.tf = self.input['t_end']		
		self.traj = traj
		self.coords = traj.atomC
		self.time = traj.nFrames
		self.n_atomes = traj.nAtoms
		self.q = self.input['charge']
		base = os.path.basename(traj.filename)
		self.fprefix = os.path.splitext(base)[0]
		d = self.input['direction']
		if d == 'X' or d == 'x':
			self.d = 0
		elif d == 'Y' or d == 'y':
			self.d = 1
		else:
			self.d = 2
		# box size did not change for all frames 
		self.bz = traj.boxSize[0][self.d]
		
	def method1(self):
		
		count_ions = np.zeros((self.time, 1))
		tmp = 0
		
		bar = FillingCirclesBar('Processing', max=self.time-1)
		for t in range(1, self.time):
			counter = 0
			for n in range(self.n_atomes):
				if self.coords[t-1][n][self.d] < self.bz/2. and self.coords[t-1][n][self.d] >= (self.bz/2. - 5.):
					if self.coords[t][n][self.d] > self.bz/2. and self.coords[t][n][self.d] <= (self.bz/2. + 5.):
						counter += 1
				if self.coords[t-1][n][self.d] > self.bz/2. and self.coords[t-1][n][self.d] <= (self.bz/2. + 5.):
					if self.coords[t][n][self.d] < self.bz/2. and self.coords[t][n][self.d] >= (self.bz/2. - 5.):
						counter -= 1
				
			tmp += counter
			count_ions[t] = tmp
			
			bar.next()
		bar.finish()
		
		return count_ions
		
	def func(self, x, k):
		
		return x * k 
		
		
	def out_put(self, taskname = 'ionic'):
		"""output raw data"""
		crossIons = np.squeeze(self.method1())
		raw_data = self.fprefix + '_' + taskname + '.dat'
		raw_image = self.fprefix + '_' + taskname + '.png'
		with open(raw_data, 'w') as f:
			for i,j in enumerate(crossIons):
				f.write('{}\t{}\n'.format(i, j))
				
		# ploting
		X = np.arange(0, self.tf-self.ti, 1)
		Y = crossIons[self.ti:self.tf] - crossIons[self.ti]
		popt, pcov = curve_fit(self.func, X, Y)
		fig = pl.figure()
		pl.plot(X, Y,'b', label=self.fprefix)
		pl.plot(X, self.func(X, *popt),'g-', label='fit: k=%5.3f' % tuple(popt))
		pl.xlabel('time(ps)')
		pl.ylabel('Number of ion cross the plane')
		pl.legend()
		fig.savefig(raw_image, dpi=fig.dpi)
		
				
		# print conductivity		
		c = 6.242 * 10 ** 18
		k = popt
		I = self.q * k / (c * 10 ** -12)
		J = I / (self.bz * 10 ** -10) ** 2
		sigma = J / self.E
		print( '           ' + self.fprefix + ' ion conductivity is %f'%sigma)
		
		
		
		
		
	
		
