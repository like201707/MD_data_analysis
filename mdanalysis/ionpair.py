#!/usr/bin/python
import numpy as np
import distance
import os
from progress.bar import FillingCirclesBar

class Ionpair_distribution(distance.Distance):
	"""
	calculate the distribution of the traget atom 
	around the central atom
	"""
	
	def __init__(self, traj1, traj2, cutoff):
		
		self.cutoff = cutoff
		base1 = os.path.basename(traj1.filename)
		base2 = os.path.basename(traj2.filename)
		self.fprefix1 = os.path.splitext(base1)[0]
		self.fprefix2 = os.path.splitext(base2)[0]
		super(Ionpair_distribution, self).__init__(traj1, traj2)
		
	def counter(self):
		"""
		return a 2-d array, 1st dimention is number of frames
		2nd dimention is the number of traget atom for each central atom
		"""
		print('      -----calculatiing distance between two atoms-----')
		cter = np.empty([self.Frames, self.n_atoms1], dtype=np.int)
		dist = self.distance_all()
		print('      -------calculatiing ion pair distribution-------')
		bar = FillingCirclesBar('Processing', max=self.Frames)

		for i in range(self.Frames):
			for j in range(self.n_atoms1):
				c = 0
				for k in range(self.n_atoms2):
					if dist[i][j][k] < self.cutoff:
						c += 1
				cter[i][j] = c
				
			bar.next()
		bar.finish()
			
		return cter
		
	def distribution(self, num=5):
		"""
		calculate the distribution
		"""
		
		cter = self.counter()
		prb = dict()
		
		for key in range(num):
			c = 0
			for i in range(self.Frames):
				for j in range(self.n_atoms1):
					if cter[i][j] == key:
						c += 1
			
			value = c / float(self.n_atoms1 * self.Frames)
			prb[str(key)] = value
		
		return prb
		
	def out_put(self, taskname = 'ionpair'):
		"""output raw data"""
		prb = self.distribution()
		raw_data = self.fprefix1 + '_' + self.fprefix2 + '_' + taskname + '.dat'
		with open(raw_data, 'w') as f:
			for key, value in prb.iteritems():
				f.write('{}\t{}\n'.format(key, value))
		
	
				
		