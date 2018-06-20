#!/usr/bin/python
# encoding: utf-8

import numpy as np

def output_welcome():
	"""print welcome information"""
	MD = """
			_|      _|  _|_|_|
			_|_|  _|_|  _|    _|
			_|  _|  _|  _|    _|
			_|      _|  _|    _|
			_|      _|  _|_|_|	          ¿ⓧ_ⓧﮌ ¿ⓧ_ⓧﮌ ¿ⓧ_ⓧﮌ
	"""
	print(MD)
	print('-' * 76)
	print('\n-------- A Tool to do MD data analysis--------\n')
	print('----------Loading System Information-----------\n')
	
def output_system_info(filename, n_atoms, n_frames, boxSize):
	"""print system information"""
	print('\nSystem information:\n')
	print('        System trajectory file name:	\t{:s}'.format(filename))
	print('        Number of atoms in the box:	\t{:d}'.format(n_atoms))
	print('        Total number of Frames: 	\t{:d}'.format(n_frames))
	print('        First Frame Box size information: \t{:f} {:f} {:f}\n'.format(boxSize[0][0], boxSize[0][1], boxSize[0][2]))
	
def output_task(name):
	"""print task information"""
	print('\nTask information:\n')
	if name is 'diffusion':
		print('        Task name: calculating diffusion coefficient')
		
	if name is 'ionpair':
		print('        Task name: calculating ion pair distribution')
		
	if name is 'ionic':
		print('        Task name: calculating ionic condctivity')
		
	if name is 'wXYZ':
		print('        Task name: output traget types atom XYZfile')
		
def output_end(t_start, t_end):
		"""print total running time"""
		print('\n' + '-' * 76)
		print ('\nTotal elapsed time = {:.2f} seconds.'.format(t_end - t_start))
		end = """
	_|_|_|      _|_|    _|      _|  _|_|_|_|
	_|    _|  _|    _|  _|_|    _|  _|
	_|    _|  _|    _|  _|  _|  _|  _|_|_|
	_|    _|  _|    _|  _|    _|_|  _|
	_|_|_|      _|_|    _|      _|  _|_|_|_|       (◠﹏◠)  (◠﹏◠)  (◠﹏◠)
					
			 """
		print(end) 