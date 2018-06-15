#!/usr/bin/python

import argparse
import time
import XYZReader
import lmpsr 
import diffusion
import ionpair
import ionic
import util
import wXYZ

def get_parser():
	
	parser = argparse.ArgumentParser(description='A tool to do MD data analysis')
	parser.add_argument('input1', type=str, nargs='?',help='XYZ trajectory of system, or Lammps trajectory of system')
	parser.add_argument('input2', type=str, nargs='?',help='XYZ trajectory of system.(it is optional, take it only necessary) ')
	parser.add_argument('-t','--task', default='XYZreader', type=str, help='type of task: XYZReader, diffusion, ionpair, ionic, wXYZ (default: XYZReader)')
	parser.add_argument('-dt','--dtmax', default=1000, type=int, help='set maximum delta t (default: 1000)')
	parser.add_argument('-c','--cutoff', default=4., type=float, help='set the cutoff value (default: 4.)')
	parser.add_argument('-tp','--types', default='13', type=str, help='set which types need to output as XYZ file (default: "13")')

	return parser
	
def command_line_runner():
	parser = get_parser()
	args = vars(parser.parse_args())
	
	if not args['input1']:
		parser.print_help()
		return
		
	if not args['input2']:
		pass
		
	if args['task']:
		t_start = time.clock()
		util.output_welcome()
		try:
			reader = XYZReader.XYZReader(args['input1'])
			util.output_system_info(reader.filename, reader.nAtoms, reader.nFrames, reader.boxSize)
		except ValueError:
			reader = lmpsr.LMPS_Reader(args['input1'])
			util.output_system_info(reader.filename, reader.nAtoms, reader.nFrames, reader.boxSize)
			
		try:
			reader2 = XYZReader.XYZReader(args['input2'])
			util.output_system_info(reader2.filename, reader2.nAtoms, reader2.nFrames, reader2.boxSize)
		except TypeError:
			pass
			
		if 'diffusion' in args['task']:
			util.output_task('diffusion')
			task = diffusion.Diffusion_coefficent(reader, args['dt_max'])
			task.out_put()
			
		if 'ionpair' in args['task']:
			util.output_task('ionpair')
			task = ionpair.Ionpair_distribution(reader, reader2, args['cutoff'])
			task.out_put()
			
		if 'ionic' in args['task']:
			util.output_task('ionic')
			task = ionic.Ionic_Conductivity(reader, 'ionic.in')
			task.out_put()
			
		if 'wXYZ' in args['task']:
			util.output_task('wXYZ')
			task = wXYZ.WXYZ(reader, list(map(int, args['types'].split(' '))))
			task.out_put()
			
	t_end = time.clock()
	util.output_end(t_start, t_end)
			
		
if __name__ == '__main__':
	command_line_runner()
