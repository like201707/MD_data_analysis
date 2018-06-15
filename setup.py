#!/usr/bin/python

from setuptools import setup, find_packages
import mdanalysis
import os

setup(
	name='mdanalysis',
	version='0.1',
	description='A tool use to do molecular dynamic data analysis',
	url='http://github.com/like201707/mdanalysis',
	author='Ke Li',
	author_email='like201125@gmail.com',
	license='MIT',
	packages=find_packages(),
	entry_points={
		'console_scripts': [
			'mdanalysis = mdanalysis.mdanalysis:command_line_runner',
		]
	},
	install_requires=[
		'numpy',
		'progress',
		'scipy',
		'matplotlib',
		'yaml',
	],
)