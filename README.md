# MD Data Analysis Tool

<img width="630" alt="screen shot 2018-06-20 at 10 09 20 pm" src="https://user-images.githubusercontent.com/29931166/41695874-a9d3c352-74d6-11e8-9a7d-aa93f91e55e1.png">

## install

$ git clone https://github.com/like201707/MD_data_analysis

$ cd MD_data_analysis

$ python setup.py install

## Running

usage: mdanalysis [-h] [-t TASK] [-dt DTMAX] [-c CUTOFF] [-tp TYPES]
                  [input1] [input2]

A tool to do MD data analysis

positional arguments:
  input1                XYZ trajectory of system, or Lammps trajectory of
                        system
  input2                XYZ trajectory of system.(it is optional, take it only
                        necessary)

optional arguments:
  -h, --help            show this help message and exit
  
  -t TASK, --task TASK  type of task: XYZReader, diffusion, ionpair, ionic,
                        wXYZ (default: XYZReader)
			
  -dt DTMAX, --dtmax DTMAX
                        set maximum delta t for task diffusion (default: 1000)
			
  -c CUTOFF, --cutoff CUTOFF
                        set the cutoff value for task ionpair (default: 4.)
			
  -tp TYPES, --types TYPES
                        set which types need to output as XYZ file for task
                        wXYZ (default: "13")

## One Example Run

		_|      _|  _|_|_|
		_|_|  _|_|  _|    _|
		_|  _|  _|  _|    _|
		_|      _|  _|    _|
		_|      _|  _|_|_|	          ¿ⓧ_ⓧﮌ ¿ⓧ_ⓧﮌ ¿ⓧ_ⓧﮌ

----------------------------------------------------------------------------

          -------- A Tool to do MD data analysis--------
Processing ◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉ 100%

 System information:

        System trajectory file name:	nvt-efield.lammpstrj
        Number of atoms in the box:		6102
        Total number of Frames: 		30001
        First Frame Box size information: 	42.850510 42.850510 42.850510


Task information:

        Task name: output traget types atom XYZfile
Processing ◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉◉ 100%

----------------------------------------------------------------------------

Total elapsed time = 1767.49 seconds.

	_|_|_|      _|_|    _|      _|  _|_|_|_|
	_|    _|  _|    _|  _|_|    _|  _|
	_|    _|  _|    _|  _|  _|  _|  _|_|_|
	_|    _|  _|    _|  _|    _|_|  _|
	_|_|_|      _|_|    _|      _|  _|_|_|_|       (◠﹏◠)  (◠﹏◠)  (◠﹏◠)



