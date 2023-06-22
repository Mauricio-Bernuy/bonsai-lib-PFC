## @namespace bonsai 
#  The bonsai module provides a set of wrapping functions to the Bonsai tree code (https://github.com/treecode/Bonsai)

"""
Bonsai.py module

Used to easily hadle .tipsy file format.

Offers nbody.txt/.tipsy conversion, galaxy manipulation, plotting and video creation.

Author: Michael M. Folkerts
E-Mail: mmfolkerts@gmail.com
Project: UC San Diego Physics 241, Winter 2014, Prof. J. Kuti
"""

import os
from subprocess import call

##\short	Path to Bonsai binary
##\details 	Default path, assuming bonsai_phys241 and Bonsai share parent folders
BONSAI_BIN = "../Bonsai/runtime/bonsai2_slowdust"

def run_tipsy(tipsy_file,snap_prefix,T,dt, dSnap, eps, bonsai_bin=None, mpi_n=0,mpi_log_file="mpiout.log"):
	"""
	Runs Bonsai with initial conditions defined by tipsy file

	@param[in]	tipsy_file		containing initial conditions
	@param[in]	snap_prefix		path prefix for snapshot files (time will be appended)
	@param[in]	T				total simulation time
	@param[in]	dt				internal time step
	@param[in]	dSnap			interval at which snapshot files are generated
	@param[in]	bonsai_bin		path to bonsai exe
	@param[in]	mpi_n			specifies the number of mpi processes (0 = mpi not used)
	@param[in]	mpi_log_file	single log file for mpi output (when mpi_n > 0)

	@returns None
	"""
	run_mode('infile',tipsy_file,snap_prefix,T,dt, dSnap, eps, bonsai_bin, mpi_n,mpi_log_file)



def run_mode(mode,nPart_or_file,snap_prefix,T,dt, dSnap, eps, bonsai_bin, mpi_n,mpi_log_file,direct):
	"""
	Run Bonsai in mode "plummer", "sphere" or "infile"

	This is an internal function, use the other interfaces instead.

	@param[in]	mode			"plummer" or "sphere" or "infile"
	@param[in]	nPart_or_file	number of particles per mpi process, or path to tipsy file for "infile" mode
	@param[in]	snap_prefix		path prefix for snapshot (tipsy) files (simulation time will be appended)
	@param[in]	T				total simulation time
	@param[in]	dt				internal time step
	@param[in] 	dSnap			interval at which snapshot files are generated
	@param[in]	bonsai_bin		path to bonsai exe
	@param[in]	mpi_n			specifies the number of mpi processes (0 = mpi not used)
	@param[in]  mpi_log_file	single log file for mpi output (when mpi_n > 0)
	@param[in]  direct			enables direct (N^2) computation

	@returns None
	@sa run_tipsy(), run_plummer(), run_sphere()
	""" 
	if mode != 'plummer' and mode != 'sphere' and mode != 'infile':
		raise Exception("Error: model '%s' is not known." % mode)

	if bonsai_bin is None:
		#use default
		bonsai_bin = BONSAI_BIN

	log = False

	if mpi_n > 0:
		#run mpi
		#mpirun -n 2 --output-filename mpiout.txt ./bonsai2_slowdust -i model3_child_compact.tipsy -T1000 --logfile logfile.txt
		if call(['mpirun','-n',str(mpi_n),
				 '--output-filename',mpi_log_file,bonsai_bin,
				 '--'+mode,str(nPart_or_file),
				 '--snapname',snap_prefix,'--snapiter',str(dSnap),
				 '-T',str(T),'-dt',str(dt),
				 '--eps',str(eps),'--direct' if direct else '',
				]):
			return "Error"
		else:
			return "Done"

	else:
		#single GPU mode
		print(f"{bonsai_bin} {'--log' if log else ''} {'--direct' if direct else ''} --{mode} {str(nPart_or_file)} --snapname {snap_prefix} --snapiter {str(dSnap)} -T {str(T)} -dt {str(dt)} --eps {str(eps)}")
		if call([bonsai_bin,'--log' if log else '','--direct' if direct else '','--'+mode,str(nPart_or_file),'--snapname',snap_prefix,'--snapiter',str(dSnap),'-T',str(T),'-dt',str(dt),'--eps',str(eps),]):
			return "Error"
		else:
			return "Done"


def run_plummer(nParticles,snap_prefix,T=2,dt=0.0625, dSnap = 0.0625, eps=0.05, bonsai_bin = None, mpi_n = 0, mpi_log_file = "mpiout.log",direct=None):
	"""
	Run a Bonsai's built in plummer model

	@param[in]	nParticles		number of particles (per mpi process)
	@param[in]	snap_prefix		path prefix for snapshot files (time will be appended)
	@param[in]	T				total simulation time
	@param[in]	dt				internal time step
	@param[in]	dSnap			interval at which snapshot files are generated
	@param[in]	bonsai_bin		path to bonsai exe
	@param[in]	mpi_n			specifies the number of mpi processes (0 = mpi not used)
	@param[in]	mpi_log_file	single log file for mpi output (when mpi_n > 0)
	@param[in] direct			enables direct (N^2) computation

	@returns	None
	"""
	run_mode("plummer",nParticles,snap_prefix,T,dt, dSnap, eps, bonsai_bin, mpi_n,mpi_log_file,direct)


def run_sphere(nParticles,snap_prefix,T=2,dt=0.0625, dSnap = 0.0625, eps=0.05, bonsai_bin = None, mpi_n = 0, mpi_log_file = "mpiout.log",direct=None):
	"""
	Run a Bonsai's built in plummer model

	@param[in]	nParticles		number of particles (per mpi process)
	@param[in]	snap_prefix		path prefix for snapshot files (time will be appended)
	@param[in]	T				total simulation time
	@param[in]	dt				internal time step
	@param[in]	dSnap			interval at which snapshot files are generated
	@param[in]	bonsai_bin		path to bonsai exe
	@param[in]	mpi_n			specifies the number of mpi processes (0 = mpi not used)
	@param[in]	mpi_log_file	single log file for mpi output (when mpi_n > 0)
	@param[in] direct			enables direct (N^2) computation

	@returns	None
	"""
	run_mode("sphere",nParticles,snap_prefix,T,dt, dSnap, eps, bonsai_bin, mpi_n,mpi_log_file,direct)

import json

def parse_save_log(logfile="gpuLog.log-1-0", outfile="output"):
	"""
	Parse the resulting gpuLog file

	@param[in]	logfile			filename

	@returns	total_timings 	dictionary with data
	"""

	file = open(logfile)
	lines = file.readlines()
	file.close()

	total_timings = {
		'Sorting':0, # sorting (ms)
		'Data-reordering':0, # moving (ms)
		'Predict':0, # integrated prediction based on current values (ms)
		'Direct_gravity':0, # replaces construction and interaction (ms)
		'Tree-construction':0, # time to build whole tree (ms)
		'setActiveGrpsFunc':0, # ? (ms)
		'Memory':0, # allocation of memory (ms)
		'Grav:':0, #  Time spent to compute approx gravity (s)
		'Build:':0, # Time spent in constructing the tree (incl sorting, making groups, etc.) (s)
		'tPredCor:':0, # predict + correct + energy (s)
		'Correct':0, # time spent correcting predicted after simulation step (ms)
		'Energy':0, # energy error calculation (ms)
		'TOTAL:':0, # Time spent between the start of 'iterate' and the final time-step  (very first step is not accounted) (s)
	}
	# TREE REBUILD  counted in sort+data reordering+tree construction

	count = 0
	for line in lines:
		count += 1
		line = line.replace('\t', ' ')
		line = line.replace('\n', '')
		tokens = line.split(" ")
		if len(tokens) < 4:
			continue
		if tokens[3] in total_timings:
			if tokens[3] == 'TOTAL:':
				# print(tokens)
				total_timings['Grav:'] += float(tokens[7])
				total_timings['Build:'] += float(tokens[15])
				total_timings['tPredCor:'] += float(tokens[28])
				total_timings[tokens[3]] += float(tokens[4])
			else:
				total_timings[tokens[3]] += float(tokens[4]) * 0.001

	# print(json.dumps(total_timings, indent=4))
	if outfile: 
		counter = 0
		filename = outfile + "_{}.json"
		while os.path.isfile(filename.format(counter)):
			counter += 1
		filename = outfile + "_{}"
		outfile = filename.format(counter)


		outfile = open(outfile + ".json","w+")
		outfile.write(json.dumps(total_timings, indent=4))
		outfile.close()
	
	return total_timings