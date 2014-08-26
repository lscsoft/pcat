#!/usr/bin/env python
# encoding: utf-8

'''

This program is to be used to create databases to be used with GMM.py and
PCA.py. A wide variety of inputs and outputs is supported, depending on the
given options. For usage, run with -h.

Daniele Trifiro`
brethil@phy.olemiss.edu

''' 
from numpy import linalg

from utilities_PCAT import *


global frame_width
frame_width = 78

def __usage__():
	"""
		Usage
	"""
	
	print "Usage:\t database.py [options] file1 file2 file3 ....\n"
	
	print "\tThis program is to be used to prepare data for use in GMM.py and\n\
	PCA.py.\n\
	Input type can be supplied through options. Output is a pickled file usable\n\
	by GMM.py and PCA.py. Output name depends on the arguments:"
	
	print "\t\t\t> PCAT.py lists\t ---> \"merged_database.data\""
	print "\t\t\t> Generic files\t\t ---> \"database.data\""
	print "\t\t\t> Power Spectra lists\t ---> \"PSD_database.data\""
	print "\t\t\t> Matrix file\t\t ---> \"matrix_database.data\""
	
	print ""
	print "#"*frame_width
	
	print "Options:"
	print "-"*frame_width
	print "\t--merge\n\
		Merge two or more databases into a single file.\n\
		Databases are output from or from PCAT.py (or finder.py)."
		
	print "\t--psd\n\
		Creates a database from a list of files (pickled arrays)\n\
		power spectra. Power spectra can be computed though fourier.py.\n\
		File name should be in the form\n\
		'IFO-frame_type-channel_name-GPSSTART-GPSEND.data.extension', in order\n\
		to correctly extract the GPS start and GPS end from the file name."
	
	print "\t--matrix\n\t\tCreates a database using the input matrix. Rows are"
	print "\t\tobservations, columns are variables. Matrix file can be plain text\n\
		or pickled (numpy array)."
	
	print "\t--list\n\
		Creates a database from a list of generic files, each\n\
		containing one observation. Input files should be in pickled\n\
		format, use --ascii for ASCII input data.\n\
		for ASCII input data"
	
	print "\t--ascii\n\t\tASCII input data."
	print "-"*frame_width
	
	print "\t--output output, -o output\n\t\t Set custom output file name."
	print "\t--silent\n\t\tSilent, do not print output on std outoput."


def __check_options_and_args__():
	global PICKLED
	PICKLED = True
	global INPUT, CUSTOM_OUTPUT, OUTPUT_NAME
	CUSTOM_OUTPUT = False
	global SILENT
	SILENT = False
	if ( len( sys.argv[1:] ) == 0 ):
		print "No arguments given."
		__usage__()
		sys.exit(0)
	else:
		try:
			opts, args = getopt.getopt( sys.argv[1:], "ho:s", [ 'help', 'ascii',
			 											'output=', 'merge',
														'psd', 'matrix',
														'list', 'silent'] )
		except getopt.error, msg:
			print msg
			sys.exit(1)
		for o, a in opts:
			if o in ( '--ascii' ):
				PICKLED = False
			elif o in ( '-h', '--help' ):
				__usage__()
				sys.exit(0)
			elif o in ( '-o', '--output' ):
				CUSTOM_OUTPUT = True
				OUTPUT_NAME = a
			elif o in ( '--merge' ):
				INPUT = 'merge'
			elif o in ( '--matrix' ):
				INPUT = 'matrix'
			elif o in ( '--list' ):
				INPUT = 'list'
			elif o in ( '--psd' ):
				INPUT = 'fft'
			elif o in ( '--silent', '-s' ):
				SILENT = True
			else:
				assert False, "Unknown option."
	if not ( any( flag in o for flag in [ '--psd', '--matrix',
	 									'--list', '--merge'] for o in opts) ):
		print "Input options have to be specified.\n\tRun create_database.py with"
		print "\tno options or \"-h\" or \"--help\" for a list of possible options." 
		sys.exit(1)
	return args	


def load_database(data):
	''' 
		Loads multiple finder.py databases and merges them into a single database
		
		Arguments: 
			- data (list)
				A list containing the filenames or paths of the files to be loaded.
		
	'''
	database = list()
	n_files = len(data)
	if not SILENT:
		bar = progressBar(minValue = 0, maxValue = n_files-1 if n_files > 1 else 1, totalWidth = frame_width/2 )
		bar.__call__(0)
	for index, element in enumerate(data):
		if not SILENT:
			bar.__call__(index)
		f = open(element, "rb")
		database.extend( pickle.load(f) )
		f.close()
	if not SILENT:
		print ""
	return database


def load_fft(data, pickled):
	'''
		Loads a list of files, given a file list and returns list
		of Spike() istances, that can be pickled and loaded and used
		into GMM.py and PCA.py
		
		Arguments:
			- data (list)
				List of file names or paths to the files to be loaded
			- pickled (bool)
				If true, input files are pickled (numpy arrays) and
				loaded through cPickle, otherwise they are plain text
				and loaded through np.loadtxt
				
		
	'''
	
	files_n = len(data)
	database = list()
	if not SILENT:
		bar = progressBar( minValue = 0, maxValue = files_n-1 if files_n > 1 else 1, totalWidth = frame_width/2 )
		bar.__call__(0)
	for index, element in enumerate(data):
		if pickled:
			f = open(element, "rb")
			tmp = pickle.load(f)
			f.close()
		else:
			tmp = np.loadtxt(element)
		
		# This part is used to extract the name of the segment from the file name
		# the if is needed to correctly hand cases when files are in the same
		# folder and when they are not (i.e. path to file is given)
		if ( "/" in element ):		# Different folder, Path to the file is given
			( start, end ) = element.split("/")[-1].split('.')[0].split('_')[-1].split('-')
		else:					# Same folder, only file name is given
			( start, end ) = element.split('.')[0].split('_')[-1].split('-')
		
		
		normalization = 1.0
		spike = Spike( int(start), int(end), 0, normalization )
		spike.waveform = list()
		spike.waveform.extend(tmp)
		database.append(spike)
		if not SILENT:
			bar.__call__(index)
	print ""
	return database


def load_matrix(data, pickled):
	'''
		Given a list of matrix files, load each matrix interpreting
		rows as observations. Return a database to be used with
		GMM.py and PCA.py
		
		Arguments:
			- data (list)
				List of paths or names to the matrix files to be loaded
			- pickled (bool)
				Pickled input.
				If pickled, data is loaded through cPickle,
				otherwise np.loadtxt is used.
	'''
	database = []
	for element in data:
		if not pickled:
			observations = np.loadtxt(data).tolist()
		else:
			f = open(data, "rb")
			observations = pickle.load(f).tolist()
			f.close()
		normalization = 1.0
		for element in observations:
			spike = Spike(0, 0, 0, normalization)
			spike.waveform = np.array(element)
			database.append(spike)
	return database


def load_list(data, pickled):
	'''
		Given a list of files, each of these is loaded as
		a single observation. A database usable with GMM.py and PCA.py
		is returned.
		
		Arguments:
			- data (list)
				List of paths or names to the files to be loaded
			- pickled (bool)
				Pickled input files.
				If pickled, data is loaded through cPickle,
				otherwise np.loadtxt is used.
	'''
	database = list()
	files_n = len(data)
	if not SILENT:
		bar = progressBar( minValue = 0, maxValue = files_n-1 if files_n > 1 else 1, totalWidth = frame_width/2 )
		bar.__call__(0)
	for index, element in enumerate(data):
		bar.__call__(index)
		if pickled:
			f = open(element, "rb")
			tmp = pickle.load(f)
			f.close()
		else:
			tmp = np.loadtxt(element)
		spike = Spike(0, 0, 0, 1.0)
		spike.waveform = np.array(tmp)
		database.append(spike)
	return database



def __main__():
	args = __check_options_and_args__()
	if not SILENT:
		print "Loading data and creating database..."
	if ( INPUT == 'merge'):
		database = load_database(args)
		output_name = "merged_database.data"
	elif ( INPUT == 'matrix'):
		database = load_matrix(args, PICKLED)
		output_name = "matrix_database.data"
	elif ( INPUT == 'fft'):
		database = load_fft(args, PICKLED)
		output_name = "psd_database.data"
	elif ( INPUT == 'list'):
		database = load_list(args, PICKLED)
		output_name = "database.data"
	print ""
	if CUSTOM_OUTPUT:
		output_name = OUTPUT_NAME
	if not SILENT:
		print "Saving "+output_name+"..."
	output = open( output_name, "wb" )
	pickle.dump( database, output )
	output.close()
	if not SILENT:
		print "Done!"


if __name__ == '__main__':
	__main__()
	