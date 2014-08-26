#!/usr/bin/env python
# encoding: utf-8
'''
	A fourth order butterworth band-pass filter.
	
	Run with -h for usage

 Daniele Trifiro´
 brethil@phy.olemiss.edu


'''
import sys
import getopt

import numpy as np
import cPickle as pickle

from scipy import signal
from data_conditioning import butterworth_band_pass

def __usage__():
	print "Usage: pickled_butter -H f_high -L f_low -S f_sampling file"
	print "\tTakes f_high (high cutoff freq), f_low (low cutoff freq),"
	print "\tsampling frequency and a plain text file and applies"
	print "\ta 4th order (back and forth) butterworth filter."
	print 
	print "Options:"
	print "\t-H, --high f_high\n\
		High Cutoff Frequency"
	print "\t-L, --low f_low\n\
		Low Cutoff Frequency"
	print "\t-S, --sampling f_sampling\n\
		Sampling Frequency"


def __check_options_and_args__():
	global band_high
	global band_low
	global f_sampling
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hH:L:S:" , [ "help", "high=", "low=", "sampling=" ] )
	except getopt.error, msg:
		print msg
		sys.exit(1)	
	if ( len(args) == 0 ):
		print "No arguments given. Exiting."
		__usage__()
		sys.exit(1)
	for o, a in opts:
		if o in ( '-h', '--help' ):
			__usage__()
			sys.exit(1)
		elif o in ('-H', '--high'):
			band_high = float(a)
		elif o in ('-L', '--low'):
			band_low = float(a)
		elif o in ('-S', '--sampling'):
			f_sampling = float(a)	
		else:
			assert False, "Unknown Option."
	if not ( any( flag in o for flag in [ '-L', '--low' ] for o in opts ) and\
	 		any( flag in o for flag in [ '-H', '--high' ] for o in opts ) and\
	 		any( flag in o for flag in [ '-S', '--sampling', 's' ] for o in opts) ):
		print "Sampling frequency, low and high cutoff frequency have"
		print "to be supplied. Exiting.\n"
		__usage__()
		sys.exit(1)
	
	return args
 



def __main__():
	global band_low, band_high
	global f_sampling
	args = __check_options_and_args__()
	file_name = args[0]
	if ( len(args) > 1 ):
		print "Warning: Only first file is processed.\
			Use xargs -i for multiple file processing."
	print file_name
	print "Parameters:\n\t f_sampling = {0}\n\t f_high = {1}\n\t f_low = {2}\n\t file_name = {3}".format(f_sampling, band_high, band_low, file_name)
	data =  np.loadtxt( file_name )
	filtered_data = butterworth_band_pass(data, 4, band_low, band_high, f_sampling)
	pickle.dump( filtered_data, open(file_name+".filtered", "wb") )
	
	
if __name__ == '__main__':
	__main__()
