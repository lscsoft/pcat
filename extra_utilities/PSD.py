#!/usr/bin/env python
# encoding: utf-8
'''

Computes the PSDs of the 
input files (either pickled or plain text) and returns
the computed PSDs in a pickle format, with output file
name 'input_name.psd'.

Daniele Trifiro`
brethil@phy.olemiss.edu
'''

from matplotlib.mlab import psd

from data_conditioning import compute_psd

from utilities_PCAT import *


def __usage__():
	print "Usage:\t fft.py -s sampl_freq file1 file2 file3 ....\n"
	print "\t Computes PSD at the desired resolution and returns pickled \".fft\" files."
	print "\t Output is file1.fft, file2.fft, ..."
	print "\tOutput is a properly normalized PSD, performed through matplotlib.mlab.psd."
	print "\t The density values are scaled by the so that the density is in units of Hz^-1"
	print "\nOptions:"
	print "\t-s sampl_freq, --sampling\n\t\tSampling frequency of the input files."
	print "\t-p, --plot\n\t\tSaves a plot of the Fourier Transform (file1.png, file2.png, ...)"
	print "\t--ascii\n\t\tUse this flag for plain text input data."
	print "\t--out output_size, -o output_size\n\t\tSets the frequency resolution of the psd\n\
		sampling_frequency/output_size gives the frequency resolution."
	print "\t--overlap overlap\n\t\tPercentage of overlap between segments when calculating\n\
	 	Fourier Transform, e.g. \"--overlap 50\"."


def __check_options_and_args__():
	global PLOT
	PLOT = False
	global PICKLED
	PICKLED = True
	global CUSTOM_OUT
	CUSTOM_OUT = False
	global output_size
	global sampling_freq
	global FOURIER_OVERLAP
	FOURIER_OVERLAP = 0.5
	
	
	if ( len( sys.argv[1:] ) == 0 ):
		print "No arguments given."
		__usage__()
		sys.exit(1)
	else:
		try:
			opts, args = getopt.getopt(sys.argv[1:], "s:po:", [ 'plot', 'ascii',
			 							'out=', 'sampling=', '--overlap=' ] )
		except getopt.error, msg:
			print msg
			sys.exit(1)
		for o, a in opts:
			if o in ( '-p', '--plot' ):
				PLOT = True
			elif o in ( '--ascii' ):
				PICKLED = False
			elif o in ( '-s', '--sampling' ):
				sampling_freq = float(a)
			elif o in ( '--out', '-o') :
				CUSTOM_OUT = True
				output_size = int(a)
			elif o in ( '--overlap' ):
					FOURIER_OVERLAP = float(a) / 100.0
			else:
				assert False, "Unknown option."
	if not ( any( flag in o for flag in ['-o', '--out'] for o in opts ) ):
		print "Output file size has to be supplied. Quitting."
		sys.exit()
	if not ( any( flag in o for flag in ['-s', '--sampling'] for o in opts ) ):
		print "Sampling frequency has to be supplied. Quitting."
		sys.exit()
	
	return args	


def plot_psd(frequencies, psd, file_name):
	'''
		Input a frequency array, a PSD array
		and plots the PSD vs frequency in log-log scale
		and saves the plot in 'file_name.png'.
	'''
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(frequencies, psd, linewidth=0.5, color='r' )
	ax.set_xlabel( "Frequency [Hz]" )
	ax.set_ylabel( "Power/Hz [Counts^2/Hz]" )
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim( (1, np.max(frequencies)) )
	ax.xaxis.grid( True, 'minor' )
	ax.yaxis.grid( True, 'minor' )
	ax.xaxis.grid( True, 'major')
	ax.yaxis.grid( True, 'major')
	ax.set_title( "Power Spectral Density" )
	ax.grid()
	fig.savefig(file_name+"_PSD.png", dpi=300)
	plt.close('all')


def __main__():
	global PLOT, pickled, n_files, CUSTOM_OUT, output_size, warning
	args = __check_options_and_args__()
	
	n_files = len(args)
	print "{0}".format(n_files)+" input file"+'s' if n_files > 1 else ''+" Processing..."
	
	input_data = load_data(args, pickled = PICKLED)
	if ( n_files > 1 ):
		bar = progressBar(minValue = 0, maxValue = n_files-1, totalWidth = 40 )
		bar(0)
	if PLOT:
		import matplotlib
		matplotlib.use('agg')
		import matplotlib.pyplot as plt
		 
	for index, element in enumerate(input_data):
		resolution = sampling_freq/(output_size)
		freqs, PSD_estimate = compute_psd(element, resolution, sampling_freq, FOURIER_OVERLAP)
		if ( PLOT ):
			plot_transforms(freqs, PSD_estimate, args[index])
		f = open( args[index]+".psd", "wb")
		pickle.dump( PSD_estimate, f )
		f.close()
		if ( n_files > 1 ):
			bar(index)


if __name__ == '__main__':
	start = time.time()
	__main__()
	print "Total execution time:\t{0:.1f}".format( time.time()-start )