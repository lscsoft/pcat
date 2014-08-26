#!/usr/bin/env python
# encoding: utf-8
# Daniele Trifiro`
# brethil@phy.olemiss.edu

from utilities_PCAT import *


def __usage__():
	print "Usage:\tspikes_time_series.py [options] file1 file2 file3 ...\n\
		file1, file2, file3,... should be finder.py databases.\n\
		Plot all the transients in the database."
	
	print "\nOptions:"
	
	print "\t--output path/to/folder/\n\
		Outputs time series to 'path/to/folder/'."
	print "\t--sampling sampling_frequency\n\
		Sampling frequency"


def __check_options_and_args__():
	global output, sampling_frequency
	sampling_frequency = 0
	output = "./"
	if ( len(sys.argv[1:]) == 0 ):
		print "No arguments given."
		__usage__()
		sys.exit(1)
	else:
		try:
			opts, args = getopt.getopt(sys.argv[1:], "h", [ 'output=', 'help', 'sampling=' ])
		except getopt.GetOptError, err:
			print str(err)
			sys.exit(1)
		for o, a in opts:
			if o in ( '--output' ):
				output = a
			elif o in ( '-h', '--help'):
				__usage__()
				sys.exit(1)
			elif o in ( '--sampling'):
				sampling_frequency = float(a)
			else:
				assert False, "Unknown option."
	return args


def __main__():
	args = __check_options_and_args__()
	spikes_list = list()
	for element in args:
		f = open( element , "rb" )
		spikes_list.extend( pickle.load(f) )
		f.close()
	spikes_waveforms = list()
	number_spikes = len(spikes_list)
	waveform_length = len(spikes_list[0].waveform)
	if ( number_spikes > 1 ):
		maxValue = number_spikes-1
	elif ( number_spikes == 1 ):
		maxValue = 1
	elif ( number_spikes == 0 ):
		print "No spikes in the selected file."
		sys.exit(0)
	bar = progressBar(minValue = 0, maxValue = maxValue-1, totalWidth = 50)
	bar.__call__(0)
	fig = plt.figure()
	if sampling_frequency:
		x_axis = [ i/sampling_frequency for i in range( 0, waveform_length ) ]
	for index, spike in enumerate(spikes_list):
		fig.clf()
		ax = fig.add_subplot(111)
		ax.set_title( "Peak at: "+str(spike.peak_GPS) )
		ax.grid( ls = '--' )
		ax.set_autoscalex_on(False)
		if sampling_frequency:
			plt.xlim( ( 0, waveform_length/sampling_frequency ) )
			ax.plot( x_axis, spike.waveform/spike.norm)
			plt.xlabel("Time [s]")
		else:
			plt.xlim( ( 0, waveform_length ) )
			ax.plot(spike.waveform/spike.norm)
		fig.savefig(output+str(spike.peak_GPS)+".pdf")
		bar.__call__(index)
	plt.close('all')
	if maxValue == 1:
		bar.__call__(1)
	print


if __name__=='__main__':
	__main__()
