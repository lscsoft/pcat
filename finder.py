#!/usr/bin/env python
# encoding: utf-8
'''
Daniele Trifiro`
brethil@phy.olemiss.edu

finder.py contains definitions and functions for a simple glitch finder.
Glitches are found by looking for the points in the time series which 
go above a certain threshold which is a multiple of the standard deviation
of the segment.

Can be used standalone. For usage: run with -h.


''' 


from numpy import linalg

from utilities_PCAT import *

from data_conditioning import median_mean_average_psd

def usage():
	print "Usage:\t finder.py -t threshold -w width --sampling sampl_freq\n\
		[--ascii, -a] [-i]\n\
		[--remove_seconds seconds] [--energy] [--timing]\n\
		[-o FILE, --output file] file1 file2 file3 ..."
	
	print "\n\tTakes data in pickled or plain text files (file1, file2, ...).\n\
	Returns a pickled (or plain text) list of Spike() class objects (see spike_class.py)\n\
	containing transients and their properties, e.g., GPS peak time, peak absolute value."
	
	print "\n\tOptions:"
	print "\t-w width, --width width\n\
		Number of sampling points defining the transient."
	
	print "\t --threshold N, -t N\n\
		Sets the trigger threshold at N times the\n\
		standard deviation of the analyzed time."
	
	print "\t--sampling\n\
		Sampling frequency of the channel."
	
	print "\t --remove_seconds seconds\n\
		Removes time (in seconds) from the\n\
		beginning and the end of the segment to avoid FFT artifacts.\n\
		(Default = 0.5 seconds)."
	
	print "\t--time_resolution\n\
		Minimum distance (in points) between 2 transients\n\
		to be considered distinct. Must be smaller or equal to width.\n\
		(Default = 500)."
	
	print "\t--output FILE, -o FILE\n\
		Sets a custom output file.\n\
		Default is 't-n_w-points.list', where 't' is the\n\
		threshold specified by the -t option and points is\n\
		the width specified by the -w option."
	
	print "\t-i\n\
		Use this flag for ASCII input data."
	
	print "\t--ascii, -a\n\
		Saves output as plain text (ASCII)."
	
	print "\t--timing\n\
		Prints timing information."
	
	print "\t--energy\n\
		Set the normalization constant to the transient's energy."
	

def check_options_and_args():
	"""
		This function parses the arguments and option of the program.
	"""
	global args
	global opts
	global OUTPUT
	# Variables Definitions
	global UNPICKLED_INPUT
	UNPICKLED_INPUT = False
	global CUSTOM_WIDTH
	CUSTOM_WIDTH = False
	
	global max_width
	global N
	
	global CUSTOM_THRESHOLD
	CUSTOM_THRESHOLD = False
		
	global CUSTOM_RSECONDS
	CUSTOM_RSECONDS = False
	
	global removed_points, threshold
	
	global CUSTOM_RESOLUTION
	CUSTOM_RESOLUTION = False
	
	global removed_seconds
	global removed_points
	
	
	global time_resolution
	CUSTOM_RESOLUTION = False
		
	global CUSTOM_OUTPUT
	CUSTOM_OUTPUT = False
		
	global NOTPICKLED
	NOTPICKLED = False
	
	global TIMING
	TIMING = False
		
	global ENERGY
	ENERGY = False
	global sampling_frequency
	
	if ( len(sys.argv[1:]) == 0 ):
		print "No arguments given."
		usage()
		sys.exit(1)
	else:
		try:
			opts, args = getopt.getopt(sys.argv[1:], "hw:t:c:o:i", ["help", "width=",
			 						"threshold=", "remove_seconds=", "cutting_width=",
			 						"output=", "ascii", "timing", "time_resolution=",
			 						'sampling=', "energy"])
		except getopt.GetoptError, err:
			print str(err)
			sys.exit(1)
		for o, a in opts:
			if o in ( '-h', '--help'):
				print "Help:\t",
				usage()
				sys.exit()	
			elif o in ( '-w', '--width' ):
				CUSTOM_WIDTH = True
				max_width = int(a)
			elif o in ( '-t', '--threshold' ):
				CUSTOM_THRESHOLD = True
				N = float(a)
			elif o in ( '--output', '-o' ):
				CUSTOM_OUTPUT = True
				OUTPUT = str(a)	
			elif o == ( '-i' ):
				UNPICKLED_INPUT = True
			elif o in ( '--remove_seconds' ):
				CUSTOM_RSECONDS = True
				removed_seconds = float(a)
			elif o in ( '--time_resolution' ):
				CUSTOM_RESOLUTION = True
				time_resolution = int(a)
			elif o in ( '--ascii' ):
				NOTPICKLED = True
			elif o == ( '--timing' ):
				TIMING = True
			elif o in ( '--sampling' ):
				sampling_frequency = float(a)
			elif o in ( '--energy' ):
				ENERGY = True
			else:
				assert False, "Unknown option."
	if not ( any( flag in o for flag in ['-w', '--width'] for o in opts ) 
			and any( flag in o for flag in ['-t', '--threshold'] for o in opts)
			and any( '--sampling' in o for o in opts) ):
		print "The options:\n-w or --width sampling_width\n-t or --threshold N"
		print "--sampling sampl_freq\n*have* to be supplied."
		#usage()
		sys.exit(1)
	
	
	# This part sets some parameters defaults if these aren't supplied as
	# arguments:
	
	if not CUSTOM_RSECONDS:
		removed_seconds = 0.5
	if not CUSTOM_RESOLUTION:
		time_resolution = 500
	# Number of points removed from the tails of the data.
	removed_points = int( np.ceil(removed_seconds*sampling_frequency) )
	if ( removed_points == 0 ):
		# removed_points will be 0 only when removed_seconds is 0,
		# in this case set removed_points to 1 point. 
		# we are not losing information anyway.
		removed_points = 1
	if not CUSTOM_WIDTH:
		max_width = 1000
	if not CUSTOM_OUTPUT:
		OUTPUT = "t-"+str(N)+"_w-"+str(max_width)
	return args


def save_segment(list, output):
	'''
	This function saves the spike information in a plain text file.
	The saved information is:
		Segment_start (GPSTIME)
		Segment_end (GPSTIME)
		Spike_start (point index)
		Spike_end (point index)
		Spike_GPS_peak (GPSTIME)
	
	The GPS time for the Spike_* variables can be obtained by dividing the 
	point index by 'sampling_frequency' and adding it
	to Segment_start. This yields the GPS time.
	'''
	
	counter = 0
	f = open(output+".txt", "w")
	f.write("# Number\t Segment_start\tSegment_end\tSpike_start\tSpike_end\tSpike_GPS_peak\tWaveform\n")
	for spike in list:
		if ( spike.width > 0 ): 
			f.write(str(counter+1)+"\t"+str(spike))
			try:
				for point in spike.waveform:
					f.write(str(point/spike.norm)+" ")
				f.write('\n')
			except AttributeError:
				f.write("\n")
				print "Error writing waveform number", counter+1
			counter += 1
	f.close()


def save_pickled(spike_list, output):
	'''
	This function saves the spike information in a pickled file.
	We are saving a list of the Spike() class istances.
	The saved information is:
		Segment_start (GPSTIME)
		Segment_end (GPSTIME)
		Spike_start (point index)
		Spike_end (point index)
		Spike_peak (point index)
		
	The GPS times for the Spike_* variables can be obtained by dividing the
	points number by 'sampling_frequency' and adding it
	to Segment_start. This yields the GPS time.
	
	'''	
	
	counter = 0
	for spike in spike_list:
		counter +=1
	f = open(output+".list", "wb")
	pickle.dump(spike_list, f)
	f.close()


def print_parameters():
	'''
		This function prints the parameters of the run.
	'''
	print "#"*frame_width
	print "Parameters for this run:\t ("+str(len(args))+" files)\n"
	print "\t Defining points:\t", max_width
	print "\t Excluded seconds:\t{0:.1f} ({1} points)".format(removed_seconds, removed_points)
	print "\t Threshold:\t"+str(N)+" sigma"
	global OUTPUT
	if NOTPICKLED:
		OUTPUT = str(OUTPUT)+".txt"
		print "\t Output:\t"+str(OUTPUT)+".txt"
		print "\t Output is plain text.\n"
	else:
		print "\t Output:\t"+str(OUTPUT)+".list"
		print "\t Output is pickled.\n"


def find_spikes_algorithm(data, removed_points, f_sampl, threshold, time_resolution, data_name, spike_width):
	'''
		This function searches for the spikes in the data segments
		and saves parameters for each found spike in the attributes of the
		"Spike" class. Spikes are saved in a list which is what is returned by the
		function.
		The data_name variable is used to retrieve the start and end GPS times, assuming
		they are in the form same form as
			L-R-L1:PSL-FSS_FAST_MON_OUT_DQ_GPSSTART-GPSEND.data.something
			IFO-FRAME_TYPE-CHANNEL_GPSSTART-GPSEND.data.something
		Arguments:
		- data (array)
			Time series
		- removed_points (integer)
			The number of points that should be excluded from the analysis
			from the start and the beginning of the time series
		- f_sampl (float)
			Sampling frequency
		- time_resolution (integer)
			Number of point between two triggers (points above threshold)
			to be considered separate transients.
		- data_name (string)
			Name of the file containing 'data'
		- spike_width (integer)
			Number of points at which the found transients should be sampled.
	Output:
		- spikes (list)
			A list fo Spike() class istances
	'''
		
	spikes = []
	"""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	( start, end ) = (data_name.split("/")[-1]).split('.')[0].split('_')[-1].split('-')
	t = np.linspace(int(start), int(end), len(data[removed_points:-removed_points]))
	ax.plot(t, data[removed_points:-removed_points])
	ax.plot(t, threshold*np.ones(len(data[removed_points:-removed_points])))
	ax.plot(t, -threshold*np.ones(len(data[removed_points:-removed_points])))
	fig.savefig("{0}-{1}.svg".format(start, end))
	plt.close()
	"""
	
	to_analyze = data[removed_points:-removed_points]
	
	last_spike_index = 0		# Last spike index is the last index of the 
								# spike above threshold
	max_spike_index = 0
	max_spike_value = 0
	HAS_SPIKE = False
	
	( start, end ) = (data_name.split("/")[-1]).split('.')[0].split('_')[-1].split('-')
	
	# Choose NFFT to have a frequency resolution of 1 Hz or better
	
	if (spike_width < f_sampl):
		window = np.hanning(spike_width)
		window_norm = (window**2).sum()/spike_width
		freqs, psd = median_mean_average_psd(to_analyze, int(f_sampl), f_sampl)
	else:
		window = np.hanning(f_sampl)
		window_norm = (window**2).sum()/float(f_sampl)
		freqs, psd = median_mean_average_psd(to_analyze, spike_width, f_sampl)
	
	delta_t = 1.0/f_sampl
	
	segment_length = len(to_analyze)
	for index, point in enumerate(to_analyze):
		if (abs(point) > threshold):
			if not HAS_SPIKE:
				first_spike_index = index
				HAS_SPIKE = True
				max_spike_index = index
				max_spike_value = abs(point)
			elif (abs(point) > max_spike_value):
				max_spike_value = abs(point)
				max_spike_index = index
			
			last_spike_index = index
		
		elif (index-max_spike_index > time_resolution ) and HAS_SPIKE:
			# We're now outside the search range for the spike: save the 
			# current spike and start over.
						
			# Discard the glitch IF there's less than 2 points above threshold
			if ( (last_spike_index-first_spike_index) < 2 ):
				HAS_SPIKE = False
				max_spike_value = 0
				max_spike_index = 0
				continue
			
			# Save waveform
			waveform = data[max_spike_index+removed_points-spike_width/2:max_spike_index+removed_points+spike_width/2] 
			# Instantiate Spike object
			first_index, last_index = first_spike_index+removed_points, last_spike_index+removed_points
			max_index = max_spike_index+removed_points
			peak_GPS = (int(start)+ ( max_index / f_sampl ))
			spike = Spike(first_index, last_index,
							max_index, max_spike_value,
							peak_GPS, int(start), int(end),
							waveform, f_sampl)
			
			# The squared SNR per unit frequency for a signal g(t) is defined as
			#	SNR^2(f) = 2 * |g(f)|^2/Pxx(f)
			# Factor of two beause the numerator should be g(f)*g_conj(f) + g_conj(f)*g(f)
			# where g(f) is the Fourier transform of g(t) and Pxx is the 
			# detector spectrum.
			# Thus the total SNR:
			#	SNR^2 = 4*\int_0^\infty |g(f)|^2/Pxx(f) df
			# Since g(f) is symmetric around f  (time series is real)/
			
			# Factor of two in psd because rfft is one sided.
			spike.psd = 2 * 1.0/window_norm * 1.0 *  np.abs(delta_t*np.fft.rfft(spike.waveform*window, n=int(f_sampl)))**2
			spike.psd[0] /= 2.0
			spike.fft_freq = freqs
			
			spike.segment_psd = psd
			
			# We don't need the factor of 4 in front of the integral because both spike.psd and psd
			# are one-sided and are correctly normalized
			spike.SNR = np.sqrt( 4* (np.array(spike.waveform)**2).sum() * 2 * f_sampl )
			
			# Check spike polarity
			if (spike.waveform[np.argmax(np.abs(spike.waveform))] > 0):
				spike.polarity = 1
			else:
				spike.waveform *= -1
				spike.polarity = -1
			
			# Save Spike object
			spikes.append(spike)
			
			HAS_SPIKE = False
			max_spike_value = 0
			max_spike_index = 0
	
	# Save last spike
	if HAS_SPIKE:
		# Save waveform
		waveform = data[max_spike_index+removed_points-spike_width/2:max_spike_index+removed_points+spike_width/2] 
		# Instantiate Spike object
		first_index, last_index = first_spike_index+removed_points, last_spike_index+removed_points
		max_index = max_spike_index+removed_points
		peak_GPS = (int(start)+ ( max_index / f_sampl ))
		spike = Spike(first_index, last_index,
						max_index, max_spike_value,
						peak_GPS, int(start), int(end),
						waveform, f_sampl)
		
		# The squared SNR per unit frequency for a signal g(t) is defined as
		#	SNR^2(f) = 2 * |g(f)|^2/Pxx(f)
		# Factor of two beause the numerator should be g(f)*g_conj(f) + g_conj(f)*g(f)
		# where g(f) is the Fourier transform of g(t) and Pxx is the 
		# detector spectrum.
		# Thus the total SNR:
		#	SNR^2 = 4*\int_0^\infty |g(f)|^2/Pxx(f) df
		# Since g(f) is symmetric around f  (time series is real)/
		
		# Factor of two in psd because rfft is one sided.
		spike.psd = 2 * 1.0/window_norm * 1.0 *  np.abs(delta_t*np.fft.rfft(spike.waveform*window, n=int(f_sampl)))**2
		spike.psd[0] /= 2.0
		spike.fft_freq = freqs
		
		spike.segment_psd = psd
		
		# We don't need the factor of 4 in front of the integral because both spike.psd and psd
		# are one-sided and are correctly normalized
		spike.SNR = np.sqrt( (np.array(spike.waveform)**2).sum() * 2 * f_sampl )
		
		# Check spike polarity
		if (spike.waveform[np.argmax(np.abs(spike.waveform))] > 0):
			spike.polarity = 1
		else:
			spike.waveform *= -1
			spike.polarity = -1
		
		# Save Spike object
		spikes.append(spike)
		
		HAS_SPIKE = False
		max_spike_value = 0
		max_spike_index = 0
		
	return spikes

def find_spikes(data, metadata, threshold, spike_width, time_resolution, removed_seconds, f_sampl, normalization=None):
	'''
		Load all the files in the 'file_list' list and searches for spikes, using the
		find_spikes_algorithm.
		
		
		Arguments:
			data:
					time series to analyze (numpy array)
			metadata:
					contains a string with channel and start-end time for the segment.
					This is used to save the GPS peak time for the found transients
			time_resolution:
					Minimum distance between two neighouring spikes to be considered different events
			removed_seconds:
					The number of seconds to exclude from the beginning and the end of the segments, in
					order to avoid fourier transform ringing artifacts.
			f_sampl:
					Sampling frequency
			normalization:
			 		Sets how the found spikes are normalized.
					Can be one in "energy", "amplitude" or *None*, default is none.
		
					If energy, spikes are normalized to unit energy, if amplitude, spikes are normalized to
					unit maximum amplitude, if None, spikes are not normalized.
		
	'''
	spikes_list = []
	spikes_number = 0
	# Calculate the standard deviation, excluding the the first and last
	# 'removed_seconds' seconds, in order to exclude ringing artifacts
	# due to the fourier transforms.
	
	removed_points = int(removed_seconds*f_sampl)
		
	sigma = np.std(data[removed_points:-removed_points])
	found_spikes = find_spikes_algorithm(data, removed_points, f_sampl, threshold*sigma,\
	 								time_resolution, metadata, spike_width)
	
	if ( len(found_spikes) != 0 ):
		spikes_number += len(found_spikes)
		spikes_list.extend( found_spikes )
	
	for spike in found_spikes:
		if ( normalization == "energy" ):
			spike.norm = linalg.norm(spike.waveform)
		elif ( normalization == "amplitude" ):
			spike.norm = np.max(np.abs(spike.waveform))
		elif ( normalization == None ):
			spike.norm = 1.0
		else:
			print "Other normalizations have yet to be implemented" 
			assert False
			
	return spikes_list


def main():
	args = check_options_and_args()
	
	spikes_list = list()
	spikes_number = 0 
	print "#"*frame_width+"\n\tLoading and processing data:\n"
	
	spikes_list = find_spikes(args, N, max_width, time_resolution, removed_seconds, sampling_frequency)
	# Save data:
	if NOTPICKLED:
		save_segment(spikes_list, OUTPUT)
	else:
		save_pickled(spikes_list, OUTPUT)



if __name__ == '__main__':
	main()
