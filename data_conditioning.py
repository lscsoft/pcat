#!/usr/bin/env python
# encoding: utf-8
"""
Created by Daniele Trifirò on 2013-07-29.
brethil@phy.olemiss.edu

data_conditioning.py contains definitions and functions to condition data.
Can be used standalone. For usage: run with -h.

Functions (and variables) defined here are imported into PCAT.py and used
for data conditioning, though it can be run as a standalone program to
condition a single (or multiple files).

Contains:
	Global variables (used in PCAT.py):
		ANALYSIS_FREQUENCY (float)
			The frequency at which data is resampled to when
			whitening the data with whiten().
		HIGH_PASS (boolean)
			Whether to apply an high pass filter in addition to
			whitening (butterworth filter).
		HIGH_PASS_ORDER (integer)
			Order of the high pass filter.
		HIGH_PASS_CUTOFF (float)
			High pass filter cutoff frequency.
		WHITENING_REMOVED_SECONDS_FACTOR  (float)
			This is used to set the number of seconds replaced
			with zeros to suppress high pass filter transients. 
	
	Functions:
		- high_pass_coefficients()
		- high_pass_filter()
		- low_pass_filter
		- whiten()
		- butterworth_band_pass()
		- median_mean_average_psd()
		- compute_psd()
		- median_bias_factor()
	
"""

from scipy import signal, interpolate
from utilities_PCAT import *


#################### PARAMETER DEFINITIONS ####################

######DEPRECATED
# Maximum sampling frequency. All channels above this frequency
# are downsampled to ANALYSIS_FREQUENCY
#global ANALYSIS_FREQUENCY
#ANALYSIS_FREQUENCY = 4096.0


# High pass and low pass Butterworth filter orders.
HIGH_PASS_ORDER = 4

# High pass filter cutoff frequency
HIGH_PASS_CUTOFF = 40.0


###############################################################
################# END OF PARAMETER DEFINITIONS ################
################################################################

##### DISCRETE FOURIER TRANSFORM DEFINITIONS (whitening) #####
# Direct Fourier transform of a signal s_n:
#  h_k = dt \sum_{n=0}^{N-1} s_n exp(2 pi i n k/N)
# Inverse Fourier transform:
# s_n = df \sum_{k=0}^{N-1} h_k exp(-2 pi i n k/N)
# Using numpy's built in DFT routines, these become:
# h = dt * np.fft.rfft(s)
# s = df * np.fft.irfft(h)
#
# dt is the sampling period, df=1/(N*dt), N is the number of points
# of s: N=len(s)

def usage():
	'''
		Usage
	'''
	print "Usage:\t data_conditioning.py -s sampl_freq  [options] file1 file2 ....\n"
	print "\t\tCondition file1, file2, file3 by applying a withening filter.\n\
		Input files are resampled to %i Hz and highpassed (%i Hz).\n\
		Use --noresample and --nohighpass to avoid this or change\n\
		the highpass cutoff frequency with --highpasscutoff.\n" % (ANALYSIS_FREQUENCY, HIGH_PASS_CUTOFF)
	
	
	print "\t Options:"
	
	print "\t\t--sampling sampl_freq, -s sampl_freq\n\
			Sampling frequency for the given files"
	
	print "\t\t[--excluded_seconds seconds]\n\
			Number of seconds that are excluded from the beginning\n\
			and the start of the segment.\n\
			(Default "+str(WHITENING_REMOVED_SECONDS_PERCENTAGE*100)+"% of the total length of input)."
	
	print "\t\t[--pickled, -p]\n\
			Pickled input files.\n\
			These are numpy arrays or python lists containing the time series."
			
	print "\t\t[--noresample]\n\
			Do not resample input data."
	
	print "\t\t[--nohighpass]\n\
			Do not apply high pass filter to input data."
			
	print "\t\t[--highpasscutoff]\n\
			High pass filter cutoff frequency (default %.1f Hz)." % HIGH_PASS_CUTOFF
	
	print "\t\t[--remove_seconds] (this has no effect if using --nohighpass)\n\
			If downloading data with download_frames.py and using the --pad option\n\
			use --remove_seconds paddding_seconds.\n\
			This is used to suppress high pass filter ringing artifacts at the beginning\n\
			and the end of the segment before whitening.\n"
	


def check_options_and_args():
	global PLOT
	PLOT = False
	global PICKLED
	PICKLED = False
	global CUSTOM_OUT
	CUSTOM_OUT = False
	global excluded_seconds
	global SAMPLING_FREQUENCY
	global RESAMPLE, HIGH_PASS, HIGH_PASS_CUTOFF
	
	excluded_seconds = -1
	
	RESAMPLE = True
	HIGH_PASS = True
	
	if ( len( sys.argv[1:] ) == 0 ):
		print "No arguments given. Run with -h or --help for usage."
		print "Sampling frequency for the input files has to be supplied through -s or --sampling."
		sys.exit(1)
	else:
		try:
			opts, args = getopt.getopt(sys.argv[1:], "s:ph",  [ 'help', 'pickled' ,'sampling=','noresample',\
			 												'nohighpass', 'highpasscutoff=', 'remove_seconds='] )
		except getopt.error, msg:
			print msg
			sys.exit(1)
		for o, a in opts:
			if o in ( [ '-h', '--help' ] ):
				usage()
				sys.exit()
			elif o in ( [ '-s', '--sampling' ] ):
				SAMPLING_FREQUENCY = float(a)
			elif o in ( [ '--pickled', '-p' ] ):
				PICKLED = True
			elif ( o == '--noresample'):
				RESAMPLE = False
			elif ( o == '--nohighpass' ):
				HIGH_PASS = False
			elif ( o == '--highpasscutoff' ):
				HIGH_PASS_CUTOFF = float(a)
			elif ( o == '--remove_seconds' ):
				excluded_seconds = float(a)
			else:
				assert False, "Unknown option."
	if not ( any( flag in o for flag in ['-s'] for o in opts ) ):
		print "Sampling frequency has to be supplied. Quitting."
		sys.exit()
	
	
	
	# If number of removed seconds has not been supplied, 
	if (excluded_seconds > 0):
		excluded_seconds /= WHITENING_REMOVED_SECONDS_FACTOR	
	return args


def high_pass_coefficients(frequencies, cutoff, f_sampl, order):
	'''
		Compute frequency response at the frequencies given by
		'frequencies' for a high pass butterworth filter with cutoff
		frequency 'cutoff' and order 'order'.
		
		Arguments:
			- frequencies (list or array)
				List of frequencies at which the frequency response
				should be computed at.
			- cutoff (float)
				Cutoff frequency
			- f_sampl (float)
				Sampling frequency
			- order (integer)
				Order of the filter
		
		Output:
			- h (array)
				Contains the (complex) frequency response of the
				filter at the frequencies given by 'frequencies'.
	'''
	nyquist_frequency = f_sampl/2.0
	
	# Compute the coefficients:
	b, a = signal.butter(order, cutoff/nyquist_frequency, btype='highpass')
		
	frequencies = (frequencies/nyquist_frequency)*np.pi
	w, h = signal.freqz(b, a, worN=frequencies)
	
	return h


def high_pass_filter(time_series, cutoff, f_sampl, order):
	'''
		Butterworth high pass filter.
		
		Arguments:
			- time_series (array)
			- cutoff (float)
				Cutoff frequency
			- f_sampl (float)
				Sampling frequency
			- order (integer)
				Order of the filter
		
		Output:
			- high_passed_time_series (array)
	'''
	nyquist_frequency = f_sampl/2.0
	
	# Compute the coefficients:
	b,a = signal.butter(order, cutoff/nyquist_frequency, btype='highpass')
	
	# Apply zero phase filter using filtfilt (forward-backward linear phase 
	# filter)
	high_passed_time_series = signal.filtfilt(b, a, time_series)
	
	return high_passed_time_series


def low_pass_filter(time_series, cutoff, f_sampl, order):
	'''
		Butterworth low pass filter.
		
		Arguments:
			- time_series (array)
			- cutoff (float)
				Cutoff frequency
			- f_sampl (float)
				Sampling frequency
			- order (integer)
				Order of the filter
		
		Output:
			- high_passed_time_series (array)
	'''
	nyquist_frequency = f_sampl/2.0
	
	# Compute the coefficients:
	b,a = signal.butter(order, cutoff/nyquist_frequency, btype='low')
	
	filtered_time_series = signal.filtfilt(b, a, time_series)
	
	return filtered_time_series


def whiten(time_series, excluded_seconds, f_sampl, resample_freq, highpass=True, highpass_cutoff=HIGH_PASS_CUTOFF, resample=True):
	'''
		Whitens the input time series.
		
		If resample=True, data is resampled to ANALYSIS_FREQUENCY (global variable). 
		
		The time series is conditioned by multiplying the coefficients
		of the Discrete Fourier Transform by the whitening coefficients
		and (if requested) by the high pass coefficients
		
		Whitening is performed by multiplying each frequency bin in
		the frequency spectrum of the original time series by the square
		root of the inverse of the Power Spectral Density.
		The whitening coefficients are thus obtained simply computing
		the square root of the inverse (and a multiplying constant).
		
		The PSD is estimated with a frequency resolution of 1 Hz through
		the median-mean average algorithm. See the FINDCHIRP Paper (qc/0509116)
		
		Outputs the whitened (and high-passed, if requested) time series.
		
		Arguments:
			- time_series (array)
				Time series to be whitened
			- excluded_seconds (int)
				Number of seconds to be replaced with zeros at the beginning
				and start of the time series to remove filtering artifacts.
				This should be the same number used with download_frames.py's
				--pad option.
			- f_sampl (float)
				Sampling frequency for the input time series.
			- highpass (boolean, optional, default=True)
				True for high-pass filtering.
			- resample (boolean, optional, default=True)
				Downsample the data to ANALYSIS_FREQUENCY if True.
		
		Output: 
			- whitened_time_series (numpy array)
				The whitened (and high-passed, if requested) 
				time series.
		
	'''
	
	# If requested and if necessary, resample the signal to ANALYSIS_FREQUENCY
	# (defined at the top of this document)
	
	if ( resample and (f_sampl > resample_freq) ):
		downsample_factor = int(f_sampl/resample_freq)
		time_series = decimate(time_series, downsample_factor)
		
		# Update sampling frequency to the new value
		f_sampl = resample_freq
		
	
	time_series_duration = len(time_series)
	x_axis = np.arange(time_series_duration)
	high_passed_time_series = high_pass_filter(time_series, HIGH_PASS_CUTOFF, f_sampl, HIGH_PASS_ORDER)
	
	# FIXME: Having 0.05 is horrible. Try to find a decent way to do this
	removed_points = int(0.05*excluded_seconds*f_sampl)
	# Suppress high pass filter transients by replacing slices at the beginning 
	# and the end of the segment with zeros
	high_passed_time_series[0:removed_points] = 0.0
	high_passed_time_series[-removed_points:] = 0.0
	
	
	# Take the FFT of the input time series,
	# accounting for the correct normalization
	delta_t = (1./f_sampl)
	transform = delta_t * np.fft.rfft(high_passed_time_series)
	transform_len = len(transform)
	delta_f = 1./(delta_t*len(transform))
	
	# Set the frequency resolution of the PSD (1 Hz)
	bin_hz = 1.0
	fft_length = int(f_sampl/bin_hz + 0.5)
	# Compute PSD using the median mean average PSD Algorithm
	# Cut the first and last "excluded_seconds" to avoid ringing artifacts
	frequencies, PSD_estimate = median_mean_average_psd( high_passed_time_series[removed_points:-removed_points], fft_length, f_sampl )
	# Create frequency array to be used for interpolation
	frequencies_to_interpolate = rfftfreq( time_series_duration, d=delta_t )
	
	# Interpolate the estimated PSD to a higher resolution, in order to have a
	# spectrum of the same resolution of the transformed time series
	interpolated_PSD = np.interp( frequencies_to_interpolate, frequencies, PSD_estimate )
	
	# The whitening coefficients are given by the inverse of square root of the
	# PSD
	# The choice of the coefficients changes the units (and the normalization)
	# of the whitened time series:
	#   - if using sqrt(2*delta_f/PSD), then whitened time series is in units 
	#     of Hz
	
	#coefficients = np.sqrt((2.0*delta_f)/interpolated_PSD)	
	
	#   - if using sqrt(2/PSD), the whitened time series is in units of sqrt(Hz)
	coefficients = np.sqrt((2.0)/interpolated_PSD)
	
	

	
	
	# Multiply coefficients by the high pass filter coefficients
	# to highpass the time series if requested.
	if highpass:
		# Replace all bins with frequencies less than highpass cutoff with zeros
		x_tmp = np.arange(len(coefficients))
		coefficients = np.where(frequencies_to_interpolate < highpass_cutoff, 0.0, coefficients)
		
		"""
		# High pass filter frequency domain coefficients
		highpass_coeffs = high_pass_coefficients(frequencies_to_interpolate, highpass_cutoff, f_sampl, HIGH_PASS_ORDER)
		
		coefficients *= highpass_coeffs
		"""
	
	# The whitened time series is given by taking the inverse of the
	# transformed time series spectrum where each fequency bin is
	# weighted by the computed coefficients
	# Multiply by delta_f to obtain the correctly normalized result. 
	whitened_time_series = np.real( 2 * delta_f * np.fft.irfft(transform*coefficients) )
	
	# The returned whitened time series normalization and units depend on the 
	# choice of the whitened coefficients, see the above definition for the 
	# coefficients.
	return whitened_time_series


def butterworth_band_pass(time_series, order, f_high, f_low, f_sampl):
	'''
		Butterworth band pass filter.
		
		Arguments:
			- time_series (array)
				The time series to be filtered.
			- order (integer)
				Order of the filter.
			- f_high, f_low (float)
				High and low cutoff frequencies.
			- f_sampl (float)
				Sampling frequency.
		
		Output:
		 	filtered_time_series (array)
				The filtered time series.
		
	'''
	nyquist_frequency = f_sampl/2.0
	[b,a] = signal.butter(order, [frequency/(nyquist_frequency) for frequency in (f_high, f_low)], btype='band');
	filtered_time_series = signal.filtfilt(b, a, time_series)
	return filtered_time_series


def median_mean_average_psd(time_series, segment_length, f_sampl):
	''' 
		PSD estimation through the median-mean-average algorithm.
		
		The median mean average algorithm splits the input time series
		into overlapping sub-segments.
		These sub-segments are split into two sets, for even and odd indexes.
		
		After being windowed (hanning window), a PSD is computed for each
		of these sub-segments. The median of both (even and odd) sets is taken.
		The PSD estimate is the weighted average (corrected for a bias factor)
		of these two medians.
		
		The reason for the splitting in odd and even sets is to have the least possible
		correlation when estimating the PSD: taking the median of even and odd segments
		separately reduces the correlation (as even segments will not share data 
		with odd segments), in the following weighted average the correlations 
		won't matter as much. This is makes the algorithm less biased by loud glitches.
		
		The PSDs for each small subsegment is computed through matpotlib.mlab.psd
		
		More info about the median-mean-average algorithm can be found 
 		in the FINDCHIRP Paper (qc/0509116)
			
		Input:
			- time_series (array)
				Time series to compute the PSD on
			- segment_length (integer)
				Sets the frequency resolution of the PSD, for faster computation
				this should be a power of two (512, 1024, 2048, ...)
			- f_sampl (float)
				Sampling frequency
		Output:
			(freqs, PSD) (tuple)
				- freqs (array)
					Contains the frequencies at which PSD was computed.
				- PSD (array)
					Power amplitudes, one-sided (only positive frequencies are returned)
					
					If integrating over this PSD, the correct integrated value is:
						2 * PSD.sum()
			
		
	'''
	
	signal_length = len(time_series)
	even_psd_list = []
	odd_psd_list = []
	step = segment_length
	
	if ( 2*signal_length//step > 1):
		Ns =  2*signal_length//step - 1
	else:
		Ns = 1
	Ns_odd =  (Ns+1) // 2 
	Ns_even = (Ns-1) // 2 
	
	
	
	if ( Ns_odd % 2 == 0 ):
		Ns_odd -= 1
	elif ( Ns_even > 1 ):
		Ns_even -= 1
	
	current_bin = 0
	
	# Create window and compute window_normalization, which is used
	# to correctly normalize the periodogram
	window = np.hanning(step)
	window_normalization = np.sum(window**2)/float(step)
	
	# Create frequencies array to return at the end
	delta_t = 1./f_sampl
	delta_f = f_sampl/float(step)
	freqs = rfftfreq(step, d=delta_t)
	
	for i in range(Ns):
		ioff = i/2
		current_bin += step/2
		
		to_transform = time_series[current_bin:current_bin+step]
		
		# Exclude the last segment if its size is less than step
		if ( len(to_transform) < step):
			continue
		
		# Compute periodogram and apply correct normalization:
		# Since we used rfft, we have to multiply
		# the result by 2 because rfft returns a one-sided version 
		# (only positive frequencies) of the FFT and multiply by delta_f
		# to obtain the correct units.
		Pxx = np.abs( delta_t * np.fft.rfft(to_transform*window) )**2
		# Apply correct normalization (from FINDCHIRP paper)
		Pxx *= ( ( 2.0 * delta_f )/(window_normalization) ) 
		# Fix normalization in the DC bin
		Pxx[0] /= 2.0
		
		# Units are now counts^2 Hz^-1
		
		# Assign PSD to the correct bin	
		if ( i % 2 == 0 ):
			if ( ioff >= Ns_odd ):
				continue
			odd_psd_list.append( Pxx )
			
		else:
			if ( ioff >= Ns_even ):
				continue
			even_psd_list.append( Pxx )
			
	
	
	even_psd_median = np.median( np.array(even_psd_list), axis=0 )
	odd_psd_median = np.median( np.array(odd_psd_list), axis=0 )
	
	
	# If there's only one segment, then the PSD estimate is just the median of
	# that segment, otherwise it is the weighted average
	if ( Ns == 1 ):
		PSD_estimate = odd_psd_median
	else:
		even_weight = Ns_even/median_bias_factor(Ns_even)
		odd_weight =  Ns_odd/median_bias_factor(Ns_odd)
		PSD_estimate = even_psd_median*even_weight + odd_psd_median*odd_weight
		PSD_estimate /= float( Ns_even + Ns_odd)
	

	return freqs, PSD_estimate

def median_mean_average_energy(time_series, step):
	"""
		Does the same thing as the median-mean-average PSD, except this is done in the time
		domain to estimate the average energy content of a segment 'step' points long.
	"""
	
	signal_length = len(time_series)
	if ( 2*signal_length//step > 1):
		Ns =  2*signal_length//step - 1
	else:
		Ns = 1
	
	Ns_odd =  (Ns+1) // 2 
	Ns_even = (Ns-1) // 2 
	
	if ( Ns_odd % 2 == 0 ):
		Ns_odd -= 1
	elif ( Ns_even > 1 ):
		Ns_even -= 1
	
	energies_even = []
	energies_odd = []
	
	current_bin = 0
	for i in range(Ns):
		ioff = i/2
		current_bin += step/2
		
		current_segment = np.array(time_series[current_bin:current_bin+step])
		current_energy = (current_segment**2).sum()
		
		# Exclude the last segment if its size is less than step
		if ( current_segment.size < step):
			continue
			
		# Assign energy to the correct list
		if ( i % 2 == 0 ):
			if ( ioff >= Ns_odd ):
				continue
			energies_odd.append( current_energy )
		
		else:
			if ( ioff >= Ns_even ):
				continue
			energies_even.append( current_energy )
			
		
		
	even_energies_median = np.median( np.array(energies_even), axis=0 )
	odd_energies_median = np.median( np.array(energies_odd), axis=0 )
	
	# If there's only one segment, then the energy is just the median of
	# that segment, otherwise it is the weighted average
	if ( Ns == 1 ):
		energy_estimate = odd_energies_median
	else:
		even_weight = Ns_even/median_bias_factor(Ns_even)
		odd_weight =  Ns_odd/median_bias_factor(Ns_odd)
		energy_estimate = even_energies_median*even_weight + odd_energies_median*odd_weight
		energy_estimate /= float( Ns_even + Ns_odd)
	
	return energy_estimate

def compute_psd(time_series, resolution, f_sampl, overlap):
	'''
	Computes the PSD of time series using matplotlib.mlab.psd
	mlab.psd computes a PSD through Welch's method.
	
	Arguments:
		- time_series (array)
			Time series 
		- resolution (float)
		 	Frequency resolution of the power spectral density in Hz.
		- f_sampl (float)
			Sampling frequency for the given time series.
	 	- overlap (float)
			Overlap between segments when computing the PSD.
			overlap should be in the interval [0,1[
			0 is no overlap, 1 is 100% overlap
	
	Output:
		(freqs, Pxx) (tuple)
			- freqs
				Frequencies at which the PSD was computed
			- Pxx
				Power amplitudes
	
	
	'''
	NFFT = int(f_sampl/resolution)
	
	(PSD, freqs) = matplotlib.mlab.psd( time_series, Fs=f_sampl, NFFT=NFFT, noverlap=int(NFFT*overlap), scale_by_freq=True)
	
	return freqs, PSD


def median_bias_factor(n):
	'''
	 MEDIANBIASFACTOR - Compute bias factor for median estimate of mean.
		usage:
		alpha = median_bias_factor(n)
		n       Scalar.  Number of samples used to estimate median.  Must be a
				positive, odd integer.  
		alpha   Scalar.  Factor by which median must be divided to give an
				unbiased estimate of the mean, assuming an exponential
				distribution.    
	'''
	if ( n < 0 ):
		print "n must be an integer ( n is %i )" % n
		assert(False)
	alpha = 0.0
	toggle = 1.0
	for i in range( 1, n+1 ):
		alpha += toggle/float(i)
		toggle = - toggle
	return alpha




def main():
	global excluded_seconds
	args = check_options_and_args()
	
	
	# Create a progress bar
	PROGRESS = False
	if ( len(args) > 1 ):
		PROGRESS = True
	print "Processing..."
	if PROGRESS:
		bar = progressBar( minValue = 0, maxValue = len(args)-1, totalWidth = 60 )
		bar(0)
	
	
	# Data processing
	for index, element in enumerate(args):
		# Load data, using pickle if PICKLED=True
		if PICKLED:
			f = open(element, "rb")
			data = pickle.load(f)
			f.close()
		else:
			data = np.loadtxt(element)
		
		(conditioned_data, PSD) = whiten(data, SAMPLING_FREQUENCY, highpass=HIGH_PASS,\
		 							highpass_cutoff=HIGH_PASS_CUTOFF, resample=RESAMPLE)
		
		# Save data
		f = open(element + ".whitened", "wb")
		pickle.dump( conditioned_data, f)
		f.close()
		
		# Update the progress bar
		if PROGRESS:
			bar(index)
			
		
	print "\nDone!"



if __name__ == '__main__':
	main()

