#!/usr/bin/env python
# encoding: utf-8
""" PCAT.py

First version by Daniele Trifirò on 2013-08-17.
brethil@phy.olemiss.edu

Run with -h or --help for usage.

This is a wrapper the the other functions in the package.

PCAT Pipeline:
	1) Retrieve data 
		Read from the frame files (and save if requested) time series
		for the given time interval.
		
		Uses retrieve_timeseries() from utilities_PCAT.py
		
	2) Condition data 
		a conditioning_function() is defined based on the required
		processing, and data is prepared for analysis by calling
		conditioning_function().
		Uses various functions from data_conditioning.py
		
		Data conditioning and retriveval is sped up by launching multiple
		processes, defined by the PARALLEL_PROCESSES variable (see below).
		
	3) Create database
		Create a database from the previously processed data.
		
		Uses create_data_matrix(...), create_data_matrix_from_psd(...)
		from utilities_PCAT.py.
		
	4) Perform PCA
		Perform PCA on the database
		PCA(...) from PCA.py
		
	5) Cluster the scores matrix from PCA.py
		gaussian_mixture(...) from GMM.py
		
		Other clustering algorithms can be easily implemented
		by changing gaussian_mixture() with another function.
		For more see GMM.py
	
	6) Plot scatterplots (with image maps), time series
		scatterplot(...), spike_time_series(...), plot_psds(...)
		from GMM.py
	
	7) Print URL to analysis results.
		
"""

import sys
from time import asctime, localtime

# Number of parallel processes to run when conditioning data
global PARALLEL_PROCESSES
PARALLEL_PROCESSES = 4

# Maximum number of principal components scores shown in the triangle plot
# since the number of subplots in the triangle plot is n(n+1)/2, just having
# 15 results in 120 subplots in a single figure, which results in a segfault
# most of the time.
TRIANGLE_COMPONENTS = 12
# Max number of principal components plotted in the principal comonent summary 
# plot. (All the requested principal components will be plotted individually
# in the Principal_components subfolder )
MAX_ALL_PRINCIPAL_COMPONENTS = 20
####################################
# Other parameters
####################################
# Order of band-pass butterworth filter
BUTTERWORTH_ORDER = 4
# Number of max plotted principal components (scores scatterplots and waveforms)
MAX_PLOTTED_COMPONENTS = 5


################################################################################
# Default Values for parameters:
################################################################################

#####################
# Download options:
# Chunk size in seconds:
global segment_size
segment_size = 8
# chunk padding in seconds
download_overlap_seconds = 4
# Chunk padding in percent (1 = 100%)
download_overlap = 0.5

#####################

#####################
# PSD Options:
#####################
# Overlap between segments when computing the PSD
psd_overlap = 0.5

##############################
# GMM default parameters
################################
# Maximum number of clusters found by GMM
max_clusters = 10



# Load modules if help not requested or there are no arguments.
# Print help if no arguments or '-h'


from glob import glob

from utilities_PCAT import *

from download_frames import retrieve_timeseries
from data_conditioning import *

from finder import find_spikes
from PCA import PCA, create_data_matrix, eigensystem, matrix_whiten
from GMM import gaussian_mixture, scatterplot, color_clusters, spike_time_series, matched_filtering_test, correlation_test
from GMM import print_cluster_info, calculate_types, plot_psds, configure_subplot_time, configure_subplot_freq

def usage():
	'''
		Usage
	'''
	"""print "Usage:\n\tPCAT.py (--time || --frequency) (--start start_time --end end_time || --list times_list)\n\
	 --frame frame_type -I IFO -c channel_name [--size segment_size]  \n\
	 [--filter] [--low low_frequency --high high_frequency]\n\
	 [--energy] [--whiten] [--nohighpass] [-v variables_number]\n\
	 [--padding_seconds padding || --padding_percent padding_perc]\n\
	 [--save_timeseries] [--reconstruct]"
	"""
	print '\033[1m' + "Usage:" + '\033[0m'
	print "\tPCAT.py --IFO IFO --frame frame -c channel\\"
	print "--start start_time --end_end_time OR --list list_of_times\\"
	
	print "\t" + '\033[1m'+ "Time Domain options:" + '\033[0m'
	print "\t\t-t threshold [-v variable_n]"
	
	print '\033[1m'+ "\tFrequency Domain options:" + '\033[0m'
	print "\t\t[--size segment_size] [-v variable_n]"
	
	print "\nThe given time interval is processed and analyzed.\n"
	print "In time domain analysis the (whitened) time series is scanned"
	print "for transients, on which PCA is performed and the resulting"
	print "Principal Component Scores clustered with GMM."
	print "In frequency domain analysis PSDs for 'segment_size' seconds"
	print "of data are computed, PCA is performed and the resulting"
	print "Principal Component Scores are clustered with GMM.\n"
	
	print "#"*80
	
	
	print '\033[1m' + "Necessary Arguments:" + '\033[0m'
	
	print "  --time or --frequency"
	print "\tType of analysis being performed, either in time or frequency."
	print "\tTime analysis also requires the following arguments"
	print "\t'-t threshold' and '--whiten' or '--filter' ."
	print "\tFrequency does not require any extra arguments."
	
	print "  -I IFO, --IFO IFO"
	print "\tSpecify IFO, can be 'H' for Hanford or 'L' for Livingston."
	
	print "\t -c channel_name, --channel channel_name"
	print "\t\tChannel name, e.g.: L1:PSL-ISS_PDA_OUT_DQ."
	
	print "  --frame frame_type"
	print "\tFrame type, 'R' should be fine most times."
	print "\tWhen in doubt, check before running using ligo_data_find."
	
	print "\n" + '\033[1m' + "Interval to Analyze:" + '\033[0m'
	print " One of the following two:"
	print "   --start start_time --end end_time"
	print "\tStart AND end time expressed in GPS time, e.g. 1043107216"
	
	print "   --list times_list"
	print "\tA file containing a list of times to perform the analysis on."
	print "\tThe list should have two tab-separated columns, with"
	print "\tstart time on the first column and end time on the second"
	print "\tcolumn, in GPS time."
	
	print '\033[1m' + "Time Domain options:" + '\033[0m'
	print "   -t threshold, --threshold threshold"
	print "\tTrigger threshold in units of the standard deviation, e.g. 5.5."
	print "   --energy"
	print "\tNormalize identified transients to unit L2 norm"
	print "\tinstead of unit maximum amplitude (default)."
	
	
	print "\n" + '\033[1m' + "Data conditioning (time domain):" + '\033[0m'
	print "  --whitening whitens the data and applies an high-pass filter,"
	
	print "  --filter applies a band-pass filter. See below for more details."
	print "\n  *WARNING* --whitening downsamples the data to 4096Hz before"
	print "\tconditioning and analyzing. This means that the -v should be tuned"
	print "\taccordingly: e.g. choosing 1024 means that the number"
	print "\tof sampled seconds is 0.25s."
	print "\n\tThis can be changed through --resample new_frequency (or --noresample)."
	print "\n\tCondtioned files are saved in binary format, and can be used"
	print "\twith python/numpy using numpy's 'load()' function."
	
	print "   --filter "
	print "\tWith this option a fourth order butterworth"
	print "\tband-pass filter is applied."
	print "\tHigh and low cutoff frequencies have to be supplied"
	print "\tthrough --high and --low."
	
	print "   --high high_frequency, --low low_frequency\t (with --filter)"
	print "\tHigh and low cutoff frequency for band-pass filter (integers)"
	print "\t(Butterworth 4th Order)."
	
	print "   --whiten"
	print "\tWhiten and high-pass input data."
	print "\tWhitening is performed using the inverse PSD method,"
	print "\tthe PSD is computed using the median-mean-average algorithm."
	print "\tDefault overlap between neighbouring segments is 50%."
	
	print "   --highpasscutoff cutoff_frequency"
	print "\tCutoff frequency for the high pass filter (with --whiten)."
	
	print "   --nohighpass    (with --whiten)"
	print "\tDo not apply high-pass filter to input data."
	
	print "  "+"-"*int(74)
			
	print '\033[1m' + "Frequency domain options:" + '\033[0m'
	print "   --low low_frequency, --high high_frequency"
	print "\tOnly perform PCA on the frequency interval from "
	print "\tlow_frequency to high_frequency (integers) with a resolution"
	print "\tof 1 Hz."
	print "\t--low and --high, replace -v. The number of variables used"
	print "\tis given by high_frequency-low_frequency, use small value"
	print "\tfor the PCA algorithm to run faster."
	
	print '\033[1m' + "Optional arguments:" + '\033[0m'
	
	print "   --size segment_size\n\tSize in seconds of the chunks in which data is split."
	print "\tto be analyzed, must be larger than 1 (smaller is faster)."
	print "\tDefault is 8 sec for time analysis, 60 sec for frequency analysis."
	
	print "   --padding_seconds padding"
	print "\tPad the segment to be analyzed with extra seconds at the start and at the beginning,"
	print "\tin order to avoid filter artifacts."
	print "\tDefault 50% overlap in time domain analysis (e.g. 4 seconds overlap with"
	print "\t8 seconds long segments) and 0% overlap in frequency domain analysis."
		
	print "   --padding_percentage padding"
	print "\tSame as above, though expressed in percentage of segment_size"

	print "   --components components_number"
	print "\tNumber of components to be used when clustering in the"
	print "\tPrincipal Components space, integer, default is 40."
	print "\tThis also sets the number of principal components used"
	print "\tto reconstruct the glitches when plotting the time series"
	print "\tin the time-domain analysis when --reconstruct is used."
	print "\tThis might have to be tweaked for optimal results."
	print "\t --variance explained_variance"
	print "\tSimilar to the above, fixes the percentage of variance to include"
	print "\tin the analysis, hence fixing the number of principal components used."
	print "\tthis has to be a float between 0 and 1 (100%)."
	print "\tWhen --components is set to 0, explained variance is set to the default"
	print "\tvalue of 75%."
	
	print "   -m number, --maxclusters number"
	print "\tSpecifies the maximum number"
	print "\tof clusters. Integer, default is 10."
	
	print "   -v variables_n, --variables variables_n"
	print "\tIf performing time-domain analysis, this sets the time resolution"
	print "\ti.e. the temporal length of the sampled transients.\n"
	print "\tThe length of the sampled transients in seconds is"
	print "\t\tvariables_n/sampling_frequency"
	print "\twith sampling_frequency being 4096 Hz if --whiten"
	print "\tor the input value if --resample.\n"
	
	print "\tIf in frequency-domain, this sets frequency resolution of "
	print "\tthe computed spectrum: resolution equals to"
	print "\t\tnyquist_frequency/variables_n. "
	print "\t(where the Nyquist frequency is half of the sampling frequency."
	print "\tWhen performing frequency-domain analysis, for faster computation,"
	print "\tvariables_n should be a power of two."
	print "\tDefault values are 512 for time-domain and 2048 for frequency-domain."
	print "\tThe higher this value, the slower the PCA decomposition."
	
	print "   --reconstruct"
	print "\tIf set, then glitches' time series are reconstructed using components_number "
	print "\tprincipal components (set with --components or --variance)."
	print "\tOnly using the first few principal components will reduce noise in the time series"
	print "\tand make the 'true' shape of the glitch more clear."
	
	print "   --save_timeseries"
	print "\tSave raw time series of the interval being analyzed."
	print "\tThe saved files (binary) can be used with python/numpy,"
	print "\tusing numpy's 'load()'"
	
	print "   --noplot"
	print "\tDo not plot transients/PSDs (makes run faster)"
	
	print "   --silent"
	print "\tDo not display progress bars"
	
	print "  "+"-"*int(74) + "\n"





################################################################################

##############################
# Helper routines:           #
##############################
	
def check_options_and_args(argv):
	global components_number, psd_overlap, max_clusters, segment_size, download_overlap
	global LIST, CUSTOM_OUT, FILTER , WHITEN, HIGH_PASS, HIGH_PASS_CUTOFF, SAVE_TIMESERIES, NORESAMPLE, RESAMPLE
	global HIGH_PASS
	
	global variables, channel, IFO, sampling
	global frame_type, variables, normalization
	global low, high, output_name, threshold, time_resolution, components_number
	
	global ANALYSIS, ANALYSIS_FREQUENCY, CLEAN, RECONSTRUCT, NOPLOT, SILENT
	global AUTOCHOOSE_COMPONENTS, VARIANCE_PERCENTAGE
	
	LIST = False
	CUSTOM_OUT = False
	FILTER = False
	WHITEN = False
	SAVE_TIMESERIES = False
	RECONSTRUCT = False
	components_number = 40
	AUTOCHOOSE_COMPONENTS = False
	
	SILENT = False
	NOPLOT = False
	CLEAN = False
	NORESAMPLE = False
	# Boolean: apply High Pass filter

	HIGH_PASS = True
	
	# Decide if data is resampled to ANALYSIS_FREQUENCY (if not supplied through
	# --resample, then it's set to 4096.0 after the argument parsing, see below)
	RESAMPLE = True
	SILENT = False
	# Booleans 
	
	
	global start_time, end_time, times_list, times, total, download_overlap_seconds
	global glitchgram_start, glitchgram_end
	glitchgram_start, glitchgram_end = None, None
	
	low, high = None, None
	
	normalization = "amplitude"
	
	if len(argv[1:]) == 0:
		print "No arguments."
		usage()
		sys.exit(1)
	try:
		opts, args = getopt.getopt(argv[1:], "hc:I:v:t:m:", [ "help", 'threshold=', 'channel=', 'IFO=',\
														'variables=', 'frame=', 'start=',\
		 												'end=', 'padding_seconds=', "padding_percentage=", 'psd_overlap=',\
		 												'high=', 'low=', 'list=', 'components=', 'time', 'frequency',\
														'filter', 'maxclusters=', 'whiten', 'size=', 'nohighpass', 'resample=',\
														'highpasscutoff=', 'energy', "save_timeseries", "clean", "noresample",\
														"reconstruct", 'noplot', 'silent', "glitchgram_start=", "glitchgram_end=",\
														"variance="])
	except getopt.error, msg:
		print msg
		sys.exit(1)
	
	# option processing
	for option, value in opts:
		if option in ( "-h", "--help" ):
			usage()
			sys.exit(1)
		elif option in ( "-v", "--variables" ):
			variables = int(value)
		elif option in ( "-t", "--threshold" ):
			threshold = float(value)
		elif option in ( "--resample" ):
			ANALYSIS_FREQUENCY = float(value)
		elif option in ( "--noresample" ):
			NORESAMPLE = True
		elif option in ( "-c", "--channel" ):
			channel = value
		elif option in ( "--energy" ):
			normalization = "energy"
		elif option in ( "-I", "--IFO" ):
			IFO = value
		elif option in ( "--start" ):
			start_time = int(value)
		elif option in ( "--end" ):
			end_time = int(value)
		elif option in ( "--padding_percentage" ):
			download_overlap = float(value)/100.0
		elif option in ( "--padding_seconds" ):
			download_overlap_seconds = int(value)
		elif option in ( "--maxclusters", "-m" ):
			max_clusters = int(value)
		elif option in ( "--psd_overlap" ):
			psd_overlap = float(value/100.0)
		elif option in ( "--frame" ):
			frame_type = value
		elif option in ( "--components" ):
			components_number = int(value)
		elif ( option == "--high" ):
			high = int(value)
		elif ( option == "--low" ):
			low = int(value)
		elif ( option == "--size" ):
			segment_size = int(value)
		elif option in ( "--list" ):
			LIST = True
			times_list = value
		elif option in ( '--time' ):
			ANALYSIS = 'time'
		elif option in ( '--frequency' ):
			ANALYSIS = 'frequency'
		elif option in ( '--filter' ):
			FILTER = True
		elif option in ( '--whiten' ):
			WHITEN = True
		elif option in ( "--nohighpass" ):
			HIGH_PASS = False
		elif option in ( "--highpasscutoff" ):
			HIGH_PASS_CUTOFF = float(value)
		elif option in ( "--save_timeseries" ):
			SAVE_TIMESERIES = True
		elif option == "--clean":
			CLEAN = True
		elif option == "--reconstruct":
			RECONSTRUCT = True
		elif option == "--noplot":
			NOPLOT = True
		elif option == "--silent":
			SILENT = True	
		elif option == "--glitchgram_start":
			glitchgram_start = int(value)
		elif option == "--glitchgram_end":
			glitchgram_end = int(value)
		elif option == "--variance":
			AUTOCHOOSE_COMPONENTS = True
			VARIANCE_PERCENTAGE = float(value)
			components_number = 0
		else:
			print "Unknown option."
			sys.exit()
		
		
	######################################################
	# Checks arguments
	####################################################
	if not ( any( '--list' in o for o in opts) ):
		if ( any( '--start' in o for o in opts) and any('--end' in o for o in opts) ):
			# To add?: check validity of start time and end time
			if ( start_time >= end_time ):
				print "Check start and end time. Quitting."
				sys.exit()
	
	if not ( any( "--components" in o for o in opts)):
		if not (any( '--variance' in o for o in opts) ):
			print "Either --components (number of principal components used for clustering)"
			print "or --variance (percentage of explained variance used) have to be"
			print "supplied. Quitting."
			sys.exit()
		elif (any('--variance' in o for o in opts)):
			if (VARIANCE_PERCENTAGE <= 0) or (VARIANCE_PERCENTAGE>1.0):
				print "Variance percentage has to be in the ]0,1] interval."
				print "Quitting."
				sys.exit()
	elif ( any( "-components" in o for o in opts)):
		if (components_number > variables):
			print "The number of principal components used for clustering has to be lower"
			print "or equal than the number of variables. Quitting."
			sys.exit()
	if (components_number == 0) and not AUTOCHOOSE_COMPONENTS:
		# If the number of supplied components is zero, fall back to 75% 
		# explained variance.
		AUTOCHOOSE_COMPONENTS = True
		VARIANCE_PERCENTAGE = 0.75
				
	if not ( any( flag in o for flag in [ "--start", "--end", "--list"] for o in opts )):
		print "Start and end GPS time or a list of GPS times have to be supplied. Quitting."
		sys.exit()
	elif not any('--list' in o for o in opts):
		if not ( any( "--start" in o for o in opts) and any("--end" in o for o in opts) ):
			print "Both --start and --end have to be supplied. Quitting."
			sys.exit()
	elif ( any( flag in o for flag in ['--start', '--end'] for o in opts) and any('--list' in o for o in opts )):
			print "Choose one between --list or --start and --end. Quitting."
			sys.exit()	
	
	if not any( flag in o for flag in [ "--channel", "-c"] for o in opts ):
		print "Channel name has to be supplied through '--channel' or '-c'. Quitting."
		sys.exit()
	if not ( any( flag in o for flag in ['--IFO', "-I"] for o in opts ) ):
		print "IFO ('L', 'H' or 'H1_R') has to be supplied"
		print "through the --IFO or -I flags. Quitting."
		sys.exit()
		
	if not ( any( flag in o for flag in ['--time', "--frequency"] for o in opts ) ):
		print "Analysis type (time domain or frequency domain) has to be supplied"
		print "through the --time or --frequency flags. Quitting."
		sys.exit()
	
	if not any( flag in o for flag in ['--frame'] for o in opts ):
		print "Frame type has to be supplied (usually 'R' for Livingston or 'H1_R' for Hanford)."
		print "Check with ligo_data_find when in doubt. Quitting."
		sys.exit()
	
	if ( ANALYSIS == "time") and ( FILTER ):
		if not ( any( flag in o for flag in ['--high', '--low'] for o in opts ) ):
			print "Filter parameters (low and high cutoff frequencies) have to be supplied. Quitting."
			sys.exit()
	
	if ( ANALYSIS == "time") :
		if not ( any( flag in o for flag in ['-t', '--threshold'] for o in opts ) ):
			print "Threshold (in units of standard deviation) has to be supplied through the -t"
			print "or --threshold option. Quitting."
			sys.exit()
		if ( any( "--whiten" in o for o in opts) and any('--filter' in o for o in opts) ):
			print "Please choose between --whiten and --filter."
			print "Quitting."
			sys.exit()
	
	############################################################################
	# End of data checking
	############################################################################
	
	# Create a list with the time intervals to be analyzed
	# times = [ [start1, end1], [start2, end2], ... ]
	times = []
	if LIST:
		f = open( times_list, "r")
		tmp = ( f.read().strip() ).split("\n")
		f.close()
		for element in tmp:
			times.append( [int(a) for a in element.split()] )
	else:
		times.append( [start_time, end_time] )
	
	
	# Retrieve sampling frequency for the given channel by retrieving a small 
	# portion of the data. Sampling frequency is accessed through the 'fs' key
	# of the 'data' dictionary (see download_frames.py for the definition
	# of retrieve_timeseries)
	start = times[0][0]
	try:
		data = retrieve_timeseries(start, start+32, channel, IFO, frame_type)
	except RuntimeError:
		assert False, "Error retrieving data. Check channel name, frame type and IFO."
	
	sampling = data['fs']
	if NORESAMPLE:
		ANALYSIS_FREQUENCY = sampling
	############################################################################
	
	
	
	# This part deals with the "bands" analysis of PCAT and is used
	# to have the correct definitions when running PCAT on fequency bands
	# Append "_bands" to ANALYSIS.
	if ( ( ANALYSIS == "frequency" ) and any( flag in o for flag in ['--high', '--low'] for o in opts ) ):
		if ( any("--high" in o for o in opts) and any("--low" in o for o in opts) ):
			ANALYSIS = "frequency_bands"
		else:
			print "If bands analysis is being performed, both low and high frequencies\n\
	have to be supplied through --low and --high. Quitting."
			sys.exit()
	
	# Compute the overlap in seconds between neighouring segments
	# if the download overlap has not been supplied through the
	# --download_overlap_seconds option.
	if ( "--padding_percentage" in argv ):
		download_overlap_seconds = int(segment_size*download_overlap)
	elif ( "--padding_seconds" in argv):
		pass
	elif ( "--size" in argv ):
		download_overlap_seconds = int(segment_size*download_overlap)
		
	# If variable input has not been supplied, set them to default values:
	if not ( any( flag in o for flag in [ '--variables', "-v" ] for o in opts ) ):
		if ( ANALYSIS == "time" ):
			variables = 512
		elif ( ANALYSIS == "frequency"):
			variables = 2048
	
	if ( ANALYSIS == "frequency_bands" ):
		# If band analysis, compute the spectrum at 1Hz resolution, then 
		# only analyze the the needed band
		# Variables number (given resolution is) [(f_s/2)/res]+1
		variables = (sampling//2)//1 +1
	
	# Set ANALYSIS_FREQUENCY to default (4096.0) if no other options have been
	# supplied
	if ( "--resample" not in argv) and not NORESAMPLE:
		ANALYSIS_FREQUENCY = 4096.0
		
	if ( "frequency" in ANALYSIS):
		if not ( any("--size" in o for o in opts)):
			segment_size = 60
		if ( any("--padd" not in o for o in opts) ):
			download_overlap = 0.0
			download_overlap_seconds = 0
			
	# Count the total number of analyzed seconds
	total = 0
	for element in times:
		total += int(element[1])-int(element[0])
		
	# Set default value for find_spikes()'s time_resolution
	time_resolution = int(variables/2)
	



def get_server_url():
	"""
		This retrieves the hostname on the server this program is being run on
		in order to give an URL for the output of PCAT.
	"""
	hostname = getstatusoutput( "echo $(hostname)| sed \"s/-[^.]*/-jobs/\"" )[1]
	if ( ".edu" ) in hostname:
		server = "http://"+hostname+"/"
	else:
		server = "http://"+getstatusoutput("echo $(hostname) | sed \"s/-[^.]*/-jobs/\" | xargs -i echo {}.$(dnsdomainname)")[1]+"/"
	
	# The following line is used to give the correct URL when running PCAT 
	# through condor
	# If the string "node" is in the URL, then remove the "node*" part and
	# substitute it with "ldas-jobs"
	if ("node") in server:
		tmp = join(server.split("/")[3:], ".")
		if ( "cit" in server ):
			server = "http://ldas-jobs.ligo.caltech.edu/" + tmp
		elif ("ligo-la" in server):
			server = "http://ldas-jobs.ligo-la.caltech.edu/" + tmp
		elif ("ligo-wa" in server):
			server = "http://ldas-jobs.ligo-wa.caltech.edu/" + tmp
		else:
			server = "Error parsing url: \t"  + tmp
	return server


def print_parameters():
	'''
	This function prints the parameters of the run.
	'''
	global times, total, units, variables
	global ANALYSIS_FREQUENCY
	print "#"*int(0.8*frame_width)
	
	print "Parameters for this run:"
	
	
	if ( "time" in ANALYSIS ):
		print "Time Domain Analysis:"
	elif ( "frequency" in ANALYSIS  ):
		print "Frequency Domain Analysis:"
	print "\t - General:\t"
	if LIST:
		print "\t\t List of times:\t\t\t", times_list 
	else:
		print "\t\t Start time:\t\t\t", start_time, "("+str(getstatusoutput("tconvert "+str(start_time) )[1]) + ")\n\
		 End time:\t\t\t", end_time, "("+str(getstatusoutput("tconvert "+str(end_time) )[1]) + ")"
	units = "%i seconds"
	global dividend
	dividend = 1.0
	if ( total > 60.0 ):
		dividend *= 60.0
		units = "minutes" 
		if ( total/dividend > 60.0 ):
			dividend *= 60.0
			units = "hours" 
			if ( total/dividend > 24.0 ):
				dividend *= 24.0
				units = "days"
	print "\t\t Total analyzed time:\t\t{0:.2f} {1}".format(total/dividend, units)
	print "\t\t Channel name:\t\t\t", channel
	print "\t\t IFO:\t\t\t\t", IFO
	print "\t\t Data Type:\t\t\t", frame_type
	print "\t\t Chunk size:\t\t\t", segment_size, "seconds"
	print "\t\t Padding:\t\t\t", download_overlap_seconds, "seconds"
	print "\t\t (Total size: ", segment_size+2*download_overlap_seconds, "seconds)"
	print "\t\t Sampling frequency:\t\t", sampling
	if (RESAMPLE and ( "time" in ANALYSIS) and ( sampling > ANALYSIS_FREQUENCY) ):
		print "\t\t (Downsampled to %.1f for analysis)" % ANALYSIS_FREQUENCY
	print ""
	if ( "time" in ANALYSIS ):
		if FILTER:
			print "\t - Butterworth filter:"
			print "\t\t Low Cutoff Frequency:\t\t", low
			print "\t\t High Cutoff Frequency:\t\t", high
			print 
		if WHITEN:
			print "\t\t Whitening:\t\t\tON"
			if HIGH_PASS:
				print "\t\t High Pass Cutoff:\t\t", HIGH_PASS_CUTOFF
			else:
				pass
		else:
			print "\t\t Whitening:\t\t\tOFF\n"
			if HIGH_PASS:
				print "\t\t High Pass Cutoff:\t\t", HIGH_PASS_CUTOFF
		print "\t - Trigger parameters:\n\
		 Number of variables:\t\t", variables
		if ( sampling > ANALYSIS_FREQUENCY ):
			seconds_per_trigger = float(variables)/ANALYSIS_FREQUENCY
		else:
			seconds_per_trigger = float(variables)/sampling
		print "\t\t (%.3f s per trigger)" % seconds_per_trigger
		print "\t\t Threshold:\t\t\t", threshold
		if RECONSTRUCT and not AUTOCHOOSE_COMPONENTS:
			print "\t\t Identified transients plots:\tReconstructed using the"
			print "\t\t\t\t\t\tfirst {0} principal components".format(components_number)
		elif AUTOCHOOSE_COMPONENTS and RECONSTRUCT:
			print "\t\t Identified transients plots:\tReconstructed using the"
			print "\t\t\t\t\t\tusing principal component accounting"
			print "\t\t\t\t\t\tfor {0:.1%} of the variance".format(VARIANCE_PERCENTAGE)
		else:
			print "\t\t Identified transients plots:\t Raw time series"
	elif ( "frequency" in ANALYSIS ):
		print "\t- Frequency Domain:\n"
		# The 'resolution' variable is used when computing the PSD, and is not 
		# the actual resolution of the PSD, but yields the correct result
		# (I know, this is complicated)
		global resolution
		resolution = sampling/(2*(variables-1))
		if ( "bands" in ANALYSIS ):
			# In bands analysis the spectrum is computed with a 1Hz resolution.
			print "\t\tPSD resolution:\t\t\t%.2f Hz" % (resolution)
			print "\t\tBand:\t\t\t\t%i to %i Hz (%i points)" % (low, high, int(high-low+1))
		else:
			print "\t\tPSD resolution:\t\t\t%.2f Hz (%i points)" % ( resolution, variables )
	print 
	print "\t - PCA and GMM:"
	if AUTOCHOOSE_COMPONENTS:
		print "\t\tPercentage of total"
		print "\t\tvariance explained:\t\t{0:.2%}".format(VARIANCE_PERCENTAGE)
	else:
		print "\t\tComponents Number:\t\t", components_number
	print "\t\t Maximum Clusters:\t\t", max_clusters
	
	global CLEAN
	if CLEAN:
		print "\tRemoving existing data folders before running pipeline."
	
	if SILENT:
		print "\tSILENT run (no progress bars)"
	print "\n\t- Log file:\n\t  {0}\n".format(log_name)
	
	print "#"*int(0.8*frame_width)
	
	# end of print_parameters()



def pipeline(args):
	global sampling
	check_options_and_args(args)
	
	###########################################################################
	# General setup 			 											  #
	###########################################################################
	
	# Server definition
	server = get_server_url()
	username = os.path.basename(os.path.expanduser("~"))
	original_wd = os.path.abspath(".")
	
	##################################################################
	# Working directory is ~/PCAT/Data
	# Processing directory is ~/PCAT/Data/channel_name/interval/
	##################################################################
	# If a list of times is supplied, create a folder with the same name
	# as the list. If start and end time are supplied then the folder name
	# is 'start_time-end_time'.
	# All data for the run is saved to this folder.
	# The final database and plots are saved in the public_html path.
	# (See further down)
	
	
	# Set up data directory
	data_directory = os.path.expanduser("~/PCAT/Data/")
	
	# Set up processing directory	
	if LIST:
		# Extract the name from the path
		# This is probably a stupid way to do this,
		# but it seems to work
		if ( "/" in times_list):
			name = times_list.split("/")[-1]
			if ( name == "/" ):
				name = times_list.split("/")[-2]
		else:
			name = times_list
		processing_directory = data_directory + channel + "/" + name + "/"
	else:
		global start_time
		global end_time
		processing_directory = data_directory + channel + "/" + "%i-%i/" % ( start_time, end_time ) 
	
	
	# Clean processing directories if requested
	if (CLEAN and os.path.exists(processing_directory)):
		from shutil import rmtree
		try:
			rmtree(processing_directory)
		except:
			print "Error cleaning directory:"
			print "\t" + processing_directory
	
	# Set up time series download directory (if requested)
	if SAVE_TIMESERIES:
		download_directory = processing_directory + "raw/"
	
	# Set up log file for error output and/or debug
	global log_name 
	log_name = processing_directory + "Analysis-{0}".format(ANALYSIS)
	
		
	if ("time" in ANALYSIS ):
		if WHITEN:
			log_name += "_whitening"
			if HIGH_PASS:
				log_name += "_highpassed_" + str(HIGH_PASS_CUTOFF)
		elif FILTER:
			log_name += "_filter_%i-%i" % (low, high)
	else:
		global resolution
		resolution = sampling/(2*(variables-1))
		log_name += "-{0:.2f}Hz".format(resolution)
	
	log_name += ".log"
	
	
	# Output folders setup
	start_date, end_date = getstatusoutput( "tconvert " + str(np.min(times)) )[1], getstatusoutput( "tconvert " + str(np.max(times)) )[1]
	if ( ANALYSIS == "time" ):
		OUTPUT = "time_PCAT/"
		OUTPUT += channel + "/"
		
		# Set up database name and output path
		database_name = "t-%.1f_w-%i.list" % ( threshold, variables) 
		
		if LIST:
			if ( "/" in times_list):
				name = times_list.split("/")[-1]
				if ( name == "/" ):
					name = times_list.split("/")[-2]
			else:
				name = times_list
			if ( "." in name ):
				name = name.split(".")[0]
			OUTPUT += name + "/" 
		else:
			OUTPUT += (start_date + "_" + end_date).replace(" ", "_") + "/"
			
		if FILTER:
			OUTPUT += "{0}_{1}".format(low, high)
		elif WHITEN:
			OUTPUT += "whitened_"
			if HIGH_PASS:
				OUTPUT += "highpassed_{0}_".format(HIGH_PASS_CUTOFF)
		OUTPUT += "{0}_".format(normalization)
		OUTPUT += "t-{0}_w-{1}/".format(threshold, variables)
	elif ( "frequency" in ANALYSIS ):
		OUTPUT = "frequency_PCAT/"
		OUTPUT += channel + "/"
		
		if ( "bands" in ANALYSIS):
			database_name = "PSD_database-band_%i-%i-1Hz.list" % ( low, high )
		else:
			database_name = "PSD_database-%.1f_Hz.list" % resolution
			
		if LIST:
			if ( "/" in times_list):
				name = times_list.split("/")[-1]
				if ( name == "/" ):
					name = times_list.split("/")[-2]
			else:
				name = times_list
			OUTPUT += name + "/" 
		else:
			OUTPUT += (start_date + "_" + end_date + "/").replace(" ", "_")
		
		if ( "bands" in ANALYSIS):
			OUTPUT += "bands_%i-%i/" % (low, high)
		else:	
			OUTPUT += "%.1f-Hz_resolution/" % resolution
		
	else:
		assert False, "DIFF ANALYSIS NOT IMPLEMENTED"
		
	
	
	# Create output folders
	output_dir = os.path.expanduser("~/public_html/" + OUTPUT)
	try:
		os.makedirs( output_dir )
	except:
		pass
	
	# Make sure 'output_dir' is empty, if not, delete all contents
	folder_contents = os.listdir(output_dir)
	if not ( folder_contents == [] ):
		from shutil import rmtree
		for element in folder_contents:
			try:
				rmtree(output_dir + element)
			except:
				os.unlink(output_dir + element)
	
	# Create results folder URL address
	results_URL = "{0}~{1}/{2}".format(server, username, OUTPUT)
	# Print parameters parsed by check_options_and_args()
	print_parameters()
	
	print "Saved URL:\n  {0}\n".format(results_URL + "parameters.txt")
	
	# Open new file to save the parameters for the run in 
	f = open(output_dir + "parameters.txt", "w")
	# Set stdout to the opened file
	sys.stdout = f
	# Print parameters: these will be written to the open file 'f'
	# due to 'sys.stdout = f'
	print_parameters()
	# 'args' is just a list, use join to join the list into a string
	# separated by a single space (s.join(iterable) returns a string
	# with the elements in args separated by the content of 's')
	original_command = " ".join(args)
	f.write("\n\n Original command:\n{0}".format(original_command))
	f.close()
	# Return to the original stout:
	sys.stdout = sys.__stdout__	
	#print "\tParameters for the run saved: 'parameters.txt'"
	

	
	# Create working folders (if not already present)
	try:
		os.makedirs( processing_directory )
	except:
		pass
	try:
		os.makedirs( download_directory )
	except:
		pass
	
	# Create log file	
	log = open(log_name, "w")
	log.write("Log start:\t " + asctime(localtime()) + "\n")
	log.write("-"*30)
	log.write("\n")
	if SAVE_TIMESERIES:
		log.write("Saving timeseries.\n")
	log.close()
	
	#TODO: it would probably be better to have a function which does all of the setup done above, instead of having it over 80+lines
	
		
		
	##################################################################
	# Data processing/conditioning definitions                       #
	##################################################################
	# When doing time-domain analysis, time series are first retrieved 
	# and processed and then scanned for transients.
	# Found transients (Spike() objects) are stored in a python list.
	# If frequency-domain analysis, time series are again retrieved and
	# processed (a PSD is computed).
	# The computed PSDs are again stored in a python list as numpy
	# arrays.
	
	# In both cases the python list is converted to a matrix on which
	# PCA is performed and then the GMM algorithm is used to cluster
	# the observations (transients/spectra).
	###################################################
	
	
	# First define conditioning_function(time_series), which,
	# depending on the kind of analysis requested
	# either band-passes and/or whitens or computes
	# the PSD of the supplied time interval.
	
	# The files returned by retrieve_timeseries are python dictionaries 
	# containing the time series with keys 'waveform', 'dt', and 'fs'
	# (see download_frames.py)
	
	if ( "time" in ANALYSIS ):
		# Define band-pass or whitening filter)
		if FILTER:
			conditioned_folder = str(low) + "_" + str(high)+"/"
			conditioning_function = lambda x: butterworth_band_pass(x['waveform'], order=BUTTERWORTH_ORDER,\
			 							f_high=high, f_low=low, f_sampl=x['fs'])
		elif WHITEN:
			if HIGH_PASS:
				conditioned_folder = "high_passed_%1.f-whitened/" % HIGH_PASS_CUTOFF
			else:
				conditioned_folder = "whitened/"
			conditioning_function = lambda x: whiten(x['waveform'], download_overlap_seconds, f_sampl=x['fs'], resample_freq=ANALYSIS_FREQUENCY, \
									highpass=HIGH_PASS, highpass_cutoff=HIGH_PASS_CUTOFF,\
									resample=RESAMPLE)
		else:
			# If no kind of processing is required return the unprocessed time
			# series.
			conditioned_folder = "raw/"
			if HIGH_PASS:
				conditioning_function = lambda x: high_pass_filter(x['waveform'], HIGH_PASS_CUTOFF, x['fs'], 4)
			else:
				conditioning_function = lambda x: x['waveform']
	elif ( "frequency" in ANALYSIS ):
		# The conditioning function simply has to compute the PSD of the input
		# segment
		conditioned_folder = "PSDs-%1.f_Hz/" % resolution
		def conditioning_function(x):
			global resolution
			freqs, PSD = compute_psd(x['waveform'], resolution=resolution, f_sampl=x['fs'], overlap=psd_overlap)
			return PSD
	try:
		os.makedirs( processing_directory + conditioned_folder )
	except:
		pass
	
	# Define the worker function, depending on the type of analyis
	# which was requested.
	# The worker function is called through parmap() (defined in 
	# utilities_PCAT.py) to speed up computation time using multiple processes.
	# - If doing a time domain analysis, the time series is analyzed as soon as
	#   it is retrieved from the frame files. The found transients are stored
	#   as Spike() objects in a python list.
	# - If doing a frequency domain analysis, the PSDs are computed from the
	#   time series as soon as they are retrieved from the frame files. They
	#   then are saved to files, which are later used to create a database for 
	#   the analysis.
	
	# Workfunction() takes in both case a tuple ('arguments') with the start
	# and end time expressed in GPStime.
	
	if ( "time" in ANALYSIS):
		# The workfunction finds transients in the given time series
		# and returns a list of spike objects. These are later used
		# for analysis
		def workfunction(arguments):
			""" arguments is a tuple, unpack it to use """
			conditioning_function, start, end = arguments
			
			out_name = "{0}-{1}-{2}_{3}-{4}.data".format(IFO, frame_type, channel,
																	start, end)
			
			out_file = processing_directory + conditioned_folder + out_name + ".conditioned"
			
			# Check if conditioned file already exists.
			# If it exists, load it, else use conditioning_function() on 
			# time_series
			downsample_factor = int(sampling/ANALYSIS_FREQUENCY)
			if os.path.isfile(out_file):
				try:	# Try loading files
					log = open(log_name, "a")
					f = open(out_file, "rb")
					conditioned = np.load(f)
					# Print part of the path to the log file
					log.write( "FILE EXISTS - LOADED:\t'./{0}'\n".format(join(out_file.split("/")[-5:], "/")))
					f.close()
					log.close()
				except Exception, error: # Try retrieving the time series and
										 # conditioning the file once again
					# If there is an exception, try loading
					# the time series from the frame files
					# and conditioning it
					try:
						time_series = retrieve_timeseries(start, end, channel, IFO, frame_type)
					except Exception, error:
						log = open(log_name, "a")
						log.write("!!! ERROR OPENING SEGMENT {0}-{1}:\t{2}\n".format(start, end, str(error)))
						log.close()
						return []
					conditioned = conditioning_function(time_series)
					
				
			else:	# First time processing data/conditioned file does not 
					# exist.
					# Retrieve the time series and condition/save it.
				try:
					time_series = retrieve_timeseries(start, end, channel, IFO, frame_type)
				except Exception, error:
					log = open(log_name, "a")
					log.write("!!! ERROR OPENING SEGMENT {0}-{1}:\t{2}\n".format(start, end, str(error)))
					log.close()
					return []
				
				if SAVE_TIMESERIES:
					try:
						f = open(download_directory + out_name, "wb")
						np.save(f, time_series['waveform'])
						f.close()
					except Exception, error:
						log = open(log_name, "a")
						log.write("ERROR SAVING:\t{0}\n".format(out_name))
						log.close()
				
				# Condition the retrieved time series
				conditioned = conditioning_function(time_series)
				del time_series
				
				with open(out_file, "wb") as f:
					np.save(f, conditioned)
			
			# Search the conditioned time series for transients
			if (WHITEN and RESAMPLE and (sampling > ANALYSIS_FREQUENCY)):
				found_spikes = find_spikes(conditioned, out_name, threshold, variables,
											time_resolution=time_resolution,
											removed_seconds=download_overlap_seconds, 
											f_sampl=ANALYSIS_FREQUENCY,
											normalization=normalization)
				
			else:
				found_spikes = find_spikes(conditioned, out_name, threshold, variables,
											time_resolution=time_resolution,
											removed_seconds=download_overlap_seconds, 
											f_sampl=sampling,
											normalization=normalization)
			del conditioned
			# Update progress bar
			if not SILENT:
				global progress, progress_index
				progress_index += 1
				print "\r\t", sprog.next(), "\r",
				progress(progress_index)
			return found_spikes
	elif ( "frequency" in ANALYSIS):
		# In this case the workfunction gets the timeseries, computes the PSD,
		# saves it to a file and returns the path to the saved file.
		# The paths to these files are stored in PSDs and then used to create
		# a database for analysis.
		def workfunction(arguments):
			conditioning_function, start, end = arguments
			out_name = "{0}-{1}-{2}_{3}-{4}.data".format(IFO, frame_type, channel,
																	start, end)
			out_file = processing_directory + conditioned_folder + out_name + ".conditioned"
			# Check if conditioned file already exists.
			# If it exists, load it, else condition
			if os.path.isfile(out_file):
				try:
					f = open(out_file, "rb")
					PSD = np.load(f)
					log = open(log_name, "a")
					log.write( "FILE EXISTS - LOADED:\t'~/{0}'\n".format(join(out_file.split("/")[-5:], "/")))
					log.close()
					f.close()
				except:
					try:
						time_series = retrieve_timeseries(start, end, channel, IFO, frame_type)
					except Exception, error:
						log = open(log_name, "a")
						log.write("!!! ERROR OPENING SEGMENT {0}-{1}:\t{2}\n".format(start, end, str(error)))
						log.close()
					PSD = conditioning_function(time_series)
					del time_series
			else:   
				try:
					time_series = retrieve_timeseries(start, end, channel, IFO, frame_type)
				except Exception, error:
					log = open(log_name, "a")
					log.write("!!! ERROR OPENING SEGMENT {0}-{1}:\t{2}\n".format(start, end, str(error)))
					log.close()
					return []
				if SAVE_TIMESERIES:
					try:
						f = open(download_directory + out_name)
						np.save(f, time_series['waveform'])
						f.close()
					except:
						log = open(log_name, "a")
						log.write("ERROR SAVING:\t'{0}'\n".format(out_name))
						log.close() 
				
				PSD = conditioning_function(time_series)
				del time_series
				try:	# Save the PSD, for later use in database
					f = open(out_file, "w")
					np.save(f, PSD)
					f.close()
				except:
					log = open(log_name, "a")
					log.write( "ERROR SAVING:\t'{0}'\n".format(out_file))
					log.close()
			# Update progress bar
			if not SILENT:
				global progress, progress_index
				progress_index += 1
				print "\r\t", sprog.next(), "\r",
				progress(progress_index)
			return out_file
		
	
	
	
	##########################################################################
	# Actual data retrieval/processing conditioning                          #
	##########################################################################
	###
	# Prepare the segments, store the segment start and end time into the 
	# 'segments' list. Data is processed by looping overt this list.
	###
	
	# Split the given time interval(s) into shorter length overlapping 
	# subsegments, store the start and end time of the subsegments
	# in segments as a list of tuples.
	segments = []
	for index, interval in enumerate(times):
		start, end = interval[0], interval[1]
		new_start = int(start) - download_overlap_seconds
		new_end = new_start + segment_size + 2*download_overlap_seconds
		segments.append((new_start, new_end))
		while (new_end < int(end) + 2*download_overlap_seconds):
			new_start = new_end - 2*download_overlap_seconds
			new_end = new_start + segment_size + 2*download_overlap_seconds
			segments.append((new_start, new_end))
	
	# Create directories for conditioned files
	os.chdir(processing_directory)
	try:
		os.makedirs(folder + processing_directory + channel)
		if SAVE_TIMESERIES:
			os.makedirs(download_directory + processing_directory + channel + "/" + conditioned_folder)
	except:
		pass
	try:
		os.chdir(processing_directory)
	except:
		pass
	
	####
	# Actual data processing
	###
	# Call workfunction() on segment's elements, parallelizing using 
	# PARALLEL_PROCESSES processes (this number is defined at the top of this 
	# file).
	# parmap() is defined in utilities_PCAT.py and uses the multiprocessing
	# module.
	if not SILENT:
		global progress, progress_index
		progress_index = 0
		bar_len = (len(segments)-1)//PARALLEL_PROCESSES
		progress = progressBar(minValue = 0, maxValue=bar_len, totalWidth = 40 )
	
	print "\tDownloading and processing..."
	if not SILENT:
		progress(0)
		print "\t", sprog.next(), "\r",
		sys.stdout.flush()
		print " "*(frame_width-3), "\r",
	
	
	"""
	# Old parallel code, better to avoid using this as it suppresses all
	# errors when one of its workers fail
	worker = lambda x: workfunction((conditioning_function, x[0], x[1]))
	tmp_result = parmap(worker, segments, nprocs=PARALLEL_PROCESSES)
	"""
	# Worker function for parallel processing
	def worker(in_list, out_q):
		out_arr = []
		for segment in in_list:
			out_arr.append(workfunction((conditioning_function, segment[0], segment[1])))
		out_q.put(out_arr)
	
	# Each process will get 'chunksize' segments and a queue to put his out
	# dict into
	out_q = multiprocessing.Queue()
	chunksize = int(np.ceil(len(segments)/float(PARALLEL_PROCESSES)))
	procs = []
	
	for i in range(PARALLEL_PROCESSES):
		p = multiprocessing.Process(
				target=worker,
				args=(segments[chunksize * i:chunksize * (i + 1)],
						out_q))
		procs.append(p)
		p.start()
	
	# Collect all results into a single array
	tmp_result = []
	for i in range(PARALLEL_PROCESSES):
		tmp_result.extend(out_q.get())
		
	# Wait for all worker processes to finish
	for p in procs:
		p.join()
		
		
	# Complete progress bar
	if not SILENT:
		progress(bar_len)
		print "\t[ O ]"
	
	if ("time" in ANALYSIS):
		results = []
		for item in tmp_result:
			results.extend(item)
		print "\tFound {0} transients.".format(len(results))
		del tmp_result
	else:
		results = tmp_result
	

	"""
	# The following block is used for sequential (non-parallel) processing
	results = []
	for segment in segments:
		tmp_results.extend(worker(segment))
	
	"""

	
	# Data retrieval and processing is complete.
	# Output from the pipeline, either a list of Spike
	# objects or a list of paths to files, is stored in 'result'
	# for later use.
	
	log = open(log_name, "a")
	log.write("-"*30)
	log.write("\nLog end ("+ str(asctime(localtime())) +  ")\n")
	log.close()
	
	log = open(log_name, "r")
	txt = log.read()
	errors = txt.count("ERROR OPENING SEGMENT")
	if ( errors != 0 ):
		print "\tWarning: {0} segment{1} not processed, check log file for more information.".format(errors, "s" if errors > 1 else "")
		#print "Log file:\t {0}".format(log_name.split("\n")[:-3])
	log.close()	
	###########################################################################
	# End of data conditioning.                                               #
	###########################################################################
		
	# If there are 0 transients, write a warning to the output directory and 
	# return the URL to the warning
	if ( len(results) == 0 ):
		with open(output_dir + "no_transients.txt", "w") as f:
			print>>f, "\n No transients found. Quitting."
		return results_URL + "no_transients.txt"
	
	###########################################################################
	# Create database with conditioned data:
	# Search for glitches through the trigger finder if
	# time-domain analysis is requested, otherwise
	# create a PSD database for frequency-domain analysis.
	###########################################################################
	
	# Create data_matrix. PCA will be perfomed on this matrix.
	if ( ANALYSIS == "time" ):
		# Create data_matrix
		data_matrix = create_data_matrix(results, ANALYSIS)
		data_list = results
	elif ( "frequency" in ANALYSIS ):
		data_list, data_matrix = create_data_matrix_from_psds(results, ANALYSIS, sampling, low, high)
	else:
		assert False, "DIFF ANALYSIS NOT IMPLEMENTED"
		
	# Change working directory to 'output_dir'
	os.chdir( output_dir )
		
	print "\tDone!"
	print "#"*int(0.8*frame_width)
	
	
	############################################################
	# PCA and GMM											   #
	############################################################
	# If resampling time series, set sampling frequency
	# to ANALYSIS_FREQUENCY (to which time series have been downsampled)
	if (RESAMPLE and (sampling > ANALYSIS_FREQUENCY) ) and ("time" in ANALYSIS) :
		sampling = ANALYSIS_FREQUENCY
	print "\tPerforming PCA and clustering..."
	observations, samples = data_matrix.shape
	print "\t%ix%i data matrix: %i observations of %i variables " % ( observations, samples, observations, samples )
	
	# Get score matrix and principal components, PCA()
	# Columns means and standard deviations are stored in means
	# stds should be a numpy array of ones, unless matrix_whiten(..., std=True)
	# in PCA()
	global components_number
	if AUTOCHOOSE_COMPONENTS:
		score_matrix, principal_components, means, stds, eigenvalues = PCA(data_matrix, components_number=components_number, variance=VARIANCE_PERCENTAGE )
	else:
		score_matrix, principal_components, means, stds, eigenvalues = PCA(data_matrix, components_number=components_number)
	# Save pickled Principal Components 
	f = open("Principal_components.dat", "wb")
	pickle.dump(principal_components, f)
	f.close()
	print "\tSaved 'Principal_components.dat'"
	
	explained_variance = np.cumsum(np.abs(eigenvalues))/np.sum(np.abs(eigenvalues))
	
	
	# If automatically chosing components update parameters dump to include
	# number of PCs used when clustering	
	if AUTOCHOOSE_COMPONENTS:
		components_number = max(1,np.argmax(np.where(explained_variance<VARIANCE_PERCENTAGE, explained_variance,0)))
		with open("parameters.txt", "r") as f:
			text = f.read()
			text = text.replace("PCA and GMM:", "PCA and GMM:\n\t\tPrincipal components used:\t{0}".format(components_number))
			text = text.replace("of the variance", "of the variance ({0} principal components).".format(components_number))
			
		with open("parameters.txt", "w") as f:
			f.write(text)
		
	print "\n\tClustering using the first {0} principal components, accounting for {1:.1%} of the total variance".format(components_number, explained_variance[components_number])
	print "\t  Maximum clusters: {0}".format(max_clusters)
	
	# Cluster the score matrix using gaussian_mixture().
	# Only use the first 'components_number' components, and standardize input 
	# data (set zero mean and unit variance) using matrix_whiten to improve
	# clustering.
	
	reduced_whitened_scores, tmp_means, tmp_stds = matrix_whiten(score_matrix[:, :components_number], std=True)
	labels = gaussian_mixture(reduced_whitened_scores, upper_bound=max_clusters, SILENT=SILENT)
	
	
	# Print information about found clusters:
	cluster_number = len( np.unique(labels) )
	colored_clusters = color_clusters( score_matrix, labels )
	
	# Save the type the glitch belongs to in the "type" attribute
	for index, spike in enumerate(data_list):
		spike.type = labels[index]
	
	print_cluster_info(colored_clusters)	
	
	if ("frequency" in ANALYSIS):
		correlation_test(data_list, labels, ANALYSIS)
	elif (ANALYSIS == "time"):
		matched_filtering_test(data_list, labels, ANALYSIS)
	elif ( ANALYSIS == "generic"):
		pass
	else:
		assert False, "Wrong analysis type"
	
	# Save scatterplot with image maps:
	# images are saved to a subfolder, "Scatterplot_images", defined in 
	# scatterplot()
	output = "Scatterplot-%i_clusters-" % cluster_number
	for x in range( 1, MAX_PLOTTED_COMPONENTS+1 ):
		for y in range ( 1, MAX_PLOTTED_COMPONENTS+1 ):
			if ( x != y ) & ( x < y ):
				scatterplot( score_matrix, data_list, colored_clusters, labels, x, y, output+"%i_vs_%i" % (x, y), ANALYSIS )
		
	print "\n\tSaved scatterplots (html imagemaps and png)."
	
	# Save triangle plots with components_number (if not greater than
	# TRIANGLE_COMPONENTS) scatterplots
	triangle_components = TRIANGLE_COMPONENTS if components_number> TRIANGLE_COMPONENTS else components_number
	fig = plt.figure(figsize=(8*triangle_components,6*triangle_components), dpi=100)
	axs = []
	# n is the plot number in pyplot's grid system
	n = 1
	# If using coordinates x and y, then the number of the subplot being used 
	# is = MAX_PLOTTED_COMPONENTS*(y-1) + x
	# Things get a bit confusing here: x is the column index, y is the row index
	x = 1
	with warnings.catch_warnings():
		for y in range(1, triangle_components+1):
			while (x <= y):
				if ( x == y ):
					axs.append( fig.add_subplot(triangle_components, triangle_components, n) )
					warnings.filterwarnings( "ignore", category=UserWarning )
					axs[-1].set_title("Principal component {0} scores distribution".format(x))
					axs[-1].hist(score_matrix[:,x-1], bins=np.sqrt(observations), log=True, histtype="stepfilled", alpha=0.8)
				else:
					# Share axis with the plot above the current
					if (n > triangle_components):
						axs.append( fig.add_subplot(triangle_components, triangle_components, n, sharex=axs[(len(axs)-(y-1))]))
					else:
						axs.append( fig.add_subplot(triangle_components, triangle_components, n))
					plot_legend = []
					for index, element in enumerate(colored_clusters):
						tmp = np.array(element)
						axs[-1].plot( tmp[:,x-1] , tmp[:,y-1], markers_and_colors[index], label = str(index), markersize = 5 )
						plot_legend.append( str(index+1) )
					axs[-1].label_outer()
					axs[-1].legend(plot_legend, bbox_to_anchor=(0, 0, 1.12, 1), loc = 'best', markerscale = 2, numpoints = 1)
					if ( x == 1 ):
						axs[-1].set_ylabel("PC score: {0}".format(y))
					if ( y == triangle_components):
						axs[-1].set_xlabel("PC score: {0}".format(x))
				axs[-1].grid(linestyle = '--')
				n += 1
				x += 1
			# Reset x
			x = 1
			n = y*triangle_components + 1
	
		fig.savefig("All scatterplots.png", dpi=fig.get_dpi(), bbox_inches='tight', pad_inches=0.2)
		plt.close('all')
		del fig
		del axs
	print "\tSaved triangle plot ('All scatterplots.png')"	
	
	
	# Compute an "average transient" for each type.
	# The average transient is obtained by averaging all the glitches
	# in a certain type in the principal component space, by inverting
	# the PCA one obtains a representation in the original space.
	if ("bands" in ANALYSIS):
		calculate_types(data_list, colored_clusters, score_matrix, principal_components, means, stds, labels, ANALYSIS, sampling, low, high)
	else:
		calculate_types(data_list, colored_clusters, score_matrix, principal_components, means, stds, labels, ANALYSIS, sampling)
	
	
	plotted_components = MAX_ALL_PRINCIPAL_COMPONENTS
	if ("time" in ANALYSIS):
		fig_all, ax_all = configure_subplot_time(MAX_ALL_PRINCIPAL_COMPONENTS)
	elif ("frequency" in ANALYSIS):
		fig_all, ax_all = configure_subplot_freq(MAX_ALL_PRINCIPAL_COMPONENTS)
		for element in ax_all:
			element.autoscale(True, axis="both", tight=True)
	
	
	os.makedirs("Principal_components")
	# Plot the first  MAX_ALL_PRINCIPAL_COMPONENTS principal components:
	for i in range(0, MAX_ALL_PRINCIPAL_COMPONENTS):
		fig = plt.figure(figsize=(12, 6), dpi=100)
		ax = fig.add_subplot(111)
		
		ax.grid(which="both")
		
		if ( "time" in ANALYSIS ):
			time_axis = np.linspace(0, len(principal_components)/ANALYSIS_FREQUENCY, len(principal_components))
			# Convert to ms
			time_axis *= 1000
			time_axis -= time_axis[time_axis.size//2]
			
			ax.plot(time_axis, principal_components.transpose()[i], 'r', linewidth = 0.4 )
			ax.set_xlim((np.min(time_axis), np.max(time_axis)))
			ax.set_xlabel("Time [ms]")
			ax_all[i].set_ylabel("")
			ax_all[i].plot(time_axis, principal_components.transpose()[i], 'r', linewidth = 0.4 )
			ax_all[i].set_xlim(min(time_axis), max(time_axis))
			ax_all[i].set_autoscalex_on(False)
		elif ( "frequency" in ANALYSIS ):
			#ax.set_yscale("log")
			ax.set_xscale("log")
			ax.autoscale(True, axis="both", tight=True)
			if ("bands" in ANALYSIS):
				freq_array = np.linspace(low, high, num=len(principal_components.transpose()[i]) )
				ax.set_xlim((np.min(freq_array), np.max(freq_array)) )
				ax_all[i].set_xlim((np.min(freq_array), np.max(freq_array)) )
			else:
				freq_array = np.linspace(0, sampling/2.0, num=len(principal_components.transpose()[i]) )
				ax.set_xlim((np.min(freq_array), np.max(freq_array)) )
				ax.set_autoscalex_on(False)
				ax_all[i].set_xlim((10, np.max(freq_array)) )
				ax_all[i].set_autoscalex_on(False)
			ax.plot(freq_array, (np.power(10, principal_components.transpose()[i])), 'r', linewidth = 0.4 )
			ax_all[i].plot(freq_array, np.power(10, principal_components.transpose()[i]), 'r', linewidth = 0.4 )
			ax_all[i].set_ylabel("")
		
		ax_all[i].set_title("Principal Component: " + str(i+1) )
			
		ax.set_title( "Principal Component: %i" % (i+1) )
		fig.savefig("Principal_components/Principal_component-%i.png" % (i+1), dpi=DPI, bbox_inches='tight', pad_inches=0.2)
		plt.close(fig)
	fig_all.savefig("All_principal_components.pdf", bbox_inches='tight', pad_inches=0.2)
	plt.close('all')
	del ax_all
	del fig_all
	del ax
	del fig
	
	print "\tSaved Principal components plots."


	# Save information about the analyzed interval (intervals and total time)
	f = open("Analyzed_interval.txt", "w")
	f.write("#Start\t\tEnd\tDuration [s]\n")
	for interval in times:
		start_tmp, end_tmp = interval[0], interval[1]
		f.write("{0}\t{1}\t{2}\n".format(start_tmp, end_tmp, end_tmp-start_tmp) )
	f.write("Total analyzed time:\t{0:.1f} {1}".format(total/dividend, units))
	f.close()
	
	# Plot a glitchgram for the given time interval if
	# doing time-domain analysis.
	# If COLOR_ENERGY is True, glitches in the glitchgram are colored
	# according to energy, otherwise they are colored by type
	if "time" in ANALYSIS:
		start_time = times[0][0]
		end_time = times[-1][1]
		global glitchgram_start, glitchgram_end
		if glitchgram_start and glitchgram_end:
			plot_glitchgram(data_list, times, glitchgram_start, glitchgram_end, HIGH_PASS_CUTOFF, sampling, labels)
			for segments in times:
				plot_glitchgram(data_list, [segments], segments[0], segments[1], HIGH_PASS_CUTOFF, sampling, labels, name="Glitchgram_{0}-{1}".format(segments[0], segments[1]))
		else:
			plot_glitchgram(data_list, times, start_time, end_time, HIGH_PASS_CUTOFF, sampling, labels)
	
	# If time-domain analysis, create a folder, with subfolders for each
	# type, to save the transient's time series
	if ( "time" in ANALYSIS ) and not NOPLOT:
		print "\tPlotting time series..."
		for i in range(cluster_number):
			try:
				os.makedirs( "time_series/Type_"+str(i+1) )
			except:
				pass
		spike_time_series(data_list, (score_matrix, principal_components, means, stds), components_number, labels, sampling, RECONSTRUCT, SILENT)
	elif ( "frequency" in ANALYSIS ) and not NOPLOT:
		print "\tPlotting PSDs..."
		# Plot a list of PSDs with the relative times:
		for i in range(cluster_number):
			try:
				os.makedirs( "PSDs/Type_"+str(i+1) )
			except:
				pass
		if ("bands" in ANALYSIS):
			plot_psds(data_list, (score_matrix, principal_components, means, stds), components_number, labels, sampling, ANALYSIS, low, high, RECONSTRUCT, SILENT)
		else:
			plot_psds(data_list, (score_matrix, principal_components, means, stds), components_number, labels, sampling, ANALYSIS, RECONSTRUCT, SILENT)
	print ""
	
	# Dump the database to a pickle file
	pickle_dump(data_list, database_name)
	print "\tSaved {0}".format(database_name)
	
	# Analysis finished. Print output URL	
	print "#"*int(0.8*frame_width)
	print "\n\tResults at:"
	print "\t" + results_URL
	
	# Return to original working directory
	os.chdir(original_wd)
	

	return results_URL

def main():
	pipeline(sys.argv)
	
if __name__ == "__main__":
	main()
