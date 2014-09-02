#!/usr/bin/env python
# encoding: utf-8
# Daniele Trifiro`
# brethil@phy.olemiss.edu
'''
This program is used to perform PCA and cluster in the principal component
space data from finder.py

It's important to note that, even if this is written to use GMM, with a few
changes, every clustering algorithm (currently 'gaussian_mixture()')
can be implemented with few lines of code.

The clustering algorithm should output a set of labels for the input data. 
Input data is a matrix with the observations on the rows, output is a list
of the same length of the number of the observations. 
Each element of this list corresponds to the cluster number of which the
observation was assigned to (ranging from 0 to number of clusters-1).

If we had 5 observations, then:

	labels[0] = 1
	labels[1] = 0
	labels[2] = 1
	labels[3] = 1
	labels[4] = 0

in this case the observations 0, 2 and 3 are assigned to cluster 1, and
observations 1 and 4 are assigned to cluster' type 0, with a total of
2 clusters.

A function for plotting 3-D scatterplots is supplied (three_plot), this can be
used by simply un-commenting a line (search for 3DPLOT).

'''

SMALL_COMPONENTS_BOOL = False
SMALL_COMPONENTS = 50


import warnings

from utilities_PCAT import *
from PCA import standardize, eigensystem, PCA, load_data, matrix_whiten
import matplotlib.mlab
from matplotlib.image import NonUniformImage


import sklearn.mixture as mix

# mplot3d is used for the 3D plots.
#from mpl_toolkits.mplot3d import axes3d

def __usage__():
	'''
		Usage help.
	'''
	print "Usage:\tGMM.py [--time || --frequency] (--sampling, -s) sampl_freq [-m gauss_comp] [-p princ_comp]"
	print "\t\t [--low low, --high high] [--plots number] [--repr_components number]" #[--model model_file (not working)]"
	print "\t\t file1 file2 file3 ..."
	
	print "\n\tGMM performs PCA and clusters the data using a"
	print "\tGaussian Mixture Model."
	print "#"*frame_width
	print "\n\tOptions:"
	print "\t--time, --frequency\n\
		Type of analysis being performed, either in time\n\
		or in frequency.\n\
	 	--time should be used for databases output from finder.py,\n\
	 	and also plots the time series of the transients in the database\n\
		\"--frequency\" should be used for PSDs\n\
		transform databases.\n\
		For frequency bands analysis the syntax is the following:\n\
		\tGMM.py --frequency --low low_freq --high high_freq database.data"
		
	
	print "\t--log\n\
		Takes base-10 logarithm of the input data. This improves clustering\n\
		when data spans several order of magnitudes."
	
	print "\t--sampling sampl_freq, -s sampl_freq\n\t\t Specify a sampling frequency for the data's channel"
	
	#print "\t[--compare file]\n\t\t PCA is performed on the difference between observations\n\t\t\
	#		and input model (plain-text)."
		
	print "\t-m number, --maxclusters number\n\t\tSpecifies the maximum number\n\
		of clusters. (Default = 10 )."
	print "\t-p number, --components\n\t\tNumber of principal components to be used when clustering\n\t\t (default is 40).\n\
		One can guess this number by simply looking at the explained variance\n\
		graph by using PCA.py with the --extras option (explained-variance.pdf).\n\
		Rule of thumb: choose a number of components which accounts for about 70%\n\
		of the total variance."
	#print "\t-S, --standardize\n\t\tStandardizes the data in standard deviation,\n\t\t\
	#			so that each variable has unit variance (default is NOT standardized)."
	
	print "\t--plots number\n\t\tPlots scatterplots up to the first 'number' components."
	
	print "\t--repr_components number\n\t\tSets the number of components with which the representative\n\
		transients are computed."
	
	print "\t-r\n\t\t Interactive mode, chose wich clusters to remove from the database."
	print "\t--noplot\n\t\t Do not plot time series/PSDs, use this to save time if not\n\
		interested in time series and/or using -r."
	"""
	print "\t--model model_file\n\t\t Specify a model file to use, this model is used\n\
		in the frequency domain analysis. To plot the differences against an 'ideal' spectrum."
	"""
	
	print "\t--silent\n\t\tDo not print progress bars."
	#----print "\t--list\n\t\t Saves GPS start and end times for each type's spectra."


def __check_options_and_args__():
	global MAX_CLUSTERS, max_clusters
	MAX_CLUSTERS = False
	max_clusters = 10
	
	global MAX_COMPONENTS
	MAX_COMPONENTS = 4	
	global principal_components_number
	principal_components_number = 8
	
	global marker
	marker = '+'
	
	global SILENT
	SILENT = False
	
	global REMOVE
	REMOVE = False
	
	global MODEL
	MODEL = False
	global MODEL_FILE
	MODEL_FILE = ""
	
	global SAMPLING
	
	global PRINT_LIST
	PRINT_LIST = False
	
	global ANALYSIS
	
	global PLOT
	PLOT = True
	
	global low, high
	#	low, high = None, None
	
	global COMPONENTS
	COMPONENTS = 4
	try:
		opts, args = getopt.getopt(sys.argv[1:],"hm:s:p:r", ["help", 'maxclusters=', 'list',\
		 							'standardize', 'marker=', 'model=', "sampling=", "plots=",\
									'time', 'frequency', 'low=', 'high=', 'log', 'bands', "repr_components=", "noplot", "silent"] )
	except getopt.error, msg:
		print msg
		sys.exit(1)
	if ( len(sys.argv[1:]) == 0 and not (any( flag in o for flag in [ '--help', '-h'] for o in opts)) ) :
		print "No input files."
		print "Use GMM.py -h for usage." 
		print "Example of usage:"
		print "\t GMM.py --time -s 32768 file.list"
		print "\t GMM.py --frequency -s 32768 file.list\n"
		print "\t GMM.py [--log] matrix_database.data\n"
		print "\t GMM.py --freq --low --high matrix_database.data\n"
		
		sys.exit(1)
	for o, a in opts:
		if o in ( '-h', '--help' ):
			__usage__()
			sys.exit(1)
		elif o in ( '-m', 'maxclusters' ):
			max_clusters = int(a)+1
		elif o in ( '--marker' ):
			marker = a
		elif o in ( '-p' ):
			principal_components_number = int(a)
		elif o in  ( '-r' ):
			REMOVE = True
		elif o in ( '--model' ):
			MODEL = True
			MODEL_FILE = a
		elif o in ( '-s', '--sampling' ):
			SAMPLING = float(a)
		elif o in ( '--list' ):
			PRINT_LIST = True
		elif o in ( '--low' ):
			low = int(a)
		elif o in ( '--high' ):
			high = int(a)
		elif o in ( '--time' ):
			ANALYSIS = 'time'
		elif o in ( '--frequency' ):
			ANALYSIS = 'frequency'
		elif o in ( '--plots' ):
			MAX_COMPONENTS = int(a)+1
		elif ( o == "--log"):
		 	# Pass, "_log" will later be appended to ANALYSIS
			pass
		elif o in ( '--repr_components' ):
			SMALL_COMPONENTS_BOOL = True
			SMALL_COMPONENTS_BOOL = int(a)
		elif (o == "--noplot"):
			PLOT = False
		elif (o == "--silent"):
			SILENT = True
		else:
			assert False, "Unknown Option."
	
	
	if not ( any( flag in o for flag in ['--time', "--frequency"] for o in opts ) ):
		#print "Analysis type (time domain or frequency domain) has to be
		# supplied"
		#print "through the --time or --frequency flags. Quitting."
		#sys.exit()
		ANALYSIS = "generic"
		if ( any( "--log" in o for o in opts) ):
			ANALYSIS += "_log"
	
	if any( flag in o for flag in ['--low'] for o in opts) and\
	 	any(flag in o for flag in ["--high"] for o in opts):
		ANALYSIS += "_bands"	
	
	elif not ( any( flag in o for flag in ['-s', "--sampling"] for o in opts )):
		
		if not any( flag in o for flag in ['--low', '--high'] for o in opts ):
			print "Sampling frequency has to be supplied. Quitting."
			sys.exit()
		elif not ( any( flag in o for flag in ['--low'] for o in opts ) and ( any( flag in o for flag in ["--high"] for o in opts ))):
			print "Both --low and --high have to be supplied if performing bands analysis. Quitting."
			sys.exit()
		
	return args
 


def gaussian_mixture(matrix, upper_bound, SILENT=False):
	'''
	This function clusters the input matrix using the GMM algorithm (gaussian mixture model)
	The number of clusters is found by running the algorithm for n_components = 2 to upper_bound
	and chosing the model which minimized the BIC.
	
	Returns the labels for each observation.
	'''
	if ( len(matrix) < upper_bound+1 ):
		print "\n\tWARNING: Not enough samples (less than the minimum %i) to run GMM." % (upper_bound)
		print "\t Only one cluster is returned.\n"
		return [0]*len(matrix)
	# Create progress bar
	if not SILENT:
		progress = progressBar(minValue = 0, maxValue = 4*upper_bound-1, totalWidth = 40 )
		progress(0)
	
	j = 0
	lowest_bic = np.infty
	bic = []
	n_components_range = range (1, upper_bound+1)
	cv_types = ['spherical', 'tied', 'diag', 'full']
	for cv_type in cv_types:
		for n in n_components_range	:
			gmm = mix.GMM(n_components = n, covariance_type = cv_type)
			gmm.fit(matrix)
			bic.append( gmm.bic(matrix) )
			if bic[-1] < lowest_bic:
				lowest_bic = bic[-1]
				best_gmm = gmm
			if not SILENT:
				progress(j)
			j += 1
		if not SILENT:
			progress(j)
			j += 1
		
	best_gmm.fit(matrix)
	res = best_gmm.predict(matrix)
	
	# Print an empty line to avoid printing on the progress bar. 
	if not SILENT:
		print ""
	
	return res


"""
def gaussian_mixture(matrix, upper_bound, weights):
	'''
	This function clusters the input matrix using the DPGMM algorithm (Dirichlet Process Gaussian Mixture Model)
	Returns the labels for each observation.
	'''
	if ( len(matrix) < upper_bound ):
		print "\n\tWARNING: Not enough samples (less than the minimum %i) to run GMM." % (upper_bound-1)
		print "\t Only one cluster is returned.\n"
		return [0]*len(matrix)
	
	gmm = mix.DPGMM(n_components = upper_bound)
	gmm.weights_ = weights
	gmm.fit(matrix)
	res = gmm.predict(matrix)
	
	return res
"""

def color_clusters(score_matrix, labels):
	'''
	This function "colors" the cluster, that is, for each different cluster in 'labels',
	a list with all the observations corresponding to that cluster is created.
	
	Returns a list, where each item is the list of observations corresponding to a cluster.
	'''
	cluster_number = len( np.unique(labels) )
	if ( cluster_number == 1 ):
		return [score_matrix.tolist()]
	
	colored_list = [ list() for i in range(0, cluster_number)]
	for index, spike in enumerate(score_matrix):
		colored_list[ labels[index] ].append(spike)
	return colored_list


def three_plot(colored_clusters, x, y, z, output):
	'''
	Plots 3-D scatterplot of colored_clusters, using components
	x, y and z, output is the output file name.
	'''
	global DPI
	global colors
	fig = plt.figure()
	ax = fig.add_subplot(111, projection = '3d')
	ax.set_xlabel( "Principal component score: "+str(x) )
	ax.set_ylabel( "Principal component score: "+str(y) )
	ax.set_zlabel( "Principal component score: "+str(z) )
	plotlabels = list()
	ax.set_title("{0} clusters and noise".format(len(colored_clusters)))
	for index, element in enumerate(colored_clusters):
		tmp = np.matrix(element)
		ax.scatter(tmp[:,0].tolist(), tmp[:,1].tolist(),tmp[:,2].tolist(),\
				 			c = str(colors[index]), marker = "." , markersize = 1)
		plotlabels.append(str(index+1))
	ax.legend(plotlabels, loc = 'upper right', markerscale = 5)
	ax.grid( linestyle = '--' )
	plt.show()
	fig.savefig(output)
	print "Written:\t"+output+"."


def print_cluster_info(colored_clusters):
	print "\t\tCluster\t\t\t Size" 
	print "\t\t" + "#"*45
	# Output to file:
	f = open( "Types_detail.txt", "w")
	f.write("\t\tCluster\t\t Length\n\t\t"+"#"*45+"\n")
	all_transients = [len(element) for element in colored_clusters]
	total = float(np.sum(all_transients))
	for index, element in enumerate(all_transients):
		percent = (element/total)*100
		output = "\t\tType {0:2d} \t\t {1:3d}\t ({2:4.1f}%)".format(index+1, element, percent)
		print output
		f.write(output + "\n")
	




def set_axes_frequency(freq_array, logy=True):
	fig = plt.figure(figsize=(12, 6), dpi=300)
	ax = fig.add_subplot(111)
	
	ax.set_xscale( 'log' )
	
	ax.set_xlabel( "Frequency [Hz]" )
	ax.set_ylabel( "Power Spectral Density [Counts^2/Hz]" )
	
	ax.xaxis.grid( True, 'minor' )
	ax.xaxis.grid( True, 'major' )
	if ( logy == True ):
		ax.set_yscale( 'log' )
		ax.yaxis.grid( True, 'minor' )
		ax.yaxis.grid( True, 'major' )
	
	ax.set_xbound( lower=10, upper=np.max(freq_array) )
	ax.set_autoscalex_on(False)
	return fig, ax


def set_axes_frequency_band(freq_array, waveform_len,logy=True):
	fig = plt.figure(figsize=(12, 6), dpi=300)
	ax = fig.add_subplot(111)
	ax.set_xscale( 'log' )
	ax.set_xlabel( "Frequency [Hz]" )
	ax.set_ylabel( "Power Spectral Density [Counts^2/Hz]" )
	
	if ( logy == True ):
		ax.set_yscale( 'log' )
	
	ax.yaxis.grid( which="both")
	ax.xaxis.grid( which="both")
	ax._autoscaleXon = False
	ax.set_xbound( lower=np.min(freq_array), upper=np.max(freq_array) )
	return fig, ax


def set_axes_time(time_array):
	fig = plt.figure(figsize=(12, 6), dpi=300)
	ax = fig.add_subplot(111)
	ax.set_xlabel( " Time [ms]" )
	ax.set_ylabel( "Amplitude [Counts]" )
	ax.grid(which="both")#, axis="both")
	
	x_min, x_max = np.min(time_array), np.max(time_array)
		
	ax.set_xbound( lower=x_min, upper=x_max )
	ax.set_autoscalex_on(False)
	return fig, ax

def configure_subplot_time(subplot_number):
	"Configure 'subplot_number' subplots, arranged in a column"
	fig = plt.figure(figsize=(12, 6*subplot_number), dpi=300)
	ax = []
	ax.append( fig.add_subplot(subplot_number, 1, 1) )

	for index in range(subplot_number):
		if (index != 0):
			ax.append( fig.add_subplot(subplot_number, 1, index+1, sharex=ax[index-1]) )
		ax[index].grid(which='both')#, axis='both')
		ax[index].set_title("Type " + str(index+1) )
		ax[index].set_ylabel("Amplitude [counts]")
		ax[index].set_xlabel("Time [ms]")
		
	
	return fig, ax

def configure_subplot_freq(subplot_number, logy=True):
	"Configure 'subplot_number' subplots, arranged in a column"
	fig = plt.figure(figsize=(12, 6*subplot_number), dpi=300)
	ax = []
	ax.append( fig.add_subplot(subplot_number, 1, 1) )
	for index in range(subplot_number):
		if (index != 0):
			ax.append( fig.add_subplot(subplot_number, 1, index+1, sharex=ax[index-1]) )
		ax[index].set_title("Type " + str(index+1) )
		ax[index].set_ylabel("Power Spectral Density [Counts^2/Hz]")
		ax[index].set_xlabel("Frequency [Hz]")
		ax[index].set_xscale('log')
		ax[index].xaxis.grid( True, 'minor' )
		ax[index].xaxis.grid( True, 'major' )
		if ( logy == True ):
			ax[index].set_yscale( 'log' )
			ax[index].yaxis.grid( True, 'minor' )
			ax[index].yaxis.grid( True, 'major' )
		
	
	return fig, ax


def calculate_types(database, clusters, score_matrix, principal_components, means, stds, labels, ANALYSIS, f_sampl, low=None, high=None):
	'''
		For each cluster compute an average observation from all observations in that cluster using a median.
		The median is used for outlier rejection.
		
		For each average observation, perform inverse PCA (multiplying by the transpose of the 
		principal components matrix and then adding the mean for each observation)
		
		After the average types have been computed, these are plotted.
	
	'''
	
	DPI = 100
	cluster_number = len( clusters )
	"""
	# The following commented code gives out weird forms for the representative transient due to the adding of means at the end
	# but might be a good way to get a good estimate of the "true" transient form
	if not SMALL_COMPONENTS_BOOL:
		rows, columns = np.shape(score_matrix)
	else:
		rows, columns = np.shape(score_matrix[:,0:SMALL_COMPONENTS])
	cluster_matrices = []
	cluster_medians = []
	# For each cluster in 'clusters' compute the median and save it in 
	# cluster_medians
	for cluster in clusters:
		if len(cluster) > 1:
			cluster_matrices.append( np.array(cluster) )
			
			# SMALL_COMPONENTS_BOOL is defined at the top of the document.
			# If TRUE, only the first SMALL_COMPONENTS principal component 
			# scores
			# are used to reconstruct
			# the average observations. This should remove noise, as the higher 
			# order principal components
			# mainly consist in noise.
			
			if not SMALL_COMPONENTS_BOOL:
				cluster_medians.append( np.median(cluster_matrices[-1], axis=0) )
			else:
				cluster_medians.append( np.median(cluster_matrices[-1], axis=0)[0:SMALL_COMPONENTS] )
		else:
			cluster_medians.append(cluster[0])
	
	# Invert PCA: multiply by the transpose of the principal components matrix
	# and add means
	average_observations = []
	for median in cluster_medians:
		average_observations.append( np.dot( median, principal_components.transpose() ) )
	average_observation_matrix = np.array(average_observations)
	average_observation_matrix *= stds
	average_observation_matrix += means
	"""
	
	time_domain_clusters = []
	for i in range(cluster_number):
		time_domain_clusters.append([])
	
	for index, spike in enumerate(database):
		time_domain_clusters[labels[index]].append(spike.waveform)
	
	cluster_medians = []
	for cluster in time_domain_clusters:
		cluster = np.array(cluster)
		cluster_medians.append(np.median(cluster, axis=0))
		
		
	average_observation_matrix = np.array(cluster_medians)
	
	
	
	'''
			THIS PART IS TO BE REVIEWD
			IS THIS NECESSARY AT ALL?
	####################################################
	# Mean max amplitude is calculated for
	# every transient in the cluster, and used as peak amplitude
	# for the representative transient.
	# This should be used only when normalizing all spikes to unit amplitude
	if ( "time" in ANALYSIS ):
		peaks = [ list() for element in np.unique(labels) ]
		for index, spike in enumerate(database):
			if ( len(spike.waveform) > 0 ):
				peaks[ labels[index] ].append( spike.norm )
		peaks_means = [ np.mean(element) for element in peaks ]
		for index, element in enumerate(datamatrix):
			element *= peaks_means[index]
	####################################################
	'''
	
	observations, waveform_len = np.shape(average_observation_matrix)
	
	# Initialize axes for summary plot
	if ( "frequency" in ANALYSIS ):
		if ("bands" in ANALYSIS):	
			freq_array = np.linspace(low, high, waveform_len)
		else:
			#freq_array = rfftfreq( 2*(waveform_len-1), d=1./f_sampl )
			freq_array, tmp = matplotlib.mlab.psd
		fig_all, ax_all = configure_subplot_freq(len(average_observation_matrix)+1)
		if ("bands" in ANALYSIS):
			for element in ax_all:
				element.autoscale(True, "both", tight=True)
		else:
			for element in ax_all:
				element.set_xbound( lower=10, upper=np.max(freq_array) )
				element.autoscale(True, "y", tight=True)
	elif ( "time" in ANALYSIS ):
		time_axis =  (np.array(range(waveform_len))/f_sampl)*1000.0
		max_index = waveform_len//2
		time_axis -= time_axis[max_index] 
		
		fig_all, ax_all = configure_subplot_time(len(average_observation_matrix))
	elif ( "generic" in ANALYSIS ):
		fig_all = plt.figure()
		ax_all = fig_all.add_subplot(111)
		ax_all.grid()
		if ( 'log' in ANALYSIS ):
			ax_all.set_yscale('log')
	else:
		assert False, "Fatal error with analysis type. Quitting."
	plotlabels = []
	
	# Default line marker is a continous line, switch to
	# dotted line if there are more than 7 types
	marker = "-"
	with warnings.catch_warnings():
		warnings.filterwarnings( "ignore", category=UserWarning )
		if not "generic" in ANALYSIS:
			ax_all[-1].set_title("Summary")
		for index, element in enumerate(average_observation_matrix):
			percent = (len(clusters[index])/float(len(database)))*100.0
			if ( "frequency" in ANALYSIS ):
				
				if ( "bands" in ANALYSIS):
					fig, ax = set_axes_frequency_band(freq_array, waveform_len)
					ax.set_xticks(np.logspace(np.log10(low), np.log10(high), num=10))
					ax.set_xticklabels([ "%.2f" % el for el in ax.get_xticks()])
				else:
					fig, ax = set_axes_frequency(freq_array)
				
				# Only plot frequencies above 10 Hz:
				# only choose indexes corresponding to freq_array>10
				if not ( "bands" in ANALYSIS):
					ax.plot( freq_array, np.power(10, element), "r-", linewidth = 0.4)
					ax_all[index].plot( freq_array, np.power(10, element), "r-", linewidth = 0.4)
					ax_all[index].set_title("Type {0:d}: {1:d} of {2:d} observations ({3:.1f}%)".format(index+1, len(clusters[index]), len(database), percent) )
					ax_all[-1].plot( freq_array, np.power(10, element), markers_and_colors[index][1]+marker, linewidth=0.4)
				elif ("bands" in ANALYSIS):
					ax.plot(freq_array, np.power(10, element), "r-", linewidth = 0.4)
					ax_all[index].plot( freq_array, np.power(10, element), "r-", linewidth = 0.4)
					ax_all[index].set_title("Type {0:d}: {1:d} of {2:d} observations ({3:.1f}%)".format(index+1, len(clusters[index]), len(database), percent) )
					ax_all[index].set_xticks( np.logspace(np.log10(low), np.log10(high), num=10 ))
					ax_all[index].set_xticklabels([ "%.2f" % el for el in ax.get_xticks()])
					ax_all[-1].plot( freq_array, np.power(10, element), markers_and_colors[index][1]+marker, linewidth=0.4)
				# Change marker: there only are 7 colors
				if (index > 5) and (marker == "-"):
					marker = "-."
				ax.autoscale(True, "both", tight=True)
				ax_all[index].axis("tight")
			elif ( "time" in ANALYSIS ):
				fig, ax = set_axes_time(time_axis)
				ax.plot(time_axis, element, 'b-', linewidth = 0.4 )
				ax_all[index].plot(time_axis, element, "b-", linewidth = 0.4)
				ax_all[index].autoscale(True, "both", tight=True)
				ax_all[index].set_title("Type {0:d}: {1:d} of {2:d} observations ({3:.1f}%)".format(index+1, len(clusters[index]), len(database), percent) )
			elif ( "generic" in ANALYSIS ):
				fig = plt.figure()
				ax = fig.add_subplot(111)
				ax.grid()
				if ( "log" in ANALYSIS ):
					ax.set_yscale('log')
					ax.plot(np.power(10, element), 'r-', linewidth = 0.4)
					ax_all.plot(np.power(10, element), markers_and_colors[index]+'-', linewidth = 0.4)
				else:
					ax.plot(element, 'r-', linewidth = 0.4)
					ax_all.plot(element, markers_and_colors[index]+'-', linewidth = 0.4)
				
			output = str(cluster_number) + "-clusters_#" + str(index+1) + ".pdf"
			plotlabels.append( "{0:d} ({1:.1f}%)".format(index+1, percent) )
			
			ax.set_title("Type {0:d}: {1:d} of {2:d} observations ({3:.1f}%)".format(index+1, len(clusters[index]), len(database), percent) )
			if ( "frequency" in ANALYSIS):
				plt.autoscale(True, axis="y", tight=True)
			
			fig.savefig(output, dpi = DPI)
	
	
	if ( "time" in ANALYSIS):
		x_min, x_max = np.min(time_axis), np.max(time_axis)
		ax_all[0].set_xbound( lower=x_min, upper=x_max )
		ax_all[0].set_autoscalex_on(False)
	else:
		if ("frequency") in ANALYSIS and ("bands" in ANALYSIS):
			x_min, x_max = np.min(freq_array), np.max(freq_array)
			ax_all[0].set_xlim( (x_min, x_max) )
			ax_all[0].set_autoscalex_on(False)
			ax_all[0].autoscale(enable=True, axis="y", tight=True)
			ax_all[-1].set_xticks( np.logspace(np.log10(low), np.log10(high), num=10) )
			ax_all[-1].set_xticklabels([ "%.2f" % el for el in ax.get_xticks()])
		elif ("frequency" in ANALYSIS):
			ax_all[0].set_xbound( lower=10, upper=np.max(freq_array) )
			ax_all[-1].axis("tight")
		box = ax_all[-1].get_position()
		ax_all[-1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax_all[-1].legend( plotlabels, loc="best", bbox_to_anchor=(1.25, 1))
		ax_all[-1].axis("tight")
	
	print "\tSaving types details..."
	fig_all.savefig( str(cluster_number) +" - All_types.pdf", dpi = DPI, bbox_inches='tight', pad_inches=0.2)
	print "\tSaved \"All_types.pdf\"."
	plt.close('all')	
	print_lists(database, labels, cluster_number, ANALYSIS)
	
	# If time analysis, plot a spectrogram for each of the average observations
	if ( "time" in ANALYSIS ):
		if not ( os.path.exists( "spectrograms" ) ):
			os.mkdir( "spectrograms" )
		print "\n\tPlotting spectrograms..."
		for index, element in enumerate(average_observation_matrix):
			with warnings.catch_warnings():
				warnings.filterwarnings( "ignore", category=UserWarning )
				
				create_spectrogram(element, f_sampl, "Type_"+str(index+1)+"_Spectrogram", "Type " + str(index+1)) 
		if cluster_number > 1:
			print "\tSaved: Type_[{0}-{1}]_Spectrogram.pdf".format(1, len(average_observation_matrix))
		else:
			print "\tSaved: Type_1_Spectrogram.pdf"
		
	plt.close('all')


def remove_clusters(database, labels):
	OUTPUT = "database_new"
	cluster_number = len( np.unique(labels) )
	keyboard_input = raw_input( "Insert the cluster number to be removed: ('q' or 'quit' to exit)\n" )
	clusters_to_remove = list()
	while ( keyboard_input != ( 'q', 'quit' ) ):
		if ( keyboard_input in str(range( 1, cluster_number+1) ) ):
			var = int(keyboard_input)-1
			clusters_to_remove.append(var)
			keyboard_input = raw_input( "Any other clusters? ('q' or 'quit' to exit)\n" )
			if ( keyboard_input in ( 'q', 'quit') ):
				break
		else:
			if ( keyboard_input in ( 'q', 'quit') ):
				break
			else:
				keyboard_input = raw_input( "Input has to be a number in the range 1 to "+str(cluster_number)+". Try again:\n" )
	indexes_to_remove = list()
	new_database = []
	removed_transients = 0
	for index, observation in enumerate(database):
		if labels[index] not in clusters_to_remove:
			new_database.append(observation)
		else:
			removed_transients += 1	
	
	print removed_transients, "transients removed."
	print "Saving "+OUTPUT+".list"
	f = open( OUTPUT+".list", "wb" )
	pickle.dump(new_database, f)
	f.close()
	print "You now may re-run GMM.py on "+OUTPUT+".list."


def print_lists(database, labels, cluster_number, ANALYSIS):
	if not ( os.path.exists( "Types_detail" ) ):
		os.mkdir( "Types_detail" )
	
	for index in range(0, cluster_number):
		OUTPUT = "Type_" + str(index+1) + ".txt"
		f = open( "Types_detail/" + OUTPUT, "w")
		
		if ( 'frequency' in ANALYSIS ):
			for j_index, observation in enumerate(database):
				if ( labels[j_index] == index ):
					f.write( str(observation.start)+"\t"+str(observation.end)+"\n" )
		elif ( 'time' in ANALYSIS ):
			for j_index, observation in enumerate(database):
				if ( labels[j_index] == index ):
					f.write( str(observation.peak_GPS)+"\n")
		f.close()
		#print "\tSaved: "+OUTPUT+"."


def create_spectrogram(data, sampling, output, title=""):
	"""
		Creates a spectrogram of the input data, and saves the spectrogram
		in "spectrograms/filename.pdf" where 'filename' is the string 'output'.
		Resolution is the spectrogram's time resolution in milliseconds
	"""
	fig = plt.figure( figsize=(12, 6), dpi=100 )
	ax = fig.add_subplot(111)
	
	# Set FFT to len(data) to obtain the maximum resolution in frequency
	NFFT = len(data)//8
	# Compute a spectrogram using matplotlib.mlab.specgram
	pxx, freq, time = matplotlib.mlab.specgram( data, NFFT=NFFT, Fs=sampling, noverlap=NFFT//4, pad_to=int(sampling) )
	
	halfbin_time = (time[1] - time[0]) / 2.0
	halfbin_freq = (freq[1] - freq[0]) / 2.0
	# Center time and frequency bins
	time -= halfbin_time
	freq -= halfbin_freq
	
	
	# Change time axis to be in milliseconds and centered on the
	# transient's maximum (which at vector's midpoint)
	half_time = int(len(time)/2.0)
	time = (time - time[half_time])*1000
	
	# Plot the spectrogram	
	cax = ax.contourf(time, freq, 10*np.log10(pxx), 50)
	# Plot the colorbar
	cbar = plt.colorbar(cax, ax=ax, orientation="vertical", shrink=0.7, fraction=0.05, pad=0.01, spacing="uniform")
	cbar.set_label("Energy (dB)")
	
	# Values for interpolation karg are *None*, 'none', 'nearest', 'bilinear',
	#'bicubic', 'spline16', 'spline36', 'hanning', 'hamming',
	#'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian',
	#'bessel', 'mitchell', 'sinc', 'lanczos'
	
	#plt.autoscale()
	#title += ", Time resolution: %i ms  (%i points) " % (RESOLUTION, NFFT)
	title += " , Spectrogram - {0:.1f}Hz Frequency Resolution".format(sampling/NFFT)
	ax.set_title( title )
	
	ax.set_xlabel( "Time [ms]" )
	ax.set_ylabel( "Frequency [Hz]" )
	ax.yaxis.grid(which="both")
	ax.xaxis.grid(which="major")
	ax.set_yscale('log')
	#ax.set_xlim( (np.min(time), np.max(time)) )
	#plt.minorticks_on()
	#ax.tick_params(axis='y', which='both', labelleft='on', left='on')
	
	plt.autoscale(True, tight=True, axis="both")
	ax.set_ylim( (10, np.max(freq)) )
	#plt.clim(0,1)
	
	
	plt.savefig("spectrograms/" + output + ".png", bbox_inches='tight', pad_inches=0.2, dpi=100)
	
	plt.close('all')


def spike_time_series(database, PCA_info, components_number, labels, f_sampl, RECONSTRUCT=False, SILENT=False):
	'''
		Create a time series for each element contained in database.
		- database is a list of Spike() istances
		- PCA_info is a tuple with: (score_matrix, principal_components, means, stds )
			where means and stds are arrays with the column means and column standard
			deviations of the original data matrix (np.array([spike.waveform for spike in database])).
			This is only used if RECONSTRUCT=True (see below)
		- components_number is the number of principal components to be used when reconstructing the glitch.
			In order to plot the reconstructed glitches one has to use RECONSTRUCT=True, otherwise the raw
			time series are plotted.
		- labels contains the information about which cluster each glitch belongs to (gaussian_mixture() output)
		- f_sampl sampling frequency for the given glitches' time series
		- RECONSTRUCT (bool), if True, then glitches' time series are reconstructed using components_number 
		  principal components. Only using the first few principal components will reduce noise in the time series
		  and make the "true" shape of the glitch more clear.
		
		The reconstructed time series are also saved in the spike.reconstructed attribute of spike, with the number
		of components used saved in spike.PCs_used

	'''
	spikes_number = len(database)
	waveform_length = len(database[0].waveform)
	
	# If RECONSTRUCT is true, generate an array of 'reconstructed' using the 
	# first few principal components
	if RECONSTRUCT:
		reconstructed = []
		# Unpack the tuple
		score_matrix, principal_components, means, sigmas = PCA_info
		# Replace all the coefficients for principal components with index 
		# greater than "components_number" with 0
		score_matrix[:, components_number:] = 0.0
		# Invert the transformation and append to reconstructed list
		new_data = np.dot(score_matrix, principal_components.transpose())
		new_data *= sigmas
		new_data += means
		for index, spike in enumerate(database):
			# Before the PCA is computed, all waveforms are normalized
			# by divind them by the norm attribute. To correctly invert
			# the PCA, one has to multiply by norm.
			tmp = new_data[index]*database[index].norm
			reconstructed.append(tmp)
			spike.reconstructed = tmp
			spike.PCs_used = components_number
	
	
	# Create a progress bar
	if not SILENT:
		if (spikes_number > 1):
			progress = progressBar(minValue = 0, maxValue = spikes_number-1, totalWidth = 40 )
			progress(0)
	
	if f_sampl:
		x_axis = np.array([ i/f_sampl for i in range( 0, waveform_length ) ])
		# Multiply by 1000 to show time axis in milliseconds
		x_axis *= 1000.0
		
		# Center the x axis on the transients' peak
		max_index = np.argmax(np.abs(database[0].waveform))
		x_axis -= x_axis[max_index] 
		x_min, x_max = np.min(x_axis), np.max(x_axis)
	
	# Start plotting
	for index, spike in enumerate(database):
		labels_list = []
		if RECONSTRUCT:
			fig = plt.figure(figsize=(12,12), dpi=300)
			ax = fig.add_subplot(211)
			ax1 = fig.add_subplot(212)
			ax1.grid(ls = '--')
		else:
			fig = plt.figure(figsize=(12,6), dpi=300)
			ax = fig.add_subplot(111)
		
		ax.set_title( "Peak at: " + str(spike.peak_GPS) )
		ax.grid( ls = '--' )
		ax.set_autoscalex_on(False)
		
		
		if f_sampl:
			ax.set_xlim( ( x_min, x_max ) )
			ax.plot( x_axis, spike.waveform, "b", label="Raw")
			labels_list.append("Raw time series")
			ax.set_xlabel("Time [ms]")
			ax.set_ylabel("Amplitude [counts] ")
			if RECONSTRUCT:
				ax.plot( x_axis, reconstructed[index], 'r', label="Reconstructed - {0} PCs".format(components_number))
				labels_list.append( "Reconstructed - {0} PCs".format(components_number) )
				ax.legend(labels_list, loc = 'best', markerscale = 2, numpoints = 1)
				ax1.plot( x_axis, spike.waveform, "b", label="Raw 2")
				ax1.set_xlim( ( x_min, x_max ) )
				ax1.set_xlabel("Time [ms]")
				ax1.set_ylabel("Amplitude [counts] ")
				ax1.legend(labels_list, loc='best', markerscale = 2, numpoints =1)	
		else:
			plt.xlim( ( 0, waveform_length ) )
			ax.plot(spike.waveform, label="Raw")
			if RECONSTRUCT:
				ax.plot( reconstructed[index], "r", label="Reconstructed - {0} PCs".format(components_number) )
				labels_list.append( "Reconstructed - {0} PCs".format(components_number) )
				ax.legend(labels_list, loc = 'best', markerscale = 2, numpoints = 1)
				ax1.plot( x_axis, spike.waveform, "b", label="Raw" )
				ax1.legend(labels_list, loc="best", markerscale=2, numpoints=1)
			
		
		fig.savefig( "time_series/Type_%i/%.3f.pdf" % (labels[index]+1, spike.peak_GPS), bbox_inches='tight', pad_inches=0.2)
		plt.close('all')
		del fig
		
		if not SILENT:
			if ( spikes_number > 1 ):
				progress(index+1)
	del labels_list


def plot_psds(database, PCA_info, components_number, labels, f_sampl, ANALYSIS="frequency", low=None, high=None, RECONSTRUCT=False, SILENT=False):
	'''
		Plot all the PSDs in database, putting each in a different folder according to
		the labels in 'labels'.
		
		f_sampl is the sampling frequency for the given PSDs.
		
		If "bands" in analysis, rescale the axis to be between 'low' and 'high' (global variables).
		
		The reconstructed time series are also saved in the spike.reconstructed attribute of spike, with the number
		of components used saved in spike.PCs_used
	'''
	psd_number = len(database)
	waveform_length = len(database[0].waveform)
	
	# If RECONSTRUCT is true, generate an array of 'reconstructed' using the 
	# first few principal components
	if RECONSTRUCT:
		reconstructed = []
		# Unpack the tuple
		score_matrix, principal_components, means, stds = PCA_info
		# Replace all the coefficients for principal components with index 
		# greater than "components_number" with 0
		score_matrix[:, components_number:] = 0.0
		# Invert the transformation
		new_data = np.dot(score_matrix, principal_components.transpose())
		new_data *= stds
		new_data += means
		# Replace time series in the database with the new reconstructed
		# time series
		for index, spike in enumerate(database):
			tmp = new_data[index]*database[index].norm
			reconstructed.append(tmp)
			spike.reconstructed = tmp
			spike.PCs_used = components_number
	
	
	# Create a progress bar
	if not SILENT:
		if ( psd_number > 1):
			progress = progressBar(minValue = 0, maxValue = psd_number-1, totalWidth = 40 )
			progress(0)
	
	# Initialize axes
	if ( "bands" in ANALYSIS ):
		freq_array = np.linspace(low, high, waveform_length)
	else:
		freq_array = rfftfreq( 2*(waveform_length-1), d=1./f_sampl )
	
	# Start plotting
	for index, spike in enumerate(database):
		labels_list = []
		if ( "bands" in ANALYSIS):
			fig, ax = set_axes_frequency_band(freq_array, waveform_length)
		else:
			fig, ax = set_axes_frequency(freq_array)
		ax.set_title( "PSD: GPS %i to %i " % ( spike.segment_start, spike.segment_end) )
		ax.grid( ls = '--', which="both" )
		ax.set_autoscalex_on(False)
		ax.set_autoscaley_on(True)
		if RECONSTRUCT:
			if ( "bands" in ANALYSIS):
				ax.plot( freq_array, np.power(10, reconstructed[index]), label="Reconstructed - {0} PCs".format(components_number) )
			else:
				ax.plot( freq_array, np.power(10,reconstructed[index]), label="Reconstructed - {0} PCs".format(components_number) )
			labels_list.append( "Reconstructed - {0} PCs".format(components_number) )
			ax.legend(labels_list, loc = 'best', markerscale = 2, numpoints = 1)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlabel("Frequency [Hz]")
		ax.set_ylabel("Power Spectral Density [Counts^2/Hz)]")
		if ( "bands" in ANALYSIS):
			ax.plot( freq_array, np.power(10, spike.waveform), 'r-', linewidth = 0.4 )
			ax.set_xticks( np.logspace(np.log10(low), np.log10(high), num=10))
			ax.set_xticklabels([ "%.2f" % el for el in ax.get_xticks()])
		else:
			ax.plot( freq_array, np.power(10, spike.waveform), 'r-', linewidth = 0.4 )
		fig.savefig( "PSDs/Type_%i/%i-%i.png" % (labels[index]+1, spike.segment_start, spike.segment_end), bbox_inches='tight', pad_inches=0.2)
		plt.close('all')
		del fig
		if not SILENT:
			if ( psd_number > 1 ):
				progress(index)
	del labels_list


def scatterplot(score_matrix, spike_database, colored_clusters_list, labels, x, y, output, ANALYSIS):
	'''
		This function plots a scatterplot of the columns x vs y of the
		score matrix, also generating an image map linked to the corresponding
		observation e.g.:
			observation x links to 'time_series/Type_y/event_gps_time.pdf'
		This way each point in the scatterplot can be clicked to obtain the 
		time series of the clicked observation.
		
	'''
	DPI = 100
	
	
	fig = plt.figure(dpi=DPI, edgecolor='k', frameon=True)#, figsize=(8.15, 6))
	ax = fig.add_subplot(111)
	
	ax.set_xlabel( "Principal component score: " + str(x) )
	ax.set_ylabel( "Principal component score: " + str(y) )
	ax.grid(linestyle = '--')
		
	plotlabels = []
	cluster_number = len( np.unique(labels) )
	ax.set_title( "%i clusters" % cluster_number  )
	
	# Plot points, saving info to 'xs','ys' and 'info' to
	# create the image map
	xs = []
	ys = []
	info_list = []
	for index, element in enumerate( np.array(score_matrix) ):
		xs.append( element[x-1] )
		ys.append( element[y-1] )
		if ( ANALYSIS == "time"):
			info_list.append( ( labels[index]+1, spike_database[index].peak_GPS ) )
		elif ( "frequency" in ANALYSIS):
			info = "%s-%s" % (spike_database[index].segment_start, spike_database[index].segment_end )
			info_list.append( ( labels[index]+1, info) )
		elif ( "generic" in ANALYSIS ):
			pass
		else:
			print "Analysis type '" + ANALYSIS + "' not supported. Quitting."
			sys.exit()
	
	for index, element in enumerate(colored_clusters_list):
		tmp = np.array(element)
		ax.plot( tmp[:,x-1] , tmp[:,y-1], markers_and_colors[index], label = str(index), markersize = 5 )
		plotlabels.append( str(index+1) )
	
	# Create a legend
	ax.legend(plotlabels, bbox_to_anchor=(0, 0, 1.12, 1), loc = 'best', markerscale = 2, numpoints = 1)
	
	# If analysis is generic, no additional steps are needed
	# The image map code does not need to be run. Return.
	if ( "generic" in ANALYSIS ):
		fig.savefig(output + ".png", dpi=fig.get_dpi())
		plt.close('all')
		return
	
	###
	# Finished plotting. Now creating the image map (only for 'time' and
	# 'frequency' ANALYSIS)
	
	
	###
	# Create an array with the x and y coordinates, saved from the previous
	# plots
	xys = zip(xs, ys)
	dpi = fig.get_dpi()
	
	height = fig.get_figheight() * dpi
	
	ixs = [0]*len(xs)
	iys = [0]*len(ys)
	i = 0
	for x1, y1 in xys:
		ixs[i], iys[i] = ax.transData.transform_point( [x1, y1] )
		i += 1
	
	icoords = zip(ixs, iys)
	
	# The minimal 'template' to generate an image map.
	tmpl = """
	<html><head><title> Scatterplot Imagemap - PC scores {0}vs{1}</title></head><body>
	<img src="%s.png" usemap="#points" border="0">
	<map name="points">%s</map>
	</body></html>""".format(x, y)
		
	# Get the correct name for the working directory, starting from public_html
	start = 0
	split_cwd = os.getcwd().split("/")
	for index, element in enumerate(split_cwd):
		if (element == "public_html"):
			start = index+1
			break
	directory = join(split_cwd[start:], "/") + "/"
	
	if ( "time" in ANALYSIS):
		fmt = "<area shape='circle' coords='%f,%f,2' href='time_series/Type_%i/%0.3f.pdf' title='GPS %.3f - Type %i'>"
	else:
		fmt = "<area shape='circle' coords='%f,%f,2' href='PSDs/Type_%i/%s.png' title='%s - Type %i'>"
		
	# need to do height - y for the image-map
	fmts = [fmt % (ix, height-iy, x, y, y, x) for (ix, iy), (x, y) in zip(icoords, info_list) ]
		
	fig.savefig(output + ".png", dpi=fig.get_dpi() )
	plt.close('all')
	
	print >> open(output + ".html", 'w'), tmpl % (output, "\n".join(fmts))
	#print "\tWritten: " + output + " (html and png)"


def chisquare_test(database, labels):
	"""
		Chi square test for the clusters
		
		Takes as input a list of (colored) clusters, labels for the (full) database
		and tests for the goodness of the clustering using a chisquare-based
		test.
	"""
	
	# Create a database for 
	cluster_number = len(np.unique(labels))
	colored_database = [[] for i in range(cluster_number)]
	for index, spike in enumerate(database):
		colored_database[labels[index]].append(spike)
		
		
	# Compute a median for each cluster, used a representative
	# time series for the cluster
	representatives = []
	for cluster in colored_database:
		median = np.median([spike.waveform for spike in cluster], axis=0 )
		representatives.append(median)
	
	# Compute differences between representative and glitches
	# for each cluster
	cluster_differences = [ [] for i in range(cluster_number)]
	for index, cluster in enumerate(colored_database):
		diffs = [representatives[index]-spike.waveform for spike in cluster]
		cluster_differences[index] = diffs
	
	# Compute chisquares, which are simply the sum squares of the difference
	# vectors
	cluster_chi_squares = [ [] for i in range(cluster_number)]
	for index, diffs in enumerate(cluster_differences):
		chi_squares = [(zeta**2).sum() for zeta in diffs]
		cluster_chi_squares[index].extend(chi_squares)
		
	
	# Plot chi squares:
	fig = plt.figure(figsize=(12, 6*cluster_number), dpi=300)
	ax = []
	for index, chi_squares in enumerate(cluster_chi_squares):
		glitch_indexes = range(1, len(chi_squares)+1)
		
		
		tmp = fig.add_subplot(cluster_number, 1, index+1)
		if (len(chi_squares) > 1):
			tmp.plot(glitch_indexes, np.array(chi_squares)/len(chi_squares))
			tmp.scatter(glitch_indexes, np.array(chi_squares)/len(chi_squares))
		else:
			tmp.set_title("Cluster #{0} (1 element)".format(index+1))
			tmp.plot(range(10), range(10))
			tmp.scatter([0], [0])
		
		
		
		tmp.set_title("Cluster #{0}".format(index+1))
		tmp.grid(which="both")
		
		minor_ticks = glitch_indexes
		for i, item in enumerate(minor_ticks):
			if ( item % 5 == 0):
				minor_ticks.pop(i)
		major_ticks = [i for i in range(0, len(chi_squares),5)]
		tmp.set_xticks(major_ticks)
		tmp.set_xticks(minor_ticks, minor=True)
		
		tmp.set_xlim(1, len(chi_squares)+1)
		tmp.set_ylim(0,1.05*np.max(chi_squares)/len(chi_squares))
		tmp.set_xlabel("Glitch Number")
		tmp.set_ylabel("\chi^2/ndf")
		ax.append(tmp)
		
	ax[0].set_title("Chi Squares/ndf: Cluster #{0}".format(index+1))
	fig.savefig("Chisquare.pdf", dpi = DPI, bbox_inches='tight', pad_inches=0.2)
	plt.close('all')
	
	print "\tSaved 'Chisquare.pdf'."
	
	return


####################################################

def main():
	global means, observations, samples, marker
	global STANDARDIZE, PRINT_LIST
	
	args = __check_options_and_args__()
	
	matrix, spike_database = load_data(args, ANALYSIS)
	observations, samples = matrix.shape
	print "Data matrix is %ix%i, %i observations of %i variables" % (observations, samples, observations, samples)
	
	time0 = time.time()
	
	# PCA create a plot of the explained variance vs. principal components 
	# number.
	# The default number of plotted components is 40
	score_matrix, principal_components, means, stds, eigenvalues = PCA(matrix, components_number=40)
	
	print "PCA timing:\t %.2f s" % float(time.time()-time0)
	time1 = time.time()
	
	# gaussian_mixture is called with 
	# 'score_matrix[:,:principal_components_number]'
	# as argument to cluster only using the first 'principal_components_number' 
	# components, as specified by the '-p' option.
	
	print "Clustering using the first %i principal components..." % principal_components_number
	reduced_score_matrix = score_matrix[:,:principal_components_number]
	
	mat, tmp, tmp1 = matrix_whiten(reduced_score_matrix, std=True)
	labels = gaussian_mixture(mat, upper_bound=max_clusters)
	
	cluster_number = len( np.unique(labels) )
	
	
	print "GMM timing:\t %.1f s" % float( time.time() - time1 ) 
	print "GMM algorithm found %i cluster" % cluster_number + \
						('s' if ( cluster_number > 1 or cluster_number == 0) else '')+"."
	
	colored_clusters = color_clusters( score_matrix, labels )
	print_cluster_info(colored_clusters)
	
	
	# DEVELOPMENT TEST 03/25/2014
	# Testing the goodness of the clutering with a
	# chisquare test on the glitches in each cluster
	with warnings.catch_warnings():
		warnings.filterwarnings( "ignore", category=UserWarning )
		chisquare_test(spike_database, labels)
	
	
	'''
	# 3D Plots
	# Uncomment for 3D scatterplots
	x = 1
	y = 2
	z = 3
	three_plot(colored_clusters, x, y, z, "Colored_scatter-"+str(clusternumber)+
				"_clusters"+str(x)+"_vs_"+str(y)+"_vs_"+str(z)+".pdf" )
	'''
	
	output = "Colored_scatter-"+str(cluster_number)+"_clusters_"
	for x in range( 1, MAX_COMPONENTS ):
		for y in range ( 1, MAX_COMPONENTS ):
			if ( x != y ) & ( x < y ):
				scatterplot( score_matrix, spike_database, colored_clusters, labels, x, y, output+str(x)+"_vs_"+str(y), ANALYSIS )
			
		
	
	if ( "bands" in ANALYSIS ):
		calculate_types(spike_database, colored_clusters, score_matrix, principal_components, means, stds, labels, ANALYSIS, SAMPLING, low, high)
	else:
		calculate_types(spike_database, colored_clusters, score_matrix, principal_components, means, stds, labels, ANALYSIS, SAMPLING)
	
	if PLOT:
		if ( "generic" not in ANALYSIS ):
			# If time-domain analysis, create a folder containing the time 
			# series for each
			# of the spikes containing a folder for each cluster
			if ( "time" in ANALYSIS ):
				print "\n\tPlotting time series..."
				if not ( os.path.exists( "time_series" ) ):
					os.mkdir( "time_series" )
				for i in range(cluster_number):
					if not ( os.path.exists( "time_series/Type_"+str(i+1) ) ):
						os.mkdir( "time_series/Type_"+str(i+1) )
				
				spike_time_series(spike_database, (score_matrix, principal_components, means, stds), components_number, labels, SAMPLING, SILEN)
				
				
			elif ( "frequency" in ANALYSIS ):
				print "\n\tPlotting PSDs..."
				# Plot a list of PSDs with the relative times:
				for i in range(cluster_number):
					try:
						os.makedirs( "PSDs/Type_"+str(i+1) )
					except:
						pass
				if ( "bands" in ANALYSIS):
					plot_psds(spike_database, (score_matrix, principal_components, means, stds), components_number, labels, SAMPLING, ANALYSIS, low, high)
				else:
					plot_psds(spike_database, (score_matrix, principal_components, means, stds), components_number, labels, SAMPLING, ANALYSIS)
				
	
	# Cluster removal
	if REMOVE:		
		print "Entering interactive mode..."
		remove_clusters(spike_database, labels)
	print "\n\t\tFinished!"



	

if __name__ == "__main__":
	start = time.time()
	main()
	endtime = float(time.time()-start)
	print "Total Execution: {0:.1f} s".format(endtime if endtime > 0 else endtime/60.0 )
