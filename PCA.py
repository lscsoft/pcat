#!/usr/bin/env python
# encoding: utf-8


# uncomment for 3D plots
#from mpl_toolkits.mplot3d import axes3d

from utilities_PCAT import *
STANDARDIZE_COLUMNS = False

global DPI
DPI = 300
global frame_width
frame_width = 80


def __usage__():
	print "Usage:   PCA.py --time || --frequency [-s] [--extras number] "
	print "\t\t [...] file1 file2 file3"
	print ""
	print "#"*frame_width
	print "Options:"
	print "\t--time, --frequency\n\t\t Type of analysis being performed, either in time"
	print "\t\t or in frequency.\n\t\t For most (generic) databases one can use the "
	print "\t\t \"--time\" option, \"--frequency\" should be used for fourier"
	print "\t\t transform databases or when the input spans several orders"
	print "\t\t of magnitude. Base 10 logarithm is applied in this case before"
	print "\t\t performing PCA."
	print "\t[--compare file]\n\t\t PCA is performed on the difference between observations"
	print "\t\t and input model (plain-text)."
	print "\t[-s file] or [--scores file]\n\t\t Reads scores matrix form a pickled file."
	#print "\t-o file, --output file\n\t\t Use 'file' as output, output will be 'file.svg'."
	print "\t[-D, --display]\n\t\t Shows the plot GUI."
	print "\t[--startspike start --endspike end]\n\t\t Plots ony the spikes between 'start' and 'end'."
	#print "\t--standardize\n\t\t Sets unit variance for every variable (columns of the data matrix)."
	#print "\t\t (Disabled by default)"
	print "\t\tand a plot of the explained variance vs number of components."
	print "\t[--plot number]\n\t\t Number of principal components scores to be plotted"
	print "\t\t in the scatterplots, default up to the 3rd component."


def __check_options_and_args__():
	'''
	This function parses options and arguments.
	'''
	# Definitions:
	global PICKLED_SCORES 
	PICKLED_SCORES = False
	
	global CUSTOM_COMPONENTS
	CUSTOM_COMPONENTS = False
	global CUSTOM_RANGE
	CUSTOM_RANGE = False
	
	global DISPLAY
	DISPLAY = False
	global SILENT
	SILENT = False
	
	global extras_n
	extras_n = 40
	global components_number
	components_number = 4
	
	global plot_number
	plot_number = 3
	global marker
	CUSTOM_MARKER = False
	
	global ANALYSIS, MODEL
	compare = False
	
	global HISTOGRAM, bins
	HISTOGRAM = False
	
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hsx:y:o:Dm:", ["help", "scores", "xaxis=",\
		 			"yaxis=", "output=", "display", "startspike=", "endspike=", "suppress",\
		 			'extras=', 'marker=', 'plot=', 'time', 'frequency', 'compare=',\
					'histogram='])
	except getopt.error, msg:
		print msg
		sys.exit(1)
	
	if ( len(sys.argv[1:]) == 0 and not (any( flag in o for flag in [ '--help', '-h'] for o in opts)) ) :
		print "No arguments."
		__usage__()
		sys.exit(0)
	
	for o, a in opts:
		if o in ( '-h', '--help' ):
			__usage__()
			sys.exit(1)
		elif o in ('--startspike', '--endspike'):
			CUSTOM_RANGE = True
			if o == '--startspike':
				global startspike
				startspike = int(a)
			else:
				global endspike
				endspike = int(a)
		elif o in ( '-s' , '--scores'):
			PICKLED_SCORES = True
		elif o in ( '-o', '--output' ):
			CUSTOM_OUTPUT = True
			output = str(a)
		elif o in ( '-D', '--display' ):
			DISPLAY = True
		elif o == ( '--extras' ):
			extras_n = int(a)
		elif o in ( '-m', '--marker'):
			CUSTOM_MARKER = True
			marker = a
		elif o in ( '--plot' ):
			components_number = int(a)+1
		elif o in ( '--time' ):
			ANALYSIS = 'time'
		elif o in ( '--frequency' ):
			ANALYSIS = 'frequency'
		elif o in ( '--compare' ):
			compare = True
			MODEL = a
		elif o in ( '--histogram' ):
			HISTOGRAM = True
			bins = int(a)
		else:
			assert False, "Unknown option."
	if compare:
		ANALYSIS = ANALYSIS+"_diff"
	if not ( any( flag in o for flag in [ '--time', '--frequency'] for o in opts) ):
		print "Type of analysis being performed, either \"--time\" or \"--frequency\","
		print "has to be supplied. Exiting."
		sys.exit(1)
	if not CUSTOM_MARKER:
		marker = '+'
	return args


def scatterplot(matrix, x, y, startspike=0, endspike=0, marker = '+'):
	'''
	This function generates a scatterplot of the input matrix, using columns x-1 and y-1
	If DISPLAY == True then the plot is shown in the GUI.
	'startspike' and 'endspike' represents the 'slice' of observations (rows) that will be
	plotted
	'''
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.grid( True, linestyle = '--' )
	ax.set_xlabel( "Principal component score: "+str(x) )
	ax.set_ylabel( "Principal component score: "+str(y) )
	if ( startspike != 0 and endspike != 0):
		ax.plot(matrix[startspike-1:endspike, x-1], matrix[startspike-1:endspike, y-1], 'k'+marker, markersize = 5 ) 
	else:
		ax.plot(matrix[:,x-1], matrix[:,y-1], 'k'+marker, markersize = 5)
	if DISPLAY:
		print "Displaying plot..."
		plt.show()
	plt.close('all')
	return fig


def standardize(matrix):
	'''
	This function standardizes the input matrix. 
	Column mean is removed, so to have zero-mean columns.
	
	Returns a tuple with the standardized matrix and removed means.
	
	'''
	# Verified 05/11/2013
	
	means = np.mean(matrix, axis=0)
	for index, element in enumerate( matrix.transpose() ):
	# Subtract the mean from each column.
		element -= means[index]
	return matrix, means

def matrix_whiten(data_matrix, std=False):
	"""
	This function whitens the input matrix:
	Each column of the matrix is transformed to have
	zero mean, and, if std=True, unit variance.
	"""
	means = np.mean( data_matrix, axis=0 )
	stds = np.std( data_matrix, axis=0 )
	# Loop over columns
	for index, column in enumerate( data_matrix.transpose() ):
		if std:
			column /= stds[index]
		column -= means[index]
	
	if not std:
		stds = np.ones(len(means))
		
	return data_matrix, means, stds


def eigensystem(matrix):
	"""
	Returns eigenvalues and eigenvectors of the input matrix (symmetric)
	sorted by decreasing eigenvalue absolute value
	Also fix the eigenvector's phase. 
	Convention: positive first nonzero component. This fixes the phase.
	
	
	"""
	
	# Verified 05/11/2013
	
	eigenvalues, eigenvectors = np.linalg.eigh( matrix )
	# Sort eigenvalues by order of decreasing absolute value
	idx = eigenvalues.argsort()
	eigenvalues = eigenvalues[idx[::-1]]
	eigenvectors = eigenvectors[:,idx[::-1]]
	index = 0
	# Sets phase, so the sign is fixed.
	for eigenvector in eigenvectors.transpose():
		if ( eigenvector[0] < 0 ):
			eigenvector *= -1
		elif ( eigenvector[0] == 0 ):
			# If the first component is zero, check the next one
			# until a non-zero component is found
			try:
				while ( eigenvector[index] == 0 ):
					index += 1
				if ( eigenvector[index] < 0 ):
					eigenvector *= -1
			except IndexError:
				pass
	
	return eigenvalues, eigenvectors


def PCA(matrix, components_number=50):
	'''
	This function performs PCA on the input matrix
	
	A couple plots are also generated:
		- explained variance in function of the components number
		- Latent roots (eigenvalues) values in function of their index
    
	
	Returns:
	(score_matrix, principal_components, column_means, std_dev, eigenvalues)
	'''
	
	# matrix_whiten removes means and if std=True divides by the standard 
	# deviation
	# output matrix has columns with zero means (if std=True then it has also
	# unit variance in the columns)
	matrix, means, stds = matrix_whiten(matrix, std=STANDARDIZE_COLUMNS)
	
	covariance_matrix = np.cov( matrix.transpose() )
	eigenvalues, principal_components = eigensystem( covariance_matrix )
    
    
	if ( components_number != 0 ):	
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_title("Explained variance")
		ax.set_xlabel("Component number")
		ax.set_ylabel("Variance explained")
		if components_number < 50:
			components_number = 50
			ax.set_xticks(range(0, 51, 5))
		ax.grid( True, linestyle = '--' )
		total = np.sum(eigenvalues)
		variances = np.cumsum(eigenvalues[:components_number])
		# Plot twice, once for the line, once for the markers
		ax.plot( range(1, components_number+1), np.divide( variances, float(total) ), "b-" )
		ax.plot( range(1, components_number+1), np.divide( variances, float(total) ), "bo", markersize=3)
		plt.ylim( (0, 1) )
		fig.savefig("explained-variance.png", dpi=DPI)
		print "\tSaved explained-variance.png"
		plt.close()
		fig = plt.figure()
		ax = fig.add_subplot(111)
		#ax.hist(np.abs(eigenvalues), np.sqrt(len(eigenvalues)), log=True)
		ax.plot(range(components_number), eigenvalues[:components_number])
		ax.plot(range(components_number), eigenvalues[:components_number], "bo", markersize=3)
		ax.plot(range(components_number), np.ones(components_number))
		ax.grid(which="both")
		ax.set_title("Latent Roots - {0} above threshold (1.0)".format((eigenvalues>1.0).tolist().count(True)))
		ax.set_xlabel("Index")
		ax.set_ylabel("Value")
		fig.savefig("Latent_root_distribution.png", bbox_inches='tight', pad_inches=0.2)
		print "\tSaved Latent_root_distribution.png"
		plt.close()
	scores_matrix = np.dot( matrix, principal_components )
		
	return np.array(scores_matrix), np.array(principal_components), np.array(means), np.array(stds), np.array(eigenvalues)


def histogram(data, component_N, bins):
	'''
	Plots an histogram of the component_N column of the data matrix
	'data'. 'bins' is the number of bins you want.
	This is useful if data is the scores matrix, to be able to see
	the distribution of the values of the k-th principal component 
	score.
	'''
	for i in range(0, component_N):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.grid(ls = '--')
		ax.set_title("Histogram for principal component "+str(i+1))
		ax.hist(data[:,i], bins)
		output = "Histogram_component_"+str(i+1)+".svg"
		fig.savefig(output, dpi = DPI)
		print "Written:\t{0}.".format(output)	


def load_data(file_list, pickled, ANALYSIS="time"):
	'''Loads the database (pickled) and returns a matrix with the observations
	in the rows and the loaded list of observations.
	
	ANALYSIS sets how data is loaded, depending on the input options, and
	whether a model is given (--compare):
	- ANALYSIS = 'time'				---> Loads data from list of Spike() class objects
	 									using the 'waveform' attribute and
	 									normalizing by the 'norm' attribute
	- ANALYSIS = 'time_diff'		---> Loads the difference between the model and
	 									the observation, normalized by the 'norm'
										attribute.
	- ANALYSIS = 'frequency'		---> Loads the data from a list of Spike() class objects
	 									'waveform' attribute and takes base-10 logarithm.
	- ANALYSIS = 'time_diff'		---> Loads the difference between the model and
										the base-10 logarithm of the observation.			
	- ANALYSIS = 'generic_log'		---> Loads the data and takes base-10 logarithm.
	- ANALYSIS = 'generic'		---> Loads the data and takes base-10 logarithm.
	'''
	
	database = list()
	if pickled:
		for element in file_list:
			f = open(element, "rb")
			database.extend( pickle.load(f) )
			f.close()
	else:
		for element in file_list:
			database.extend( np.loadtxt(element) )
	# Retrieve the waveforms, normalizing them by spike.norm
	# If called with the --frequency option, logarithm is applied to
	# the normalized waveform.
	# If called with the --frequency_compare option, PCA is performed on the 
	# difference
	# between the input model and the obsevations.
	waveforms = list()
	if ( "diff" in ANALYSIS ):
		model_waveform = np.loadtxt(MODEL)
	for observation in database:
			if ( ANALYSIS == 'time' ):
				waveforms.append( (observation.waveform)/(observation.norm) )
			elif ( ANALYSIS == 'time_diff' ):
				waveforms.append( np.abs( model-( (observation.waveform)/(observation.norm) ) ) )
			elif ( ANALYSIS == 'frequency' ):
				waveforms.append( np.log( (observation.waveform)/(observation.norm) )/np.log(10) )
			elif ( ANALYSIS == 'frequency_diff' ):
				waveforms.append( np.log(model)/np.log(10)-np.log( (observation.waveform)/(observation.norm) )/np.log(10) )
			elif ( ANALYSIS == "generic_log"):
				waveforms.append( np.log(observation.waveform)/np.log(10) )
			elif ( ANALYSIS == "generic"):
				waveforms.append( observation.waveform )
	data = np.array(waveforms)
	return data, database



def generate_plots(scores, components=3):
	''' Plots scatterplots for the given score matrix, if option is \
	given, histograms are also generated'''
	if CUSTOM_RANGE:
		for x in range(1, components):
			for y in range (1, components):
				if ( x != y ) & ( x < y ):
					output = "Scatterplot_"+str(x)+"vs"+str(y)+".pdf"
					fig = scatterplot(scores, x, y, startspike = startspike,\
					 				endspike = endspike, marker=marker)
					fig.savefig(output)
					print "Saved "+str(output)
	else:
		for x in range(1,components):
			for y in range (1, components):
				if ( x != y ) & ( x < y ):
					output = "Scatterplot_"+str(x)+"vs"+str(y)+".pdf"
					fig = scatterplot(scores, x, y, marker=marker)
					fig.savefig(output)
					print "Saved "+str(output)
	if ( HISTOGRAM ):
		for x in range(1, components):
			histogram(scores, x, bins)
	




def main():
	args = __check_options_and_args__()
	if ( len(args) == 0 ) and not ( PICKLED_SCORES ):
		print "No input files."
		sys.exit(0)
		
	# Load and analyze data
	if not ( PICKLED_SCORES ):
		data, observations = load_data(args, ANALYSIS)
		rows, columns = np.shape(data)
		print "Input is a", str(rows)+"x"+str(columns), "matrix."
		print "Performing PCA..."
		scores, principal_components, means, stds, eigenvalues = PCA(data, extras_n)
		
		# Dump score matrix and principal components matrix
		f = open("score_matrix.data", "wb")
		pickle.dump(scores, f )
		f.close()
		f = open("principal_components_matrix.data", "wb")
		pickle.dump(principal_components, f)
		f.close()
		print "Saved scores matrix and principal components matrix as:"
		print "score_matrix.data\t\t(pickled)"
		print "principal_components_matrix.data\t\t(pickled)"
		
	else:
		# A pickled score_matrix has been supplied, we just have to load the 
        # data.
		f = open("score_matrix.data", "rb")
		scores = pickle.load( f )
		f.close()
		print "Input is a", str(np.shape(scores)[0])+"x"+str(np.shape(scores)[1]), "scores matrix"
	
	# Generate plots:
	generate_plots(scores, components_number)




if __name__ == '__main__':
	total = time.time()
	main()
	print "Total execution time: ", time.time()-total
