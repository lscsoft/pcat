#!/usr/bin/env python
# encoding: utf-8
"""
utilities_PCAT.py

Contains various PCAT definitions and functions.


Daniele Trifirò
brethil@phy.olemiss.edu

2013-08-17

"""


DPI = 300

import os, sys
import getopt
import time
from string import join
from itertools import cycle
from commands import getstatusoutput
import warnings


import numpy as np
import cPickle as pickle

import matplotlib
import matplotlib.colors
matplotlib.use('agg')
import matplotlib.pyplot as plt


from scipy.signal import cheby1, firwin, lfilter, resample

from spike_class import *

"""def get_terminal_size():
	'''
		Gets and returns terminal weidth and height as a tuple
	'''
	env = os.environ
	def ioctl_GWINSZ(fd):
		try:
			import fcntl, termios, struct
			cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
			'1234'))
		except:
			return None
		return cr
	cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
	if not cr:
		try:
			fd = os.open(os.ctermid(), os.O_RDONLY)
			cr = ioctl_GWINSZ(fd)
			os.close(fd)
		except:
			pass
	if not cr:
		try:
			cr = (env['LINES'], env['COLUMNS'])
		except:
			cr = (25, 80)
	return int(cr[1]), int(cr[0])


frame_width, height = get_terminal_size()
"""
frame_width = 79

# Markers and colors, supports up to 21 different types, 
# add more if needed, in the form 'markercolor'
# '+r' means '+' marker and red color ('r')
#
# If more types are needed, simply add more elements to the
# array defined below, the allowed markers and colors are
# listed in the matplotlib.pyplot documentation

global sprog
s = ["[ | ]", "[ / ]", "[ - ]", "[ \\ ]"]
sprog = cycle(s)


global marker_and_colors
#markers_and_colors = ['+k', '+r', '+g', '+b', '+c', '+m', '+y',\
#					'.k', '.r', '.g', '.b', '.c', '.m', '.y']
# Use points (first line), crosses (second line), squares (third line)
# colors are: black, red, green, blue, cyan, magenta , yellow
markers_and_colors = ['.k', '.r', '.g', '.b', '.c', '.m', '.y',\
					'+k', '+r', '+g', '+', '+c', '+m', '+y',\
					'sk', 'sr', 'sg', 'sb', 'sc', 'sm', 'sy']
#TODO: This will have to be changed: a quick way to do this is the following
# t = np.linspace(0,1,things_to_plot)
# cmap = plt.get_cmap("jet")   # or an other colormap
# colors = [cmap[t[i]] for i in range(things_to_plot) ]
# the way markers are plotted should also be changed in all of the plotting
# functions (only one marker is now needed!!)

class progressBar:
	""" Creates a text-based progress bar. Call the object with the `print'
		command to see the progress bar, which looks something like this:
		
		[=======>        22%                  ]
		
		You may specify the progress bar's width, min and max values on init.
		
		
		Example of usage
		
		N = 100000
		bar = progressBar(minValue = 0, maxValue = N, totalWidth = 60)
		for i in range(0,N):
			for j in range(0,N/100):
				j = i**2
			bar(i+1)
		print "\n"
		'''
	"""
		
	def __init__(self, minValue = 0, maxValue = 100, totalWidth=80):
		self.progBar = "[]"   # This holds the progress bar string
		self.min = minValue
		self.max = maxValue
		self.span = maxValue - minValue
		self.width = totalWidth
		self.amount = 0       # When amount == max, we are 100% done
		self.updateAmount(0)  # Build progress bar string
	
	
	def updateAmount(self, newAmount = 0):
		""" Update the progress bar with the new amount (with min and max
			values set at initialization; if it is over or under, it takes the
			min or max value as a default. 
		"""
		if newAmount < self.min:
			newAmount = self.min
		if newAmount > self.max:
			newAmount = self.max
		self.amount = newAmount
			
		# Figure out the new percent done, round to an integer
		diffFromMin = float(self.amount - self.min)
		percentDone = (diffFromMin / float(self.span)) * 100.0
		percentDone = int(round(percentDone))
			
		# Figure out how many hash bars the percentage should be
		allFull = self.width - 2
		numHashes = (percentDone / 100.0) * allFull
		numHashes = int(round(numHashes))
			
		# Build a progress bar with an arrow of equal signs; special cases for
		# empty and full
		if numHashes == 0:
			self.progBar = "\t\t[>%s]" % (' '*(allFull-1))
		elif numHashes == allFull:
			self.progBar = "\t\t[%s]" % ('='*allFull)
		else:
			self.progBar = "\t\t[%s>%s]" % ('='*(numHashes-1),' '*(allFull-numHashes))
			
		# figure out where to put the percentage, roughly centered
		percentPlace = (len(self.progBar) / 2) - len(str(percentDone))
		percentString = str(percentDone) + "%"
			
		# slice the percentage into the bar
		self.progBar = ''.join( [self.progBar[0:percentPlace], percentString, self.progBar[percentPlace+len(percentString):]] )
	
		
	def __str__(self):
		return str(self.progBar)
	
	
	def __call__(self, value):
		""" Updates the amount, and writes to stdout. Prints a carriage return
			first, so it will overwrite the current line in stdout.
		"""
		print '\r',
		self.updateAmount(value)
		sys.stdout.write(str(self))
		sys.stdout.flush()
		print '\r',
	


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
	return server


def decimate(x, q, n=None, ftype='iir', axis=-1):
	'''
	This code is copied from scipy.signal
		
	Downsample the signal by using a filter.
	By default, an order 8 Chebyshev type I filter is used.  A 30 point FIR
	filter with hamming window is used if `ftype` is 'fir'.
		
		Parameters
	----------
	x : ndarray
		The signal to be downsampled, as an N-dimensional array.
	q : int
		The downsampling factor.
	n : int, optional
		The order of the filter (1 less than the length for 'fir').
	ftype : str {'iir', 'fir'}, optional
		The type of the lowpass filter.
	axis : int, optional
		The axis along which to decimate.
		Returns
	-------
	y : ndarray
		The down-sampled signal.
		See also
		--------
	scipy.resample
	'''
		
	if not isinstance(q, int):
		raise TypeError("q must be an integer")
		
	if n is None:
		if ftype == 'fir':
			n = 30
		else:
			n = 8
			
	if ftype == 'fir':
		b = firwin(n + 1, 1. / q, window='hamming')
		a = 1.
	else:
		b, a = cheby1(n, 0.05, 0.8 / q)
		
	y = lfilter(b, a, x, axis=axis)
	sl = [slice(None)] * y.ndim
	sl[axis] = slice(None, None, q)
	return y[sl]


def pickle_dump(data, file_name):
	'''
		Saves 'data' as a binary file named 'file_name'
		The binary file is pickled.
	'''
	f = open(file_name, "wb")
	pickle.dump(data, f)
	f.close()


def load_data(args, pickled):
	'''
		Loads given list of file names. 'pickled' is boolean,
		if True, input data is supposed to be pickled, and loaded
		through the cPickle module, otherwise input is supposed
		plain text and np.loadtxt is used.
		
		Returns a list with the loaded data.
	'''
	input_data = list()
	if pickled:
		for index, element in enumerate(args):
			f = open(element, "wb")
			input_data.append( pickle.load(f) )
			f.close()
	else:
		for index, element in enumerate(args):
			input_data.append( np.loadtxt(element) )
	return input_data


def load_time_series(file_name, pickled=False):
	''' DEPRECATED
		Given a file name, loads the file, using numpy's loadtxt
		if pickled=False, using cPickle if pickled=True.
		
		Returns the time series.
	'''
			
	if pickled:
		f = open(file_name, "rb")
		time_series = pickle.load(f)
		f.close()
	else:
		time_series = np.loadtxt(file_name)
	return time_series


def create_data_matrix(data_list, ANALYSIS, model_waveform=None):
	''' 
		Takes as input a list of Spike instances and type of analysis being performed:
		'time', 'frequency', 'time_diff', 'frequency_diff'
		
		The last two require a Model file to be supplied (array)
		
		Retrieve the waveforms, normalizing them by spike.norm
		If ANALYSIS = "frequency", base 10 logarithm of the input data is taken.
		
		Returns a matrix with the observations in the rows and the variables in the columns.
		
	'''
	waveforms = []
	for observation in data_list:
			if ( ANALYSIS == 'time' ):
				waveforms.append( (observation.waveform)/(observation.norm) )
			elif ( ANALYSIS == 'time_diff' ):
				waveforms.append( np.abs( model_waveform-( (observation.waveform)/(observation.norm) ) ) )
			elif ( ANALYSIS == 'frequency' ):
				waveforms.append( np.log10( (observation.waveform)/(observation.norm) ))
			elif ( ANALYSIS == 'frequency_diff' ):
				waveforms.append( np.log10(model_waveform)-np.log10( (observation.waveform)/(observation.norm) ) )
	data_matrix = np.array(waveforms)
	return data_matrix
	

	
def create_data_matrix_from_psds(file_name_list, ANALYSIS, f_sampl, low=None, high=None):
	'''
		Creates a list and a matrix from the input file names list
		
		ANALYSIS is a string. If it contains 'bands', only a slice of the loaded
		input file is taken as input:
			x[low:high]
		
		This is used for band analysis in frequency-domain PCAT.
		
	'''
	database = []
	waveforms = []
	
	if ("bands" in ANALYSIS):
		tmp = open(file_name_list[0])
		tmp1 = np.load(tmp)
		waveform_length = tmp1.size
		tmp.close()
		freq_array = rfftfreq( 2*(waveform_length-1), d=1./f_sampl )
		# We can now select the correct values from the spectra
		# by choosing only the values defined by the mask
		mask = (freq_array>=low) & (freq_array<=high)
	for index, element in enumerate(file_name_list):
		f = open(element, "rb")
		tmp = np.log10(np.load(f))
		f.close()
		# Cut the spectra to the defined interval if required:
		if ( "bands" in ANALYSIS ):
			tmp = tmp[mask]
		waveforms.append(tmp)
		
		# This part is used to extract the name of the segment from the file 
		# name the if is needed to correctly hand cases when
		# files are in the same folder and when they are not (i.e. path to file 
		# is given)
		if ( "/" in element ):		# Different folder, Path to the file is 
									#given
			( start, end ) = element.split("/")[-1].split('.')[0].split('_')[-1].split('-')
			start = int(start)
			end = int(end)
		else:					# Same folder, only file name is given
			( start, end ) = element.split('.')[0].split('_')[-1].split('-')
			start = int(start)
			end = int(end)
		
		normalization = 1.0
		spike = Spike(int(start), int(end),
						 0, normalization,
						0, int(start), int(end),
						tmp, f_sampl)
		database.append(spike)
		
	return database, np.array(waveforms)



def nearest_power_of_two(number):
	"Returns the nearest power of two less than 'number'"
	i=1
	while(i*2 < number):
		i*=2
	return i

def tformat(x, start_time):
	"""tick formatter definition"""
	global units
	x -= start_time
	if units == "h":
		x /= 3600.0
	elif units == "m":
		x /= 60.0
	elif units == "d":
		x /= 24.0
	tick = "{0:.1f}{1}".format(x, units)
	return tick


def plot_glitchgram(data, times, start_time, end_time, highpass_cutoff, f_sampl, labels, name="Glitchgram"):
	"""
	Plot a glitchgram of all the glitches in 'data' (a list of Spike() istances)
	
	times is the listof locked times
	
	start_time and end_time are the earliest and the latest GPS times
	
	"""
	DPI = 100
	fig = plt.figure(figsize=(12,3*6), dpi=DPI)
	plt.subplots_adjust(left=0.10, right=0.95, top=0.97, bottom=0.05)
	ax = fig.add_subplot(311, axisbg="gray", alpha=0.05)
	ax2 = fig.add_subplot(312, axisbg="gray", alpha=0.05)
	ax3 = fig.add_subplot(313)#, axisbg="gray", alpha=0.05)
	ax.set_xlabel("Time since GPS time {0}".format(start_time))
	ax.set_ylabel("Frequency [Hz]")
	ax.set_yscale('log')
	global units
	if (end_time-start_time)>4*3600:
		units = "h" # hours
	elif (end_time-start_time)>5*60:
		units = "m" # minutes
	elif (end_time-start_time>3*86400):
		units = "d" # days
	
	time_axis = []
	peak_frequencies = []
	color = []
	color_type = []
	SNRs = []
	
	glitch_number = len(data)
	
	for index, spike in enumerate(data):
		time_axis.append(spike.peak_GPS)
		
		(PSD, freqs) =  spike.psd, spike.fft_freq
	
		central_freq = (np.sum(PSD*freqs))/PSD.sum()
		peak_frequency = freqs[np.argmax(PSD)]
		spike.peak_frequency = peak_frequency
		spike.central_freq = central_freq
		
		peak_frequencies.append(peak_frequency)
		SNRs.append(spike.SNR)
		
	cluster_number = len(np.unique(labels))
	
	loudest_event_index = np.argmax(SNRs)
	ax.set_title("Glitchgram - Loudest at GPS {0}, SNR={1:.1g}, Peak Frequency={2:.1f}Hz, Type={3}".format(data[loudest_event_index].peak_GPS, np.max(SNRs), peak_frequencies[loudest_event_index], labels[loudest_event_index]+1))
	
	# Define ticks and define ax (glitchgram with color according to SNR)
	x_ticks = [start_time]
	x_ticks_labels = ["0"]
	interval = end_time-start_time
	step = max((interval)//15, 1)
	
	range_end = end_time
	for i in range(start_time+step, range_end-step, step):
		x_ticks.append(i)
		x_ticks_labels.append(tformat(i, start_time))
	x_ticks.append(end_time)
	x_ticks_labels.append(tformat(end_time, start_time))
	
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_ticks_labels, fontsize=12)
	
	
	# Color in green the part of the plot corresponding to locked times
	# (defined in the 'times' variable)
	locked_times_plots = []
	locked_times_plots2 = []
	locked_times_plots3 = []
	
	for interval in times:
		# Background color is grey (defined above), locked segments are shown 
		# in white. For green: #008000, or any other color HTML code or any
		# other format accepted my matplotlib
		locked_times_plots.append(matplotlib.collections.BrokenBarHCollection.span_where(
			interval, ymin=highpass_cutoff, ymax=f_sampl, where=np.array(peak_frequencies)>0, facecolor='#FFFFFF', alpha=1))
		
		locked_times_plots2.append(matplotlib.collections.BrokenBarHCollection.span_where(
			interval, ymin=highpass_cutoff, ymax=f_sampl, where=np.array(peak_frequencies)>0, facecolor='#FFFFFF', alpha=1))
			
		"""
		locked_times_plots3.append(matplotlib.collections.BrokenBarHCollection.span_where(
			interval, ymin=0.0, ymax=10.0*np.max(SNRs), where=np.array(SNRs)>0, facecolor='#FFFFFF', alpha=1))
		"""
		
		
		ax.add_collection(locked_times_plots[-1])
		ax2.add_collection(locked_times_plots2[-1])
		#ax3.add_collection(locked_times_plots3[-1])
		
	
	
	# Setup the third plot, the Glitch SNR Distribution
	plot_labels_SNR = []
	try:
		for i in range(cluster_number):
			mask = np.where(labels == i, True, False)
			time_axis_tmp = np.array(time_axis)[mask]
			y_axis_tmp = np.array(SNRs)[mask]
			ax3.plot(time_axis_tmp, y_axis_tmp, markers_and_colors[i], label=str(i+1), markersize=5)
			plot_labels_SNR.append(str(i+1))
		color_legend = ax3.legend(plot_labels_SNR, bbox_to_anchor=(1.12, 1), numpoints = 1)
	except:
		ax3.plot(time_axis, SNRs, "ro", markersize=5)
	
	# PLOT STAR FOR LOUDEST EVENT IN PLOT 3
	ax3.plot(time_axis[loudest_event_index], np.max(SNRs), markers_and_colors[labels[loudest_event_index]][1]+"*", markersize=10)
	ax3.yaxis.grid(which="both")
	ax3.xaxis.grid(which="major")
	ax3.set_title("SNR distribution")
	ax3.set_xlabel("Time since GPS time {0}".format(start_time))
	ax3.autoscale(True, axis="x", tight=True)
	ax3.set_yscale('log')
	ax3.autoscale(True, axis="y", tight=True)
	
	if (np.min(SNRs) >= 1 ):
		y_lim_min = 1
	else:
		y_lim_min = np.min(SNRs)
	ax3.set_ylim( (y_lim_min, np.max(SNRs)*10) )
	ax3.set_xticks(x_ticks)
	ax3.set_xlim((start_time, end_time))
	plt.xticks(x_ticks, x_ticks_labels, fontsize=12)
	ax3.set_ylabel("SNR")	
	# Resize ax3 to have the same dimensions as ax, which is smaller
	# due to the colorbar
	box = ax3.get_position()
	ax3.set_position([box.x0, box.y0, box.width*0.94, box.height])
	ax3.set_xticks(x_ticks)
	ax3.set_xticklabels(x_ticks_labels, fontsize=12)
	
	# Plot the scatterplot are colored and sized according to SNR
	
	sizes = np.log10(SNRs)
	sizes -= np.min(sizes)
	sizes /= np.max(sizes)
	
	# Now sizes is normalized in [0,1], default marker size for scatter
	# is 20, we want all the points to be between 15 and 25:
	sizes *= 20
	sizes += 3
	
	# Plot the scatterplot for the first panel and colorbar
	cax = ax.scatter(time_axis, peak_frequencies, c=SNRs, s=sizes, norm=matplotlib.colors.LogNorm(), cmap=matplotlib.cm.jet, edgecolor="none")
	cbar = plt.colorbar(cax, ax=ax, orientation="vertical", shrink=0.7, fraction=0.05, pad=0.01, spacing="proportional")

	# Fix colorbar labels
	cbar.set_label("SNR", fontsize=10)
	cmin,cmax = cbar.get_clim()
	ticks = np.logspace(np.log10(cmin), np.log10(cmax), 10)
	ticklabels = ['%.1g' % t for t in ticks]
	cbar.set_ticks(ticks)
	cbar.set_ticklabels(ticklabels)
	
	# PLOT STAR FOR LOUDEST EVENT IN PLOT 1
	ax.plot([time_axis[loudest_event_index]], [peak_frequencies[loudest_event_index]], "r*", markersize=10)
	# Plot the last panel, in which the glitches are colored according to Type
	plot_labels = []
	
	try:
		for i in range(cluster_number):
			mask = np.where(labels == i, True, False)
			time_axis_tmp = np.array(time_axis)[mask]
			y_axis_tmp = np.array(peak_frequencies)[mask]
			ax2.plot(time_axis_tmp, y_axis_tmp, markers_and_colors[i], label=str(i+1))
			plot_labels.append(str(i+1))
		box = ax2.get_position()
		ax2.set_position([box.x0, box.y0, box.width*0.94, box.height])
		color_legend = ax2.legend(plot_labels, bbox_to_anchor=(1.12, 1), numpoints = 1)
	except:
		ax2.plot(time_axis, peak_frequencies, "ro")
	
	# PLOT STAR AS LOUDEST EVENT FOR PLOT 2
	ax2.plot(time_axis[loudest_event_index], peak_frequencies[loudest_event_index], markers_and_colors[labels[loudest_event_index]][1]+"*", markersize=10)
	
    #ax2.add_artist(color_legend)
	
	ax2.set_xlim((start_time, end_time))
	ax2.autoscale(False, axis="both")
	ax2.set_ylim((highpass_cutoff, int(1.50*(f_sampl//2)) ))
	ax2.set_yscale('log')
	ax2.set_xlabel("Time since GPS time {0}".format(start_time))
	ax2.set_ylabel("Peak Frequency [Hz]")
	ax2.yaxis.grid(which="both")
	ax2.xaxis.grid(which="major")
	ax2.set_xticks(x_ticks)
	ax2.set_xticklabels(x_ticks_labels, fontsize=12)
	
		
	
	ax.set_xlim((start_time, end_time))
	ax.set_ylim((highpass_cutoff, int(1.50*(f_sampl//2)) ))
	ax.grid(which="both")
	
	
	# Plot points, saving info to 'xs','ys' and 'info' to
	# create the image map
	info_list = []
	for index, element in enumerate(data):
		info_list.append( ( labels[index]+1, data[index].peak_GPS ) )
	
	###
	# Create an array with the x and y coordinates of the points
	xys = zip(time_axis, peak_frequencies)
	xys1 = zip(time_axis, SNRs)
	xys2 = zip(time_axis, peak_frequencies)
	dpi = fig.get_dpi()
	
	
	height = fig.get_figheight() * dpi
	
	ixs = [0]*glitch_number
	iys = [0]*glitch_number
	ixs1 = [0]*glitch_number
	iys1 = [0]*glitch_number
	ixs2 = [0]*glitch_number
	iys2 = [0]*glitch_number
	
	# Get coordinates for the points in the image for each plot
	i = 0
	for x, y in xys:
		ixs[i], iys[i] = ax.transData.transform_point( [x, y] )
		i+=1
	i=0
	for x, y in xys2:
		ixs2[i], iys2[i] = ax2.transData.transform_point( [x, y] )
		i += 1
	i=0
	for x, y in xys1:
		ixs1[i], iys1[i] = ax3.transData.transform_point( [x, y] )
		i += 1	
	icoords = zip(ixs, iys)
	icoords1 = zip(ixs1, iys1)
	icoords2 = zip(ixs2, iys2)
	
	# The minimal 'template' to generate an image map.
	tmpl = """
	<html><head><title>Glitchgram</title></head><body>
	<img src="%s.png" usemap="#points" border="0">
	<map name="points">%s</map>
	</body></html>"""
	
	fmt = "<area shape='circle' coords='%f,%f,3' href='time_series/Type_%i/%0.3f.pdf' title='GPS %0.2f - Type %i ' >"
		
	# need to do height - y for the image-map
	fmts = [fmt % (ix, height-iy, x, y, y, x) for (ix, iy), (x, y) in zip(icoords, info_list) ]
	fmts1 = [fmt % (ix, height-iy, x, y, y, x) for (ix, iy), (x, y) in zip(icoords1, info_list) ]	
	fmts2 = [fmt % (ix, height-iy, x, y, y, x) for (ix, iy), (x, y) in zip(icoords2, info_list) ]	
	
	fig.savefig("{0}.png".format(name), dpi=fig.get_dpi()) # bbox_inches='tight', 
	print "\tSaved: {0}.html".format(name)
	plt.close('all')
	
	f = open("{0}.html".format(name), "w")
	print >> f, tmpl % ("{0}".format(name), "\n".join(fmts+fmts1+fmts2))
	f.close()
	


####### Functions used to parallelize 
# Can't use a simple multiprocessing Pool,
# it fails with a picklingerror
import multiprocessing

def spawn(f):
	def fun(q_in,q_out):
		while True:
			i,x = q_in.get()
			if i == None:
				break
			q_out.put((i,f(x)))
	return fun


def parmap(f, X, nprocs = multiprocessing.cpu_count()):
	q_in   = multiprocessing.Queue(1)
	q_out  = multiprocessing.Queue()
	proc = [multiprocessing.Process(target=spawn(f),args=(q_in,q_out)) for _ in range(nprocs)]
	for p in proc:
		p.daemon = True
		p.start()
		
	sent = [q_in.put((i,x)) for i,x in enumerate(X)]
	[q_in.put((None,None)) for _ in range(nprocs)]
	res = [q_out.get() for _ in range(len(sent))]	
	[p.join() for p in proc]
	return [x for i,x in sorted(res)]
	

##### Stolen from Numpy 1.8, returns the frequencies vector
##### associated to the fourier transform obtained through
#### np.fft.rfft

def rfftfreq(n, d=1.0):
	'''
		-------From Numpy 1.8
		Return the Discrete Fourier Transform sample frequencies
		(for usage with rfft, irfft).
		The returned float array `f` contains the frequency bin centers in cycles
		per unit of the sample spacing (with zero at the start). For instance, if
		the sample spacing is in seconds, then the frequency unit is cycles/second.
		
		Given a window length `n` and a sample spacing `d`::
		
		f = [0, 1, ..., n/2-1, n/2] / (d*n) if n is even
		f = [0, 1, ..., (n-1)/2-1, (n-1)/2] / (d*n) if n is odd
		
		Unlike `fftfreq` (but like `scipy.fftpack.rfftfreq`)
		the Nyquist frequency component is considered to be positive.
		
		Parameters
		----------
		n : int
			Window length.
		d : scalar, optional
		Sample spacing (inverse of the sampling rate). Defaults to 1.
		
		Returns
		-------
		f : ndarray
		Array of length ``n//2 + 1`` containing the sample frequencies.
		
		Examples
		--------
		>>> signal = np.array([-2, 8, 6, 4, 1, 0, 3, 5, -3, 4], dtype=float)
		>>> fourier = np.fft.rfft(signal)
		>>> n = signal.size
		>>> sample_rate = 100
		>>> freq = np.fft.fftfreq(n, d=1./sample_rate)
		>>> freq
		array([ 0., 10., 20., 30., 40., -50., -40., -30., -20., -10.])
		>>> freq = np.fft.rfftfreq(n, d=1./sample_rate)
		>>> freq
		array([ 0., 10., 20., 30., 40., 50.])
			
	'''
	val = 1.0/(n*d)
	N = n//2 + 1
	results = np.arange(0, N, dtype=int)
	return results * val



