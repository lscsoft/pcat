#!/usr/bin/env python
# encoding: utf-8
'''
Daniele Trifiro`
brethil@phy.olemiss.edu

download_frames.py is used do retrieve time series from frame files (.gwf).
pyl
Can be used standalone. For usage: run with -h.


Contains:
	- download_interval()
	- retrieve_timeseries()
'''

import os
import sys

import getopt

from glue import lal
import pylal.frutils


from utilities_PCAT import progressBar

def usage():
	print "Usage:\n\tdownload_data.py -s start_time -e end_time -c channel -I IFO -f frame_type "
	print "\t\t\t[--size segment_size] [--pad padding_seconds]"
	
	print "\n\tRetrieve and save (in the current folder) time series for the given"
	print "\ttime interval/channel."
	print "\tOutput name is 'IFO-frame_type-channelname-starttime-endtime.data'"
	print "\tOutput type can be used with python/numpy:"
	print "\t\t> import numpy as np"
	print "\t\t> f = open('filename.data', 'rb')"
	print "\t\t> time_series = np.load(f)"
	print "\t\t> f.close()"
		
	print ""
	
	print "\tThe segments are split into 'segment_size' chunks\n\
	(Default segment_size = 60 seconds)."
	
	print "Options:"
	
	print "\t -s time, --start time\n\
		Start time in GPSTIME"
	
	print "\t -e time, --end time\n\
		End time in GPSTIME"
	
	print "\t -c channel, --channel channel\n\
		Channel name."
	
	print "\t -I ifo, --IFO ifo\n\
		IFO 'L' for LLO or 'H' for LHO)"
	
	print "\t -f type, --frame type\n\
		Data type. Usually 'R' for Livingston, 'H1_R' for Hanford.\n\
		If not sure check with ligo_data_find."
	
	print "\t --size segment_size\n\
	 	If time interval is bigger than segment_size,\n\
		split the given times into chunks long 'segment_size' seconds.\n\
		(Default = 60)"	
	
	print "\t --pad padding_seconds\n\
		Overlap time in (integer) seconds of consecutive segments.\n\
	 	(Default = 0)."
	

def check_options_and_args():
	global startime, endtime
	global IFO, datatype, channel
	global segment_size, overlap
	overlap = 0
	segment_size = 60
	datatype = 'R'
	
	if ( len(sys.argv[1:]) == 0 ):
		print "No arguments given."
		usage()
		sys.exit(1)
	else:
		try:
			opts, args = getopt.getopt(sys.argv[1:], "hs:d:e:c:I:f:",\
			 							["help", "start=", "end=",\
			 							 "channel=", 'IFO=', 'size=', 'pad=', 'frame='])
		except getopt.GetoptError, err:
			print str(err)
			sys.exit(1)
		for o, a in opts:
			if o in ( '-h', '--help'):
				usage()
				exit()
			elif o in ( '-s', '--start'):
				startime = int(a)
			elif o in ( '-e', '--end'):	
				endtime = int(a)
			elif o in ( '-c', '--channel'):
				channel = a
			elif o in ( '-I', '--IFO'):
				IFO = a
			elif o in ( '-f', '--frame', '-d'):
				datatype = a
			elif o in ( '--size' ):
				segment_size = int(a)
			elif o in ( '--pad' ):
				overlap = int(a)
			else:
				assert False, "Unknown Option"
		# Check for all the needed options
		if not ( any( flag in o for flag in ['-s', '--start'] for o in opts )\
				and any( flag in o for flag in ['-e', '--end'] for o in opts)\
				and any( flag in o for flag in ['-c','--channel'] for o in opts)\
				and any( flag in o for flag in ['-I','--IFO'] for o in opts) \
				and any( flag in o for flag in ['-d', '--frame', '-f'] for o in opts) ):
			print "Incorrect usage. Try again supplying:"
			print "-s or --start\t GPS start time"
			print "-e or --end\t GPS end time"
			print "-c or --channel\t Channel name"
			print "-I or --IFO\t IFO"
			print "-f o --frame\t Frame Type"
			#usage()
			exit()


def download_interval(start_time, end_time, step, overlap, channel_name, IFO, frame_type ):
	'''
		Read frame files and returns time series of the given interval/channel.
		
		The segments splitted into chunks 'step+2*overlap' seconds long
		and overlapping with neighouring segments by 'overlap'. 
		
		Downloaded files are plain text and saved in the form:
			'IFO-frame_type-channel_name_start_time_end_time.data'
		
		If file name already exists, the segment is not downloaded.
		
		Arguments:
			- start_time (integer)
					Start time in GPS time
			- end_time (integer)
					End time in GPS time
			- step (integer)
					Segment size in seconds
			- overlap (integer)
					Overlap between segments in seconds
			- channel_name (string)
					Name of the channel e.g. L1:LSC-DARM_ERR
			- IFO (string)
					IFO: 'L' for Livingston, 'H' for Hanford
			- frame_type (string)
					Usually 'R' for Livingston and 'H' for hanford.
			 		Use ligo_data_find to find out the frame type.
		 
		Returns a list of names of the downloaded files,
		
		
	'''
	
	start1 = int(start_time) - overlap
	end1 = start1 + step + 2*overlap
	error_count = 0
	downloaded_list = []
	while ( end1 < int(end_time)+2*overlap ):
		start1 = end1-2*overlap
		end1 = start1+step+2*overlap
		try:
			time_series = retrieve_timeseries(start1, end1, channel_name, IFO, frame_type)['waveform']
		except Exception, error:
			print "Error retrieving time series ({0})".format(str(error))
			print "Start, end:\t{0}-{1}".format(start1, end1)
			print "IFO:\t\t{0}".format(IFO)
			print "Frame Type:\t{0}".format(frame_type)
			print "Channel\t{0}".format(channel_name)
			exit()
		
		out_name = "{0}-{1}-{2}_{3}{4}.data".format(IFO, frame_type, channel_name, start1, end1)
		try:
			f = open(out_name, "wb")
			np.save(out_name)
			f.close()
			downloaded_list.append(out_name)
		except Exception, error:
			print "Error saving time series '{0}', error: {1}".format(out_name, error)
			
	return downloaded_list



def retrieve_timeseries(start_time, end_time, channel_name, IFO, frame_type):
	"""
		Read the frame files associated to the input parameters and returns
		a the time series for the given time interval.
		
		Arguments:
			- start_time (integer)
					Start time in GPS time
			- end_time (integer)
					End time in GPS time
			- step (integer)
					Segment size in seconds
			- overlap (integer)
					Overlap between segments in seconds
			- channel_name (string)
					Name of the channel e.g. L1:LSC-DARM_ERR
			- IFO (string)
					IFO: 'L' for Livingston, 'H' for Hanford
			- frame_type (string)
					Usually 'R' for Livingston and 'H' for hanford.
			 		Use ligo_data_find to find out the frame type.
		
		Output:
			- time_series, dictionary with attributes:
				@'waveform' time_series (numpy array)
				@ 'fs' sampling frequency (float)
				@'dt' sampling period, 1/fs (float)
	"""
	
	with open("/home/brethil/classification_paper/caches/MDC3.cache") as cachefile:
		cache = lal.Cache.fromfile(cachefile)
	d = pylal.frutils.FrameCache(cache)
	
	data = d.fetch(channel_name, start_time, end_time)
	
	time_series = {
		'waveform': data.A,
		'dt'      : data.metadata.dt,
		'fs'      : 1.0/data.metadata.dt,
	}
	return time_series

def main():
	# Check arguments:
	check_options_and_args()
	
	# Download data from startime to endtime, splitting into 'segment_size' size
	# segments (in seconds) overlapping by 'overlap'
	download_interval(startime, endtime, segment_size, overlap, channel, IFO, datatype)
	print ""

if __name__ == '__main__':
	main()