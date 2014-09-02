# encoding: utf-8
'''
Daniele Trifiro`
brethil@phy.olemiss.edu

Definition of the Spike Class:
The Spike class is used to store the information about a spike, such as amplitude, 
	 peak position, start and end, so the original data can be retrieved:
		spike=Spike() 	
		spike.start			->	Starting point
		spike.end			->	Ending point
		spike.peak			->	Peak point
		spike.norm			->	Normalization constant
		spike.width			->	Width (seconds)
		spike.waveform		->	Waveform
		spike.segment_start	->	Segment start in GPS time
		spike.segment_end	->	Segment end in GPS time
		spike.peak_GPS		->	Peak in GPS time
		spike.sampling		->	Sampling Frequency
		spike.SNR			->	SNR
		spike.type			->	Type
		spike.peak_frequency->	Peak frequency
	
	import inspect
	variables = [i for i in dir(t) if not inspect.ismethod(i)
	
'''

import matplotlib.mlab
import numpy as np

class Spike:
	
	def __init__(self, start, end, peak, norm, peak_GPS, segment_start, segment_end, waveform, sampling_frequency):
		self.start = start
		self.end = end
		self.peak = peak
		self.norm = norm
		self.width = end-start
		self.peak_GPS = peak_GPS
		self.segment_start = segment_start
		self.segment_end = segment_end
		self.waveform = waveform
		self.len = len(self.waveform)
		self.sampling = sampling_frequency
		self.duration = self.width/self.sampling
			
	def __str__(self):
		'''
		String reprentation of the spike (segment_start, segment_end, spike_start, 
		spike_end, spike_peak_GPS).
		'''
		return str(self.segment_start)+"\t"+str(self.segment_end)+"\t"+str(self.start)+"\t"\
				+str(self.end)+"\t"+str(self.peak_GPS)+"\t"
	