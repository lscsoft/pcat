# encoding: utf-8
'''
Daniele Trifiro`
brethil@phy.olemiss.edu

Definition of the Spike Class:
The Spike class is used to store the information about a spike, such as amplitude, 
	 peak position, start and end, so the original data can be retrieved:
	spike=Spike()
	
	For time domain analysis the attributes are:
		spike.start				->	Starting point
		spike.end				->	Ending point
		spike.peak				->	Peak point
		spike.norm				->	Normalization constant
		spike.duration			->	Duration (in s)
		spike.waveform			->	Waveform
		spike.segment_start		->	Segment start in GPS time
		spike.segment_end		->	Segment end in GPS time
		spike.peak_GPS			->	Peak in GPS time
		spike.sampling			->	Sampling Frequency
		spike.SNR				->	SNR
		spike.type				->	Type
		spike.peak_frequency	->	Peak frequency
		spike.central_freq		->	Central Frequency
		spike.fft				->	Fourier Transform of the glitch (dimensional)
		spike.psd				->	Power Spectral Density (PSD) of the glitch (np.abs(fft)**2)
		spike.fft_freq			->	Fourier Transform frequencies
		spike.segment_psd		->	PSD of the (whitened) segment from which the glitch was extracted.
		
	
	For frequency domain analysis:
		spike.segment_start	->	Segment start in GPS time
		spike.segment_end	->	Segment end in GPS time
		spike.waveform		->	PSD for the segment
		spike.sampling		->	Sampling frequency for the original data
	
	
	To get a list of all the defined attributes of an istance "t":
	import inspect
	variables = [i for i in dir(t) if not inspect.ismethod(i)]
	
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
		self.duration = (end-start)*sampling_frequency
        
	def __str__(self):
		'''
		String reprentation of the spike (segment_start, segment_end, spike_start, 
		spike_end, spike_peak_GPS).
		'''
		return str(self.segment_start)+"\t"+str(self.segment_end)+"\t"+str(self.start)+"\t"\
				+str(self.end)+"\t"+str(self.peak_GPS)+"\t"
	