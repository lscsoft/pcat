# encoding: utf-8
'''
Daniele Trifiro`
brethil@phy.olemiss.edu

Definition of the Spike Class:
The Spike class is used to store the information about a spike, such as amplitude, 
	 peak position, start and end, so the original data can be retrieved:
		spike=Spike() 	
		spike.start			->	Spike starting point
		spike.end			->	Spike ending point
		spike.peak			->	Spike peak point
		spike.norm			->	Spike normalization constant
		spike.width			->	Spike Width (seconds)
		spike.waveform		->	Waveform
		spike.segment_start	->	Segment start expressed in GPS time
		spike.segment_end	->	Segment end expressed in GPS time
		spike.peak_GPS		->	Spike peak expressed in GPS time
		spike.sampling		->	Sampling Frequency for the channel being used
		spike.SNR			->	SNR of the glitch
	
	import inspect
	variables = [i for i in dir(t) if not inspect.ismethod(i)
	
'''

import matplotlib.mlab
import numpy as np

class Spike:
	
	def __init__(self, start, end, peak, norm, peak_GPS, segment_start, segment_end, waveform):
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
			
	def __str__(self):
		'''
		String reprentation of the spike (segment_start, segment_end, spike_start, 
		spike_end, spike_peak_GPS).
		'''
		return str(self.segment_start)+"\t"+str(self.segment_end)+"\t"+str(self.start)+"\t"\
				+str(self.end)+"\t"+str(self.peak_GPS)+"\t"
	
	def compute_SNR(self, f_sampl, original_segment_length):
		"""
			Given a PSD, computes the SNR of the glitch compared to the given PSD and stores the result
			in the SNR attribute
		"""
		self.sampling = f_sampl
		delta_f = f_sampl/float(self.len)
		delta_t = 1.0/f_sampl
		
		"""window = np.hanning(self.len)
		window_normalization = np.sum(window**2)/float(self.len)
		
		whitened_spike_PSD = np.abs( (1.0/f_sampl) * np.fft.rfft(self.waveform*window) )**2
		# Correct normalization
		whitened_spike_PSD *= ( 2.0 * (f_sampl/(self.len))/(window_normalization) )
		"""
		
		# Compute the PSD of the glitch, applying the correct normalization
		#whitened_spike_PSD = np.abs(delta_t*np.fft.rfft(self.waveform))**2
		
		# whitened_spike_PSD is simply d(SNR^2)/df, we can obtain the SNR by 
		# integrating, remembering to include a factor of 2 because the 
		# computed PSD is one-sided (has only positive frequenciese) because 
		# rfft is being used.
		
		
		#squared_SNR = 
		self.SNR = np.sqrt(squared_SNR)
		

	