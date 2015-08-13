#!/usr/bin/env python
'''
[varargout] = FFTex1(f,N,A,theta,sr,maxf,varargin)
 
 
DESCRIPTION 
------------------------------------------------------------------------| 
The purpose of this function is to illustrate the basic principles of the
Fast Fourier Transform (FFT). The function takes as input the 
frequency (f) of a sine wave and returns a plot of a sine wave of this 
frequency. Additional parameters can also be specified (see below). 

 
INPUTS 
------------------------------------------------------------------------| 
f: frequency
N: number of seconds
A: Amplitude of sinusoid
theta: phase of sinusoid (radians)
sr: sampling rate (Hz)
maxf: max frequency to display in plot (Hz)

 
OUTPUTS 
------------------------------------------------------------------------| 
A couple of plots, along with optional outputs:
varargout{1} ts: time series structure specifying sampling rate, 
                 time points, and sinusoid amplitude.
varargout{2} fs: Fourier series structure specifying spectral resolution, 
                 frequency bin centers, and Fourier coefficients.

  
NOTES 
------------------------------------------------------------------------| 
Examples of usage:
FFtex1 by itself plots a 1 Hz sinusoid over a single cycle of 1 second.
[ts,fs] = FFTex1(1,2,1); returns a 1 Hz sinusoid plotted over 2 cycles
                         with an amplitude of 2.
 
  
Written 08/13/2015 
By Sam Thorpe 
'''


# # Module Imports
# -----------------------------------------------------|
from sys import argv
from sampy.common import keyboard
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt


# # Defs
# -----------------------------------------------------|
def FFTex1(*arg):
	''' main function which builds/plots sinusoid and spectrum'''

	
	# # Parse Inputs 
	# -----------------------------------------------------|
	f, N, A, theta, sr, maxf = arg


	# # compute time series and spectrum
	# -----------------------------------------------------|
	t = np.r_[0:N:1./sr]
	y = A*np.sin(2.*np.pi*f*t + theta)
	yF = np.fft.fft(y)/len(y)
	df = sr/len(yF)
	freqs = np.r_[0:sr:df] 
	yS = np.abs(yF)
	yT = np.angle(yF)


	# # Plots
	# -----------------------------------------------------|

	# # plot spectrum
	plt.figure(1)
	plt.bar(freqs,yS)
	plt.rc('font', **{'family':'normal','weight':'bold','size':16})
	plt.xlabel('Frequency (Hz)',fontweight='bold',fontsize=16)
	plt.ylabel('Amplitude',fontweight='bold',fontsize=16)
	plt.title('%i Hz sinusoid: Amplitude Spectrum' %f ,fontweight='bold',fontsize=18)
	plt.grid(True)
	plt.axis([0, maxf, 0, A])


	# # plot time series
	plt.figure(2)
	plt.rc('font', **{'family':'normal','weight':'bold','size':16})
	plt.axhline(0,0,1,linewidth=3,linestyle='--',color=[0.65, 0.65, 0.65])
	ytmp = A*np.sin(2.*np.pi*f*t)
	plt.plot(t,ytmp,linewidth=2,linestyle='--',color=[0.65, 0.65, 0.65]) 
	plt.plot(t,y,linewidth=3,color='b')
	plt.xlabel('Time (sec)',fontweight='bold',fontsize=16)
	plt.ylabel('Amplitude',fontweight='bold',fontsize=16)
	plt.title('%i Hz sinusoid: Time Series' %f ,fontweight='bold',fontsize=18)
	plt.grid('on');
	plt.axis([0, N, -A, A])
	plt.show()


	# # Parse Varargout
	# -----------------------------------------------------|
	'''
	if nargout 
	    ts.sr = sr; ts.t = t; ts.y = y;
	    varargout{1} = ts;     
	end;
	if nargout>1  
	    fs.df = df; fs.freqs = freqs; fs.yF = yF;
	    varargout{2} = fs;     
	end;
	'''



# # Test Suite
# -----------------------------------------------------|
if __name__ == '__main__':

	f = 1.0                       #| frequency
	N = 1.0                       #| N seconds
	A = 1.0                       #| Amplitude
	theta = np.pi/2                   #| Phase
	sr = 1000.0                   #| sampling rate
	maxf = 10.0                   #| max frequency to plot
	FFTex1(f,N,A,theta,sr,maxf)
	


#                                END ALL 
# # ----------------------------------------------------------------------| 
#  
