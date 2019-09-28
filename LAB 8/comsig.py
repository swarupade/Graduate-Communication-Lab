# Module: comsig.py 
# Class and method definitions for COMmunication Signals 
import numpy as np
import copy
class sigWave:
    """ Class for 'waveform' signals """
    type = 'waveform'
    def __init__(self, sig, Fs=8000, t0=0, signame='Signal'):
        """
        sig: real or complex-valued waveform samples
        Fs: sampling rate (default 8000 samples/sec)
        t0: start time of waveform in seconds (default 0)
        """
        self.sig = np.asanyarray(sig)
        self.Fs = Fs
        self.t0 = t0
        self.Nsamp = len(self.sig)
        self.tlen = self.Nsamp/float(self.Fs)
        self.tend = self.t0 + (self.Nsamp-1)/float(self.Fs)
        self.signame = signame

    def __len__(self):
        return self.Nsamp    # Returns length in samples
    def get_Fs(self):
        return self.Fs       # Returns sampling rate
    def get_t0(self):
        return self.t0       # Returns start time
    def set_t0(self, t0):
        self.t0 = t0         # Set new start time
        self.tend = self.t0 + (self.Nsamp-1)/float(self.Fs)

    def timeAxis(self):       # Generate time axis
        return self.t0 + np.arange(self.Nsamp)/float(self.Fs)
    def signal(self):         # Return the waveform
        return self.sig
    def copy(self):           # Make a copy of a sigWave object
        return copy.deepcopy(self)
    def scale(self, factor):  # Make a scaled copy of a sigWave object
        return sigWave(factor*self.sig, self.Fs, self.t0)

class sigSequ:
    """ Class for ’sequence’ (DT) signals """
    type = 'sequence'
    def __init__(self, sig, FB = 100, n0 = 0):
        """
        sig: real or complex-valued sequence values
        FB: symbol (or Baud) rate (default 100 Baud)
        n0: start index of sequence (default 0)
        """
        self.sig = np.asanyarray(sig)
        self.FB = FB
        self.n0 = n0

    # Properties
    def __len__(self):
        return len(self.sig)
    def get_FB(self):
        return self.FB
    def get_n0(self):
        return self.n0

    # Methods+    def indexAxis(self):
        return self.n0 + np.arange(len(self.sig))
    def signal(self):
        return self.sig
    def copy(self):
        return copy.deepcopy(self)
    def scale_offset(self, a, b = 0):
        """ x[n]_out = a*x[n]_in + b """
        return sigSequ(a*self.sig+b, self.FB, self.n0)

