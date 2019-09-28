# File: wavfun.py
# Functions for reading and writing 16-bit PCM wav-files in Python.
import numpy as np
import wave
import struct

def wavread(fname):
    """
    Read N-channel 16-bit PCM wav-file
    >>>>> Fs, rt = wavread(fname) <<<<<
    where  fname   file name of wav-file
           Fs      sample rate of wav-file
           rt      data read from wav-file, N channels,
                   Nsamples data samples per channel, 
                   (Nsamples x N) numpy array of type np.float32,
                   data samples normalized to range -1 ... +1
    """
    fh = wave.open(fname,'rb')
    (nchannels, sampwidth, framerate, nframes, comptype, compname) = fh.getparams()
    if sampwidth == 2:
        frames = fh.readframes(nframes * nchannels)
        dn = struct.unpack_from('%dh' % nframes*nchannels, frames)
        if nchannels > 1:
            out = np.array([dn[i::nchannels] for i in range(nchannels)])/float(2**15)
        else:
            out = np.array(dn)/float(2**15)
    else:
        print('not a 16 bit wav-file')
        out = [0]
    fh.close()
    return (out,framerate)

def wavwrite(fname, Fs, xt):
    """
    Write N-channel 16-bit PCM wav-file
    >>>>> wavwrite(fname, Fs, xt) <<<<<
    where  fname   file name of wav-file
           Fs      sample rate of wav-file
           xt      data to be written to wav-file
                   (Nsamples x N) numpy array of floats
                   normalized to range -1...+1
                   N channels, Nsamples data samples per channel
    """
    # convert to np.int16 data type
    xt = np.array((2**15-1)*xt, np.int16)
    sio_wav.write(fname, Fs, xt)
