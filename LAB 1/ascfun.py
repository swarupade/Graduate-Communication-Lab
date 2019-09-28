from pylab import *
import numpy as np
def asc2bin(txt, bits=8): 
    """ 
    ASCII text to serial binary conversion 
    >>>>> dn = asc2bin(txt, bits) <<<<< 
    where txt         input text string 
          abs(bits)         bits per char, default=8 
          bits > 0          LSB first parallel to serial 
          bits < 0          MSB first parallel to serial dn binary output sequence 
    """
    txtnum = array([ord(c) for c in txt]) # int array 
    if bits > 0: # Neg powers of 2, increasing exp  
        p2 = np.power(2.0,arange(0,-bits,-1)) 
    else: # Neg powers of 2, decreasing exp 
        p2 = np.power(2.0,1+arange(bits,0)) 
    B = array(mod(array(floor(outer(txtnum,p2)),int),2),int8) 
         # Rows of B are bits of chars 
    dn = reshape(B,B.size) 
    return dn       # Serial binary output

def bin2asc(dn, bits=8, flg=1):
	"""
	Serial binary to ASCII text conversion
	>>>> txt = bin2asc(dn, bits, flg) <<<<
        where dn 		binary input sequence
	abs(bits)               bits per character
        bits > 0                LSB first parallel to serial
        bits < 0                MSB first parallel to serial
	flg != 0		limit range to [0...127]
	txt 			output text string
	"""
	N = int(dn.size/bits)
	C = reshape(dn,(N,bits)).astype(int)
	
	if bits > 0:
		p2 = np.power(2,arange(0,bits)).T 	 # making an array of increasing exponential of 2
	else:
		p2 = np.power(2,arange(bits,0,-1)-1).T 		# making an array of increasing exponential of 2

	num4txt = dot(C,p2)				# making an array of ASCII code from the binary coded characters
	txt = ''.join([chr(c) for c in num4txt])	# making a list of txt character and then joining it to make the required text
	return txt	