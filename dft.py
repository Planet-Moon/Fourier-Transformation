# 
# Discrete Fourier transform (Python)
# by Project Nayuki, 2017. Public domain.
# https://www.nayuki.io/page/how-to-implement-the-discrete-fourier-transform
# 


# 
# Computes the discrete Fourier transform (DFT) of the given complex vector.
# 'input' is a sequence of numbers (integer, float, or complex).
# Returns a list of complex numbers as output, having the same length.
# 
import cmath
def compute_dft_complex(input):
	n = len(input)
	output = []
	for k in range(n):  # For each output element
		s = complex(0)
		for t in range(n):  # For each input element
			angle = 2j * cmath.pi * t * k / n
			s += input[t] * cmath.exp(-angle)
		output.append(s)
	return output


# 
# (Alternate implementation using only real numbers.)
# Computes the discrete Fourier transform (DFT) of the given complex vector.
# 'inreal' and 'inimag' are each a sequence of n floating-point numbers.
# Returns a tuple of two lists of floats - outreal and outimag, each of length n.
# 
import math
def compute_dft_real_pair(inreal, inimag, N):
	assert len(inreal) == len(inimag)
	n = len(inreal)
	outreal = []
	outimag = []
	for k in range(N):  # For each output element
		sumreal = 0.0
		sumimag = 0.0
		for t in range(n):  # For each input element
			angle = 2 * math.pi * t * k / (N * 2)
			sumreal +=  inreal[t] * math.cos(angle) + inimag[t] * math.sin(angle)
			sumimag += -inreal[t] * math.sin(angle) + inimag[t] * math.cos(angle)
		outreal.append(sumreal)
		outimag.append(sumimag)
	return (outreal, outimag)

def gen_sinus(amp: float, freq: float, phase: float, time: float, Fs: float):
    sinusList = []
    timeList= []
    samples = math.floor(time * Fs)
    currentTime = 0.0
    for i in range(samples):
        timeList.append(currentTime)
        sinusList.append(amp * math.sin(2 * math.pi * freq * currentTime + phase))
        currentTime += 1/Fs
    return sinusList, timeList

def gen_rect(amp: float, freq: float, phase: float, time: float, Fs: float):
    rectList = []
    timeList = []
    samples = math.floor(time * Fs)
    currentTime = 0.0
    for i in range(samples):
        timeList.append(currentTime)
        rectList.append(amp * round(math.sin(2 * math.pi * freq * currentTime + phase)))
        currentTime += 1/Fs
    return rectList, timeList
