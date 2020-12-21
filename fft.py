import numpy as np
import math
from dft import gen_rect, gen_sinus
import matplotlib.pyplot as plt
import cmath

def fft_complex(f: list):
    n = len(f)
    if n == 1:
        return f
    else:
        even_index = f[::2]
        odd_index = f[1::2]
        g = fft_complex(even_index)
        u = fft_complex(odd_index)
        c = np.zeros(int(n), dtype=complex)
        for k in range(int(n/2)):
            angle = 2.0 * cmath.pi * 1j * k / n
            c[k] = g[k] + u[k] * cmath.exp(-angle)
            c[int(k + n/2)] = g[k] - u[k] * cmath.exp(-angle)
        return c


def fft(f_real: list, f_imag: list):
    assert len(f_real) == len(f_imag)
    n = len(f_real)
    if n == 1:
        return f_real, f_imag
    else:
        even_index_real = f_real[::2]
        even_index_imag = f_imag[::2]
        odd_index_real = f_real[1::2]
        odd_index_imag = f_imag[1::2]
        g_real, g_imag = fft(even_index_real, even_index_imag)
        u_real, u_imag = fft(odd_index_real, odd_index_imag)
        c_real = np.zeros(int(n))
        c_imag = np.zeros(int(n))
        k_factor = 2 * math.pi / n
        for k in range(0, int(n/2)):
            index = [k, int(k + (n/2))]
            angle =  k * k_factor
            c_real[index[0]] = g_real[k] + (u_real[k] * math.cos(angle) + u_imag[k] * math.sin(angle))
            c_imag[index[0]] = g_imag[k] + (u_imag[k] * math.cos(angle) - u_real[k] * math.sin(angle))
            c_real[index[1]] = g_real[k] - (u_real[k] * math.cos(angle) + u_imag[k] * math.sin(angle))
            c_imag[index[1]] = g_imag[k] - (u_imag[k] * math.cos(angle) - u_real[k] * math.sin(angle))
        return c_real, c_imag
    
    
def plot_spectrum(time, signal, complex_enable: bool):
    
    if complex_enable:
        fft_complex_result = fft_complex(signal)
        fft_abs = np.zeros(len(fft_complex_result))
        fft_arg = np.zeros(len(fft_complex_result))
        for i in range(len(fft_abs)):
            fft_abs[i] = abs(fft_complex_result[i])
            fft_arg[i] = cmath.phase(fft_complex_result[i])
    else:
        zero_list = np.zeros(len(time_vec))
        fft_result_real, fft_result_imag = fft(signal, zero_list)
        fft_abs = np.zeros(len(fft_result_real))
        fft_arg = np.zeros(len(fft_result_real))
        for i in range(len(fft_abs)):
            fft_abs[i] = math.sqrt(fft_result_real[i]**2 + fft_result_imag[i]**2)
            fft_arg[i] = math.atan2(fft_result_imag[i], fft_result_real[i])
    
    plt.figure()
    plt.plot(time_vec, signal.real)
    plt.title('Signal')
    plt.ylabel('')
    plt.xlabel('t')
    
    fftFreqList = np.linspace(0, Fs, len(fft_abs))
    plt.figure()
    plt.subplot(211)
    plt.plot(fftFreqList, 10.*np.log10(fft_abs))
    plt.title('sinusTest spectrum')
    plt.ylabel('Abs [dB]')
    plt.xlabel('f')
    plt.subplot(212)
    plt.plot(fftFreqList, fft_arg)
    plt.title('sinusTest phase')
    plt.ylabel('Arg')
    plt.xlabel('f')
    plt.show()
    pass

amp = 2.0
freq = 100 # Hz
phase = 0 # Pi
time = 2 # s
Fs = 2048 # Hz sample frequency
time_vec = np.arange(0, time, 1/Fs)
signal_list = {}
signal_list["rect_vec"] = np.array(gen_rect(amp=amp, freq=freq, phase=phase, time=time_vec))
signal_list["sinus_vec"] = np.array(gen_sinus(amp=amp, freq=freq, phase=phase, time=time_vec), dtype=complex)
signal_list["noise"] = np.random.normal(0,1,len(time_vec))
signal = np.zeros(len(time_vec), dtype=complex)
for i in signal_list:
    signal += signal_list[i]
    
plot_spectrum(time_vec, signal, False)

time_vec = np.array([0,1,2,3,4,5,6,7], dtype=complex)
signal = np.array([0,1,2,2,1,0,-1,-2], dtype=complex)
plot_spectrum(time_vec, signal, False)
pass