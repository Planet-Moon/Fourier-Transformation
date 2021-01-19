import numpy as np
import math
from dft import gen_rect, gen_sinus, gen_sawtooth
import matplotlib.pyplot as plt
import cmath
import fftUtils

def dif_fft4(f: list):
    '''
    radix-4 dif fft
    '''
    pot_of_4 = math.log(len(f))/math.log(4)
    if pot_of_4%1 > 0:
        raise ValueError("not a number of 4**{}".format(pot_of_4))
    out = f.copy()
    svar = int(np.log(len(out))/np.log(4))
    nvar = len(out)
    kmax = int(3*((float(nvar)/4.)-1))
    k_wavenumber = np.linspace(0, kmax, kmax+1)
    twiddle = np.exp(-2j * np.pi * k_wavenumber / nvar)
    tss = 1
    krange = int(float(nvar)/4.)
    block = 1
    base = 0
    for wvar in range(0, svar):
        for hvar in range(0, block):
            for kvar in range(0, krange):
                # butterfly
                offset = int(nvar / 4)
                avar = base + kvar + (0 * offset)
                bvar = base + kvar + (1 * offset)
                cvar = base + kvar + (2 * offset)
                dvar = base + kvar + (3 * offset)
                apc = out[avar] + out[cvar]
                bpd = out[bvar] + out[dvar]
                amc = out[avar] - out[cvar]
                bmd = out[bvar] - out[dvar]
                out[avar] = apc + bpd
                if kvar == 0:
                    out[bvar] = amc - (1j*bmd)
                    out[cvar] = apc - bpd
                    out[dvar] = amc + (1j*bmd)
                else:
                    r1var = twiddle[1 * kvar * tss]
                    r2var = twiddle[2 * kvar * tss]
                    r3var = twiddle[3 * kvar * tss]
                    out[bvar] = (amc - (1j * bmd)) * r1var
                    out[cvar] = (apc - bpd) * r2var
                    out[dvar] = (amc + (1j * bmd)) * r3var
                pass
            base = base + (4*krange)
        block = block * 4
        nvar = float(nvar) / 4.
        krange = int(float(krange) / 4.)
        base = 0
        tss = int(tss * 4)
    zvar = fftUtils.digitreversal(out, 4, svar, len(out))
    return zvar

def dif_fft0 (f: list):
    '''
    radix-2 dif  
    '''
    out = f.copy()
    svar = int(np.log(len(out))/np.log(2))
    nvar = len(out)
    kmax = int((float(nvar)/2.)-1)
    k_wavenumber = np.linspace(0, kmax, kmax+1)
    twiddles = np.exp(-2j*np.pi*k_wavenumber/nvar)
    b_p = 1
    nvar_p = out.size
    twiddle_step_size = 1
    for pvar in range(0, svar):           # pass loop
        nvar_pp = int(nvar_p/2)
        base_e = 0
        for bvar in range(0, b_p):       # block loop
            base_o = base_e + nvar_pp
            for nvar in range(0, nvar_pp):   # butterfly loop
                evar =  out[base_e+nvar] + out[base_o+nvar]
                if nvar == 0:
                    ovar = out[base_e + nvar] - out[base_o + nvar]
                else:
                    twiddle_factor =  nvar * twiddle_step_size
                    ovar = (out[base_e + nvar] - out[base_o + nvar]) * twiddles[twiddle_factor]
                out[base_e + nvar] = evar
                out[base_o + nvar] = ovar
            base_e = base_e + nvar_p
        b_p = b_p * 2
        nvar_p = int(nvar_p / 2)
        twiddle_step_size = 2 * twiddle_step_size
    zvar = fftUtils.digitreversal(out, 2, svar, len(out))
    return zvar

def plot_spectrum(time, signal, Fs):
    plt.figure()
    plt.plot(time, signal.real)
    plt.title('Signal')
    plt.ylabel('')
    plt.xlabel('t')
    
    fft_complex_result = dif_fft4(signal)
    
    if(False):
        fftFreqList = np.linspace(0, Fs, len(fft_complex_result), endpoint=False)
        plt.figure()
        plt.plot(fftFreqList, 20*np.log10(abs(fft_complex_result)))
        plt.title('Raw return of fft')
        plt.ylabel('')
        plt.xlabel('Hz')
    
    fft_abs = np.zeros(len(fft_complex_result))
    fft_arg = np.zeros(len(fft_complex_result))
    for i in range(len(fft_abs)):
        fft_abs[i] = abs(fft_complex_result[i])
        fft_arg[i] = cmath.phase(fft_complex_result[i])/math.pi
    fft_abs = fft_abs[:int(len(fft_abs)/2)]
    fft_arg = fft_arg[:int(len(fft_arg)/2)]
        
    np_fft = np.fft.fft(signal)
    np_fft = np_fft[range(int(len(signal)/2))]
    np_fft_freq = np.fft.fftfreq(len(signal), 1/Fs)
    np_fft_freq = np_fft_freq[range(int(len(signal)/2))]
    
    fftFreqList = np.linspace(0, int(Fs/2), len(fft_abs), endpoint=False)
    plt.figure()
    plt.subplot(211)
    plt.plot(fftFreqList, 20*np.log10(fft_abs), label='myFFT') # meine FFT
    # plt.plot(np_fft_freq, np_fft.real, label='npFFT.real')
    # plt.plot(np_fft_freq, np_fft.imag, label='npFFT.imag')
    plt.plot(np_fft_freq, 20*np.log10(abs(np_fft)), label='abs(npFFT)') # numpy FFT
    plt.legend(loc="upper right")
    plt.title('sinusTest spectrum')
    plt.ylabel('Abs [dB]')
    plt.xlabel('f [Hz]')
    
    plt.subplot(212)
    plt.plot(fftFreqList, fft_arg, label='myFFT')  # meine FFT
    np_fft_arg = list(map(lambda x: cmath.phase(x)/(1*math.pi), np_fft))
    plt.plot(np_fft_freq, np_fft_arg, label='arg(npFFT)') # numpy FFT
    plt.legend(loc="upper right")
    plt.title('sinusTest phase')
    plt.ylabel('Arg [Pi]')
    plt.xlabel('f [Hz]')
    plt.show()
    pass

if __name__ == "__main__":
    amp = 2.0
    freq = 10 # Hz
    phase = 0 # Pi
    time_duration = 5.12 # s
    Fs = 200 # Hz sample frequency
    time_vec = np.arange(0, time_duration, 1/Fs)
    signal_list = {}
    signal_list["rect_vec"] = np.array(gen_rect(amp=amp*2, freq=freq, phase=phase, time=time_vec))
    signal_list["sinus_vec"] = np.array(gen_sinus(amp=amp, freq=freq*3, phase=phase, time=time_vec), dtype=complex)
    signal_list["sinus_vec2"] = np.array(gen_sinus(amp=amp, freq=freq*2, phase=phase, time=time_vec), dtype=complex)
    signal_list["noise"] = np.random.normal(0,0.001,len(time_vec))
    signal_list["gen_sawtooth"] = np.array(gen_sawtooth(amp=amp*2, freq=freq/5, phase=phase, time=time_vec), dtype=complex)
    signal_list["constant"] = np.zeros(len(time_vec), dtype=complex) + 0.2
    signal_list["sinus_vec3"] = np.array(gen_sinus(amp=amp, freq=4, phase=phase, time=time_vec), dtype=complex)
    signal = np.zeros(len(time_vec), dtype=complex)
    for i in signal_list:
        signal += signal_list[i]
    
    plot_spectrum(time_vec, signal, Fs)
    pass