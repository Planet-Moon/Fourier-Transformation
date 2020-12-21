#%%
import matplotlib.pyplot as plt
import numpy as np
import dft

amp = 1.0
freq = 50 # Hz
phase = 0 # s
time = 1 # s
Fs = 3000 # Hz sample frequency

sinusTest, timeTest = gen_sinus(amp=amp, freq=freq, phase=phase, time=time, Fs=Fs)
sinusTest2, timeTest2 = gen_sinus(amp=amp*2, freq=freq*2, phase=phase, time=time, Fs=Fs)
rectTest, timeTest3 = gen_rect(amp=amp, freq=freq/2, phase=phase, time=time, Fs=Fs)

noise = np.random.normal(0,1,len(sinusTest))
# 0 is the mean of the normal distribution you are choosing from
# 1 is the standard deviation of the normal distribution
# 100 is the number of elements you get in array noise
plt.figure()
plt.plot(timeTest, noise)
plt.title('noise')
plt.ylabel('noise')
plt.xlabel('t [s]')

for i in range(len(sinusTest)):
    sinusTest[i] = sinusTest[i] + sinusTest2[i] + rectTest[i] + noise[i]

# Timeseries plot
plt.figure()
plt.plot(timeTest, sinusTest)
plt.title('sinusTest')
plt.ylabel('sinus')
plt.xlabel('t [s]')

zero_list = np.zeros(len(timeTest))
dft_real, dft_imag = dft.compute_dft_real_pair(sinusTest, zero_list, 2000)
dft_abs = np.zeros(len(dft_real))
dft_arg = np.zeros(len(dft_real))
for i in range(len(dft_abs)):
    dft_abs[i] = math.sqrt(dft_real[i]**2 + dft_imag[i]**2)
    dft_arg[i] = math.atan2(dft_imag[i], dft_real[i])

# Frequencydomain plot
dftFreqList = np.linspace(0, Fs/2, len(dft_abs))
plt.figure()

plt.subplot(211)
# plt.plot(dftFreqList, 10.*np.log10(dft_abs)) # logarithmic
plt.plot(dftFreqList, dft_abs) # linear
plt.title('sinusTest spectrum')
plt.ylabel('Abs') 	
plt.xlabel('f [Hz]')

plt.subplot(212)
plt.plot(dftFreqList, dft_arg)
plt.title('sinusTest phase')
plt.ylabel('Arg')
plt.xlabel('f [Hz]')

plt.show()