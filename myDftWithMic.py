import pyaudio
import matplotlib.pyplot as plt
import numpy as np
import dft
import math



FORMAT = pyaudio.paInt16 # We use 16bit format per sample
CHANNELS = 1
RATE = 44100 # Fs [Hz]
CHUNK = 1024 # 1024bytes of data red from a buffer
RECORD_SECONDS = 0.1
WAVE_OUTPUT_FILENAME = "file.wav"

audio = pyaudio.PyAudio()

# start Recording
stream = audio.open(format=FORMAT,
                    channels=CHANNELS,
                    rate=RATE,
                    input=True)#,
                    #frames_per_buffer=CHUNK)
stream.start_stream()


f,ax = plt.subplots(2)

# Prepare the Plotting Environment with random starting values
x = np.arange(10000)
y = np.random.randn(10000)

# Plot 0 is for raw audio data
li, = ax[0].plot(x, y)
ax[0].set_xlim(0,1000)
ax[0].set_ylim(-5000,5000)
ax[0].set_title("Raw Audio Signal")
# Plot 1 is for the FFT of the audio
li2, = ax[1].plot(x, y)
ax[1].set_xlim(0,RATE/2)
ax[1].set_ylim(0,100000)
ax[1].set_title("Fast Fourier Transform")
# Show the plot, but without blocking updates
plt.pause(0.01)
plt.tight_layout()

while True:
    in_data = stream.read(CHUNK)
    audio_data = np.frombuffer(in_data, np.int16)

    li.set_xdata(np.arange(len(audio_data)))
    li.set_ydata(audio_data)

    zero_list = np.zeros(len(audio_data))
    dft_real, dft_imag = dft.compute_dft_real_pair(audio_data, zero_list, 100)
    dft_abs = np.zeros(len(dft_real))
    dft_arg = np.zeros(len(dft_real))
    for i in range(len(dft_abs)):
        dft_abs[i] = math.sqrt(dft_real[i]**2 + dft_imag[i]**2)
        dft_arg[i] = math.atan2(dft_imag[i], dft_real[i])
    dftFreqList = np.linspace(0, RATE/2, len(dft_abs))
    li2.set_xdata(dftFreqList)
    li2.set_ydata(dft_abs)
    plt.pause(0.01)
    pass

pass