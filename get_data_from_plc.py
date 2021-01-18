import pyads
import numpy as np
import fft_nonRecursive
import threading
import time
import matplotlib.pyplot as plt

time_length = 256
time_vector = [np.zeros(time_length), np.zeros(time_length, dtype=complex)]
fft_vector = [np.zeros(time_length), np.zeros(time_length, dtype=complex)]
fft_plc_vector = [np.zeros(time_length), np.zeros(time_length, dtype=complex)]
plot_new_data = True

def format_complex_array(input: list):
    array = np.zeros(int(len(input)/2), dtype=complex)
    for i in range(len(array)):
        array[i] = input[i*2] + 1j * input[i*2+1]
    return array

def plot():
    global plot_new_data
    f,ax = plt.subplots(3)
    # Prepare the Plotting Environment with random starting values
    x = np.arange(10000)
    y = np.random.randn(10000)
    # Plot 0 is for raw data
    li, = ax[0].plot(x, y)
    ax[0].set_xlim(0,256)
    ax[0].set_ylim(-2,2)
    ax[0].set_title("Raw Signal")
    # Plot 1 is for the FFT
    li2, = ax[1].plot(x, y)
    ax[1].set_xlim(0,256)
    ax[1].set_ylim(-25,50)
    ax[1].set_title("Fast Fourier Transform")
    # Plot 1 is for the FFT of the plc
    li3, = ax[2].plot(x, y)
    ax[2].set_xlim(0,256)
    ax[2].set_ylim(-25,80)
    ax[2].set_title("Fast Fourier Transform of plc")
    # Show the plot, but without blocking updates
    plt.pause(0.01)
    plt.tight_layout()
    
    while True:
        if plot_new_data:
            li.set_xdata(np.arange(len(time_vector[0])))
            li.set_ydata(time_vector[1])
            li2.set_xdata(np.arange(len(fft_vector[0])))
            li2.set_ydata(20*np.log10(abs(fft_vector[1])))
            li3.set_xdata(np.arange(len(fft_plc_vector[0])))
            li3.set_ydata(20*np.log10(abs(fft_plc_vector[1])))
            plt.pause(0.01)
            plot_new_data = False
        else:
            time.sleep(0.2)
        pass
    
def thread_plc():
    global plot_new_data
    remote_ip = '192.168.20.157'
    remote_ads = '192.168.30.202.1.1'
    plc = pyads.Connection(remote_ads, pyads.PORT_TC3PLC1, remote_ip)
    while True:
        if not plot_new_data:
            plc.open()
            # symbols_list = plc.get_all_symbols()
            time_vector_temp = plc.read_by_name("KrogstrupMBE2.FFT_Data", pyads.PLCTYPE_LREAL * (256*2))
            fft_vector_temp = plc.read_by_name("KrogstrupMBE2.FFT_Result", pyads.PLCTYPE_LREAL * (256*2))
            plc.close()
            time_vector[1] = format_complex_array(time_vector_temp)
            fft_plc_vector[1] = format_complex_array(fft_vector_temp)
            time_duration = 1.28
            Fs = 200
            time_vector[0] = np.arange(0, time_duration, 1/Fs)
            fft_vector[1] = fft_nonRecursive.dif_fft4(time_vector[1])
            fft_vector[0] = np.linspace(0, Fs, len(fft_vector[1]), endpoint=False)
            plot_new_data = True
        else:
            time.sleep(0.1)
        pass
        

def main():
    plc_thread = threading.Thread(target=thread_plc)
    plc_thread.start()
    plot()
    pass

if __name__ == '__main__':
    main()