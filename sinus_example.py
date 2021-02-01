import matplotlib.pyplot as plt
import math


def main():
    f = 2
    T = 1/f
    CycleTime = 0.005
    N_cycles = int(T / CycleTime)
    
    signal = []
    time_vec = []
    for i in range(N_cycles):
        time_vec.append(i*CycleTime)
        signal.append(math.sin(time_vec[i]*f*2*math.pi))
    
    plt.figure()
    plt.plot(time_vec, signal)
    plt.title('Signal')
    plt.ylabel('A')
    plt.xlabel('t')
    plt.show()
    return(0)
    
if __name__ == '__main__':
    main()