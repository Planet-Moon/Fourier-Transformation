'''
Radix-4 DIT,  Radix-4 DIF,  Radix-2 DIT,  Radix-2 DIF FFTs
John Bryan,  2017
Python 2.7.3
'''

import numpy as np
import time
import matplotlib.pyplot as plt
import warnings
import math
from scipy import signal
# np.set_printoptions(threshold = np.nan, precision = 3, suppress = 1)
warnings.filterwarnings("ignore")


def swap(xarray, i, j):
    '''
    swap
    '''
    temp = xarray[i]
    xarray[i] = xarray[j]
    xarray[j] = temp
    return None


def digitreversal(xarray, radix, log2length, length):
    '''
    digitreversal
    '''
    if log2length % 2 == 0:
        n1var = int(np.sqrt(length))       #seed table size 
    else:
        n1var = int(np.sqrt(int(length/radix)))
    # algorithm 2,  compute seed table  
    reverse = np.zeros(n1var, dtype = int)
    reverse[1] = int(length/radix)
    for jvar in range(1, radix):
        reverse[jvar] = reverse[jvar-1]+reverse[1]
        for i in range(1, int(n1var/radix)):
            reverse[radix*i] = int(reverse[i]/radix)
            for jvar in range(1, radix):
                reverse[int(radix*i)+jvar] = reverse[int(radix*i)]+reverse[jvar]
    #algorithm 1
    for i in range(0, n1var-1):
        for jvar in range(i+1, n1var):
            uvar = i+reverse[jvar]
            vvar = jvar+reverse[i]
            swap(xarray, uvar, vvar)
            if log2length % 2 == 1:
                for zvar in range(1, radix):
                    uvar = i+reverse[jvar]+(zvar*n1var)
                    vvar = jvar+reverse[i]+(zvar*n1var)
                    swap(xarray, uvar, vvar)
    return xarray


def dif_fft4(xarray, twiddle, svar):
    '''
    radix-4 dif fft
    '''
    nvar = np.power(4, svar)
    tss = 1
    krange = int(float(nvar)/4.)
    block = 1
    base = 0
    for wvar in range(0, svar):
        for hvar in range(0, block):
            for kvar in range(0, krange):
                # butterfly
                offset = int(nvar/4)
                avar = base+kvar
                bvar = base+kvar+offset
                cvar = base+kvar+(2*offset)
                dvar = base+kvar+(3*offset)
                apc = xarray[avar]+xarray[cvar]
                bpd = xarray[bvar]+xarray[dvar]
                amc = xarray[avar]-xarray[cvar]
                bmd = xarray[bvar]-xarray[dvar]
                xarray[avar] = apc+bpd
                if kvar == 0:
                    xarray[bvar] = amc-(1j*bmd)
                    xarray[cvar] = apc-bpd
                    xarray[dvar] = amc+(1j*bmd)
                else:
                    r1var = twiddle[kvar*tss]
                    r2var = twiddle[2*kvar*tss]
                    r3var = twiddle[3*kvar*tss]
                    xarray[bvar] = (amc-(1j*bmd))*r1var
                    xarray[cvar] = (apc-bpd)*r2var
                    xarray[dvar] = (amc+(1j*bmd))*r3var
            base = base+(4*krange)
        block = block*4
        nvar = float(nvar)/4.
        krange = int(float(krange)/4.)
        base = 0
        tss = int(tss*4)
    return xarray



def fft4(xarray, twiddles, svar):
    '''
    radix-4 dit fft
    '''
    nvar = 4
    tss = np.power(4, svar-1)
    krange = 1
    block = int(xarray.size/4)
    base = 0
    for wvar in range(0, svar):
        for zvar in range(0, block):
            for kvar in range(0, krange):
                # butterfly
                offset = int(nvar/4)
                avar = base+kvar
                bvar = base+kvar+offset
                cvar = base+kvar+(2*offset)
                dvar = base+kvar+(3*offset)
                if kvar == 0:
                    xbr1 = xarray[bvar]
                    xcr2 = xarray[cvar]
                    xdr3 = xarray[dvar]
                else:
                    r1var = twiddles[kvar*tss]
                    r2var = twiddles[2*kvar*tss]
                    r3var = twiddles[3*kvar*tss]
                    xbr1 = (xarray[bvar]*r1var)
                    xcr2 = (xarray[cvar]*r2var)
                    xdr3 = (xarray[dvar]*r3var)
                evar = xarray[avar]+xcr2
                fvar = xarray[avar]-xcr2
                gvar = xbr1+xdr3
                hvar = xbr1-xdr3
                j_h = 1j*hvar
                xarray[avar] = evar+gvar
                xarray[bvar] = fvar-j_h
                xarray[cvar] = -gvar+evar
                xarray[dvar] = j_h+fvar
            base = base+(4*krange)
        block = int(block/4)
        nvar = 4*nvar
        krange = 4*krange
        base = 0
        tss = int(float(tss)/4.)
    return xarray


def dif_fft0 (xarray, twiddle, log2length):
    '''
    radix-2 dif  
    '''
    b_p = 1
    nvar_p = xarray.size
    twiddle_step_size = 1
    for pvar in range(0,  log2length):           # pass loop
        nvar_pp = int(nvar_p/2)
        base_e = 0
        for bvar in range(0,  b_p):       # block loop
            base_o = base_e+nvar_pp
            for nvar in range(0,  nvar_pp):   # butterfly loop
                evar =  xarray[base_e+nvar]+xarray[base_o+nvar]
                if nvar == 0:
                    ovar = xarray[base_e+nvar]-xarray[base_o+nvar]
                else:
                    twiddle_factor =  nvar*twiddle_step_size
                    ovar = (xarray[base_e+nvar] \
                        -xarray[base_o+nvar])*twiddle[twiddle_factor]
                xarray[base_e+nvar] = evar
                xarray[base_o+nvar] = ovar
            base_e = base_e+nvar_p
        b_p = b_p*2
        nvar_p = int(nvar_p/2)
        twiddle_step_size = 2*twiddle_step_size
    return xarray


def fft2 (xarray, twiddle, svar) :
    '''
    radix-2 dit
    '''
    nvar = xarray.size
    b_p = int(nvar/2)
    nvar_p = 2
    twiddle_step_size = nvar/2
    for pvar in range(0,  svar):
        nvar_pp = int(nvar_p/2)
        base_t = 0
        for bvar in range(0,  b_p):
            base_b = base_t+nvar_pp
            for nvar in range(0,  nvar_pp):
                if nvar == 0:
                    bot = xarray[base_b+nvar]
                else:
                    twiddle_factor = int(nvar*twiddle_step_size)
                    bot = xarray[base_b+nvar]*twiddle[twiddle_factor]
                top = xarray[base_t+nvar]
                xarray[base_t+nvar] = top+bot
                xarray[base_b+nvar] = top-bot
            base_t = base_t+nvar_p
        b_p = int(b_p/2)
        nvar_p = nvar_p*2
        twiddle_step_size = twiddle_step_size/2
    return xarray


def testr4dif():
    '''
    Test and time dif radix4 w/ multiple length random sequences
    '''
    flag = 0
    i = 0
    radix = 4
    r4diftimes = np.zeros(6)
    for svar in range (2, 8):
        xarray = np.random.rand(2*np.power(4, svar)).view(np.complex128)
        xpy = np.fft.fft(xarray)
        radix = 4
        nvar = np.power(4, svar)
        kmax = int(3*((float(nvar)/4.)-1))
        k_wavenumber = np.linspace(0, kmax, kmax+1)
        twiddlefactors = np.exp(-2j*np.pi*k_wavenumber/nvar)
        tvar = time.time()
        xarray = dif_fft4(xarray, twiddlefactors, svar)
        r4diftimes[i] =  time.time()-tvar
        xarray = digitreversal(xarray, radix, svar, nvar)
        t_f = np.allclose(xarray, xpy)
        if t_f == 0:
            flag = 1
        assert(t_f)
        i = i+1
    if flag == 0:
        print ("All radix-4 dif results were correct.")
    return r4diftimes


def testr4():
    '''
    Test and time dit radix4 w/ multiple length random sequences
    '''
    flag = 0
    i = 0
    radix = 4
    r4times = np.zeros(6)
    for svar in range (2, 8):
        xarray = np.random.rand(2*np.power(4, svar)).view(np.complex128)
        xpy = np.fft.fft(xarray)
        nvar = np.power(4, svar)
        xarray = digitreversal(xarray, radix, svar, nvar)
        kmax = int(3*((float(nvar)/4.)-1))
        k_wavenumber = np.linspace(0, kmax, kmax+1)
        twiddles = np.exp(-2j*np.pi*k_wavenumber/nvar)
        tvar = time.time()
        xarray = fft4(xarray, twiddles, svar)
        r4times[i] =  time.time()-tvar
        t_f = np.allclose(xarray, xpy)
        if t_f == 0:
            flag = 1
        assert(t_f)
        i = i+1
    if flag == 0:
        print ("All radix-4 dit results were correct.")
    return r4times


def testr2dif():
    '''
    Test and time radix2 dif w/ multiple length random sequences
    '''
    flag = 0
    i = 0
    radix = 2
    r2diftimes = np.zeros(6)
    for rvar in range (2, 8):
        svar = np.power(4, rvar)
        cpy = np.random.rand(2*svar).view(np.complex_)
        t = np.linspace(0,1,len(cpy))
        sinus_sig = np.array(list(map(lambda x: math.sin(2 * np.pi * 20 * x), t.tolist())))
        sinus_sig2 = np.array(list(map(lambda x: math.sin(2 * np.pi * 1234 * x), t.tolist())))
        cpy = 0.001 * cpy + sinus_sig + sinus_sig2 * 2
        gpy = np.fft.fft(cpy)
        nvar = svar
        kmax = int((float(nvar)/2.)-1)
        k_wavenumber = np.linspace(0, kmax, kmax+1)
        twiddles = np.exp(-2j*np.pi*k_wavenumber/nvar)
        t1time = time.time()
        gvar = dif_fft0(cpy, twiddles, int(2*rvar))
        r2diftimes[i] =  time.time()-t1time
        zvar = digitreversal(gvar, radix, int(2*rvar), svar)
        
        # plt.figure()
        # plt.plot(20*np.log10(zvar[:int(len(zvar)/2)]), label="zvar")
        # plt.plot(20*np.log10(gpy[:int(len(gpy)/2)]), label="gpy")
        # plt.legend(loc="upper right")
        # plt.title('spectrum of {}'.format(rvar))
        
        t_f = np.allclose(zvar, gpy)
        if t_f == 0:
            flag = 1
        assert(t_f)
        i = i+1
    if flag == 0:
        print ("All radix-2 dif results were correct.")
    return r2diftimes


def testr2():
    '''
    Test and time radix2 dit w/ multiple length random sequences
    '''
    radix = 2
    flag = 0
    i = 0
    r2times = np.zeros(6)
    for rvar in range (2, 8):
        svar = np.power(4, rvar)
        cpy = np.random.rand(2*svar).view(np.complex_)
        gpy = np.fft.fft(cpy)
        nvar = svar
        kmax = int((float(nvar)/2.)-1)
        k_wavenumber = np.linspace(0, kmax, kmax+1)
        twiddles = np.exp(-2j*np.pi*k_wavenumber/nvar)
        zvar = digitreversal(cpy, radix, int(2*rvar), svar)
        t1time = time.time()
        garray = fft2(zvar, twiddles, int(2*rvar))
        r2times[i] =  time.time()-t1time        
        t_f = np.allclose(garray, gpy)
        if t_f == 0:
            flag = 1
        assert(t_f)
        i = i+1
    if flag == 0:
        print ("All radix-2 dit results were correct.")
    return r2times


def plot_times(tr2, trdif2, tr4, trdif4):
    '''
    plot performance
    '''
    uvector = np.zeros(6, dtype = int)
    for i in range(2, 8):
        uvector[i-2] = np.power(4, i)
    plt.figure(figsize = (7, 5))
    plt.rc("font", size = 9)
    plt.loglog(uvector, trdif2, 'o',  ms = 5,  markerfacecolor = "None",  \
               markeredgecolor = 'red',  markeredgewidth = 1,  \
               basex = 4,  basey = 10,  label = 'radix-2 DIF')
    plt.loglog(uvector, tr2, '^',  ms = 5,  markerfacecolor = "None", \
               markeredgecolor = 'green',  markeredgewidth = 1,  \
               basex = 4,  basey = 10,  label = 'radix-2 DIT')
    plt.loglog(uvector, trdif4, 'D',  ms = 5,  markerfacecolor = "None", \
               markeredgecolor = 'blue',  markeredgewidth = 1, \
                basex = 4,  basey = 10,  label = 'radix-4 DIF')
    plt.loglog(uvector, tr4, 's',  ms = 5,  markerfacecolor = "None", \
               markeredgecolor = 'black',  markeredgewidth = 1, \
               basex = 4,  basey = 10,  label = 'radix-4 DIT')
    plt.legend(loc = 2)
    plt.grid()
    plt.xlim([12, 18500])
    plt.ylim([.00004, 1])
    plt.ylabel("time (seconds)")
    plt.xlabel("sequence length")
    plt.title("Time vs Length")
    plt.savefig('tvl2.png',  bbox_inches = 'tight')
    plt.show()
    return None


def test():
    '''
    test performance
    '''
    trdif4 = testr4dif()
    tr4 = testr4()
    trdif2 = testr2dif()
    tr2 = testr2()
    plot_times(tr2, trdif2, tr4, trdif4)
    return None

if __name__ == "__main__":
    test()

