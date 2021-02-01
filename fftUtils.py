import numpy as np

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
                reverse[(radix*i)+jvar] = reverse[(radix*i)]+reverse[jvar]
                pass
    #algorithm 1
    for i in range(0, n1var-1):
        for jvar in range(i+1, n1var):
            uvar = int(i + reverse[jvar])
            vvar = int(jvar + reverse[i])
            swap(xarray, uvar, vvar)
            if log2length % 2 == 1:
                for zvar in range(1, radix):
                    uvar = int(i + reverse[jvar] + (zvar * n1var))
                    vvar = int(jvar + reverse[i] + (zvar * n1var))
                    swap(xarray, uvar, vvar)
    return xarray