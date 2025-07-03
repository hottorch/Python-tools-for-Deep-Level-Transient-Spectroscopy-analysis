# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 13:12:30 2025

@author: 81412
"""
import numpy as np

def Q_DLTS(Temp, t, I, ratio=None, pts=None, baseline=None, SNR=False):
    """
    Interpret current transient to Q-DLTS signalgnal
    
    Argument list:
        t: time in 's'
        I: Current in 'A'
        ratio: product of t2/t1, default ratio is '2'
        pts: how many data points you want in the output, 
        baseline: the value that the current transient shoud be offset to make sure it decays to 0
    
    Returns:
        A 2 columns array, col_0 is the time constant in 'Second', col_1 is the Q-DLTS signal (delta Q) in 'Coulomb'
    """
    Q_DLTS_sig = []
    print ('\n-------------------------------')
    if ratio is None:
        ratio = 2
        print('t2/t1 Ratio: 2 (default)')
    else:
        print(f't2/t1 Ratio: {ratio}')
    if pts is None:
        pts = 50
        print('Output data points: 50 (default)')
    else:
        print(f'Output data points: {pts}')
        
    #Subtract the baseline according to the mean value of the last 1% data or the specified value
    num_of_data = int(0.01*len(t))
    mean_val = np.mean(I[-num_of_data:])
    
    if baseline is None:
        baseline = mean_val
        print(f'Back ground ubtracted by: {mean_val*1E6} \u00B5A (the mean value of the last 1% data points)')
    else:
        print(f'Back ground ubtracted by: {baseline} A (User defined)')
    print ('-------------------------------')
    I = I - baseline
    
    # Reduce the output data points of the Q-DLTS spectra by choosing t1 exponentially.
    index_base = (t.size / ratio) ** (1 / pts)
    
    for i in range(pts):
        t1 = t[int(index_base**i)]
        t2 = ratio * t1
        if t2 > t[-1]:  # break if t2 is out of the bound for t
            break
        # Define a time window
        mask = (t1 <= t) & (t <= t2)
        t_sub = t[mask]
        I_sub = I[mask]
        # Find the integral of I with respect to t in this subset (trapezoidal rule)
        delta_Q = np.trapz(I_sub, t_sub)  # Convert to coulomb

        tau = (t_sub[-1] - t_sub[0]) / np.log(t_sub[-1] / t_sub[0])  # In second
        Q_DLTS_sig.append([Temp, tau, delta_Q, t1, t2,t_sub[0], t_sub[-1], index_base])
        
    return np.array(Q_DLTS_sig)  # Construct a 2 columns array