# Q-DLTS Signal Processing Script
# --------------------------------
# This script processes current transients to compute Q-DLTS (charge-based Deep Level Transient Spectroscopy) signals.
# Developed by: hottorch
# Created on: July 3, 2025

import numpy as np

def q_dlts(Temp, t, I, ratio=None, pts=None, baseline=None):
    """
    Calculate the Q-DLTS signal from a current transient.

    Parameters
    ----------
    Temp : float
        Temperature at which the measurement is made (for labeling only).
    t : np.ndarray
        Time array (in seconds).
    I : np.ndarray
        Current array (in amperes).
    ratio : float, optional
        Ratio t2/t1 for the integration window. Default is 2.
    pts : int, optional
        Number of Q-DLTS points to generate. Default is 50.
    baseline : float, optional
        Manual baseline value to subtract from I. If None, uses mean of last 1% of data.
    SNR : bool, optional
        Currently unused. Reserved for future signal-to-noise calculations.

    Returns
    -------
    np.ndarray
        2D array where each row contains:
        [Temperature, Tau (s), delta_Q (C), t1, t2, t_start, t_end, index_base]
    
    t1, t2:
        Nominal start and end times (in seconds) of the integration window, where t2 = ratio × t1.
    t_start, t_end:
        Actual start and end times from the time array closest to t1 and t2, used for integration.
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
