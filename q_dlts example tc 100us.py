# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 14:22:35 2025

@author: 81412
"""

from q_dlts import q_dlts
import numpy as np
import matplotlib.pyplot as plt  

# Example transient
t = np.linspace(0, 1e-3, 1000)  # 0 to 1 ms
I = np.exp(-t / 1e-4) * 1e-6  # Simulated transient (A)
Temp = 300  # K

result = q_dlts(Temp, t, I, ratio=None, pts=50, baseline=None)

# Extract tau and delta_Q from result
tau = result[:, 1] * 1E3        # Column 1: time constant (tau), convert from s to ms
deltaQ = result[:, 2] * 1e9     # Column 2: delta_Q, convert from C to nC

# Plot Q-DLTS spectrum
fig, ax1 = plt.subplots()
plt.scatter(tau, deltaQ, label="Q-DLTS", s=30)
plt.xscale('log')
ax1.set_xlabel(r'Time constant ($\tau: ms$)')
ax1.set_ylabel(r'$\Delta$ Q (nC)')
ax1.set_title(f'Q-DLTS Spectrum at {Temp} K')

