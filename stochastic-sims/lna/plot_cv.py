# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 08:45:43 2018

@author: aivar
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as mpl

mpl.rc('font', family='sans-serif')
mpl.rc('text', usetex=True)

# setting up print output (if required to print)
np.set_printoptions(precision = 2, suppress = True)

data_rc = spio.loadmat('lna_omprc_kt_kp.mat')
means_rc  = data_rc['mean_rr']
stds_rc   = data_rc['std_rr']
noises_rc = stds_rc[:,:,:2,:]/means_rc


data_tc = spio.loadmat('lna_tazc_kt_kp.mat')
means_tc  = data_tc['mean_tc']
stds_tc   = data_tc['std_tc']
noises_tc = stds_tc[:,:,:2,:]/means_tc

colors_0 = ['b-', 'r-', 'c-', 'g-', 'k-']
colors_1 = ['b--', 'r--', 'c--', 'g--', 'k--']

plt.figure()
plt.xscale('log')
plt.yscale('log')
for k in range(means_rc.shape[1]):
        plt.plot(means_rc[0, k, 0, :]/602.2, noises_rc[0, k, 0, :], colors_0[k], ms = 10)
        plt.plot(means_rc[0, k, 1, :]/602.2, noises_rc[0, k, 1, :], colors_1[k], ms = 10)

plt.xlim([1e-3, 25000/602.2])
plt.ylim([0.02, 3])
plt.title(r'Coefficient of variation vs mean')
plt.xlabel(r'Mean Output Steady State ([$\mu$M])')
plt.ylabel('Coefficient of variation')


plt.figure()
plt.xscale('log')
plt.yscale('log')
for k in range(means_tc.shape[1]):
        plt.plot(means_tc[0, k, 0, :]/602.2, noises_tc[0, k, 0, :], colors_0[k], ms = 10)
        plt.plot(means_tc[0, k, 1, :]/602.2, noises_tc[0, k, 1, :], colors_1[k], ms = 10)
     
plt.xlim([1e-3, 25000/602.2])
plt.ylim([0.02, 3])
plt.title(r'Coefficient of variation vs mean')
plt.xlabel(r'Mean Output Steady State ([$\mu$M])')
plt.ylabel('Coefficient of variation')

plt.show()