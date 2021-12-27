# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 08:45:43 2018

@author: aivar
"""
import matplotlib.pyplot as plt
import numpy as np
import pickle 
import time
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import rc
import scipy.io as spio
import matplotlib.pyplot as mpl

mpl.rc('font', family='sans-serif')
mpl.rc('text', usetex=True)

#rc('font', **{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

# setting up print output (if required to print)
np.set_printoptions(precision = 2, suppress = True)

data_lna = spio.loadmat('lna_data.mat')
means_lna  = data_lna['mean_rr']
stds_lna   = data_lna['std_rr']
noises_lna = stds_lna/means_lna
asp_lna  = data_lna['asp_grid']
sin_lna  = data_lna['input_p']

#colors_ = ['kd', 'cx','cx', 'co', 'co','cd','cd','bx','bx','bo','bo','bd']
#sink_tot_grid = np.array([0, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1 ])
sink_tot_grid = np.array([0, 0.1, 0.4, 0.7, 1])
asp_grid      = np.array([0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10])

len_a    = asp_grid.shape[0]
len_s    = sink_tot_grid.shape[0]

##reading the data for the in cis design
#file_name = 'stats_ompr_ol.npz'
#result_ompr = np.load(file_name)
#means_ol  = result_ompr['means']
#noises_ol = result_ompr['noises']
#del result_ompr

file_name = 'stats_ompr_cl.npz'
result_ompr = np.load(file_name)
means_cl  = result_ompr['means']
noises_cl = result_ompr['noises']
del result_ompr


# Setting up the plots:
st_plot_grid     = [0, 1, 2, 3, 4] # we will plot only these values in the krep grid
# legend:
colors_ = ['kd', 'cx', 'cd','co','bd','bo','bx']
colors_solids = ['k-', 'c:', 'c--','c-','b--','b:','b-']
line_handles, line_labels = [None] * len_s, [None] * len_s
 
plt.figure()
#plt.xscale('log')
plt.yscale('log')
means_cur = [[means_cl[asp*len_s+st, -1, -1]/602.2 for asp in range(len_a)] for st in range(len_s)] 
noises_cur = [[noises_cl[asp*len_s+st, -1, -1] for asp in range(len_a)] for st in range(len_s)]
#print np.asarray(means_cur)
#print np.asarray(noises_cur)
for k in st_plot_grid:
        plt.plot(means_cur[k], noises_cur[k], colors_solids[k], ms = 10)

plt.xlim([0, 25000])
plt.ylim([0.02, 2])
plt.title(r'Coefficient of variation vs mean')
plt.xlabel(r'Mean GFP Steady State (Number of Molecules)')
plt.ylabel('Coefficient of variation')

plt.figure()
# Setting up the plots:
st_plot_grid     = [0, 2, 4] # we will plot only these values in the krep grid
st_plot_grid_lna = [0, 1, 2, 3, 4]
#plt.xscale('log')
plt.yscale('log')
for k in st_plot_grid:
        plt.plot(means_cur[k], noises_cur[k], colors_[k], ms = 10)
for k in st_plot_grid_lna:        
        plt.plot(means_lna[0,k,:]/602.2, noises_lna[0,k,:], colors_solids[k])
plt.title(r'Coefficient of variation vs mean')
plt.xlabel(r'Mean Output Steady State ([$\mu$M])')
plt.ylabel('Coefficient of variation')
plt.xlim([0, 25000/602.2])
plt.ylim([0.02, 2])
plt.show()    
