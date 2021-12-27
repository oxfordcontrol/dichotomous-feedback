import numpy as np
import runsSSA
import os
from tqdm import tqdm

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-s", "--datapoints", dest = 'datapoints',  default = 300, type=int, 
                    help="number of saved time snaphsots of a trajectory")
parser.add_argument("-n", "--nsim", dest = "n_sim", default= 500, type=int,
                    help="number of saved points in a trajectory")
parser.add_argument("-T", "--Tlen", dest = "simulationLength", default= 300, type=int, 
		    help="simulation length in minutes")
parser.add_argument("-kt", "--kt", dest = "kt", nargs='+',  default= [6.12*1e+3],  type = float,    
		    help="Phosphorylation rate ")
parser.add_argument("-kp", "--kp", dest = "kp", nargs='+', default= [4*1e-3], type = float, 
		    help="De-phosphorylation rate ")
parser.add_argument("-ht", "--hk_tot", dest = "hk_tot", nargs='+', default= [100], type = float,   
		    help="HK copy number ")
parser.add_argument("-rt", "--rr_tot", dest = "rr_tot", nargs='+', default= [3500], type = float, 
		    help="RR copy number ")
parser.add_argument("-sg", "--s_grid", dest = "s_grid", nargs='+', default= [0, 0.1, 0.25, 0.5, 0.75, 1], type = float, 
		    help="SR copy number grid ")
parser.add_argument("-st", "--sink_tot", dest = "sink_tot", nargs='+', default= [3500], type = float, 
		    help="Sink copy number ")
parser.add_argument("-bg", "--beta_gfp", dest = "beta_gfp", nargs='+', default= [1], type = float, 
		    help="GFP production rate ")
parser.add_argument("-bs", "--beta_sr", dest = "beta_sr", nargs='+', default= [1], type = float, 
		    help="OmpRc production rate ")
parser.add_argument("-kd", "--Kd", dest = "Kd", nargs='+', default= [0.0314], type = float, 
		    help="pOmpc dissassociation constant")
parser.add_argument("-hc", "--hill_coeff", dest = "hill_coeff", nargs='+', default= [4], type = float, 
		    help="pOmpc dissassociation constant")
parser.add_argument("-dp", "--dprot", dest = "dp",  nargs='+', default= [0.0234], type = float, 
		    help="protein dilution rate ")
parser.add_argument("-as", "--asp", dest = "asp", nargs='+', default= [0, 1, 5, 10], type = float, 
		    help="Aspartate concentration ")
args = parser.parse_args()

# simulation parameters
datapoints = args.datapoints
n_sim = args.n_sim
simulationLength = args.simulationLength
# model Parameters
cell_volume =  602.2
k_t = 1/cell_volume*np.asarray(args.kt)[0]             # 1/(min \mu M)
k_p = 1/cell_volume*np.asarray(args.kp)[0]             # 1/(min \mu M)
beta_gfp    = cell_volume*np.asarray(args.beta_gfp)[0] # \mu M/min
hill_coeff  = np.asarray(args.hill_coeff)[0]           # dimensionless
Kd          = cell_volume*np.asarray(args.Kd)[0]       # \mu M
d_prot = np.asarray(args.dp)[0]                        # 1/min
hkt = np.asarray(args.hk_tot)[0]                       # number of molecules
rrt = np.asarray(args.rr_tot)[0]                       # number of molecules
st  = np.asarray(args.sink_tot)[0]                     # number of molecules
beta_sr = cell_volume*np.asarray(args.beta_sr)[0]      # number of molecules
beta_hk = hkt*d_prot                                   # number of molecules/min
beta_rr = rrt*d_prot                                   # number of molecules/min
asp = np.asarray(args.asp)                             # mM
Kda = 2                                                # mM
kap_max = 0.02                                         # 1/min
s_grid = np.asarray(args.s_grid)

ka_grid = kap_max * asp/(asp + Kda)
beta_st_grid = s_grid*st*d_prot                               # number of molecules/min
beta_sr_grid = s_grid*beta_sr
print beta_st_grid
print beta_sr_grid
out = open("../log.txt",'a')
print >>out, 'Parameters for sRNA overexpression simulations  \n'
print >>out, args
print >>out, '\n'
out.close()

# computing stats for the open loop model
xmlModel = "Tazc_ol.xml"
species_nom = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])


#parameters_nom = np.array([1.0, beta_gfp, hill_coeff, Kd, k_t, k_p, k_p, beta_rr, beta_hk, d_prot,  beta_tazc, k_a])
parameters_nom = np.array([1.0, beta_gfp, hill_coeff, Kd, k_t, k_p, k_p, beta_rr, beta_hk, d_prot, 0,0])
fileName = "tazc_ol"
print '    Preparing simulations for the open loop system'
sim_params = list([xmlModel, fileName, simulationLength, datapoints, n_sim])
runsSSA.run(parameters_nom, species_nom, ka_grid, beta_st_grid, sim_params)
	
# computing stats for the closed loop model
xmlModel = "Tazc_cl.xml"
species_nom = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

#parameters_nom = np.array([1.0, beta_gfp, hill_coeff, Kd, k_t, k_p,  k_p, beta_rr, beta_hk, d_prot,  beta_tazc, k_a])
parameters_nom = np.array([1.0, beta_gfp, hill_coeff, Kd, k_t, k_p,  k_p, beta_rr, beta_hk, d_prot, 0,0])
fileName = "tazc_cl"
print '    Preparing simulations for the closed loop system'
sim_params = list([xmlModel, fileName, simulationLength, datapoints, n_sim])
runsSSA.run(parameters_nom, species_nom, ka_grid, beta_sr_grid, sim_params)
