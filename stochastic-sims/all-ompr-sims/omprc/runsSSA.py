import numpy as np
import time, os

import cudasim
import cudasim.EulerMaruyama as EulerMaruyama
import cudasim.Gillespie as Gillespie
import cudasim.Lsoda as Lsoda
import cudasim.SBMLParser as Parser
import pickle
from tqdm import tqdm

##### parameters #####
def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]

def run(parameters_nom, species_nom,ka_grid, beta_st_grid, sim_params, timing = False):

	# Number of equally distributed datapoints to be saved
	xmlModel = sim_params[0]
	fileName = sim_params[1]
	simulationLength = sim_params[2]
	datapoints = sim_params[3]
	n_sim = sim_params[4]

	ind_start = int(np.ceil(250*np.float(datapoints)/simulationLength))
	dt = 1e-2


	# Number of repeated simulation for each set of parameters and species (ignored for ODE)
	n_sim_per_loop = 500 
	n_loops =int( np.ceil(n_sim/np.float(n_sim_per_loop)))

	# Location of the folder for saving the results
	resultFolder = "./results"

	# Location of the temporary folder for model specific CUDA codes
	# temp folder in this file's folder will be created if not specified  
	temp = None
	
	##### initialization #####
	
	# create temp folder
	if(temp == None):
	    temp = os.path.join(os.path.split(os.path.realpath(__file__))[0],"temp")
	    try:
	        os.mkdir(temp)
	    except:
	        pass
	# create result folder
	try:
	    os.mkdir(os.path.realpath(resultFolder))
	except:
	    pass

# setting up parameters for the simulations
	parameters = []
	species = []
	for ka in ka_grid:
		for beta_st in beta_st_grid:
				parameters_changed = parameters_nom+0.0
				parameters_changed[-2]  = beta_st
				parameters_changed[-1]  = ka
		       		parameters.append(parameters_changed+0.0)
				species.append(species_nom+0.0)

	
	# Type of integration
	stats_file = resultFolder + '/stats_' + fileName

	for integrationType in [ "MJP"]: #, "ODE"

	    # Name of the model
	    name = "sRNA" + "_"+integrationType

	    # create CUDA code from SBML model
	    Parser.importSBMLCUDA([xmlModel],[integrationType],ModelName=[name],method=None,outpath=temp)

	    # determining the timepoints for the output
	    timepoints = np.array(range(datapoints+1),dtype=np.float32) * simulationLength/datapoints 

	    # reading in the CUDA code
	    cudaCode = os.path.join(temp, name + ".cu")

	    # create model
	    if(integrationType == "MJP"):
	            result = np.ndarray((len(parameters), n_sim, len(timepoints), len(species[0])))
#		    print 'Running simulations'
		    for ii in tqdm(range(n_loops)):
	                    if ii < n_loops - 1:
	                        beta = n_sim_per_loop
				end_ind = beta*(ii+1)	
		            else:
	                        beta = n_sim - n_sim_per_loop*(n_loops - 1)
				end_ind = n_sim
	                    modeInstance = Gillespie.Gillespie(timepoints, cudaCode, beta=beta, dt=dt)
	                    result_cur = modeInstance.run(parameters, species, timing=timing)
			    for kk in range(result.shape[0]):
			            result[kk][beta*ii:end_ind][:][:] = result_cur[kk]
#	    print 'Writing output'
		    means_  = np.mean(result, axis = 1) 
		    stds_   = np.std(result, axis = 1)
	            noises_ = stds_/means_
		    traces_ = result[:,:20,:,:]
		    np.savez(stats_file, timepoints = timepoints, means = means_, stds = stds_,  noises = noises_, traces=traces_)
