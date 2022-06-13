The included files reproduce the results in "Learning Anisotropic Interaction Rules from Individual Trajectories in a Heterogeneous Cellular Population", DA Messenger, GE Wheeler, X. Liu, DM Bortz, arXiv:2204.14141, 2022. Example datasets can be found at https://doi.org/10.5281/zenodo.6637041

Users have the option of running the first portion of the algorithm remotely in a trivially parallel manner, or use parfor loops on a local machine. 

------------------------------------------------------------------------------

PREPARING DATA:

(1) POSITION DATA
- collect trajectories into a cell array 'Xscell' with each cell containing data from independent experiments 

- Each cell in Xscell should be an array of dimensions N x d x M for N particles in d dimensions over M (equally spaced) timepoints

- Different cells in Xscell can have different dimensions

(2) TIME POINTS
- collect timepoints into a vector 't'

(3) VELOCITY DATA 
- If possible, supply particle velocities in a cell array 'Vscell' matching the dimensions of Xscell

(4) store variables in .mat folder, which will be loaded as 'input_data' in scripts below.
 
------------------------------------------------------------------------------ 
 
DEFINING INPUTS:

(1) SINGLE-CELL MODEL LEARNING
- in 'singlecell_inputs.m': define data subsampling, velocity computation method, neighbor cell sampling, homing cell sampling, WSINDy parameters, forward simulation parameters

(2) CLASSIFICATION OF MODELS
- in 'classify_inputs.m': define validation error metrics, validation simulation parameters, halting criteria

------------------------------------------------------------------------------

RUNNING THE ALGORITHM: 

**complete PREPARING DATA and DEFINING INPUTS steps above

-------(Local Machine)

(1) Run 'WSINDy_CellCluster_all.m' with lines 3-6 replaced by appropriate paths to data and desired output folder

-------(Remote Machine with SLURM - see subdirectory 'slurm_scripts')

(1) Run parscript_all.sh (e.g. "sbatch slurm_scripts/parscript_all.sh"), which will submit in sequence the following jobs: precomp.sh, single_cell.sh, consoldata.sh, classify.sh, gentraj.sh. Note that data_dr,input_data,save_dr will need to be changed to the relevant files/directories

------------------------------------------------------------------------------

VISUALIZING TRAJECTORIES:

(1) run 'visualize_trajectories.m' to visualize learned trajectories overlapping original trajectories. Change 'species_ind' and 'traj_inds' for different species and specific trajectories.
