#!/bin/bash

#SBATCH --job-name=singlecell_000_111
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --output=/projects/dame8201/datasets/111_test/wsindy_cellcluster_%A_%a.log
#SBATCH --mail-user=dame8201@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --qos=preemptable

data_dr=/projects/dame8201/datasets/111_test/
input_data=sim08-Feb-2022_000-111_1_comb4000.mat
save_dr=${data_dr}single-cell-models/

ml purge
module load matlab/R2019b

cp /projects/dame8201/WSINDy_CellCluster/*.m $SLURM_SCRATCH/

matlab -nodesktop -nodisplay -r "data_dr='${data_dr}'; input_data='${input_data}'; load([data_dr,input_data],'Xscell','Vscell','t'); run parprecomp.m;"
