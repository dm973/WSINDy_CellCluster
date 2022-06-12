#!/bin/bash

#SBATCH --job-name=singlecell_000_111
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --output=/projects/dame8201/datasets/111_test/classify_%A_%a.log
#SBATCH --mail-user=dame8201@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --qos=preemptable
#SBATCH --partition=blanca-bortz
#SBATCH --qos=blanca-bortz  

data_dr=/projects/dame8201/datasets/111_test/
input_data=sim08-Feb-2022_000-111_1_comb4000.mat

ml purge
module load matlab/R2019b

cp ${data_dr}*.mat $SLURM_SCRATCH/
cp /projects/dame8201/WSINDy_CellCluster/*.m $SLURM_SCRATCH/

matlab -nodesktop -nodisplay -r "data_dr='${data_dr}';input_data='${input_data}';load(['${SLURM_SCRATCH}/','precomp_',input_data],'expr');single_cell_data=['${SLURM_SCRATCH}/','singlecell_',input_data];run parclassify.m"
