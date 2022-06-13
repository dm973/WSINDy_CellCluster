#!/bin/bash

#SBATCH --job-name=singlecell_000_111
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --output=/projects/dame8201/datasets/111_test/consoldata_%A_%a.log
#SBATCH --mail-user=dame8201@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --qos=preemptable
#SBATCH --partition=blanca-bortz
#SBATCH --qos=blanca-bortz  

data_dr=/projects/dame8201/datasets/111_test/
input_data=sim08-Feb-2022_000-111_1_comb4000.mat

ml purge
module load matlab/R2019b

cp /projects/dame8201/WSINDy_CellCluster/*.m $SLURM_SCRATCH/
cp ${data_dr}precomp_${input_data} $SLURM_SCRATCH/
cp ${data_dr}single-cell-models/home_cell*.mat $SLURM_SCRATCH/

matlab -nodesktop -nodisplay -r "save_dr='${data_dr}'; data_dr='${data_dr}'; input_data='${input_data}'; homecell_dr='${SLURM_SCRATCH}/'; load(['${SLURM_SCRATCH}/','precomp_',input_data]); run ${SLURM_SCRATCH}/parconsoldata.m;"
