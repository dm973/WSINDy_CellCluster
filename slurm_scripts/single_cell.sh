#!/bin/bash
#SBATCH --job-name=singlecell_000_111
##SBATCH --account=blanca-dame8201
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=00:15:00
#SBATCH --output=/projects/dame8201/datasets/111_test/single-cell-models/home_cell_%A-%a.log
#SBATCH --mail-user=dame8201@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL                                                                   
#SBATCH --qos=preemptable                                                                            
#SBATCH --array=1-100   ## set to number of homing cells

ml purge
module load matlab/R2019b

data_dr=/projects/dame8201/datasets/111_test/
input_data=sim08-Feb-2022_000-111_1_comb4000.mat
precomp_data=precomp_${input_data}
save_dr=${data_dr}single-cell-models/

cp /projects/dame8201/WSINDy_CellCluster/*.m $SLURM_SCRATCH/
cp $data_dr$precomp_data $SLURM_SCRATCH/

matlab -nodesktop -nodisplay -r "i1=${SLURM_ARRAY_TASK_ID}; load('${SLURM_SCRATCH}/${precomp_data}'); run $SLURM_SCRATCH/parsingle_cell.m; save(['${save_dr}','home_cell_',num2st\
r(i1),'.mat'],'algouttemp','Xspred','errs','valid_cells','nu_learned')"
