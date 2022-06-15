#!/bin/bash
## script suffers dependencies issues, for current use simply copy and past each of the following lines into the terminal successively

JID_precomp=$(sbatch precomp.sh)
JID_singlecell=$(sbatch --dependency=afterok:${JID_precomp##* } single_cell.sh)
JID_consoldata=$(sbatch --dependency=afterok:${JID_singlecell##* } consoldata.sh) 
JID_classify=$(sbatch --dependency=afterok:${JID_consoldata##* } classify.sh) 
JID_gentraj=$(sbatch --dependency=afterok:${JID_classify##* } gentraj.sh) 
