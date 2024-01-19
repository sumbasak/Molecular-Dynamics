#!/bin/bash 

# This script allows taking up residues and extracting their coordinates from the trajectory file 
# without the effect of PBC. The coordinates would indicate a continuous trajectory of the selected 
# residue without any abysmal jumps near the PBC.

function process_residue() { 

    local residue_number=$1 
    local index_file="index_${residue_number}SOL.ndx" 
    local output_file="${residue_number}SOL_nopbc.gro" 

    echo -e "a OW HW1 HW2 LP1 LP2 & r $residue_number 
    q" | gmx make_ndx -f npt.gro -o "$index_file" 

    echo "7 
    7" | gmx trjconv -f traj_comp.xtc -s finalmd.tpr -n "$index_file" -o "$output_file" -pbc nojump 
}

# Process residues
process_residue 12061
process_residue 1535
process_residue 12101
process_residue 3679
process_residue 11943


