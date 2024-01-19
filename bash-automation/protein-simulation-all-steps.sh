#!/bin/bash

# This script goes to a specific folder and performs a set of simulations for 1 ns. 
# Starting from protein structure download from the data bank, it goes on to simulate NVT,
# NPT, and then, finally the trajectory. All the RMSD values calculated for a range of 
# temperatures are placed in one file. 

file_path=/home/sumon/MEGA/tutorials/pradiptabash/

cd "$file_path"

# Download the pdb file from the RCSB server
wget https://files.rcsb.org/download/1VII.pdb

# Clean the water from the pdb file
grep -v HOH 1VII.pdb > 1VII_clean.pdb

# Process the file from pdb to gmx
echo -e "15" | gmx pdb2gmx -f 1VII_clean.pdb -o 1VII_processed.gro -water spce -ignh

# Get the protein inside a box
gmx editconf -f 1VII_processed.gro -o 1VII_newbox.gro -c -d 1.0 -bt cubic

# Get the protein solvated with water
gmx solvate -cp 1VII_newbox.gro -cs spc216.gro -o 1VII_solv.gro -p topol.top

# Copy the ion.mdp file from my other directory
cp /home/sumon/MEGA/tutorials/lysozymeInWater/*.mdp .

# Get ions solvated in the water-protein system
gmx grompp -f ions.mdp -c 1VII_solv.gro -p topol.top -o ions.tpr
echo -e "13" | gmx genion -s ions.tpr -o 1VII_solv_ions.gro -p topol.top -pname NA -nname CL -neutral 

# Minimize the energy of the system
gmx grompp -f minim.mdp -c 1VII_solv_ions.gro -p topol.top -o em.tpr 
gmx mdrun -v -deffnm em 

# NVT simulation
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT simulation
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# Change the temperature inside to a known phrase
sed -i 's/300/sbk/g' md.mdp

# Run the loop to copy all the files and run individual simulations
for temp in {300..320..5}; do 

    # Make a new directory
    mkdir "$temp"

    # Copy the requisite files to that new directory
    cp md.mdp npt.cpt npt.gro topol.top $temp/  

    # Navigate to that directory
    cd "$temp"/

    # Change sbk with the temperature
    sed -i 's/sbk/'"$temp"'/g' md.mdp

    # Final mdrun for that specific temperature
    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr 
    gmx mdrun -deffnm md_0_1

    # Centering the protein inside the box
    echo -e "1 0" | gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center 
    
    # Calculation of rmsd values for all the temperatures
    echo -e "4
    1" | gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns 

    cd "$file_path"
done

# Output file
output_file="all_rmsd.xvg"

# Get the time column from the first folder and create a temporary file
awk '{print $1}' "300/rmsd.xvg" > time_column.tmp

# Loop through the folders and extract RMSD columns
for temp in {300..320..5}; do
    current_file="${temp}/rmsd.xvg"
    # Use awk to extract RMSD columns starting from the second column
    awk '{print $2}' "${current_file}" > "rmsd_${temp}.tmp"
done

# Combine the time column and all RMSD columns into the final output file
paste time_column.tmp rmsd_*.tmp > "${output_file}"

# Remove temporary files
rm time_column.tmp rmsd_*.tmp
