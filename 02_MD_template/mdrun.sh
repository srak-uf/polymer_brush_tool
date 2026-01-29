#!/bin/bash
#$ -cwd
#$ -l cpu_40=1
#$ -l h_rt=8:00:00


source ~/.bashrc
module purge
module load intel/2024.0.2 intel-mpi/2021.11
source /home/4/uf02194/tate_ssd/sasakir/apl/gmx_build/bin/GMXRC.bash

rm *trr *edr *log *tpr *xtc *pdb mdout.mdp

gmx_mpi grompp -f npt.mdp -c nvt.gro -p grafted_chain_water_restraint.top -o npt.tpr -r grafted_chain_water_box.gro -maxwarn 1 
mpirun -np 36 gmx_mpi mdrun -ntomp 1 -deffnm npt -v
