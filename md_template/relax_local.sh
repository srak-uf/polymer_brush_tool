gmx grompp -f min.mdp -c grafted_chain_water_box.gro -p hardrest_grafted_chain_restraint.top -r grafted_chain_water_box.gro -o min.tpr
gmx mdrun -deffnm min -ntmpi 8 -ntomp 1 -v

gmx grompp -f nvt.mdp -c min.gro -p hardrest_grafted_chain_restraint.top -r grafted_chain_water_box.gro -o nvt.tpr
gmx mdrun -deffnm nvt -ntmpi 8 -ntomp 1 -v
