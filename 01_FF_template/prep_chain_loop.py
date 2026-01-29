import os
from ase.io import read, write
import numpy as np
import parmed as pmd

# Setup of T-MPI, OMP
t_mpi = 8
o_mpi = 1

# x-y ratio of one chain
xyratio = 0.5  # box_y / box_x default 1

flag_loop = True

# Definition of HEAD monomer
# antechamber -fi gout -fo pdb -i mpc_pcm_resp.log -o mpc.pdb -c resp -pf y -at gaff2
# antechamber -fi gout -fo ac -i mpc_pcm_resp.log -o mpc.ac -c resp -pf y -at gaff2
head_tailname = "C2"
head_omitnames = ["H24"]
head_post_tailtype = "c3"
head_termname = "H1"
head_resname = "hmp"
head_acfile = "mpc.ac"
# Number of C-C from head to tail in HEAD monomer
head_n_cc = 1

# Definition of MID monomer
mid_headname = "C11"
mid_tailname = "C2"
mid_omitnames = ["H23", "H24"]
mid_pre_headtype = "c3"
mid_post_tailtype = "c3"
n_mid_repeat_units = 26
mid_resname = "mmp"
mid_acfile = "mpc.ac"
mid_n_cc = 1

# Definition of TAIL monomer
tail_headname = "C11"
tail_omitnames = ["H23"]
tail_pre_headtype = "c3"
tail_termname = "H24"
tail_resname = "tmp"
tail_acfile = "mpc.ac"
tail_n_cc = 1

# Graft density
rho = 0.45 / 2  # chains/nm^2
nx = 1  # number of chains in x direction
ny = 2  # number of chains in y direction
box_x = np.sqrt(nx * ny / rho) * 10  # A
box_y = np.sqrt(nx * ny / rho) * 10  # A

# Polymer chain length
d_polymer = 14.9  # nm  or  None  HEAD----(MID)n----TAIL length
if d_polymer is None:
    n_cc_all = (head_n_cc+1) + n_mid_repeat_units * (mid_n_cc+1) + (tail_n_cc+1) + 2
    d_cc = 1.54  # C-C bond length in Angstrom
    d_polymer = d_cc * n_cc_all * 0.8

if flag_loop:
    n_cc_all = (head_n_cc+1) + n_mid_repeat_units * (mid_n_cc+1) + (tail_n_cc+1) + 2
    d_cc = 1.54  # C-C bond length in Angstrom
    d_height = d_cc * n_cc_all * 0.8 / 2
############################


def write_fragment(resname, headname=None, tailname=None, omitnames=None,
                   pre_headtype=None, post_tailtype=None):
    n_omits = 0
    with open(f"{resname}.chain", mode="w") as f:
        if headname is not None:
            n_omits += 1
            if pre_headtype is None:
                raise ValueError("pre_headtype must be provided if headname is provided.")
            else:
                f.write(f"HEAD_NAME  {headname}\n")
                f.write(f"PRE_HEAD_TYPE  {pre_headtype}\n")

        if tailname is not None:
            n_omits += 1
            if post_tailtype is None:
                raise ValueError("post_tailtype must be provided if tailname is provided.")
            else:
                f.write(f"TAIL_NAME  {tailname}\n")
                f.write(f"POST_TAIL_TYPE  {post_tailtype}\n")

        if omitnames is not None:
            if isinstance(omitnames, str):
                omitnames = [omitnames]
            elif not isinstance(omitnames, list):
                raise TypeError("omitnames must be a string or a list of strings.")

            if isinstance(omitnames, list):
                if len(omitnames) == n_omits:
                    pass
                else:
                    raise ValueError("Length of omitnames does not match number of head/tail names provided.")
            for omitname in omitnames:
                f.write(f"OMIT_NAME  {omitname}\n")

        f.write("CHARGE 0.000\n")


def write_tleap(head_resname, mid_resname, tail_resname, n_mid_repeat_units):
    with open("build_chain.tleap", mode="w") as f:
        f.write("source leaprc.gaff2\n")
        f.write(f"loadamberprep {head_resname}.prepi\n")
        f.write(f"loadamberprep {mid_resname}.prepi\n")
        f.write(f"loadamberprep {tail_resname}.prepi\n")
        f.write(f"chain = sequence {{{head_resname.upper()} ")
        for _ in range(n_mid_repeat_units):
            f.write(f"{mid_resname.upper()} ")
        f.write(f" {tail_resname.upper()}}}\n")
        f.write("saveAmberParm chain chain.prmtop chain.inpcrd\n")
        f.write("savepdb chain chain.pdb\n")
        f.write("quit\n")


def amber_min_with_pull(prmtop, inpcrd,
                        file_prefix,
                        head_idx,
                        tail_idx,
                        polymer_length,
                        loop_height=None):
    with open("pull_termination.restraint", mode="w") as f:
        r3 = polymer_length
        r2 = 0.8 * r3
        r4 = 1.2 * r3
        f.write("&rst\n")
        f.write(f"iat={head_idx},{tail_idx},r1=0., r2={r2}, r3={r3}, r4={r4}, rk2=1000.7, rk3=1000.7,/\n")
        f.write("&end\n")
        if loop_height is not None:
            r3_loop = loop_height
            r2_loop = 0.8 * r3_loop
            r4_loop = 1.2 * r3_loop
            mid_idx = (head_idx + tail_idx) // 2
            f.write("&rst\n")
            f.write(f"iat={mid_idx},{tail_idx},r1=0., r2={r2_loop}, r3={r3_loop}, r4={r4_loop}, rk2=1000.7, rk3=1000.7,/\n")
            f.write("&end\n")
            f.write("&rst\n")
            f.write(f"iat={head_idx},{mid_idx},r1=0., r2={r2_loop}, r3={r3_loop}, r4={r4_loop}, rk2=1000.7, rk3=1000.7,/\n")
            f.write("&end\n")

    with open(f"{file_prefix}_min_pull.in", mode="w") as f:
        f.write("Minimize\n")
        f.write("&cntrl\n")
        f.write("imin=1,\n")
        f.write("ntb=0,\n")
        f.write("ntx=1,\n")
        f.write("irest=0,\n")
        f.write("maxcyc=50000,\n")
        f.write("ncyc=1000,\n")
        f.write("ntpr=100,\n")
        f.write("ntwx=0,\n")
        f.write("cut=999.0,\n")
        f.write("nmropt=1,\n")
        f.write("/\n")
        f.write("&wt type='END' /\n")
        f.write("DISANG=pull_termination.restraint\n")

    cmd = (f"sander -O -i {file_prefix}_min_pull.in -o {file_prefix}_min_pull.out -p {prmtop} -c {inpcrd} -r {file_prefix}_min_pull.rst7")
    print(f"EXEC: {cmd}")
    os.system(cmd)
    # ambpdb
    cmd_ambpdb = (f"ambpdb -p {prmtop} -c {file_prefix}_min_pull.rst7 > {file_prefix}_min_pull.pdb")
    print(f"EXEC: {cmd_ambpdb}")
    os.system(cmd_ambpdb)


def graft_brush(chain_pdb, box_x, box_y, nx, ny):
    atoms = read(chain_pdb)
    min_z = atoms.get_positions()[:,2].min()
    max_z = atoms.get_positions()[:,2].max()
    box_z = max_z - min_z + 30.0  # A
    with open("graft_brush.in", mode="w") as f:
        f.write("tolerance 2.0\n")
        f.write("output grafted_chain.pdb\n")
        f.write("filetype pdb\n")
        f.write(f"pbc {box_x} {box_y} {box_z}\n")
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
        id = 0
        for ix in range(nx):
            for iy in range(ny):
                f.write(f"structure {chain_pdb}\n")
                f.write("  resnumbers 1\n")
                f.write("  number 1\n")
                f.write(f"  chain {alphabet[id]}\n")
                x = (ix + 0.5) * (box_x / nx)
                y = (iy + 0.5) * (box_y / ny)
                f.write(f"  fixed {x} {y} 0.0 0.0 0.0 0.0\n")
                f.write("end structure\n")
                f.write("\n")
                id += 1
    cmd = "packmol < graft_brush.in > graft_brush.out"
    print(f"EXEC: {cmd}")
    os.system(cmd)


def insert_tip3p_top(gmx_top, output_top):
    with open(gmx_top, 'r') as f:
        lines = f.readlines()

    with open(output_top, 'w') as f:
        for i, line in enumerate(lines):
            f.write(line)
            if i < len(lines) - 1 and lines[i+1].strip() == '[ moleculetype ]':
                f.write('#include "tip3p.itp"\n')
                f.write('\n')


def insert_restraint_top(gmx_top, output_top, linker_indices):
    with open(gmx_top, 'r') as f:
        lines = f.readlines()

    # [ atoms ]からindicesのリストを作成
    all_atom_indices = []
    in_atoms_section = False
    for line in lines:
        if line.strip() == '[ atoms ]':
            in_atoms_section = True
            continue
        if in_atoms_section:
            if line.strip() == '' or line.strip().startswith('['):
                break
            if line.strip().startswith(';'):
                continue
            else:
                parts = line.split()
                if len(parts) >= 1:
                    all_atom_indices.append(int(parts[0]))

    with open("hardrest_" + output_top, 'w') as f:
        for i, line in enumerate(lines):
            f.write(line)
            if i < len(lines) - 1 and lines[i+1].strip() == '[ system ]':
                f.write('\n[ position_restraints ]\n')
                f.write('; atom  type      fx      fy      fz\n')
                for idx in linker_indices:
                    if idx in all_atom_indices:
                        f.write(f' {idx}  1  1000000  1000000  1000000\n')
                f.write('\n')

    with open(output_top, 'w') as f:
        for i, line in enumerate(lines):
            f.write(line)
            if i < len(lines) - 1 and lines[i+1].strip() == '[ system ]':
                f.write('\n[ position_restraints ]\n')
                f.write('; atom  type      fx      fy      fz\n')
                for idx in linker_indices:
                    if idx in all_atom_indices:
                        f.write(f' {idx}  1  10000  10000  10000\n')
                f.write('\n')


write_fragment(head_resname, tailname=head_tailname,
               omitnames=head_omitnames,
               post_tailtype=head_post_tailtype)

write_fragment(mid_resname, headname=mid_headname, tailname=mid_tailname,
               omitnames=mid_omitnames,
               pre_headtype=mid_pre_headtype, post_tailtype=mid_post_tailtype)

write_fragment(tail_resname, headname=tail_headname,
               omitnames=tail_omitnames,
               pre_headtype=tail_pre_headtype)

print("*******************************************")
print("** Generating HEAD monomer prepi file... **")
print("*******************************************")
cmd_head = (f"prepgen -i {head_acfile} -o {head_resname}.prepi "
            f"-f prepi -m {head_resname}.chain "
            f"-rn {head_resname.upper()} -rf {head_resname}.res")
print(f"EXEC: {cmd_head}")
os.system(cmd_head)


print("******************************************")
print("** Generating MID monomer prepi file... **")
print("******************************************")
cmd_mid = (f"prepgen -i {mid_acfile} -o {mid_resname}.prepi "
           f"-f prepi -m {mid_resname}.chain "
           f"-rn {mid_resname.upper()} -rf {mid_resname}.res")
print(f"EXEC: {cmd_mid}")
os.system(cmd_mid)

print("*******************************************")
print("** Generating TAIL monomer prepi file... **")
print("*******************************************")
cmd_tail = (f"prepgen -i {tail_acfile} -o {tail_resname}.prepi "
            f"-f prepi -m {tail_resname}.chain "
            f"-rn {tail_resname.upper()} -rf {tail_resname}.res")
print(f"EXEC: {cmd_tail}")
os.system(cmd_tail)

print("*********************")
print("** Executing tleap **")
print("*********************")
write_tleap(head_resname, mid_resname, tail_resname, n_mid_repeat_units)
os.system("tleap -f build_chain.tleap")

write("chain.pdb", read("chain.pdb"))

atoms_chain = read("chain.pdb")
get_terminal_index = lambda termname, resname: [i for i, (atype, rname) in enumerate(zip(atoms_chain.arrays["atomtypes"], atoms_chain.arrays["residuenames"])) if atype == termname and rname == resname.upper()][0]

head_idx = get_terminal_index(head_termname, head_resname)
tail_idx = get_terminal_index(tail_termname, tail_resname)

print("********************************************************")
print("** Performing AMBER minimization with pull constraint **")
print("********************************************************")
if flag_loop:
    amber_min_with_pull("chain.prmtop", "chain.inpcrd",
                        "chain",
                        head_idx + 1,
                        tail_idx + 1,
                        d_polymer,
                        loop_height=d_height)
else:
    amber_min_with_pull("chain.prmtop", "chain.inpcrd",
                        "chain",
                        head_idx + 1,
                        tail_idx + 1,
                        d_polymer)

print("** Align z-axis **")
# 2.1から100まで0.1刻みでリストを作成
scan_x = np.arange(2.1, 100.1, 0.1)

print("Please type bottom atom index... ")
print("(You can check structure by vesta chain_min_pull.pdb)")
bottom_atom_index = int(input())

for x in scan_x:
    with open("align_z.in", mode="w") as f:
        f.write("tolerance 2.0\n")
        f.write("output aligned_chain.pdb\n")
        f.write("filetype pdb\n")
        f.write("add_amber_ter\n")
        f.write("structure chain_min_pull.pdb\n")
        f.write("  number 1\n")
        f.write(f"  inside box 0. 0. 0. {x} {x * xyratio} 1000.\n")
        f.write(f"  atoms {bottom_atom_index}\n")
        f.write(f"  inside box 0. 0. 1.2 {x} {x * xyratio} 1.2\n")
        f.write("  end atoms\n")
        f.write("end structure\n")
    print(f"The new x coordinate is: {x}")
    os.system("packmol < align_z.in > align_z.out")
    # align_z.outに”ERROR"が含まれていなければループを抜ける
    with open("align_z.out", mode="r") as f:
        align_z_out = f.read()
        if "ERROR" not in align_z_out:
            break

print("** Grafting polymer brushes... **")
graft_brush("aligned_chain.pdb", box_x, box_y, nx, ny)


print("*********************************")
print("** FF assign for grafted brush **")
print("*********************************")
grafted_chain_atoms = read("grafted_chain.pdb")
La = grafted_chain_atoms.cell[0, 0]
Lb = grafted_chain_atoms.cell[1, 1]
Lc = grafted_chain_atoms.cell[2, 2]

with open("grafted_chain.tleap", mode="w") as f:
    f.write("source leaprc.gaff2\n")
    f.write(f"loadamberprep {head_resname}.prepi\n")
    f.write(f"loadamberprep {mid_resname}.prepi\n")
    f.write(f"loadamberprep {tail_resname}.prepi\n")
    f.write("mol = loadpdb grafted_chain.pdb\n")
    f.write(f"set mol box {{{La} {Lb} {Lc}}}\n")
    f.write("saveAmberParm mol grafted_chain.prmtop grafted_chain.inpcrd\n")
    f.write("quit\n")

os.system("tleap -f grafted_chain.tleap")

parm = pmd.load_file("grafted_chain.prmtop", xyz="grafted_chain.inpcrd")
parm.save("grafted_chain.top", overwrite=True)
parm.save("grafted_chain.gro", overwrite=True)

# Centering
os.system("gmx editconf -f grafted_chain.gro -c -o grafted_chain_center.gro")

print("******************************************")
print("** Shifting grafted chain to the bottom **")
print("******************************************")
print("")
print("Default: Linker RESNAME, ATOMNAME / ", head_resname.upper(), head_termname)
print("Type RESNAME1, ATOMNAME1, RESNAME2, ATOMNAME2, ... or type 'y' to accept default")
grafted_chain_atoms = read("grafted_chain_center.gro")
ans = input()
linker_resname = []
linker_atomname = []
if ans.lower() != "y":
    linker_info = ans.split(",")
    for i in range(0, len(linker_info), 2):
        linker_resname.append(linker_info[i].strip())
        linker_atomname.append(linker_info[i+1].strip())
else:
    linker_resname.append(head_resname.upper())
    linker_atomname.append(head_termname)

# grafted_chain_atomsでlinker_resnameとlinker_atomnameに該当する原子のz座標を最小値に揃える
linker_indices = []
for resname, atomname in zip(linker_resname, linker_atomname):
    indices = [i for i, (aname, rname) in enumerate(zip(grafted_chain_atoms.arrays["atomtypes"], grafted_chain_atoms.arrays["residuenames"])) if aname == atomname and rname == resname]
    linker_indices.extend(indices)
min_z = grafted_chain_atoms.get_positions()[:, 2].min()
for i in linker_indices:
    grafted_chain_atoms.positions[i, 2] = min_z
linker_indices = [i + 1 for i in linker_indices]
print("Linker atoms indices:", linker_indices)

grafted_chain_atoms.positions[:, 2] -= min_z
write("grafted_chain_shifted.gro", grafted_chain_atoms)
insert_restraint_top("grafted_chain.top",
                     "grafted_chain_restraint.top",
                     linker_indices)

os.system("gmx grompp -f min_vac.mdp -p hardrest_grafted_chain_restraint.top -c grafted_chain_shifted.gro -o min_vac.tpr -r grafted_chain_shifted.gro -maxwarn 2")
os.system(f"gmx mdrun -deffnm min_vac -ntmpi {t_mpi} -ntomp {o_mpi}")

os.system("gmx grompp -f nvt_vac.mdp -p hardrest_grafted_chain_restraint.top -c min_vac.gro -o nvt_vac.tpr -r grafted_chain_shifted.gro -maxwarn 2")
os.system(f"gmx mdrun -deffnm nvt_vac -ntmpi {t_mpi} -ntomp {o_mpi}")

os.system("gmx solvate -cp nvt_vac.gro -p hardrest_grafted_chain_restraint.top -o grafted_chain_water.gro")
atoms = read("grafted_chain_water.gro")
Lc = (atoms.cell[2, 2] + 4) / 10  # nm
La = atoms.cell[0, 0]/10  # nm
Lb = atoms.cell[1, 1]/10  # nm
os.system(f"gmx editconf -f grafted_chain_water.gro -box {La} {Lb} {Lc} -o grafted_chain_water_box.gro")

insert_tip3p_top("grafted_chain_restraint.top", "grafted_chain_water_restraint.top")
insert_tip3p_top("hardrest_grafted_chain_restraint.top",
                 "hardrest_grafted_chain_water_restraint.top")

print("Necessary files for GROMACS simulation: ")
print("  - grafted_chain_water_box.gro")
print("  - grafted_chain_water_restraint.top")
print("  - hardrest_grafted_chain_water_restraint.top")
print("  - tip3p.itp")
