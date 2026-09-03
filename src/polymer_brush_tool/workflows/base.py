"""Shared workflow logic for polymer brush system preparation."""

from __future__ import annotations

import shutil
from pathlib import Path

from polymer_brush_tool import runner
from polymer_brush_tool.config import BrushConfig
from polymer_brush_tool.ff import fragments, minimize, topology
from polymer_brush_tool.structure import atoms as atom_utils
from polymer_brush_tool.structure import graft


class BrushWorkflowBase:
    """Base class shared by :class:`LinearBrushWorkflow` and :class:`LoopBrushWorkflow`.

    Parameters
    ----------
    config:
        Complete configuration object for the simulation system.
    work_dir:
        Directory where all intermediate and output files are written.
        Created if it does not exist.
    mdp_dir:
        Directory containing the GROMACS MDP template files
        (``min_vac.mdp``, ``nvt_vac.mdp``).  Defaults to the
        ``01_FF_template`` directory next to this package.
    """

    def __init__(
        self,
        config: BrushConfig,
        work_dir: str | Path = ".",
        mdp_dir: str | Path | None = None,
    ) -> None:
        self.config = config
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        if mdp_dir is None:
            # Default: 01_FF_template next to the installed package tree
            pkg_root = Path(__file__).parent.parent.parent.parent
            mdp_dir = pkg_root / "01_FF_template"
        self.mdp_dir = Path(mdp_dir)

    # ------------------------------------------------------------------
    # Step 1 – Fragment and prepi files
    # ------------------------------------------------------------------

    def step_prepgen(self) -> None:
        """Write ``.chain`` files and run ``prepgen`` for HEAD, MID, and TAIL."""
        cfg = self.config
        head, mid, tail = cfg.head, cfg.mid, cfg.tail

        print("\n=== Step 1: Generating prepi files via prepgen ===")

        # HEAD – only a tail side (bonds to MID)
        fragments.write_fragment(
            head.resname,
            tailname=head.tailname,
            omitnames=head.omitnames,
            post_tailtype=head.post_tailtype,
            work_dir=self.work_dir,
        )

        # MID – both head and tail sides
        fragments.write_fragment(
            mid.resname,
            headname=mid.headname,
            tailname=mid.tailname,
            omitnames=mid.omitnames,
            pre_headtype=mid.pre_headtype,
            post_tailtype=mid.post_tailtype,
            work_dir=self.work_dir,
        )

        # TAIL – only a head side (bonds to MID)
        fragments.write_fragment(
            tail.resname,
            headname=tail.headname,
            omitnames=tail.omitnames,
            pre_headtype=tail.pre_headtype,
            work_dir=self.work_dir,
        )

        for label, spec in [("HEAD", head), ("MID", mid), ("TAIL", tail)]:
            print(f"\n** Generating {label} monomer prepi file... **")
            runner.run(
                f"prepgen -i {spec.ac_file}"
                f" -o {spec.resname}.prepi"
                f" -f prepi"
                f" -m {spec.resname}.chain"
                f" -rn {spec.resname.upper()}"
                f" -rf {spec.resname}.res",
                cwd=self.work_dir,
            )

    # ------------------------------------------------------------------
    # Step 2 – Build chain with tleap
    # ------------------------------------------------------------------

    def step_build_chain(self) -> None:
        """Write tleap script and build the polymer chain topology."""
        cfg = self.config
        print("\n=== Step 2: Building polymer chain with tleap ===")
        fragments.write_tleap(
            cfg.head.resname,
            cfg.mid.resname,
            cfg.tail.resname,
            cfg.n_mid_repeat_units,
            work_dir=self.work_dir,
        )
        runner.run("tleap -f build_chain.tleap", cwd=self.work_dir)

        # Re-write PDB to normalise formatting (ASE round-trip)
        from ase.io import read, write
        chain_pdb = self.work_dir / "chain.pdb"
        write(str(chain_pdb), read(str(chain_pdb)))

    # ------------------------------------------------------------------
    # Step 3 – AMBER minimisation (implemented by subclass)
    # ------------------------------------------------------------------

    def step_amber_minimize(self) -> None:  # pragma: no cover
        raise NotImplementedError

    # ------------------------------------------------------------------
    # Step 4 – Z-axis alignment (interactive if index not in config)
    # ------------------------------------------------------------------

    def step_align_z(self) -> None:
        """Align the minimised chain along z using PACKMOL."""
        cfg = self.config
        print("\n=== Step 4: Z-axis alignment via PACKMOL ===")

        bottom_idx = cfg.bottom_atom_index
        if bottom_idx is None:
            print(
                "Please enter the bottom atom index "
                "(inspect chain_min_pull.pdb in VESTA, 1-based index):"
            )
            bottom_idx = int(input())

        graft.align_chain_z(
            self.work_dir / "chain_min_pull.pdb",
            bottom_idx,
            xyratio=cfg.xyratio,
            work_dir=self.work_dir,
        )

    # ------------------------------------------------------------------
    # Step 5 – Graft brush
    # ------------------------------------------------------------------

    def step_graft(self) -> None:
        """Place nx × ny chains on a 2-D grid with PACKMOL."""
        cfg = self.config
        print("\n=== Step 5: Grafting polymer brushes ===")
        graft.graft_brush(
            self.work_dir / "aligned_chain.pdb",
            cfg.box_x(),
            cfg.box_y(),
            cfg.nx,
            cfg.ny,
            work_dir=self.work_dir,
        )

    # ------------------------------------------------------------------
    # Step 6 – Force-field assignment for grafted system
    # ------------------------------------------------------------------

    def step_assign_ff_grafted(self) -> None:
        """Run tleap on the multi-chain PDB to assign GAFF2 parameters."""
        from ase.io import read
        import parmed as pmd

        cfg = self.config
        print("\n=== Step 6: Force-field assignment for grafted brush ===")

        grafted_pdb = self.work_dir / "grafted_chain.pdb"
        atoms = read(str(grafted_pdb))
        La, Lb, Lc = (
            atoms.cell[0, 0],
            atoms.cell[1, 1],
            atoms.cell[2, 2],
        )

        fragments.write_tleap_grafted(
            cfg.head.resname,
            cfg.mid.resname,
            cfg.tail.resname,
            La, Lb, Lc,
            work_dir=self.work_dir,
        )
        runner.run("tleap -f grafted_chain.tleap", cwd=self.work_dir)

        # Convert AMBER topology → GROMACS
        parm = pmd.load_file(
            str(self.work_dir / "grafted_chain.prmtop"),
            xyz=str(self.work_dir / "grafted_chain.inpcrd"),
        )
        parm.save(str(self.work_dir / "grafted_chain.top"), overwrite=True)
        parm.save(str(self.work_dir / "grafted_chain.gro"), overwrite=True)

        # Centre the box
        runner.run(
            f"gmx editconf"
            f" -f {self.work_dir / 'grafted_chain.gro'}"
            f" -c"
            f" -o {self.work_dir / 'grafted_chain_center.gro'}",
        )

    # ------------------------------------------------------------------
    # Step 7 – Shift to bottom and add position restraints
    # ------------------------------------------------------------------

    def step_position_restraints(self) -> None:
        """Shift chains to z=0, identify linker atoms, write restraint topology."""
        from ase.io import read, write

        cfg = self.config
        print("\n=== Step 7: Position restraints and chain shifting ===")

        grafted_atoms = read(str(self.work_dir / "grafted_chain_center.gro"))

        # Resolve linker atoms (interactive if not specified in config)
        if cfg.linker_atoms is not None:
            linker_specs = cfg.linker_atoms
        else:
            default_resname = cfg.head.resname.upper()
            default_atomname = cfg.head.termname
            print(
                f"Default linker: resname={default_resname}, atomname={default_atomname}"
            )
            print(
                "Enter RESNAME1,ATOMNAME1,RESNAME2,ATOMNAME2,... "
                "or press Enter / type 'y' to accept default:"
            )
            ans = input().strip()
            if ans.lower() in ("", "y"):
                linker_specs = [
                    {"resname": default_resname, "atomname": default_atomname}
                ]
            else:
                parts = [p.strip() for p in ans.split(",")]
                linker_specs = [
                    {"resname": parts[i], "atomname": parts[i + 1]}
                    for i in range(0, len(parts), 2)
                ]

        # 0-based indices
        linker_0based = atom_utils.find_linker_indices(grafted_atoms, linker_specs)
        print(f"Linker atom indices (0-based): {linker_0based}")

        # Shift all atoms so the minimum z is at 0
        min_z = grafted_atoms.get_positions()[:, 2].min()
        # Pin linker atoms exactly to z=0
        for i in linker_0based:
            grafted_atoms.positions[i, 2] = min_z
        grafted_atoms.positions[:, 2] -= min_z

        shifted_gro = self.work_dir / "grafted_chain_shifted.gro"
        write(str(shifted_gro), grafted_atoms)

        # 1-based indices for GROMACS topology
        linker_1based = [i + 1 for i in linker_0based]
        print(f"Linker atom indices (1-based, for topology): {linker_1based}")

        topology.insert_restraint_top(
            self.work_dir / "grafted_chain.top",
            self.work_dir / "grafted_chain_restraint.top",
            linker_1based,
        )

    # ------------------------------------------------------------------
    # Step 8 – Vacuum GROMACS relaxation
    # ------------------------------------------------------------------

    def step_vacuum_relax(self) -> None:
        """Run vacuum energy minimisation + NVT in GROMACS."""
        cfg = self.config
        wd = self.work_dir
        print("\n=== Step 8: Vacuum GROMACS relaxation ===")

        # Copy MDP templates if needed
        for mdp in ("min_vac.mdp", "nvt_vac.mdp"):
            src = self.mdp_dir / mdp
            dst = wd / mdp
            if not dst.exists() and src.exists():
                shutil.copy(src, dst)

        hard_top = wd / "hardrest_grafted_chain_restraint.top"
        shifted_gro = wd / "grafted_chain_shifted.gro"

        runner.run(
            f"gmx grompp"
            f" -f min_vac.mdp"
            f" -p {hard_top}"
            f" -c {shifted_gro}"
            f" -o min_vac.tpr"
            f" -r {shifted_gro}"
            f" -maxwarn 2",
            cwd=wd,
        )
        runner.run(
            f"gmx mdrun -deffnm min_vac"
            f" -ntmpi {cfg.t_mpi} -ntomp {cfg.t_omp} -v",
            cwd=wd,
        )

        runner.run(
            f"gmx grompp"
            f" -f nvt_vac.mdp"
            f" -p {hard_top}"
            f" -c min_vac.gro"
            f" -o nvt_vac.tpr"
            f" -r {shifted_gro}"
            f" -maxwarn 2",
            cwd=wd,
        )
        runner.run(
            f"gmx mdrun -deffnm nvt_vac"
            f" -ntmpi {cfg.t_mpi} -ntomp {cfg.t_omp} -v",
            cwd=wd,
        )

    # ------------------------------------------------------------------
    # Step 9 – Solvation
    # ------------------------------------------------------------------

    def step_solvate(self) -> None:
        """Solvate the system with TIP3P water."""
        wd = self.work_dir
        cfg = self.config
        print("\n=== Step 9: Solvation with TIP3P water ===")

        hard_top = wd / "hardrest_grafted_chain_restraint.top"

        runner.run(
            f"gmx solvate"
            f" -cp {wd / 'nvt_vac.gro'}"
            f" -p {hard_top}"
            f" -o {wd / 'grafted_chain_water.gro'}",
        )

        from ase.io import read
        atoms = read(str(wd / "grafted_chain_water.gro"))
        Lc = (atoms.cell[2, 2] + 4) / 10  # +4 Å headroom, convert to nm
        La = atoms.cell[0, 0] / 10
        Lb = atoms.cell[1, 1] / 10

        runner.run(
            f"gmx editconf"
            f" -f {wd / 'grafted_chain_water.gro'}"
            f" -box {La:.6f} {Lb:.6f} {Lc:.6f}"
            f" -o {wd / 'grafted_chain_water_box.gro'}",
        )

        # Patch topologies with TIP3P include
        topology.insert_tip3p_top(
            wd / "grafted_chain_restraint.top",
            wd / "grafted_chain_water_restraint.top",
        )
        topology.insert_tip3p_top(
            wd / "hardrest_grafted_chain_restraint.top",
            wd / "hardrest_grafted_chain_water_restraint.top",
        )

        print("\n=== Done! ===")
        print("Necessary files for GROMACS simulation:")
        print(f"  - {wd / 'grafted_chain_water_box.gro'}")
        print(f"  - {wd / 'grafted_chain_water_restraint.top'}")
        print(f"  - {wd / 'hardrest_grafted_chain_water_restraint.top'}")
        print("  - tip3p.itp  (copy from 02_MD_template/)")
