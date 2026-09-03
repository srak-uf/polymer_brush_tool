"""Shared workflow logic for polymer brush system preparation.

The pipeline is split into ``step_*`` methods so that a run can be resumed
or inspected between stages.  Every external tool is executed with
``cwd=self.work_dir`` and every file name written into a tool input is
relative to that directory; :attr:`work_dir` is therefore resolved to an
absolute path once, in ``__init__``.
"""

from __future__ import annotations

import shutil
from importlib import resources
from pathlib import Path

from polymer_brush_tool import runner
from polymer_brush_tool.config import BrushConfig
from polymer_brush_tool.ff import fragments, minimize, topology
from polymer_brush_tool.structure import atoms as atom_utils
from polymer_brush_tool.structure import graft, solvent

_TEMPLATE_FILES = ("min_vac.mdp", "nvt_vac.mdp", "tip3p.itp")


def package_template_dir() -> Path:
    """Directory holding the mdp/itp templates shipped inside the package."""
    return Path(resources.files("polymer_brush_tool") / "templates")


class BrushWorkflowBase:
    """Base class shared by :class:`LinearBrushWorkflow` and :class:`LoopBrushWorkflow`.

    Parameters
    ----------
    config:
        Complete configuration object; ``config.validate()`` is called.
    work_dir:
        Directory where all intermediate and output files are written.
        The antechamber ``.ac`` files named in the config must be there.
    mdp_dir:
        Directory containing ``min_vac.mdp``, ``nvt_vac.mdp`` and
        ``tip3p.itp``.  Defaults to the templates bundled with the package.
    """

    def __init__(
        self,
        config: BrushConfig,
        work_dir: str | Path = ".",
        mdp_dir: str | Path | None = None,
    ) -> None:
        config.validate()
        self.config = config
        self.work_dir = Path(work_dir).resolve()
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.mdp_dir = Path(mdp_dir).resolve() if mdp_dir else package_template_dir()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _run(self, cmd: str, **kw) -> None:
        runner.run(cmd, cwd=self.work_dir, **kw)

    def _read(self, name: str):
        from ase.io import read
        return read(str(self.work_dir / name))

    def _write(self, name: str, atoms) -> None:
        from ase.io import write
        write(str(self.work_dir / name), atoms)

    def _copy_template(self, name: str) -> None:
        dst = self.work_dir / name
        if dst.exists():
            return
        src = self.mdp_dir / name
        if not src.exists():
            src = package_template_dir() / name
        shutil.copy(src, dst)

    # ------------------------------------------------------------------
    # Step 1 – .chain files and prepgen
    # ------------------------------------------------------------------

    def step_prepgen(self) -> None:
        """Write ``.chain`` connectivity files and run ``prepgen`` for HEAD, MID, TAIL."""
        cfg = self.config
        print("\n=== Step 1: prepgen (HEAD / MID / TAIL .prepi) ===")

        for spec in (cfg.head, cfg.mid, cfg.tail):
            ac = self.work_dir / spec.ac_file
            if not ac.exists():
                raise FileNotFoundError(
                    f"{spec.resname}: AC file not found: {ac}. "
                    f"Run antechamber first and place the .ac file in the work directory."
                )

        fragments.write_fragment(
            cfg.head.resname,
            tailname=cfg.head.tailname,
            omitnames=cfg.head.omitnames,
            post_tailtype=cfg.head.post_tailtype,
            work_dir=self.work_dir,
        )
        fragments.write_fragment(
            cfg.mid.resname,
            headname=cfg.mid.headname,
            tailname=cfg.mid.tailname,
            omitnames=cfg.mid.omitnames,
            pre_headtype=cfg.mid.pre_headtype,
            post_tailtype=cfg.mid.post_tailtype,
            work_dir=self.work_dir,
        )
        fragments.write_fragment(
            cfg.tail.resname,
            headname=cfg.tail.headname,
            omitnames=cfg.tail.omitnames,
            pre_headtype=cfg.tail.pre_headtype,
            work_dir=self.work_dir,
        )

        for spec in (cfg.head, cfg.mid, cfg.tail):
            self._run(
                f"prepgen -i {spec.ac_file} -o {spec.resname}.prepi -f prepi"
                f" -m {spec.resname}.chain -rn {spec.resname.upper()}"
                f" -rf {spec.resname}.res"
            )

    # ------------------------------------------------------------------
    # Step 2 – tleap chain build
    # ------------------------------------------------------------------

    def step_build_chain(self) -> None:
        """Polymerise HEAD + MID×n + TAIL with tleap → chain.prmtop / chain.inpcrd / chain.pdb."""
        cfg = self.config
        print("\n=== Step 2: tleap chain build ===")
        fragments.write_tleap(
            cfg.head.resname, cfg.mid.resname, cfg.tail.resname,
            cfg.n_mid_repeat_units, work_dir=self.work_dir,
        )
        self._run("tleap -f build_chain.tleap")
        # ASE round-trip normalises the PDB written by tleap.
        self._write("chain.pdb", self._read("chain.pdb"))

    # ------------------------------------------------------------------
    # Step 3 – sander minimisation
    # ------------------------------------------------------------------

    def _terminal_indices_1based(self) -> tuple[int, int]:
        cfg = self.config
        atoms = self._read("chain.pdb")
        head = atom_utils.find_atom_index(atoms, cfg.head.termname, cfg.head.resname.upper())
        tail = atom_utils.find_atom_index(atoms, cfg.tail.termname, cfg.tail.resname.upper())
        return head + 1, tail + 1

    def step_amber_minimize(self) -> None:
        """Stretch the chain with restrained sander minimisation → chain_min_pull.pdb."""
        cfg = self.config
        print(f"\n=== Step 3: sander minimisation with restraints ({cfg.topology}) ===")
        head_idx, tail_idx = self._terminal_indices_1based()
        minimize.amber_min_with_pull(
            self.work_dir / "chain.prmtop",
            self.work_dir / "chain.inpcrd",
            file_prefix="chain",
            head_idx=head_idx,
            tail_idx=tail_idx,
            polymer_length=cfg.polymer_length(),
            loop_height=cfg.loop_height() if cfg.topology == "loop" else None,
            work_dir=self.work_dir,
        )

    # ------------------------------------------------------------------
    # Step 4 – z alignment
    # ------------------------------------------------------------------

    def step_align_z(self) -> None:
        """Orient the stretched chain along z with its bottom atom at the substrate."""
        cfg = self.config
        print("\n=== Step 4: z-axis alignment (PACKMOL) ===")
        bottom_idx = cfg.bottom_atom_index
        if bottom_idx is None:
            print("Enter the 1-based index of the bottom (grafting) atom.")
            print(f"  (inspect {self.work_dir / 'chain_min_pull.pdb'} in VESTA)")
            bottom_idx = int(input("> "))
        graft.align_chain_z(
            self.work_dir / "chain_min_pull.pdb",
            bottom_idx,
            xyratio=cfg.xyratio,
            work_dir=self.work_dir,
        )

    # ------------------------------------------------------------------
    # Step 5 – graft
    # ------------------------------------------------------------------

    def step_graft(self) -> None:
        """Place nx × ny chains on the grafting grid (PACKMOL) → grafted_chain.pdb."""
        cfg = self.config
        print("\n=== Step 5: graft chains on grid (PACKMOL) ===")
        graft.graft_brush(
            self.work_dir / "aligned_chain.pdb",
            cfg.box_x(), cfg.box_y(), cfg.nx, cfg.ny,
            work_dir=self.work_dir,
        )

    # ------------------------------------------------------------------
    # Step 6 – force field for grafted system, convert to GROMACS
    # ------------------------------------------------------------------

    def step_assign_ff_grafted(self) -> None:
        """tleap on the multi-chain PDB, ParmEd → .top/.gro, then centre with editconf."""
        import parmed as pmd

        cfg = self.config
        print("\n=== Step 6: force field for grafted system, AMBER → GROMACS ===")
        cell = self._read("grafted_chain.pdb").cell
        fragments.write_tleap_grafted(
            cfg.head.resname, cfg.mid.resname, cfg.tail.resname,
            cell[0, 0], cell[1, 1], cell[2, 2],
            work_dir=self.work_dir,
        )
        self._run("tleap -f grafted_chain.tleap")

        parm = pmd.load_file(
            str(self.work_dir / "grafted_chain.prmtop"),
            xyz=str(self.work_dir / "grafted_chain.inpcrd"),
        )
        parm.save(str(self.work_dir / "grafted_chain.top"), overwrite=True)
        parm.save(str(self.work_dir / "grafted_chain.gro"), overwrite=True)

        self._run("gmx editconf -f grafted_chain.gro -c -o grafted_chain_center.gro")

    # ------------------------------------------------------------------
    # Step 7 – linker atoms, shift, position restraints
    # ------------------------------------------------------------------

    def _resolve_linker_specs(self) -> list[dict]:
        cfg = self.config
        default = [{"resname": cfg.head.resname.upper(), "atomname": cfg.head.termname}]
        if cfg.linker_atoms is not None:
            specs = cfg.linker_atoms
        else:
            print(f"Default linker: {default[0]['resname']} {default[0]['atomname']}")
            print("Enter RESNAME1,ATOMNAME1,RESNAME2,ATOMNAME2,... or press Enter / 'y' for default")
            ans = input("> ").strip()
            if ans.lower() in ("", "y"):
                specs = default
            else:
                parts = [p.strip() for p in ans.split(",")]
                if len(parts) % 2:
                    raise ValueError("linker input must be RESNAME,ATOMNAME pairs")
                specs = [{"resname": parts[i], "atomname": parts[i + 1]} for i in range(0, len(parts), 2)]
        return [{"resname": s["resname"].upper(), "atomname": s["atomname"]} for s in specs]

    def step_position_restraints(self) -> None:
        """Pin linker atoms to z = linker_height, shift the system, write soft/hard restraint topologies."""
        cfg = self.config
        print("\n=== Step 7: linker atoms, z-shift, position restraints ===")
        atoms = self._read("grafted_chain_center.gro")
        specs = self._resolve_linker_specs()

        linker_0 = atom_utils.find_linker_indices(atoms, specs)
        if not linker_0:
            raise ValueError(f"No atoms matched linker specification {specs}")

        # Linkers are covalently bonded to the substrate: place them one bond
        # length above the wall (z = 0) instead of on it.
        atom_utils.pin_linkers_to_substrate(atoms, linker_0, height=cfg.linker_height)
        self._write("grafted_chain_shifted.gro", atoms)

        linker_1 = [i + 1 for i in linker_0]
        print(f"Linker atom indices (1-based): {linker_1}, placed at z = {cfg.linker_height} Å")
        topology.insert_restraint_top(
            self.work_dir / "grafted_chain.top",
            self.work_dir / "grafted_chain_restraint.top",
            linker_1,
        )

    # ------------------------------------------------------------------
    # Step 8 – vacuum relaxation
    # ------------------------------------------------------------------

    def step_vacuum_relax(self) -> None:
        """GROMACS steepest-descent + short NVT in vacuum with hard restraints."""
        cfg = self.config
        print("\n=== Step 8: vacuum relaxation (gmx) ===")
        for name in ("min_vac.mdp", "nvt_vac.mdp"):
            self._copy_template(name)

        top = "hardrest_grafted_chain_restraint.top"
        ref = "grafted_chain_shifted.gro"
        mdrun = f"gmx mdrun -ntmpi {cfg.t_mpi} -ntomp {cfg.t_omp} -v"

        self._run(f"gmx grompp -f min_vac.mdp -p {top} -c {ref} -r {ref} -o min_vac.tpr -maxwarn 2")
        self._run(f"{mdrun} -deffnm min_vac")
        self._run(f"gmx grompp -f nvt_vac.mdp -p {top} -c min_vac.gro -r {ref} -o nvt_vac.tpr -maxwarn 2")
        self._run(f"{mdrun} -deffnm nvt_vac")

    # ------------------------------------------------------------------
    # Step 9 – solvation
    # ------------------------------------------------------------------

    def step_solvate(self) -> None:
        """Solvate with TIP3P, extend the box in z, finalise both topologies."""
        cfg = self.config
        print("\n=== Step 9: solvation (gmx solvate) ===")
        hard = "hardrest_grafted_chain_restraint.top"
        soft = "grafted_chain_restraint.top"

        # gmx solvate appends "SOL N" to the topology passed with -p; drop any
        # SOL line left by an earlier run so that re-running this step is safe.
        topology.remove_molecule(self.work_dir / hard, "SOL")
        self._run(f"gmx solvate -cp nvt_vac.gro -p {hard} -o grafted_chain_water_raw.gro")
        # gmx solvate ignores the wall at z = 0 and also fills the slab below
        # the grafting atoms; delete that water and fix the SOL count.
        kept, removed = solvent.remove_water_below(
            self.work_dir / "grafted_chain_water_raw.gro",
            self.work_dir / "grafted_chain_water.gro",
            cfg.solvent_min_z,
        )
        print(f"Removed {removed} water molecules with O below z = {cfg.solvent_min_z} Å ({kept} kept)")
        topology.set_molecule_count(self.work_dir / hard, "SOL", kept)
        # The soft topology must carry the same [ molecules ] block.
        topology.copy_molecules_section(self.work_dir / hard, self.work_dir / soft)

        cell = self._read("grafted_chain_water.gro").cell
        La, Lb = cell[0, 0] / 10, cell[1, 1] / 10
        Lc = (cell[2, 2] + 4) / 10                 # +4 Å vacuum gap above the water
        # -noc: editconf re-centres the system when -box is given, which would
        # lift the linkers off z = linker_height.  Keep the coordinates as they are.
        self._run(
            f"gmx editconf -f grafted_chain_water.gro -box {La:.6f} {Lb:.6f} {Lc:.6f}"
            f" -noc -o grafted_chain_water_box.gro"
        )

        topology.insert_tip3p_top(self.work_dir / soft, self.work_dir / "grafted_chain_water_restraint.top")
        topology.insert_tip3p_top(self.work_dir / hard, self.work_dir / "hardrest_grafted_chain_water_restraint.top")
        self._copy_template("tip3p.itp")

        print("\n=== Done. Files for the MD stage (copy together with 02_MD_template/*.mdp) ===")
        for name in (
            "grafted_chain_water_box.gro",
            "grafted_chain_water_restraint.top",
            "hardrest_grafted_chain_water_restraint.top",
            "tip3p.itp",
        ):
            print(f"  - {self.work_dir / name}")

    # ------------------------------------------------------------------
    # full pipeline
    # ------------------------------------------------------------------

    STEPS = (
        "step_prepgen",
        "step_build_chain",
        "step_amber_minimize",
        "step_align_z",
        "step_graft",
        "step_assign_ff_grafted",
        "step_position_restraints",
        "step_vacuum_relax",
        "step_solvate",
    )

    def run(self, start: str | None = None) -> None:
        """Run the pipeline, optionally starting from the named step.

        Parameters
        ----------
        start:
            Name of a ``step_*`` method to resume from, e.g. ``"step_graft"``.
            Earlier steps are assumed to have left their outputs in *work_dir*.
        """
        steps = list(self.STEPS)
        if start is not None:
            if start not in steps:
                raise ValueError(f"unknown step {start!r}; choose from {steps}")
            steps = steps[steps.index(start):]
        for name in steps:
            getattr(self, name)()
