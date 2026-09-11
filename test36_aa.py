#!/usr/bin/env python3
"""Explicit-water montmorillonite pilot MD -> input for test36.py.

Pinned ClayCode CD21 unit cells + its ClayFF/SPC/IOD ion parameters. No ML
model is used to generate the all-atom trajectory. See TEST36_AA.md.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import signal
import time

import ase.io
from ase import Atoms
import numpy as np
import openmm as mm
from openmm import app, unit
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parent
COMMIT = "fb52753bbad9c21369442589ac66b16a3bdb24e2"
SOURCE = ROOT / "test36-aa/sources" / f"ClayCode-{COMMIT}"
STOP = False


def stop(signum, frame):
    global STOP
    STOP = True


def dump_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fresh(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=True)
    if any(path.iterdir()):
        raise ValueError(f"Refusing to overwrite nonempty directory: {path}")
    return path


def sections(path):
    """Read simple data sections; never execute source-repository Python."""
    result, section = {}, None
    for raw in Path(path).read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("["):
            section = line.strip("[] ")
            result.setdefault(section, [])
        elif section:
            result[section].append(line.split())
    return result


def mic(delta, box):
    return delta - np.rint(delta / box) * box


def read_uc(data, name):
    gro_path = data / f"UCS/CD21/{name}.gro"
    itp_path = data / f"UCS/CD21/{name}.itp"
    lines = gro_path.read_text().splitlines()
    count = int(lines[1])
    xyz = np.array([[float(line[a:a+8]) for a in (20, 28, 36)]
                    for line in lines[2:count+2]])
    box = np.array([float(x) for x in lines[count+2].split()])
    top = sections(itp_path)
    if box.shape != (3,) or len(top["atoms"]) != count:
        raise ValueError("Only the pinned orthorhombic CD21 cells are supported")
    for i, (line, atom) in enumerate(zip(lines[2:count+2], top["atoms"])):
        if int(atom[0]) != i+1 or line[10:15].strip() != atom[4]:
            raise ValueError("Unit-cell coordinate/topology ordering mismatch")
    return xyz, box, top


def build(args):
    # Count checks precede creation of output files.
    if args.nx < 4 or args.ny < 3 or args.nx * args.ny % 8:
        raise ValueError("Require nx>=4, ny>=3, and nx*ny divisible by 8")
    if args.layers < 2 or args.waters_per_cell < 1 or args.spacing_nm < 1.4:
        raise ValueError("Require >=2 layers, water, and spacing >=1.4 nm")
    source = args.source.resolve()
    manifest = json.loads((ROOT / "test36-aa/source_manifest.json").read_text())
    if manifest["commit"] != COMMIT:
        raise ValueError("Source manifest/implementation version mismatch")
    for relative, expected in manifest["files"].items():
        if sha(source / relative) != expected:
            raise ValueError(f"Source differs from pinned version: {relative}")
    data = source / "package/ClayCode/data/data"
    cells = {n: read_uc(data, n) for n in ("C2000", "C2013", "C2014")}
    ff = {}
    ff_paths = [data / "FF/ClayFF_Fe.ff/ffnonbonded.itp",
                data / "FF/Ions.ff/ffnonbonded.itp"]
    for path in ff_paths:
        ff.update({r[0]: r for r in sections(path)["atomtypes"]})
    # Upstream's st row labels Si as Z=16. Correct ONLY the element label;
    # its mass, charge, sigma and epsilon are left unchanged.
    ff["st"] = list(ff["st"])
    ff["st"][1] = "14"
    base_box = cells["C2000"][1]
    box = np.array([base_box[0]*args.nx, base_box[1]*args.ny,
                    args.spacing_nm*args.layers])
    rng = np.random.default_rng(args.seed)
    positions, atom_rows, bonds, layer_ids, uc_names = [], [], [], [], []
    for layer in range(args.layers):
        for ix in range(args.nx):
            for iy in range(args.ny):
                name = "C2000" if (ix+iy+layer) % 2 else (
                    "C2013" if (ix+layer) % 2 else "C2014")
                xyz, _, top = cells[name]
                xyz = xyz.copy()
                xyz[:, 2] -= 0.4  # CD21 octahedral plane in supplied gro (nm)
                xyz += [ix*base_box[0], iy*base_box[1],
                        (layer+0.5)*args.spacing_nm]
                offset = len(atom_rows)
                positions.extend(xyz)
                for r in top["atoms"]:
                    atom_rows.append(dict(type=r[1], name=r[4],
                                          charge=float(r[6]), mass=float(r[7]),
                                          residue=name, residue_id=len(uc_names)+1,
                                          atomic_number=int(ff[r[1]][1])))
                for r in top["bonds"]:
                    if int(r[2]) != 1:
                        raise ValueError("Only harmonic clay OH bonds supported")
                    bonds.append([offset+int(r[0])-1, offset+int(r[1])-1,
                                  float(r[3]), float(r[4])])
                layer_ids.extend([layer]*len(xyz))
                uc_names.append(name)
    nclay = len(atom_rows)
    clay = np.array(positions)
    # Verify mixing unit cells has NOT put substituted-O charges next to the
    # wrong metal. Coordination must agree across tile and periodic boundaries.
    mg = np.array([i for i, a in enumerate(atom_rows) if a["type"] == "mgo"])
    oxy = np.array([i for i, a in enumerate(atom_rows) if a["atomic_number"] == 8])
    distances = np.linalg.norm(mic(clay[mg, None]-clay[None, oxy], box), axis=-1)
    affected = (distances < 0.25).any(axis=0)
    assigned = np.array([atom_rows[i]["type"] in ("obos", "ohs") for i in oxy])
    if not np.array_equal(affected, assigned) or not np.all((distances < .25).sum(1) == 6):
        raise ValueError("Substitution charge types do not match Mg coordination")
    charge_clay = sum(a["charge"] for a in atom_rows)
    if abs(charge_clay + len(mg)) > 1e-4:
        raise ValueError(f"Unexpected layer charge {charge_clay}")
    # Equal Na/Ca charge equivalents, not equal ion numbers. Each interlayer
    # gets enough countercharge to neutralize one sheet's charge.
    q_per_sheet = args.nx*args.ny//2
    n_na, n_ca = q_per_sheet//2, q_per_sheet//4
    occupied = list(clay[[a["atomic_number"] != 1 for a in atom_rows]])
    centers = []
    for gap in range(args.layers):
        zmid = (gap+1)*args.spacing_nm
        for ion_type in ["Na"]*n_na + ["Ca"]*n_ca:
            for attempt in range(100000):
                candidate = np.array([rng.uniform(0, box[0]), rng.uniform(0, box[1]),
                                      zmid + rng.uniform(-.05, .05)]) % box
                if np.linalg.norm(mic(np.array(occupied)-candidate, box), axis=1).min() > .30:
                    break
            else:
                raise ValueError("Cannot place counterions without overlap")
            positions.append(candidate)
            centers.append(candidate)
            occupied.append(candidate)
            p = ff[ion_type]
            atom_rows.append(dict(type=ion_type, name=ion_type,
                                  residue=ion_type, residue_id=len(uc_names)+len(centers),
                                  charge=float(p[3]), mass=float(p[2]),
                                  atomic_number=int(p[1])))
            layer_ids.append(-1)
    nions = len(centers)
    # Rigid SPC geometry from spc.itp SETTLE distances; random placement is ONLY
    # an initial condition, never an MD training trajectory by itself.
    spc = sections(data / "FF/ClayFF_Fe.ff/spc.itp")
    doh, dhh = map(float, spc["settles"][0][2:4])
    theta = 2*np.arcsin(dhh/(2*doh))
    water_shape = np.array([[0, 0, 0], [doh, 0, 0],
                            [doh*np.cos(theta), doh*np.sin(theta), 0]])
    water_centers = []
    water_count = args.nx*args.ny*args.waters_per_cell
    # Surface basal O is at |z-z_sheet|~0.33 nm. Keep O-water >=0.24 nm
    # from a surface O, but allow its initial z closer in surface cavities.
    half_gap = (args.spacing_nm-.66)/2-.12
    for gap in range(args.layers):
        for _ in range(water_count):
            for attempt in range(200000):
                candidate = np.array([rng.uniform(0, box[0]), rng.uniform(0, box[1]),
                                      (gap+1)*args.spacing_nm+rng.uniform(-half_gap, half_gap)]) % box
                if np.linalg.norm(mic(np.array(occupied)-candidate, box), axis=1).min() < .235:
                    continue
                if water_centers and np.linalg.norm(mic(np.array(water_centers)-candidate, box), axis=1).min() < .24:
                    continue
                water = candidate + Rotation.random(random_state=rng).apply(water_shape)
                # Avoid singular H--heavy-atom contacts during minimization.
                heavy = np.array(occupied + water_centers)
                if np.linalg.norm(mic(water[1:, None]-heavy[None], box), axis=-1).min() < .13:
                    continue
                break
            else:
                raise ValueError("Water packing failed; reduce hydration or increase spacing")
            water_centers.append(candidate)
            positions.extend(water)
            for r in spc["atoms"]:
                p = ff[r[1]]
                atom_rows.append(dict(type=r[1], name=r[4], residue="SOL",
                                      residue_id=len(uc_names)+nions+len(water_centers),
                                      charge=float(r[6]), mass=float(p[2]),
                                      atomic_number=int(p[1])))
                layer_ids.append(-1)
    charge = sum(a["charge"] for a in atom_rows)
    if abs(charge) > 1e-4:
        raise ValueError(f"Non-neutral system: {charge} e")
    out = fresh(args.output)
    (out / "ClayCode.LICENSE.txt").write_text((source / "LICENSE.txt").read_text())
    # Do not wrap H atoms separately: OpenMM's rigid-water constraints need a
    # whole molecule even when its O lies beside a periodic boundary. PME and
    # our periodic clay-OH bonds correctly accept coordinates outside the box.
    positions = np.array(positions)
    # A self-contained Gromacs topology is portable; no source path includes.
    used = sorted({a["type"] for a in atom_rows})
    top = ["; Generated by test36_aa.py; ClayCode sources in provenance.json",
           "; st atomic number corrected from 16 to 14; numerical FF unchanged",
           "[ defaults ]", "1 2 yes 1.0 1.0", "", "[ atomtypes ]"]
    top += [" ".join(ff[name]) for name in used]
    top += ["", "[ moleculetype ]", "CLAY 1", "", "[ atoms ]"]
    for i, a in enumerate(atom_rows[:nclay]):
        top.append(f"{i+1} {a['type']} {a['residue_id']} {a['residue']} {a['name']} {i+1} {a['charge']:.8f} {a['mass']:.8f}")
    top += ["", "[ bonds ]"]
    top += [f"{i+1} {j+1} 1 {d:.8f} {k:.8f}" for i, j, d, k in bonds]
    for name in ("Na", "Ca"):
        top += ["", "[ moleculetype ]", f"{name} 1", "[ atoms ]",
                f"1 {name} 1 {name} {name} 1 {ff[name][3]} {ff[name][2]}"]
    top += ["", "[ moleculetype ]", "SOL 2", "[ atoms ]"]
    top += [" ".join(r) for r in spc["atoms"]]
    top += ["[ settles ]", " ".join(spc["settles"][0]), "[ exclusions ]",
            "1 2 3", "2 1 3", "3 1 2", "", "[ system ]",
            "Idealized hydrated Na/Ca cis-montmorillonite pilot", "", "[ molecules ]", "CLAY 1"]
    for a in atom_rows[nclay:nclay+nions]:
        top.append(f"{a['type']} 1")
    top.append(f"SOL {len(water_centers)}")
    (out / "system.top").write_text("\n".join(top)+"\n")
    gro = ["test36 explicit-water montmorillonite initial condition", str(len(atom_rows))]
    for i, (a, xyz) in enumerate(zip(atom_rows, positions)):
        gro.append(f"{a['residue_id']%100000:5d}{a['residue']:<5s}{a['name']:>5s}{(i+1)%100000:5d}"
                   + "".join(f"{p:8.3f}" for p in xyz))
    gro.append(" ".join(f"{p:.8f}" for p in box))
    (out / "initial.gro").write_text("\n".join(gro)+"\n")
    np.save(out / "initial_positions_nm.npy", positions)
    atoms = Atoms(numbers=[a["atomic_number"] for a in atom_rows],
                  positions=positions*10, cell=box*10, pbc=True,
                  masses=[a["mass"] for a in atom_rows])
    atoms.new_array("atom_id", np.arange(1, len(atoms)+1))
    atoms.new_array("ff_type", np.array([a["type"] for a in atom_rows]))
    atoms.new_array("layer_id", np.array(layer_ids))
    atoms.set_initial_charges([a["charge"] for a in atom_rows])
    ase.io.write(out / "initial.extxyz", atoms)
    condition = dict(temperature_k=300., clay_type="idealized cis-montmorillonite",
                     layer_charge_e_per_Si8O20_OH4=-.5,
                     basal_spacing_angstrom=args.spacing_nm*10,
                     water_model="rigid SPC", waters_per_unit_cell=args.waters_per_cell,
                     sodium_count=nions*2//3, calcium_count=nions//3,
                     sodium_charge_fraction=.5, added_salt_molar=0., ensemble="NVT",
                     layers=args.layers, lateral_repeats=[args.nx, args.ny])
    # Minimal CG mapping: eliminate ONLY water, retain every clay atom and ion.
    # This is solvent coarse-graining, not a rigid-platelet model.
    retained = atom_rows[:nclay+nions]
    species = []
    for name in sorted({a["type"] for a in retained}):
        a = next(a for a in retained if a["type"] == name)
        species.append(dict(name=name, atomic_number=a["atomic_number"], mass_amu=a["mass"]))
    mapping = dict(source_is_equilibrium=False, condition=condition,
                   description="Remove SPC water only; all clay atoms including hydroxyl H and counterions retained. Pilot, not proven equilibrium.",
                   species=species, sites=[dict(species=a["type"], indices=[i]) for i, a in enumerate(retained)])
    dump_json(out / "mapping.json", mapping)
    used_files = [source / "LICENSE.txt", *ff_paths,
                  data / "FF/ClayFF_Fe.ff/forcefield.itp",
                  data / "FF/ClayFF_Fe.ff/spc.itp"]
    used_files += [data / f"UCS/CD21/{name}.{ext}" for name in cells for ext in ("gro", "itp")]
    provenance = dict(source_repo="https://github.com/Erastova-group/ClayCode",
                      source_commit=COMMIT, source_files=[dict(path=str(p.relative_to(source)), sha256=sha(p)) for p in used_files],
                      citations=["https://doi.org/10.1021/acs.jctc.4c00987",
                                 "https://doi.org/10.1021/acs.jpcc.1c04600",
                                 "https://doi.org/10.1021/acs.jctc.2c01255"],
                      forcefield="ClayCode ClayFF_Fe (no Fe used) + rigid SPC + distributed IOD Na/Ca parameters",
                      element_label_correction="st atomic number 16 -> 14 (Si); numerical FF parameters unchanged",
                      units="nm, ps, kJ/mol, amu, elementary charge; extxyz in angstrom",
                      seed=args.seed, unit_cell_counts=dict(Counter(uc_names)),
                      condition=condition, atoms=len(atoms), clay_atoms=nclay,
                      water_molecules=len(water_centers), ion_count=nions,
                      clay_charge_e=charge_clay, net_charge_e=charge,
                      mg_coordination_verified=True, box_nm=box.tolist(),
                      clay_bonds=bonds, atom_table=atom_rows,
                      source_is_equilibrium=False,
                      limitation="Infinite periodic layers, chosen fixed spacing/composition; no edges/tactoids, no equilibrium or water/ion parameter-transfer validation.")
    dump_json(out / "provenance.json", provenance)
    print(json.dumps({k: provenance[k] for k in ("atoms", "clay_atoms", "water_molecules", "net_charge_e", "box_nm")}), flush=True)


def make_system(input_dir, cutoff):
    gro = app.GromacsGroFile(str(input_dir / "initial.gro"))
    top = app.GromacsTopFile(str(input_dir / "system.top"),
                             periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir=str(input_dir))
    system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=cutoff*unit.nanometer,
                              constraints=None, rigidWater=True, ewaldErrorTolerance=1e-5)
    # Homogeneous dispersion correction is not appropriate for this layered,
    # inhomogeneous pilot. Plain truncated LJ + PME Coulomb; record explicitly.
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            force.setUseDispersionCorrection(False)
            force.setExceptionsUsePeriodicBoundaryConditions(True)
        if isinstance(force, mm.HarmonicBondForce):
            force.setUsesPeriodicBoundaryConditions(True)
    return top, system


def validate_system(input_dir, system, metadata):
    n = metadata["atoms"]
    if system.getNumParticles() != n or system.getNumConstraints() != metadata["water_molecules"]*3:
        raise ValueError("Incorrect particle/water constraint count")
    nonbond = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    for i, atom in enumerate(metadata["atom_table"]):
        q, sig, eps = nonbond.getParticleParameters(i)
        if abs(q.value_in_unit(unit.elementary_charge)-atom["charge"]) > 1e-9:
            raise ValueError("Charge changed during topology import")
        if abs(system.getParticleMass(i).value_in_unit(unit.dalton)-atom["mass"]) > 1e-9:
            raise ValueError("Mass changed during topology import")
    exceptions = {}
    for i in range(nonbond.getNumExceptions()):
        a, b, q, s, e = nonbond.getExceptionParameters(i)
        if abs(q.value_in_unit(unit.elementary_charge**2)) > 1e-15 or abs(e.value_in_unit(unit.kilojoule_per_mole)) > 1e-15:
            raise ValueError("Unexpected scaled nonbonded exception")
        exceptions[tuple(sorted((a, b)))] = True
    expected = {tuple(sorted((a, b))) for a, b, _, _ in metadata["clay_bonds"]}
    for i in range(metadata["clay_atoms"]+metadata["ion_count"], n, 3):
        expected.update({(i, i+1), (i, i+2), (i+1, i+2)})
    if set(exceptions) != expected:
        raise ValueError("Imported nonbonded exclusions differ from clay OH + SPC intramolecular pairs")


def run(args):
    input_dir = args.input.resolve()
    meta = json.loads((input_dir / "provenance.json").read_text())
    if min(meta["box_nm"]) <= 2*args.cutoff_nm:
        raise ValueError("Cutoff must be smaller than half the shortest box length")
    if args.dt_fs > 0.5:
        raise ValueError("This pilot leaves clay OH flexible: use dt <=0.5 fs")
    if args.production_ps < args.save_every*args.dt_fs/1000:
        raise ValueError("Production must be at least one output interval")
    top, system = make_system(input_dir, args.cutoff_nm)
    validate_system(input_dir, system, meta)
    integrator = mm.LangevinMiddleIntegrator(args.temperature*unit.kelvin,
                                             args.friction_per_ps/unit.picosecond,
                                             args.dt_fs*unit.femtosecond)
    integrator.setConstraintTolerance(1e-6)
    integrator.setRandomNumberSeed(args.seed)
    platform = mm.Platform.getPlatformByName(args.platform)
    properties = {}
    if args.platform == "CPU":
        properties = {"Threads": str(args.threads), "DeterministicForces": "true"}
    elif args.platform in ("CUDA", "OpenCL"):
        properties = {"Precision": "mixed"}
    sim = app.Simulation(top.topology, system, integrator, platform, properties)
    out = args.output.resolve()
    if args.resume:
        if not out.is_dir() or not (out / "latest.chk").is_file():
            raise ValueError("--resume requires an existing run directory with latest.chk")
    else:
        out = fresh(out)
    (out / "system.xml").write_text(mm.XmlSerializer.serialize(system))
    sim.context.setPositions(np.load(input_dir / "initial_positions_nm.npy")*unit.nanometer)
    if args.resume:
        sim.loadCheckpoint(str(out / "latest.chk"))
    start = time.monotonic()
    settings = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items() if k != "func"}
    settings.update(openmm_version=mm.__version__, platform_actual=sim.context.getPlatform().getName(),
                    input_topology_sha256=sha(input_dir / "system.top"),
                    input_positions_sha256=sha(input_dir / "initial_positions_nm.npy"),
                    periodic_clay_bonds=True, dispersion_correction=False, lj_switching=False,
                    source_is_equilibrium=False, constrained_water_only=True)
    dump_json(out / "settings.json", settings)
    mapping = json.loads((input_dir / "mapping.json").read_text())
    mapping["condition"]["temperature_k"] = args.temperature
    dump_json(out / "mapping.json", mapping)
    if not args.resume:
        print(f"Minimizing {meta['atoms']} atoms on {settings['platform_actual']}", flush=True)
        initial_energy = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer, maxIterations=args.minimize_iterations)
        min_state = sim.context.getState(getPositions=True, getEnergy=True)
        print(f"Energy: {initial_energy:.3f} -> {min_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.3f} kJ/mol", flush=True)
    template = ase.io.read(input_dir / "initial.extxyz")

    def frame(state, stage):
        atoms = template.copy()
        atoms.positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        atoms.info.update(time_ps=state.getTime().value_in_unit(unit.picosecond), stage=stage,
                          source_is_equilibrium=False)
        return atoms

    if not args.resume:
        ase.io.write(out / "minimized.extxyz", frame(min_state, "minimized"))
        sim.context.setVelocitiesToTemperature(50*unit.kelvin, args.seed+1)
        sim.context.applyVelocityConstraints(1e-6)
    dof = 3*meta["atoms"]-system.getNumConstraints()-3
    warm = round(args.warmup_ps*1000/args.dt_fs)
    equil = round(args.equilibration_ps*1000/args.dt_fs)
    production = round(args.production_ps*1000/args.dt_fs)
    total = warm+equil+production
    records, frames = [], 0
    previous_step, previous_time = 0, time.monotonic()
    with (out / "thermo.csv").open("a" if args.resume else "w") as stream:
        writer = csv.DictWriter(stream, fieldnames=["step", "time_ps", "stage", "target_k", "temperature_k",
                                                      "potential_kj_mol", "kinetic_kj_mol", "max_force_kj_mol_nm", "steps_per_s"])
        if not args.resume:
            writer.writeheader()
        while sim.currentStep < total:
            step = sim.currentStep
            if step < warm:
                stage, end = "warmup", warm
                target = 50+(args.temperature-50)*min(1., (step+args.save_every)/max(1, warm))
            elif step < warm+equil:
                stage, end, target = "equilibration", warm+equil, args.temperature
            else:
                stage, end, target = "production_pilot", total, args.temperature
            integrator.setTemperature(target*unit.kelvin)
            sim.step(min(args.save_every, end-step))
            state = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
            temperature = 2*ke/(dof*unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin))
            if not np.isfinite(forces).all() or not np.isfinite(pe+ke) or temperature > 2000:
                raise RuntimeError("Nonfinite/unstable trajectory; do not use as training data")
            now = time.monotonic()
            record = dict(step=sim.currentStep, time_ps=state.getTime().value_in_unit(unit.picosecond),
                          stage=stage, target_k=target, temperature_k=temperature,
                          potential_kj_mol=pe, kinetic_kj_mol=ke,
                          max_force_kj_mol_nm=float(np.linalg.norm(forces, axis=1).max()),
                          steps_per_s=(sim.currentStep-previous_step)/(now-previous_time))
            previous_step, previous_time = sim.currentStep, now
            writer.writerow(record)
            stream.flush()
            records.append(record)
            atoms = frame(state, stage)
            if stage == "production_pilot":
                ase.io.write(out / "production.extxyz", atoms, append=True)
                frames += 1
            ase.io.write(out / "latest.extxyz", atoms)
            # Checkpoint includes RNG state on this platform; XML state is a
            # portable starting point (not bitwise/RNG-identical continuation).
            checkpoint_tmp = out / "latest.chk.tmp"
            sim.saveCheckpoint(str(checkpoint_tmp))
            checkpoint_tmp.replace(out / "latest.chk")
            sim.saveState(str(out / "latest.state.xml"))
            progress = dict(record, total_steps=total, production_frames=frames,
                            wall_seconds=now-start, source_is_equilibrium=False)
            dump_json(out / "progress.json", progress)
            print(f"{sim.currentStep}/{total} {stage} {record['time_ps']:.3f} ps T={temperature:.1f} K {record['steps_per_s']:.1f} step/s", flush=True)
            if STOP or now-start > args.time_budget_hours*3600:
                dump_json(out / "run_metrics.json", dict(progress, complete=False))
                return 75
    ase.io.write(out / "final.extxyz", atoms)
    prod = [r for r in records if r["stage"] == "production_pilot"]
    metrics = dict(progress, complete=True, atoms=meta["atoms"],
                   production_ps=args.production_ps,
                   production_temperature_mean_k=float(np.mean([r["temperature_k"] for r in prod])),
                   production_temperature_std_k=float(np.std([r["temperature_k"] for r in prod])),
                   warning="A short pilot at fixed spacing. Thermal stability is not proof of equilibrium or force-field accuracy.")
    dump_json(out / "run_metrics.json", metrics)
    print(json.dumps(metrics, indent=2), flush=True)
    return 0


def positive(text):
    value = float(text)
    if not np.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("must be positive and finite")
    return value


def analyze(args):
    input_dir, run_dir = args.input.resolve(), args.run.resolve()
    meta = json.loads((input_dir / "provenance.json").read_text())
    with (run_dir / "thermo.csv").open() as stream:
        rows = [r for r in csv.DictReader(stream) if r["stage"] == "production_pilot"]
    if not rows:
        raise ValueError("No production-pilot records")
    frames = list(ase.io.iread(run_dir / "production.extxyz"))
    if len(frames) != len(rows):
        raise ValueError("Trajectory/thermo frame counts disagree")
    nclay = meta["clay_atoms"]
    atom_types = np.array([a["type"] for a in meta["atom_table"]])
    clay_o = np.array([i for i, a in enumerate(meta["atom_table"][:nclay]) if a["atomic_number"] == 8])
    box = np.array(meta["box_nm"])*10

    def coordination(frame):
        result = {}
        for name, cutoff, expected in (("st", 2.2, 4), ("ao", 2.6, 6), ("mgo", 2.7, 6)):
            ids = np.flatnonzero(atom_types == name)
            dist = np.linalg.norm(mic(frame.positions[ids, None]-frame.positions[None, clay_o], box), axis=-1)
            counts = (dist < cutoff).sum(axis=1)
            result[name] = dict(atoms=len(ids), cutoff_angstrom=cutoff,
                                expected_coordination_fraction=float(np.mean(counts == expected)),
                                mean_coordination=float(counts.mean()),
                                nearest_shell_mean_angstrom=float(np.sort(dist, axis=1)[:, :expected].mean()))
        return result

    blocks = []
    for ids in np.array_split(np.arange(len(rows)), min(4, len(rows))):
        block = [rows[i] for i in ids]
        blocks.append(dict(first_ps=float(block[0]["time_ps"]), last_ps=float(block[-1]["time_ps"]),
                           mean_temperature_k=float(np.mean([float(r["temperature_k"]) for r in block])),
                           potential_kj_mol_per_atom=float(np.mean([float(r["potential_kj_mol"]) for r in block])/meta["atoms"])))
    water_first = nclay+meta["ion_count"]
    water = frames[-1].positions[water_first:].reshape(-1, 3, 3)
    water_oh = np.linalg.norm(mic(water[:, 1:]-water[:, :1], box), axis=-1)
    water_hh = np.linalg.norm(mic(water[:, 1]-water[:, 2], box), axis=-1)
    bonds = np.array(meta["clay_bonds"])
    clay_oh = np.linalg.norm(mic(frames[-1].positions[bonds[:, 0].astype(int)]-
                                   frames[-1].positions[bonds[:, 1].astype(int)], box), axis=-1)
    result = dict(frames=len(frames), frame_atoms=len(frames[0]),
                  first_ps=float(rows[0]["time_ps"]), last_ps=float(rows[-1]["time_ps"]),
                  production_blocks=blocks,
                  coordination_initial=coordination(ase.io.read(input_dir / "initial.extxyz")),
                  coordination_first_production=coordination(frames[0]),
                  coordination_last_production=coordination(frames[-1]),
                  last_water_oh_max_error_angstrom=float(np.max(np.abs(water_oh-1.))),
                  last_water_hh_max_error_angstrom=float(np.max(np.abs(water_hh-1.633))),
                  last_clay_oh_mean_angstrom=float(clay_oh.mean()),
                  last_clay_oh_range_angstrom=[float(clay_oh.min()), float(clay_oh.max())],
                  identity_order_preserved=all(np.array_equal(f.numbers, frames[0].numbers) and
                                               np.array_equal(f.arrays["atom_id"], frames[0].arrays["atom_id"]) for f in frames),
                  finite_coordinates=all(np.isfinite(f.positions).all() for f in frames),
                  source_is_equilibrium=False,
                  limitation="This checks numerical stability/coordination only, not equilibrium, independence, swelling pressure, or force-field accuracy.")
    dump_json(run_dir / "analysis.json", result)
    print(json.dumps(result, indent=2), flush=True)


def count(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be >=1")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("build")
    p.add_argument("--source", type=Path, default=SOURCE)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--nx", type=count, default=4)
    p.add_argument("--ny", type=count, default=4)
    p.add_argument("--layers", type=count, default=2)
    p.add_argument("--waters-per-cell", type=count, default=12)
    p.add_argument("--spacing-nm", type=positive, default=1.55)
    p.add_argument("--seed", type=count, default=3601)
    p.set_defaults(func=build)
    p = sub.add_parser("run")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--platform", choices=["CPU", "OpenCL", "CUDA", "Reference"], default="CPU")
    p.add_argument("--threads", type=count, default=4)
    p.add_argument("--dt-fs", type=positive, default=.5)
    p.add_argument("--cutoff-nm", type=positive, default=1.0)
    p.add_argument("--temperature", type=positive, default=300.)
    p.add_argument("--friction-per-ps", type=positive, default=1.)
    p.add_argument("--warmup-ps", type=positive, default=2.)
    p.add_argument("--equilibration-ps", type=positive, default=5.)
    p.add_argument("--production-ps", type=positive, default=20.)
    p.add_argument("--save-every", type=count, default=200)
    p.add_argument("--minimize-iterations", type=count, default=2000)
    p.add_argument("--time-budget-hours", type=positive, default=1.)
    p.add_argument("--seed", type=count, default=3602)
    p.add_argument("--resume", action="store_true")
    p.set_defaults(func=run)
    p = sub.add_parser("analyze")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--run", type=Path, required=True)
    p.set_defaults(func=analyze)
    args = parser.parse_args()
    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, stop)
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
