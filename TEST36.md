# test36: clay CG using the test32 denoiser

`test36.py` preserves test32's actual model/method:

- graphite `NequIP`, the same irreps, 3 convolutions, radial network and two species embeddings;
- graphite `RattleParticles(sigma_max=0.75)` and `DownselectEdges(cutoff=5)`;
- unweighted MSE against applied displacement `data.dx`, AdamW (2e-4);
- **no sigma/time/thermodynamic condition input**;
- generation with `pos -= model(current_graph) + sigma * randn`, linearly
  annealed from 1 to 0.001 A over 2900 steps, then 100 polish steps;
- generation graphs use the large 10 A cutoff, just as the actual test32 code
  does (training downselects to 5 A). The resulting mismatch is preserved and
  should be considered when interpreting generated structures.

It does not import/execute the top-level training code in `test32.py` or alter
its SiO2 checkpoint. Dataset handling, temporal splitting, CLI configuration,
checkpointing and species mapping are new. Training samples real frames with
replacement; it does not materialize 128 identical copies. Validation uses
up to four held-out frames. Correlated trajectory frames still require a
proper independent-trajectory/condition validation for scientific claims.

## What this implements from the figure

External all-atom trajectory -> user-defined CG mapping -> CG training ->
generation / experimental force evaluation -> fixed-cell Langevin MD and
partial RDF comparisons. The coordinates and effective masses of clay sites,
Na and Ca are explicit; water atoms can be omitted by the mapping.

The figure's scalar PMF, pressure, rigid platelets, condition embeddings,
explicit/screened Coulomb terms, pore hierarchy, tactoid detection and transport
coefficients are **not implemented**. These cannot be claimed from the unchanged
test32 denoiser. Species atomic numbers are labels for geometry/file export;
they do not introduce charge interactions. Clay platelets are flexible collections
of unconstrained sites in this version. There is no automatic independence or
equilibrium guarantee for generated samples.

Use one independently trained model per temperature/composition/salinity/clay
condition. `metadata.condition` records the condition; it is **not an input to
the neural network**. Do not mix multiple conditions into one training trajectory.

## Input and mapping

An explicit-water clay pilot is now available through `test36_aa.py`; see
[TEST36_AA.md](TEST36_AA.md) for its sourced structure/force field, local
trajectory and limitations. It is short MD, not verified equilibrium evidence.
Alternatively supply a periodic trajectory (e.g. extxyz with a cell, or a LAMMPS dump). Frame ordering,
atom IDs and species must stay fixed. ASE sorts LAMMPS dumps by ID; mapping
indices refer to the resulting **zero-based row positions**, not raw LAMMPS IDs.
All lengths in prepared files are angstrom, masses amu, temperatures K.

Copy and customize `test36_mapping.example.json`. It is a schema example,
not a validated smectite mapping. Define all CG species, their representative
atomic numbers and **effective masses**, and all retained sites:

```json
{"species": "clay_Al", "indices": [0, 1], "weights": [0.7, 0.3]}
```

Groups must be disjoint; omitted atoms (e.g. water) disappear. Weights must
sum to one; absent weights give the arithmetic centroid. Choose mass-weighted
weights explicitly for a center of mass. Groups are unwrapped by minimum-image
displacements about their first atom, so every group must fit within that
minimum-image neighborhood. This is not a whole-platelet unwrapping algorithm.
Do not delete all O atoms: surface O and water O must be distinguished by mapping.
Set `source_is_equilibrium` only if the source trajectory actually represents
the intended equilibrium ensemble. That flag is provenance, not a validation.

```bash
python test36.py prepare --trajectory clay.extxyz --mapping clay_mapping.json \
  --output test36-data
python test36.py train --dataset test36-data --output test36-train --device cuda
python test36.py generate --checkpoint test36-train/checkpoint.pt \
  --output test36-generated --device cuda
python test36.py md --checkpoint test36-train/checkpoint.pt \
  --output test36-md --device cuda --steps 10000 --sigma-ref 0.1
python test36.py analyze --reference test36-data --samples test36-md/md.extxyz \
  --output test36-analysis --r-max 5
```

For local verification use the existing ScoreMD Pixi Python
`../.pixi/envs/default/bin/python` and `--device cpu` or `mps`. The MD time step
defaults to 0.1 fs and damping to 0.1 ps. The thermostat uses the recorded
training condition's temperature. `sigma-ref` must be chosen explicitly:
`F = -kBT * predicted_dx / sigma_ref^2` is an approximation, not an exact
conversion for a sigma-agnostic model. MD begins at the first mapped reference
frame. The generator also begins there, exactly as test32 begins at its reference.

An optional LAMMPS export uses the existing driver (orthorhombic fixed box only):

```bash
python test36.py export --checkpoint test36-train/checkpoint.pt \
  --output test36-forcefield --sigma-ref 0.1
python test32-CGMD-repo/test32_lammps.py \
  --bundle test36-forcefield/cg_score_forcefield.pt --output-dir test36-lammps
```

The export uses mapped reference coordinates and actual CG masses. No dummy
potential energy is returned. NPT, minimization and shear-stress measurements
are unavailable. Matching structure does not calibrate CG physical times.

## Restart and output

Training saves model, optimizer, RNG states and dataset fingerprints to
`checkpoint.pt`. Re-run with `--resume` in the same output directory; the
device and training settings must match. The total `--updates` may be extended.
Generation also supports `--resume`, with its checkpoint and memory-mapped
`positions.npy`; only `generation.json:valid_frames` have been generated.
Both stages stop at step boundaries on SIGINT/SIGTERM or the time budget,
save their state, and exit **75** (incomplete), not zero.

MD saves `md.extxyz`, `thermo.csv`, `progress.json`, `latest.extxyz`, and on
completion `final.extxyz` and `run_metrics.json`. MD is not exactly restartable
in this version; its Langevin RNG/integrator state is not serialized. A time
budget/SIGTERM also yields exit 75 and leaves the last saved structure. Set
`--save-every` small enough to finish a chunk before a job walltime expires.

Analysis writes normalized partial RDFs in `rdf.json`; compare samples after
an appropriate discarded transient (the current CLI averages all selected
frames). It does not certify equilibrium or compute d001, tactoid sizes,
swelling pressure, shear stress or transport coefficients.

## Supercomputer

Copy `test36.py`, `run_test36.pbs`, `Singularity.test36.def`, and a DM2 checkout
containing `src/graphite`. A prepared dataset is independent of the local path.
The optional LAMMPS route additionally needs `test32-CGMD-repo` and its LAMMPS
dependencies; the built-in ASE MD route does not.

```bash
singularity build --fakeroot test36.sif Singularity.test36.def
qsub -P <ProjectGroup_ID> \
  -v STAGE=train,DATASET_PATH=/path/to/test36-data,OUTPUT_DIR=/path/to/test36-train \
  run_test36.pbs
# Resume an unfinished training job:
qsub -P <ProjectGroup_ID> \
  -v STAGE=train,RESUME=1,DATASET_PATH=/path/to/test36-data,OUTPUT_DIR=/path/to/test36-train \
  run_test36.pbs
qsub -P <ProjectGroup_ID> \
  -v STAGE=md,CHECKPOINT_PATH=/path/to/test36-train/checkpoint.pt,OUTPUT_DIR=/path/to/test36-md,SIGMA_REF=0.1,STEPS=10000 \
  run_test36.pbs
```

PBS defaults follow the existing repository (sg8, 1 GPU, 1 MPI rank, 12 h).
Adapt queue/module requirements to the site. `CONTAINER_RUNTIME=apptainer`
is supported. `DM2_ROOT`, `SIF_IMAGE`, `BATCH_SIZE`, `UPDATES` and
`TIME_BUDGET_HOURS` can be overridden. The script binds dataset/model read-only
and output read/write. Training is single GPU; submitting multiple MPI ranks
would launch duplicate work, not distributed training. The container build
and GPU execution must be verified on the target cluster.

## Local software verification

`../.pixi/envs/default/bin/python -m unittest test36_checks -v` uses explicitly
synthetic coordinates and checks the original test32 architecture/forward
output, periodic CG mapping, input rejection, train/resume equivalence,
test32's generation update and resumed noise sequence, a two-step MD,
partial RDF analysis, and LAMMPS bundle export. These checks do not validate
a clay force field; the target-cluster CUDA/container run is not yet tested.
