# AB′ bilayer CrI₃ atomistic LLGS model

This workspace implements the approved class-4, material-specific atomistic model described in [plan.md](plan.md): BiDB AB′ bilayer CrI₃ geometry combined with a fitted scalar reduction of the Soenen exchange model. One fixed-length classical unit spin is attached to each Cr site. Iodine coordinates remain in the immutable source record but are not dynamical variables.

The baseline includes scalar exchange, out-of-plane single-ion anisotropy, Zeeman flux density, and signed layer-resolved damping-like and field-like SOT. It excludes DMI, Kitaev and other tensor exchange, dipolar fields, temperature, defects, electronic bands, and topological-magnon claims.

## Run a simulation

Run the documented candidate open-flake example directly from MATLAB's Run
button or the Command Window:

```matlab
result = run_open_sot_example();
```

The example creates the ignored `results/` directory when needed and saves
`results/candidate.mat`. It refuses to overwrite an existing result.

For a custom simulation, start from the documented configuration and explicitly
choose damping and final time:

```matlab
config = cri3_default_config();
config.nx = 8;
config.ny = 8;
config.boundary = 'open';
config.alpha = 0.001;       % explicit choice; see qualification below
config.t_final_s = 2e-12;
config.torque.dl_T = [0.01, -0.01]; % tesla: [bottom, top]
config.torque.fl_T = [0, 0];
config.output_file = fullfile('results', 'candidate.mat');
result = cri3_run(config);
```

For custom simulations, create the output directory first if saving is requested. Existing files are not overwritten unless `config.overwrite=true`. Leaving `output_file` empty returns the result without writing a file.

`config.alpha` intentionally has no production default. Soenen *et al.* used $\alpha=0.001$ for their separate $T=0$ stability calculation, but its transfer to a particular SOT device is a modeling choice. The code therefore requires the caller to select it explicitly. Likewise, torque amplitudes are direct field-equivalent inputs in tesla; this repository does not claim a material-specific switching threshold without an approved device and drive parameterization.

For a time-dependent pulse, set `config.torque_function` to a function handle that accepts time in seconds and returns a struct with `sigma`, `dl_T`, and `fl_T` fields.

## Model and sign conventions

The unique-bond Hamiltonian is

$$
\mathcal H=\sum_{\{i,j\}}J_{ij}\mathbf s_i\cdot\mathbf s_j
-\sum_i K_i s_{i,z}^2
-\sum_i\mu_i\mathbf s_i\cdot\mathbf B_{\rm ext}.
$$

Negative $J$ favors parallel spins, positive $J$ favors antiparallel spins, and every physical bond is stored once. Energies are in joules, moments in J/T, and effective fields and torque amplitudes in tesla. `cri3_llgs_rhs.m` implements the explicit fixed-length LLGS equation `eq:ll-with-flt` from the registered `spintronic_methods.tex`, including both damping-dependent torque cross terms.

Periodic boundaries are used for bulk construction checks. Open boundaries retain edge spins and omit every bond whose endpoint lies outside the flake; they do not wrap, clamp, or introduce ghost spins.

## Main files

- `cri3_source_structure.m` — exact 16-atom BiDB CIF record and provenance.
- `cri3_build_supercell.m` — four-Cr spin basis replication.
- `cri3_exchange_table.m` — approved Soenen-to-BiDB mapping, $S^2$ conversion, and fitted interlayer correction.
- `cri3_energy_fields.m` — energy and conservative tesla-valued fields.
- `cri3_llgs_rhs.m` and `cri3_integrate_rk4.m` — explicit LLGS and fixed-step RK4.
- `cri3_run.m` — validated simulation and provenance-rich result interface.
- `run_open_sot_example.m` — no-argument runnable candidate example.
- `cri3_assess_ground_state.m` and `cri3_transverse_hessian.m` — reproducible local-minimum diagnostics. These do not prove a global ground state.

Results contain the complete source record, site and bond tables, bond SHA-256, configuration, atom-resolved trajectory, layer and honeycomb averages, $\mathbf M$, $\mathbf L$, energy components, solver metadata, units, and source/mapping provenance.

## Verification

Run the deterministic smoke suite:

```bash
matlab -batch "addpath('tests'); run_all_tests"
```

Run the longer multistart, larger-supercell, Hessian-refinement, and material-timestep validation:

```bash
matlab -batch "addpath('tests'); run_ground_state_validation"
```

The extended acceptance thresholds were fixed after the first exploratory 2×2 comparison: residual torque below 0.002 T after 50 ps for the selected starts, no negative AFM transverse mode, 2×2/3×3 minimum-curvature agreement within $10^{-4}$ meV, and final-spin agreement below $10^{-10}$ between 1 fs and 0.5 fs for the deterministic short SOT trajectory.

The short open-flake test compares the complete central-spin histories for 5×5, 7×7, and 9×9 flakes and requires the successive size difference to decrease by at least a factor of 20. It validates boundary-aware propagation for that deterministic interval; it is not a finite-device switching or edge-reconstruction result.

## Evidence and limitations

Scientific sources, exact claim-to-implementation mappings, coefficient conversions, exclusions, and known hybrid-parameterization risks are recorded in [plan.md](plan.md). The source geometry is BiDB UID `2CrI3-1-2-1_0_0_1-0_-0.69`; the exchange model is based on Maarten Soenen *et al.*, *Physical Review Materials* **7**, 024421 (2023), DOI `10.1103/PhysRevMaterials.7.024421`.

The fitted scalar model reproduces the approved collinear identity $E_{\rm AFM}-E_{\rm FM}=-0.06$ meV per source cell. It does not reproduce the source exchange tensor's general angular dependence and must not be used to claim quantitative magnon bands, topology, or microscopic fitted bond interactions.
