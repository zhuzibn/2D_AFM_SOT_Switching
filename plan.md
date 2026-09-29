# Atomistic AB′ bilayer CrI₃ development plan

## Goal and Success Criteria

Build a standalone MATLAB material-specific atomistic model of monoclinic AB′ bilayer CrI₃. Construct the crystallographic cell from the supplied Bilayer Database (BiDB) record, retain all Cr and I coordinates in the structure provenance, and evolve one fixed-length classical unit spin $\mathbf s_i$ on each Cr site. Use the explicit Soenen-to-BiDB pair mapping recorded below. For the reduced scalar Hamiltonian, use the Table S.V isotropic intralayer exchange and the Table S.IV longitudinal interlayer components, converted by $S^2$, with one declared uniform interlayer fit correction constrained by the published collinear energy difference.

The first implementation will include isotropic exchange, out-of-plane single-ion anisotropy, Zeeman field, and layer-resolved spin-orbit torque (SOT). It will not claim a self-consistent single-DFT parameterization: the BiDB geometry and the Soenen Hamiltonian were generated with different electronic-structure workflows. DMI, Kitaev/anisotropic exchange, dipolar fields, temperature, defects, electronic bands, and topological-magnon calculations are excluded from the first baseline.

Success requires:

- exact reconstruction of the 16-atom BiDB source cell and its four Cr spin sites without a runtime network dependency;
- the documented Cr-site, layer, cell, and periodic-image mapping between the BiDB cell and every adopted Soenen exchange pair;
- explicit conversion of every adopted coefficient to joules per unique unit-spin bond, including separate provenance for the fitted correction and no use of the BiDB layer-level exchange field as a bond energy;
- analytical confirmation that the connected intralayer bond graph with the mapped negative exchange coefficients is minimized by layerwise FM order when considered on its own;
- separate validation of the fitted collinear AFM–FM calibration identity and assessment of the full zero-temperature ground state using perturbed and multistart noncollinear minimization, with an exchange-spectrum or transverse-stability check where practical;
- a staged boundary-condition workflow in which periodic boundaries establish bulk coordination and term-level correctness before open boundaries are used for finite-flake spatial switching;
- term-level energy/field checks, reciprocal and duplicate-free bonds, periodic supercell invariance, open-edge coordination, finite-size convergence, fixed spin norms, and time-step convergence; and
- preservation of source identity, units, transformations, and all modeling reductions in the saved result.

## Model Classification and Equation Family

This is **class 4: a material-specific crystallographic atomistic model**. The structural source cell contains Cr$_4$I$_{12}$; the dynamical degrees of freedom are the four magnetic Cr sites per cell. Iodine coordinates remain part of the structure record and validation but are not dynamical spins. Each monolayer is a two-site honeycomb lattice; the antiferromagnetic order parameter is primarily layer-resolved rather than the honeycomb A/B label.

Use fixed-length Gilbert LLGS in the explicit Landau–Lifshitz form of Eq. `eq:ll-with-flt` in [llg_equation.tex](</mnt/c/Users/zzf-m/OneDrive/presentations/equations/sections/llg_equation.tex>). The state is $\mathbf s\in\mathbb R^{N_{\mathrm{Cr}}\times3}$ with $|\mathbf s_i|=1$. Per-site moments $\mu_i$ are in J/T, energies are stored in J, effective flux densities and field-equivalent torque amplitudes are in T, and $\gamma$ is in rad/(s T).

Inputs are a validated embedded source cell, supercell size, boundary conditions, parameter-set identifier, applied field, SOT pulse, initial state, solver settings, and output path. Outputs are atom-resolved spins, layer and honeycomb-sublattice averages, net and Néel vectors, energy components, topology data, and complete provenance. The fixed-length invariant is $|\mathbf s_i|=1$ to the declared numerical tolerance. Electronic-structure quantities remain outside the classical spin model.

## Existing Behavior and Reference Implementations

At plan approval, the active workspace contained planning documents only and had no MATLAB implementation, tests, accepted numerical baseline, or usable Git repository. The approved implementation now provides the staged source, topology, Hamiltonian, LLGS, RK4, ground-state diagnostic, simulation-interface, and test files listed below. It remains a new standalone model rather than a modification of an existing solver; because usable Git metadata was absent during implementation, no pre-implementation Git diff or history is available.

The ASTRO reference repository was inspected at `atomistic_model/2D/atomistic/atomistic_v2`, specifically `AGENTS.md`, `README.md`, `astro_default_config.m`, `astro_evaluate_field_terms.m`, and `integrate_llg.m`. It demonstrates separation of deterministic fields, explicit tesla-valued effective fields, RK4 substage field updates, post-step normalization, and saved-time checks. Its square Gd–FeCo lattice, four-neighbor topology, species mapping, coefficients, and production defaults are not transferable to CrI₃.

The governing convention is taken from [hamiltonian.tex](</mnt/c/Users/zzf-m/OneDrive/presentations/equations/sections/hamiltonian.tex>), [llg_equation.tex](</mnt/c/Users/zzf-m/OneDrive/presentations/equations/sections/llg_equation.tex>), and [material_parameters.tex](</mnt/c/Users/zzf-m/OneDrive/presentations/equations/sections/material_parameters.tex>). The ordered-pair Hamiltonian in the equation notes carries $1/2$; this plan instead stores each physical bond once and omits that factor.

## State Variables, Geometry, and Boundaries

### Adopted source cell

The authoritative construction input is BiDB UID `2CrI3-1-2-1_0_0_1-0_-0.69`, derived from C2DB monolayer `2CrI3-1` [1]. The downloadable CIF is P1-expanded, while the database reports layer group `c2/m11` (number 18) with inversion symmetry. BiDB structures are assembled from C2DB monolayers, with the monolayers frozen during the interlayer-distance scan and sliding workflow [2, Methods]. The raw static-force array is therefore not a full-geometry relaxation acceptance test.

Use the BiDB lattice vectors exactly as stored in the CIF:

$$
\mathbf a_1=(7.007944101379477,0,0)\ \text{Å},
$$

$$
\mathbf a_2=(-3.5039720506897396,6.069057620095938,0)\ \text{Å},
$$

$$
\mathbf a_3=(0,0,27.197312819559507)\ \text{Å}.
$$

The source uses periodicity in $\mathbf a_1$ and $\mathbf a_2$ only. Preserve the complete Cr$_4$I$_{12}$ fractional-coordinate table from the downloaded CIF in the durable structure artifact. The four dynamical Cr sites are:

| BiDB site | Layer | Fractional coordinates $(f_1,f_2,f_3)$ |
|---|---:|---|
| Cr1 | bottom | $(0.6666666668402268,\ 0.3333333338772982,\ 0.3320831285767547)$ |
| Cr2 | bottom | $(0.3333333330534472,\ 0.6666666661068942,\ 0.33208314549019036)$ |
| Cr3 | top | $(0.6663373607355287,\ 0.6455589525179144,\ 0.578228585828896)$ |
| Cr4 | top | $(0.333004026948749,\ 0.9788922847475102,\ 0.5782286027423317)$ |

These coordinates and vectors are direct record values [1, `Crystal structure` CIF download]. Derived checks are a Cr-plane separation of approximately $6.694495$ Å, an intralayer nearest-Cr distance of approximately $4.046038$ Å, and a relative top-layer translation of approximately $(-0.0003293,0.3122256)$ in the BiDB fractional basis, equivalent modulo a lattice vector to the UID translation near $(0,-0.69)$. These are validation values, not replacements for the stored coordinates.

### Structure-construction evidence table

| Physical detail | Adopted value or rule | Builder or interaction mapping | Qualifying reference and exact locator | Evidence status |
|---|---|---|---|---|
| Material and cell contents | AB′ Cr$_4$I$_{12}$ bilayer; 16 atoms and four Cr spin sites | Preserve all source sites; select Cr1–Cr4 for dynamics | [1, `Overview` and CIF] | Direct database record |
| Cell and boundaries | Vectors listed above; $a=b=7.00794410137948$ Å, $\gamma=120^\circ$; periodic in-plane for bulk verification, open in-plane for finite-flake spatial studies, and always open out of plane | Wrap $f_1,f_2$ only in periodic mode; in open mode omit bonds whose mapped endpoint lies outside the finite flake; never wrap $f_3$ | [1, `Crystal structure` and CIF] for the cell; [5, abstract] and [6, abstract] for DMI-sensitive finite boundaries | Direct record plus source-supported software boundary choices |
| Symmetry | Layer group `c2/m11`, number 18, inversion present; distributed CIF is P1-expanded | Validate the coordinate set against reported symmetry; do not generate missing atoms from P1 | [1, `Crystal structure`] | Direct database record |
| Site coordinates and species | Full BiDB CIF; Cr coordinates reproduced above | Stable source IDs, original fractional coordinates, Cartesian coordinates, layer, honeycomb site, and supercell translation | [1, CIF] | Direct database record |
| AB′ registry | Top layer displaced by the source-record translation encoded in the UID; derived Cr-sublattice displacement stated above | Preserve coordinates unchanged; do not substitute an idealized $(1/3,0)$ shift | [1, UID and CIF]; qualitative AB′ comparison in [4, Figure 1] | Direct record plus transparent derivation |
| Cr-plane spacing and intralayer nearest distance | $6.694495$ Å and $4.046038$ Å, respectively | Geometry regression checks computed from the source vectors and coordinates | Derived from [1, CIF] | Transparent derivation |
| Intralayer topology | Honeycomb nearest neighbors; three per Cr in the periodic bulk cell | Generate periodic candidates, then validate endpoints and coordination against the source geometry | [4, main text and Supporting Information Eqs. (1)–(2)]; distances from [1, CIF] | Source topology plus coordinate derivation |
| BiDB interlayer distance groups | 7.043–7.117 Å (4 unique bonds/cell), 8.122–8.248 Å (4), and 9.017–9.188 Å (8) under the adopted unique-bond enumeration | Geometry diagnostics for the mapped NN, NNN, and 3NN pairs; distance never determines $J$ by itself | Derived from [1, CIF] using the mapping below | Transparent derivation; frozen mapping and multiplicities |
| Soenen site labels | Supplement Figure S10 uses a $3\times3\times1$ AB′ cell with Cr1–Cr18 in the top layer and Cr19–Cr36 in the bottom layer | Use the explicit label formula, lattice transformation, site bijection, and image offsets below | [3, Supplement Figure S10] | Direct labels plus documented graph/geometry derivation; resolved |
| Pair-resolved exchange | Intralayer: Table S.V scalar $J$. Interlayer: Table S.IV $J_{zz}$ scalarized for the reduced collinear model, plus the declared fit correction | Store source pair ID, BiDB endpoints, image offset, shell, source value, $S^2$ conversion, fitted term, final value, and provenance | [3, Supplement Tables S.IV–S.V and main Table II] | Source values plus explicit fitted reduction; resolved for the baseline |
| BiDB magnetic classification | AFM | Use only as an independent sign-level cross-check | [1, `Magnetic properties`] | Direct database classification; not a bond parameter |

### Construction rules

1. Embed an exact, auditable copy of the BiDB CIF data in a structure-source function or machine-readable artifact; include the record URL and access date.
2. Convert fractional positions with $\mathbf r=f_1\mathbf a_1+f_2\mathbf a_2+f_3\mathbf a_3$ and replicate $N_x\times N_y$ cells in deterministic `(cell, layer, source site)` order.
3. Create candidate pairs using integer in-plane image offsets. Canonicalize each unordered pair once and derive reciprocal solver adjacency from that table.
4. Generate the Soenen $3\times3\times1$ labels from the formulas and canonical transformation below, then verify every published pair against its expected BiDB distance. Similar distances alone are insufficient to establish a mapping.
5. Persist the source cell, complete atomic table, Cr selection, transformations, pair keys, shell/source labels, multiplicities, coefficient conversions, units, and hashes.

### Boundary-condition development workflow

Boundary selection follows the scientific question rather than a single universal default. Begin model development with periodic in-plane boundaries because they remove physical edges and make the repeated crystal, bulk coordination, pair counting, energy per cell, and supercell invariance easiest to audit. This periodic stage is both pedagogical and a verification reference; it represents an infinite tiled bilayer, so opposite displayed edges are neighbors and it must not be interpreted as an isolated device.

After the periodic structure, bond table, Hamiltonian, fields, and dynamics pass their gates, use open in-plane boundaries for spatially resolved finite-flake switching. Open means that Cr sites at the edge remain dynamical while every exchange or other finite-range bond whose mapped endpoint falls outside the flake is omitted. Do not clamp edge spins or create zero-valued ghost spins unless a separately sourced physical pinning mechanism is introduced. Open boundaries remove periodic wraparound, but opposite edges can still interact through a texture or propagating excitation when the flake is too small; therefore require finite-size convergence of both central-region and edge-resolved observables.

DMI is excluded from the first baseline, but the boundary design must remain compatible with it. A future bond-resolved DMI implementation must retain the directed-bond convention $\mathbf D_{ji}=-\mathbf D_{ij}$ for every retained bond and omit DMI bonds that cross an open edge. Free DMI boundaries can generate physical chiral edge canting, boundary-localized twists, and edge-assisted nucleation rather than a removable numerical artifact [5, abstract; 6, abstract]. Periodic boundaries suppress those physical edges and instead permit interaction with periodic images. Edge-modified material coefficients are not determined by the present bulk parameterization; truncating bulk bonds is the declared first finite-flake approximation, not a claim of an atomically reconstructed CrI$_3$ termination.

There is no out-of-plane wrapping. Mixed boundaries, such as periodic along a nanoribbon and open across its width, require a geometry-specific scientific objective and separate verification. Neither ordinary periodic nor open boundaries are nonreflecting: if outgoing spin-wave absorption becomes necessary, a damping buffer must be introduced and validated as a distinct numerical model.

### Resolved Soenen-to-BiDB mapping

For $u,v\in\{0,1,2\}$ and honeycomb index $s\in\{0,1\}$, decode Figure S10 as

$$
n_{\rm top}=1+6u+2v+s,
\qquad
n_{\rm bottom}=19+6u+2v+s.
$$

Choose the canonical, symmetry-equivalent mapping

$$
\mathbf a_{\rm Soenen}=-\mathbf a_{2,\rm BiDB},
\qquad
\mathbf b_{\rm Soenen}=\mathbf a_{1,\rm BiDB},
$$

with top $s=0\rightarrow$ Cr4, top $s=1\rightarrow$ Cr3, bottom $s=0\rightarrow$ Cr2, and bottom $s=1\rightarrow$ Cr1. Image offsets below are integer coefficients of $(\mathbf a_{1,\rm BiDB},\mathbf a_{2,\rm BiDB})$. Other mappings that satisfy Figure S10 and all distances are global lattice rotations/reflections; selecting one is a coordinate convention, and all device axes must use the same convention.

For intralayer exchange, the adopted Table S.V values and unit-spin coefficients are:

| Source pair(s) | BiDB pair(s) and image offset | $J^{\rm source}$ (meV) | $J^{\rm unit}=S^2J^{\rm source}$ (meV) |
|---|---|---:|---:|
| 1–2 | Cr4→Cr3 $(0,0)$ | $-4.62$ | $-10.395$ |
| 2–3 | Cr3→Cr4 $(1,0)$ | $-4.20$ | $-9.450$ |
| 2–7 | Cr3→Cr4 $(0,-1)$ | $-4.64$ | $-10.440$ |
| 1–3; 4–8 | Cr4→Cr4 $(1,0)$; Cr3→Cr3 $(-1,-1)$ | $-0.68$ | $-1.530$ |
| 1–7; 2–4 | Cr4→Cr4 $(0,-1)$; Cr3→Cr3 $(1,0)$ | $-0.62$ | $-1.395$ |
| 3–7; 2–8 | Cr4→Cr4 $(-1,-1)$; Cr3→Cr3 $(0,-1)$ | $-0.62$ | $-1.395$ |

For interlayer exchange, define the fitted scalar coefficient

$$
J_b^{\rm fit}=S^2J_{zz,b}^{\rm source}+\delta J_{\rm fit},
\qquad
S=\frac32,
\qquad
\delta J_{\rm fit}=0.0075\ \mathrm{meV}
$$

for each of the 16 unique interlayer bonds per cell:

| Shell | Source pair | BiDB pair and image offset | Distance (Å) | $J_{zz}^{\rm source}$ (meV) | $J_b^{\rm fit}$ (meV) |
|---|---|---|---:|---:|---:|
| NN | 7–20 | Cr4→Cr1 $(0,1)$ | 7.116569 | $-0.30$ | $-0.6675$ |
| NN | 1–19 | Cr4→Cr2 $(0,0)$ | 7.043360 | $-0.29$ | $-0.6450$ |
| NN | 2–19 | Cr3→Cr2 $(0,0)$ | 7.115432 | $-0.28$ | $-0.6225$ |
| NN | 2–20 | Cr3→Cr1 $(0,0)$ | 7.043360 | $-0.29$ | $-0.6450$ |
| NNN | 9–21 | Cr4→Cr2 $(0,1)$ | 8.248420 | $+0.13$ | $+0.3000$ |
| NNN | 9–26 | Cr4→Cr1 $(-1,0)$ | 8.121772 | $+0.18$ | $+0.4125$ |
| NNN | 10–22 | Cr3→Cr1 $(0,1)$ | 8.248420 | $+0.13$ | $+0.3000$ |
| NNN | 10–29 | Cr3→Cr2 $(1,0)$ | 8.123762 | $+0.19$ | $+0.4350$ |
| 3NN | 9–20 | Cr4→Cr1 $(-1,1)$ | 9.186441 | $+0.12$ | $+0.2775$ |
| 3NN | 9–23 | Cr4→Cr2 $(1,1)$ | 9.132501 | $+0.06$ | $+0.1425$ |
| 3NN | 9–25 | Cr4→Cr2 $(-1,0)$ | 9.129844 | $+0.06$ | $+0.1425$ |
| 3NN | 9–28 | Cr4→Cr1 $(0,0)$ | 9.018277 | $+0.01$ | $+0.0300$ |
| 3NN | 10–23 | Cr3→Cr2 $(1,1)$ | 9.188200 | $+0.12$ | $+0.2775$ |
| 3NN | 10–24 | Cr3→Cr1 $(1,1)$ | 9.132501 | $+0.06$ | $+0.1425$ |
| 3NN | 10–26 | Cr3→Cr1 $(-1,0)$ | 9.129844 | $+0.05$ | $+0.1200$ |
| 3NN | 10–33 | Cr3→Cr2 $(0,-1)$ | 9.017380 | $+0.01$ | $+0.0300$ |

## Hamiltonian

The first baseline is a scalar reduction of the Soenen Hamiltonian. Each physical bond is stored once:

$$
\mathcal H
=\sum_{\{i,j\}\in\mathcal B}J_{ij}\,\mathbf s_i\cdot\mathbf s_j
-\sum_i K_i(\mathbf s_i\cdot\hat{\mathbf z})^2
-\sum_i\mu_i\mathbf s_i\cdot\mathbf B_{\rm ext}.
$$

Negative $J_{ij}$ favors ferromagnetic alignment, positive $J_{ij}$ favors antiferromagnetic alignment, and positive $K_i$ gives an out-of-plane easy axis. $J_{ij}$ is energy per unique bond and $K_i$ is energy per Cr site. SOT is nonconservative and is not part of $\mathcal H$.

Soenen *et al.* use a generalized exchange matrix and single-ion anisotropy for spins of magnitude $S=3/2$ [3, main Eqs. (1)–(2) and Supplement Tables S.IV–S.V]. Replacing each source spin by $\hat{\mathbf S}_i=S\mathbf s_i$ gives the unit-spin conversion

$$
J_{ij}^{\rm unit}=S^2J_{ij}^{\rm source},
\qquad
K_i=-S^2A_{zz}.
$$

The source uses a $1/2$ ordered-pair sum; the code uses each physical bond once with the same converted coefficient and no $1/2$. With $A_{zz}=-0.07$ meV, the reduced easy-axis coefficient is $K=0.1575$ meV per Cr.

For intralayer links, use the isotropic $J$ defined by the local-eigenbasis construction in Table S.V. For interlayer links, use the global $J_{zz}$ entries from Table S.IV as effective scalar coefficients for the deliberately collinear reduction, then add $\delta J_{\rm fit}$ to every interlayer bond. This does not reproduce the full angular dependence of the source exchange tensor and must not be used for quantitative magnon or topology claims. Do not enable the tabulated Kitaev, DMI, off-diagonal SIA, or other tensor components in this baseline.

The rounded source interlayer values satisfy

$$
\sum_{b\in\mathcal B_\perp}S^2J_{zz,b}^{\rm source}
=-0.09\ \mathrm{meV}.
$$

The adopted correction gives

$$
\sum_{b\in\mathcal B_\perp}J_b^{\rm fit}
=-0.09+16(0.0075)
=+0.03\ \mathrm{meV},
$$

and therefore

$$
E_{\rm AFM}-E_{\rm FM}
=-2\sum_{b\in\mathcal B_\perp}J_b^{\rm fit}
=-0.06\ \mathrm{meV}
$$

per source cell, exactly matching the calibration target in [3, main Table II]. The $0.0075$ meV correction is a fitted model coefficient, not a reported first-principles bond interaction.

## Parameters and Evidence

| Quantity | Adopted/provisional value | Units and representation | Source and status | Model use |
|---|---:|---|---|---|
| BiDB source UID | `2CrI3-1-2-1_0_0_1-0_-0.69` | Stable database identifier | [1] adopted | Geometry provenance |
| Lattice and sites | Values in the geometry section | Å and fractional coordinates | [1, CIF] adopted | Structure builder |
| Cr spin magnitude | $S=3/2$ | Source spin variable | [3, main Sec. II] adopted | $S^2=2.25$ unit-spin conversion |
| Cr moment | $\mu=3\mu_B$ | $2.7822\times10^{-23}$ J/T per Cr | [3, atomistic LLG setup] adopted | Effective field and dynamics |
| Intralayer scalar exchange | Exact Table S.V $J$ rows and mappings listed above | $S^2J^{\rm source}$ meV per unique unit-spin bond; convert with $1\ \mathrm{meV}=1.602176634\times10^{-22}$ J | [3, Supplement Table S.V] adopted reduced model | Intralayer bonds |
| Interlayer scalar exchange | Exact mapped Table S.IV $J_{zz}$ rows listed above | $S^2J_{zz}^{\rm source}+0.0075$ meV per unique unit-spin bond | [3, Supplement Table S.IV and main Table II] adopted fitted reduction | Interlayer bonds |
| Interlayer fit correction | $\delta J_{\rm fit}=+0.0075$ meV per each of 16 links | $1.2016324755\times10^{-24}$ J/link | Derived from [3, main Table II] and rounded Table S.IV; user-approved modeling choice | Enforces the collinear splitting only |
| Single-ion anisotropy | $K=-S^2A_{zz}=+0.1575$ | meV per Cr in $-Ks_z^2$ | $A_{zz}=-0.07$ meV from [3, main Table II] | Easy-axis term |
| AB′ collinear calibration target | $E_{\rm AFM}-E_{\rm FM}=-0.06$ | meV per source cell | [3, main Table II] adopted fit target | Exact reduced-model constraint |
| Gyromagnetic ratio | $\gamma=g\mu_B/\hbar$, with $g\approx2$ | rad/(s T) | [3, atomistic LLG setup] adopted | LLGS |
| Gilbert damping | No production default; caller must set `config.alpha` explicitly. Soenen used $\alpha=0.001$ for a separate $T=0$ stability calculation | Dimensionless, scalar or per-site | [3, Supplement Sec. S.IV and Eq. (1)] source value documented but not assumed transferable to every SOT run | LLGS and numerical relaxation |
| DL/FL torque amplitudes | Direct signed inputs $\tau_{{\rm DL},\ell}(t)$ and $\tau_{{\rm FL},\ell}(t)$ | T, layer resolved | Project interface choice | SOT drive |
| Temperature | $T=0$ | Deterministic | Deliberate baseline choice | No thermal field |

The Soenen AB′ interlayer entries retain their competing ferromagnetic NN and antiferromagnetic NNN/3NN signs [3, Supplement Table S.IV]. The common positive correction changes none of those signs but adjusts their nearly cancelling total to the published collinear splitting. Sivadas reports that the weak AB′ magnetic preference is sensitive to the electronic-structure setup [4]; therefore the fitted model is a calibrated reduced baseline, not a new first-principles parameter set.

## Effective Fields and Units Check

Using $\mathbf B_{\mathrm{eff},i}=-(1/\mu_i)\,\partial\mathcal H/\partial\mathbf s_i$ from [hamiltonian.tex](</mnt/c/Users/zzf-m/OneDrive/presentations/equations/sections/hamiltonian.tex>),

$$
\mathbf B_{\mathrm{eff},i}
=-\frac{1}{\mu_i}\sum_{j\in\mathcal N(i)}J_{ij}\mathbf s_j
+\frac{2K_i}{\mu_i}(\mathbf s_i\cdot\hat{\mathbf z})\hat{\mathbf z}
+\mathbf B_{\rm ext}.
$$

With $J$ and $K$ in J and $\mu$ in J/T, each ratio is in T. No factor $\mu_0$ appears because the integrated quantity is flux density $\mathbf B$, not field strength $\mathbf H$. If an input is supplied in A/m, convert it first with $\mathbf B=\mu_0\mathbf H$. Finite differences of the implemented energy must reproduce each conservative field contribution under the unit-vector convention.

## Dynamical Equation and Nonconservative Torques

For Cr site $i$ in layer $\ell(i)$, integrate Eq. `eq:ll-with-flt`:

$$
\begin{aligned}
(1+\alpha_i^2)\dot{\mathbf s}_i={}&
-\gamma_i\mathbf s_i\times\mathbf B_{\mathrm{eff},i}
-\alpha_i\gamma_i\mathbf s_i\times(\mathbf s_i\times\mathbf B_{\mathrm{eff},i})\\
&+\alpha_i\gamma_i\tau_{{\rm DL},\ell}\mathbf s_i\times\boldsymbol\sigma
-\gamma_i\tau_{{\rm DL},\ell}\mathbf s_i\times(\mathbf s_i\times\boldsymbol\sigma)\\
&+\alpha_i\gamma_i\tau_{{\rm FL},\ell}\mathbf s_i\times(\mathbf s_i\times\boldsymbol\sigma)
+\gamma_i\tau_{{\rm FL},\ell}\mathbf s_i\times\boldsymbol\sigma .
\end{aligned}
$$

Divide the complete right-hand side by $1+\alpha_i^2$. Retain both damping-dependent torque cross terms and validate the explicit equation against the implicit Gilbert residual. Use layer-resolved signed amplitudes so both uniform and staggered torque configurations are representable without changing the solver.

## Initial Conditions, Solver, and Outputs

Initialize each layer internally ferromagnetic and the two layer magnetizations antiparallel along $\pm\hat{\mathbf z}$, with a separately configurable small deterministic canting only when needed to break an exact symmetry. Use $\boldsymbol\sigma=+\hat{\mathbf y}$ for the initial SOT convention and record current, layer-normal, and polarization directions with every run.

Use deterministic fixed-step RK4, recomputing fields and torques at every substage and normalizing spins only after each complete step. Start with $\Delta t=1$ fs and require convergence against $0.5$ fs; reduce further if the largest mapped exchange field demands it. Include the exact initial and requested final states, and require `state(k)` to correspond to `time(k)`. Develop and verify the periodic in-plane configuration first; after its structure, topology, static-field, RHS, and integrator gates pass, use open boundaries as the production configuration for spatially resolved finite-flake switching. Use no RNG in the baseline production dynamics; ground-state validation may use a fixed, recorded seed solely to construct reproducible noncollinear initial states.

Save atom-resolved spins, layer magnetizations $\mathbf m_1,\mathbf m_2$, honeycomb-site averages, $\mathbf L=(\mathbf m_1-\mathbf m_2)/2$, $\mathbf M=(\mathbf m_1+\mathbf m_2)/2$, energy components, spin norms, configuration, units, source IDs, source and simulation coordinates, bond table and hash, exchange-conversion audit, solver, timestamps, and code version when available. Exploratory runs must use non-destructive output paths.

## Implementation Mapping

Treat the stages below as explicit implementation-and-test checkpoints. Do not begin a later stage until the focused checks for the preceding stage pass; a failed gate returns the work to the component responsible for that behavior.

| Stage | Implementation | Verification gate |
|---|---|---|
| 1. Structure | Source cell, Cr-site selection, layer/site labels, supercell replication, and explicit periodic/open boundary mode | Exact 16-atom cell, four Cr sites with two per layer, correct derived distances, no duplicate wrapped sites, and no wrapping in open mode |
| 2. Interaction topology | Map every adopted Soenen pair to BiDB endpoints and periodic-image offsets; construct periodic bulk and truncated open-flake unique-bond tables with reciprocal adjacency | Periodic unique bonds, reciprocal adjacency, correct 4/4/8 interlayer multiplicities, fitted bond sum, and supercell invariance; open edges have only retained in-flake bonds and no wraparound |
| 3. Static Hamiltonian | Exchange, anisotropy, and Zeeman energy and effective flux-density contributions | Unit conversions, pair counting, finite-difference field checks, the isolated intralayer FM minimum, and the restricted collinear AFM–FM calibration identity |
| 4. LLGS right-hand side | Precession and damping first, followed by layer-resolved DL and FL torque terms | Tangent dynamics $\mathbf s_i\cdot\dot{\mathbf s}_i=0$, precession direction, implicit Gilbert residual, and torque signs |
| 5. Time integrator | RK4 with fields and torques recomputed at every substage | Finite values, spin norms, exact saved-time indexing, energy conservation/relaxation limits, and timestep convergence |
| 6. Ground-state assessment | Reproducible no-drive perturbed and multistart noncollinear relaxations, with exchange-spectrum or transverse-Hessian analysis where practical | Converged energies and residual torques, layer/order-parameter classification, periodic-size checks, and no claim of globality from relaxation alone |
| 7. SOT simulations and outputs | Configuration entry point, layer-resolved pulses, boundary-tagged observables, and result provenance; periodic control before open finite-flake production runs | DL-only/FL-only behavior, layer polarity, approved switching behavior, reproducible non-destructive saved results, and open-flake size convergence |

The architectural boundary is:

```text
structure -> bond table -> energy/effective field -> LLGS RHS -> integrator -> simulation
```

`cri3_llgs_rhs.m` must not construct geometry, infer periodic images, interpret source exchange tables, or apply coefficient conversions. It receives an already verified effective flux density and explicit torque inputs. Conversely, the conservative-field evaluator does not integrate spins or embed nonconservative SOT. This separation makes a failed trajectory traceable to structure/topology, energy-to-field derivation, torque signs, integration, or simulation configuration. In particular, implementation must not jump from interaction mapping directly to LLGS: the static Hamiltonian and effective fields must pass Stage 3 before dynamical work begins.

1. **Source structure and supercell:** implement `cri3_source_structure.m` with the immutable BiDB vectors, full 16-site fractional table, and provenance, then implement `cri3_build_supercell.m` for Cr selection, layer/site labels, replication, explicit periodic/open semantics, and deterministic site ordering. Gate: recover the exact source-cell contents, four Cr sites with two per layer, the declared geometric checks, no duplicate wrapped sites in periodic mode, and no coordinate wrapping in open mode.
2. **Interaction topology and coefficients:** implement `cri3_exchange_table.m` with the resolved Soenen pair mapping, source values, $S^2$ conversion, $\delta J_{\rm fit}$, final joule-valued coefficients, per-row provenance, and boundary-aware bond retention. Gate: reproduce every endpoint/image mapping, periodic unique and reciprocal bonds, the 4/4/8 interlayer multiplicities, the fitted bond sum, and supercell-invariant bulk topology; then verify that open mode removes every crossing bond without pinning or ghost neighbors.
3. **Static Hamiltonian and effective fields:** implement `cri3_energy_fields.m` for scalar exchange, anisotropy, Zeeman energy, and tesla-valued fields without introducing time integration. Gate: pass unit-conversion and pair-counting audits, FM/AFM two-spin limits, finite-difference field checks, the isolated intralayer FM minimum, and the restricted collinear AFM–FM calibration identity. Do not interpret this gate as proof of the full model ground state.
4. **LLGS right-hand side:** implement `cri3_llgs_rhs.m` against the already verified effective-field interface, first checking precession and damping and then the DL/FL terms. Gate: verify tangent dynamics, precession direction, the explicit-versus-implicit Gilbert residual, and layer-resolved DL-only and FL-only torque directions.
5. **Time integration:** implement `cri3_integrate_rk4.m`, recomputing fields and torques at every substage and normalizing only after each complete step. Gate: verify finite states, fixed spin norms, exact saved-time alignment, zero-damping energy behavior, positive-damping relaxation, and 1 fs versus 0.5 fs convergence.
6. **Ground-state assessment:** after the static Hamiltonian, RHS, and integrator gates pass, run the no-drive perturbed and multistart noncollinear relaxation campaign defined below and, where practical, the independent exchange-spectrum or transverse-Hessian analysis. Gate: report converged energies, residual torques, order parameters, periodic-size dependence, and the remaining limit on any global-ground-state claim.
7. **Simulation interface, SOT runs, and outputs:** add one documented configuration/entry point, layer-resolved signed torque inputs, observables, provenance-rich saving, and the focused end-to-end tests. Gate: verify the approved project-specific deterministic SOT behavior and reproducible outputs.

Use arrays shaped `(NCr,3)` for spins/positions and a bond table containing `(i,j,image_offset,source_pair,shell,J_Joule,provenance)`. Keep structure construction, Hamiltonian mapping, SOT calibration, solver behavior, and optional physics as separate changes. Do not import ASTRO's square-lattice indexing or material parameters.

## Verification and Acceptance Criteria

### Structure and topology

- Reproduce all 16 BiDB fractional coordinates and both in-plane lattice vectors to stored floating-point precision; recover Cr$_4$I$_{12}$ with no duplicate wrapped sites.
- Recover four Cr sites per cell, two per layer, Cr-plane separation $6.694495$ Å, and intralayer nearest distance $4.046038$ Å within $10^{-5}$ Å.
- Reproduce the frozen periodic interlayer distance/multiplicity table: four NN bonds at 7.043–7.117 Å, four NNN bonds at 8.122–8.248 Å, and eight 3NN bonds at 9.017–9.188 Å per cell; every unordered bond key must be unique and every solver adjacency reciprocal.
- Reproduce the documented bijection between every adopted Soenen source pair and `(BiDB endpoints, integer image offset)` in the reconstructed $3\times3\times1$ cell. No exchange is assigned by cutoff alone.
- Verify that two supercell sizes have identical bulk coordination, per-cell pair multiplicities, and energy per cell for a repeated spin state.
- On a deterministic small flake, verify that open boundaries contain no wrapped bonds, every retained adjacency is reciprocal, edge coordination is reduced only by explicitly omitted crossing bonds, and interior coordination matches the periodic bulk value.
- Repeat the spatial switching calculation for at least two increasing open-flake sizes. Require the central-region trajectory and selected switching observables to converge within a tolerance declared after the first comparison, while retaining edge-resolved outputs instead of averaging away boundary physics.
- When DMI is later added, test periodic seam orientation and $\mathbf D_{ji}=-\mathbf D_{ij}$ separately from open-edge bond truncation and edge canting; do not reuse the scalar-exchange boundary test as sufficient DMI validation.

### Hamiltonian and dynamics

- Verify each source coefficient, the fixed $S^2=2.25$ conversion, the separately stored $+0.0075$ meV fit correction, meV-to-J conversion, and ordered-versus-unique pair convention in a machine-readable audit table.
- Check FM and AFM two-spin limits and finite-difference derivatives for exchange and anisotropy.
- Require the 16 uncorrected converted interlayer coefficients to sum to $-0.09$ meV/cell, the correction to sum to $+0.12$ meV/cell, and the final coefficients to sum to $+0.03$ meV/cell.
- Require direct energy evaluation to give $E_{\rm AFM}-E_{\rm FM}=-0.06$ meV per source cell to floating-point tolerance within the collinear layerwise-FM subspace. This is a calibration, sign, and pair-counting identity, not independent ground-state evidence.
- Confirm analytically that every mapped intralayer exchange coefficient is negative and that the periodic intralayer bond graph is connected. Under the declared Hamiltonian convention, this proves that the isolated intralayer exchange energy is minimized by a uniform spin direction in each layer; it does not by itself prove the ground state after competing interlayer terms are included.
- Assess the full zero-temperature model with no SOT, applied field, or thermal noise using reproducible starts that include small perturbations of both collinear FM and AFM states, prescribed canted or modulated states, and multiple fixed-seed noncollinear states. Do not infer stability from an exactly collinear trajectory because such a state can have zero torque even when it is not the global minimum.
- For every relaxation, require energy nonincrease within the converged numerical tolerance, an energy plateau, and a small maximum residual torque $\max_i|\mathbf s_i\times\mathbf B_{{\rm eff},i}|$, with tolerances declared from convergence studies rather than chosen retroactively. Compare final energies and report $|\mathbf m_1|$, $|\mathbf m_2|$, $\mathbf m_1\cdot\mathbf m_2$, $|\mathbf M|$, and $|\mathbf L|$ instead of classifying the result from the net magnetization alone.
- Repeat the ground-state assessment for increasing periodic supercells. Where practical, independently examine the Fourier-space exchange matrix or the transverse Hessian about the lowest-energy candidate to check for a lower finite-$q$ mode or unstable direction. Treat multistart LLGS relaxation alone as evidence for a low-energy local minimum, not as proof of the global ground state.
- At zero damping and without drive, check approximate energy conservation; with positive damping, check monotonic relaxation toward a local minimum.
- Check finite values, $\max_i||\mathbf s_i|-1|$, precession direction, explicit-versus-implicit Gilbert residual, and DL-only/FL-only torque directions on both layers.
- Require the 1 fs and 0.5 fs trajectories to agree within a declared observable tolerance chosen after a convergence study; do not invent the tolerance before the first deterministic run.

## References

[1] Bilayer Database (BiDB), record `2CrI3-1-2-1_0_0_1-0_-0.69`, derived from C2DB monolayer `2CrI3-1`, fields `Overview`, `Crystal structure`, `Magnetic properties`, downloadable CIF, and downloadable JSON. Direct record: [https://bidb.fysik.dtu.dk/material/2CrI3-1-2-1_0_0_1-0_-0.69](https://bidb.fysik.dtu.dk/material/2CrI3-1-2-1_0_0_1-0_-0.69). Accessed 2026-09-27.

[2] Sahar Pakdel, Asbjørn Rasmussen, Alireza Taghizadeh, Mads Kruse, Thomas Olsen, Kristian S. Thygesen, High-throughput computational stacking reveals emergent properties in natural van der Waals bilayers, *Nature Communications* **15**, 932 (2024). DOI: [10.1038/s41467-024-45003-w](https://doi.org/10.1038/s41467-024-45003-w).

[3] Maarten Soenen, Cihan Bacaksiz, Raí M. Menezes, Milorad V. Milošević, Stacking-dependent topological magnons in bilayer CrI₃, *Physical Review Materials* **7**, 024421 (2023). DOI: [10.1103/PhysRevMaterials.7.024421](https://doi.org/10.1103/PhysRevMaterials.7.024421). Local paper ID: `2D-TI-van der vaals/2023-Stacking-dependent topological magnons in bilayer CrI3-PRM-Maarten Soenen`. The reviewed summary and indexed main and supplementary source PDFs supplied the Hamiltonian, Figure S10 labels, Table S.IV $J_{zz}$ values, and Table S.V scalar values; the PDFs were used because the Markdown tables are flattened.

[4] Nikhil Sivadas, Satoshi Okamoto, Xiaodong Xu, Craig J. Fennie, Di Xiao, Stacking-Dependent Magnetism in Bilayer CrI3, *Nano Letters* **18**, 7658–7664 (2018). DOI: [10.1021/acs.nanolett.8b03321](https://doi.org/10.1021/acs.nanolett.8b03321). Local paper ID: `2D-TI-van der vaals/2018-Stacking-Dependent Magnetism in Bilayer CrI-nano letters-Nikhil Sivadas`. The reviewed combined summary and indexed main/SI records supplied the evidence.

[5] Stanislas Rohart, André Thiaville, Skyrmion confinement in ultrathin film nanostructures in the presence of Dzyaloshinskii-Moriya interaction, *Physical Review B* **88**, 184422 (2013). DOI: [10.1103/PhysRevB.88.184422](https://doi.org/10.1103/PhysRevB.88.184422). The online journal abstract supplied the evidence that DMI changes the natural magnetization behavior at finite edges; it is not used as a CrI$_3$ material-parameter source.

[6] Oleksandr V. Pylypovskyi, Artem V. Tomilo, Denis D. Sheka, Jürgen Fassbender, Denys Makarov, Boundary conditions for the Néel order parameter in a chiral antiferromagnetic slab, *Physical Review B* **103**, 134413 (2021). DOI: [10.1103/PhysRevB.103.134413](https://doi.org/10.1103/PhysRevB.103.134413). The online journal abstract supplied qualitative antiferromagnetic evidence that DMI and finite boundaries modify spatial textures; its G-type slab results are not transferred quantitatively to A-type bilayer CrI$_3$.

## Assumptions, Missing References, Risks, and Approval Questions

1. **Coordinate convention:** the chosen Soenen-to-BiDB mapping is one of four symmetry-equivalent global orientations. The device and torque axes must remain consistent with this choice.
2. **Hybrid parameterization risk and approved fidelity label:** BiDB uses frozen C2DB monolayers and its own stacking workflow, whereas Soenen uses a separately relaxed PBE-D2+$U_{\rm eff}=2.8$ eV structure. The user approved the label `BiDB geometry + fitted Soenen scalar reduction` on 2026-09-28. Use the BiDB coordinates unchanged and do not silently rescale either geometry.
3. **Approved fitted interlayer correction:** the user approved $\delta J_{\rm fit}=+0.0075$ meV/link on 2026-09-28. It is derived from one collinear observable and serves as the reduced-model fallback. It corrects the net AFM–FM splitting but is not a microscopic first-principles interaction. The model must not claim quantitative tensor anisotropy, magnon bands, topology, or general noncollinear energetics.
4. **Weak-interlayer-magnetism risk:** AB′ CrI₃ lies near a magnetic phase boundary and published signs depend on electronic-structure choices [4]. The BiDB AFM label and the fitted Soenen reduction are cross-checks, not interchangeable parameterizations.
5. **Finite-boundary risk:** open boundaries remove direct wraparound but do not guarantee independent edges in a small flake, so size convergence remains mandatory. The first open-edge model truncates bulk bonds and does not model edge reconstruction or edge-specific parameters.
6. **Approved Soenen-to-BiDB mapping:** the user explicitly approved the revised mapping on 2026-09-29. Implementation must preserve the selected global coordinate convention and every documented endpoint/image offset; any material remapping requires renewed approval.
7. **Unspecified production drive:** no material-specific damping transfer, SOT pulse amplitude, duration, or switching acceptance threshold has been approved. The implementation therefore requires damping explicitly, accepts signed layer-resolved torque amplitudes in tesla, and limits built-in checks to equation signs, numerical convergence, and deterministic open-boundary behavior rather than claiming a physical switching threshold.

There are no remaining approval questions for the implemented baseline. A physical switching study requires a separately approved damping/drive/device parameterization and acceptance observable.
