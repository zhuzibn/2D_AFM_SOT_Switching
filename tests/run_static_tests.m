function run_static_tests()
%RUN_STATIC_TESTS Deterministic Hamiltonian and effective-field gates.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
meV_to_J = 1.602176634e-22;
config = cri3_default_config();
parameters.mu_J_per_T = config.mu_J_per_T;
parameters.K_J = config.K_J;
parameters.B_external_T = [0, 0, 0];

system1 = cri3_build_supercell(1, 1, 'periodic');
bonds1 = cri3_exchange_table(system1);
fm1 = repmat([0, 0, 1], system1.n_sites, 1);
afm1 = fm1;
afm1(system1.layer == 2, 3) = -1;
[energy_fm1, ~] = cri3_energy_fields(fm1, system1, bonds1, parameters);
[energy_afm1, ~] = cri3_energy_fields(afm1, system1, bonds1, parameters);
delta_meV = (energy_afm1.total_J - energy_fm1.total_J) / meV_to_J;
assert(abs(delta_meV - (-0.06)) < 1e-11);

system2 = cri3_build_supercell(2, 2, 'periodic');
bonds2 = cri3_exchange_table(system2);
fm2 = repmat([0, 0, 1], system2.n_sites, 1);
[energy_fm2, ~] = cri3_energy_fields(fm2, system2, bonds2, parameters);
assert(abs(energy_fm2.total_J / 4 - energy_fm1.total_J) < 1e-32);

rng(1729, 'twister');
spins = randn(system2.n_sites, 3);
spins = spins ./ vecnorm(spins, 2, 2);
parameters.B_external_T = [0.013, -0.021, 0.034];
[energy0, fields0] = cri3_energy_fields(spins, system2, bonds2, parameters);
epsilon = 1e-7;
sites_to_check = [1, 6, 13];
for site = sites_to_check
    for component = 1:3
        shifted = spins;
        shifted(site, component) = shifted(site, component) + epsilon;
        energy_plus = cri3_energy_fields(shifted, system2, bonds2, parameters);
        shifted(site, component) = shifted(site, component) - 2 * epsilon;
        energy_minus = cri3_energy_fields(shifted, system2, bonds2, parameters);
        derivative = (energy_plus.total_J - energy_minus.total_J) / (2 * epsilon);
        expected = -parameters.mu_J_per_T * fields0.total_T(site, component);
        scale = max([abs(derivative), abs(expected), 1e-30]);
        assert(abs(derivative - expected) / scale < 2e-7);
    end
end
assert(isfinite(energy0.total_J));

two_site.n_sites = 2;
two_bonds.i = 1;
two_bonds.j = 2;
two_bonds.J_J = -2.3e-22;
two_bonds.count = 1;
two_parameters.mu_J_per_T = config.mu_J_per_T;
two_parameters.K_J = 0;
two_parameters.B_external_T = [0, 0, 0];
parallel = [0, 0, 1; 0, 0, 1];
antiparallel = [0, 0, 1; 0, 0, -1];
energy_parallel = cri3_energy_fields(parallel, two_site, two_bonds, two_parameters);
energy_antiparallel = cri3_energy_fields(antiparallel, two_site, two_bonds, two_parameters);
assert(energy_parallel.total_J < energy_antiparallel.total_J);
two_bonds.J_J = -two_bonds.J_J;
energy_parallel = cri3_energy_fields(parallel, two_site, two_bonds, two_parameters);
energy_antiparallel = cri3_energy_fields(antiparallel, two_site, two_bonds, two_parameters);
assert(energy_antiparallel.total_J < energy_parallel.total_J);

fprintf('Stage 3 static Hamiltonian tests passed; AFM-FM = %.12g meV/cell.\n', ...
    delta_meV);
end
