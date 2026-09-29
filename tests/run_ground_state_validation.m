function run_ground_state_validation()
%RUN_GROUND_STATE_VALIDATION Extended, reproducible scientific validation.
%   Runtime is longer than the smoke suite. Tolerances were fixed after the
%   first 2x2 exploratory comparison documented during implementation.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
meV_to_J = 1.602176634e-22;
config = cri3_default_config();
config.nx = 2;
config.ny = 2;
options.t_final_s = 5e-11;
options.dt_s = 1e-15;
options.save_every = 1000;
options.random_seeds = [17, 29];
options.hessian_epsilon = 1e-4;
assessment = cri3_assess_ground_state(config, options);

assert(strcmp(assessment.lowest_candidate_label, 'collinear-afm'));
assert(max([assessment.candidates.max_residual_torque_T]) < 0.002);
assert(max([assessment.candidates.maximum_energy_increase_J]) < 1e-29);
assert(abs(assessment.candidates(1).final_energy_J / (4 * meV_to_J) ...
    - (-78.51)) < 1e-9);
assert(assessment.lowest_transverse_hessian.negative_mode_count == 0);
assert(assessment.lowest_transverse_hessian.minimum_eigenvalue_meV > 0.314);

system2 = cri3_build_supercell(2, 2, 'periodic');
bonds2 = cri3_exchange_table(system2);
parameters.mu_J_per_T = config.mu_J_per_T;
parameters.K_J = config.K_J;
parameters.B_external_T = [0, 0, 0];
afm2 = cri3_initial_state(system2, 'afm', 0);
hessian_refined = cri3_transverse_hessian(afm2, system2, bonds2, ...
    parameters, 5e-5);
assert(abs(hessian_refined.minimum_eigenvalue_meV - ...
    assessment.lowest_transverse_hessian.minimum_eigenvalue_meV) < 1e-4);

system3 = cri3_build_supercell(3, 3, 'periodic');
bonds3 = cri3_exchange_table(system3);
afm3 = cri3_initial_state(system3, 'afm', 0);
[energy3, fields3] = cri3_energy_fields(afm3, system3, bonds3, parameters);
hessian3 = cri3_transverse_hessian(afm3, system3, bonds3, parameters, 1e-4);
assert(abs(energy3.total_J / (9 * meV_to_J) - (-78.51)) < 1e-9);
assert(max(vecnorm(cross(afm3, fields3.total_T, 2), 2, 2)) < 1e-12);
assert(hessian3.negative_mode_count == 0);
assert(abs(hessian3.minimum_eigenvalue_meV - ...
    assessment.lowest_transverse_hessian.minimum_eigenvalue_meV) < 1e-4);

config.alpha = 0.001;
config.t_final_s = 5e-14;
config.initial_canting = 1e-3;
config.torque.dl_T = [0.01, -0.01];
config.dt_s = 1e-15;
coarse = cri3_run(config);
config.dt_s = 0.5e-15;
fine = cri3_run(config);
spin_difference = max(abs(coarse.spins(:, :, end) - ...
    fine.spins(:, :, end)), [], 'all');
assert(spin_difference < 1e-10);

fprintf(['Extended ground-state validation passed; residual <0.002 T, ' ...
    '3x3 hessian min %.9g meV, dt difference %.3e.\n'], ...
    hessian3.minimum_eigenvalue_meV, spin_difference);
end
