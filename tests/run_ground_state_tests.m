function run_ground_state_tests()
%RUN_GROUND_STATE_TESTS Smoke gate for reproducible local-minimum analysis.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
config = cri3_default_config();
config.nx = 1;
config.ny = 1;
options.t_final_s = 2e-14;
options.dt_s = 1e-15;
options.save_every = 5;
options.random_seeds = 17;
options.hessian_epsilon = 2e-4;
assessment = cri3_assess_ground_state(config, options);

assert(numel(assessment.candidates) == 4);
assert(assessment.lowest_candidate_index >= 1);
assert(all(isfinite([assessment.candidates.final_energy_J])));
assert(all(isfinite([assessment.candidates.max_residual_torque_T])));
assert(numel(assessment.lowest_transverse_hessian.eigenvalues_meV) == 8);
assert(contains(assessment.warning, 'not proof'));

fprintf('Stage 6 ground-state assessment smoke test passed; lowest start %s.\n', ...
    assessment.lowest_candidate_label);
end
