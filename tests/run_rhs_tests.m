function run_rhs_tests()
%RUN_RHS_TESTS Deterministic LLGS equation and torque-sign gates.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
gamma = 1.76085963023e11;
rng(41, 'twister');
spins = randn(8, 3);
spins = spins ./ vecnorm(spins, 2, 2);
B = randn(8, 3) * 0.4;
layer = [ones(4, 1); 2 * ones(4, 1)];
dynamics.alpha = linspace(0.001, 0.02, 8)';
dynamics.gamma_rad_per_s_T = gamma;
torque.sigma = [0, 1, 0];
torque.dl_T = [0.012, -0.007];
torque.fl_T = [-0.003, 0.005];

ds = cri3_llgs_rhs(spins, B, layer, dynamics, torque);
tangent_error = max(abs(sum(spins .* ds, 2))) / max(vecnorm(ds, 2, 2));
assert(tangent_error < 1e-14);

alpha = dynamics.alpha;
dl = torque.dl_T(layer)';
fl = torque.fl_T(layer)';
sigma = repmat(torque.sigma, size(spins, 1), 1);
residual = ds + gamma * cross(spins, B, 2) ...
    - alpha(:, ones(1, 3)) .* cross(spins, ds, 2) ...
    + gamma * dl(:, ones(1, 3)) .* cross(spins, cross(spins, sigma, 2), 2) ...
    - gamma * fl(:, ones(1, 3)) .* cross(spins, sigma, 2);
residual_scale = max([vecnorm(ds, 2, 2); gamma * vecnorm(B, 2, 2)]);
assert(max(vecnorm(residual, 2, 2)) / residual_scale < 2e-15);

single_spin = [1, 0, 0];
single_B = [0, 0, 1];
single_layer = 1;
zero_damping.alpha = 0;
zero_damping.gamma_rad_per_s_T = gamma;
no_torque.sigma = [0, 1, 0];
no_torque.dl_T = [0, 0];
no_torque.fl_T = [0, 0];
precession = cri3_llgs_rhs(single_spin, single_B, single_layer, ...
    zero_damping, no_torque);
assert(norm(precession - [0, gamma, 0]) / gamma < 1e-15);

spin_z = [0, 0, 1];
zero_B = [0, 0, 0];
dl_only = no_torque;
dl_only.dl_T = [0.01, -0.02];
dl_bottom = cri3_llgs_rhs(spin_z, zero_B, 1, zero_damping, dl_only);
dl_top = cri3_llgs_rhs(spin_z, zero_B, 2, zero_damping, dl_only);
assert(norm(dl_bottom - [0, gamma * 0.01, 0]) / gamma < 1e-15);
assert(norm(dl_top - [0, -gamma * 0.02, 0]) / gamma < 1e-15);

fl_only = no_torque;
fl_only.fl_T = [0.01, -0.02];
fl_bottom = cri3_llgs_rhs(spin_z, zero_B, 1, zero_damping, fl_only);
fl_top = cri3_llgs_rhs(spin_z, zero_B, 2, zero_damping, fl_only);
assert(norm(fl_bottom - [-gamma * 0.01, 0, 0]) / gamma < 1e-15);
assert(norm(fl_top - [gamma * 0.02, 0, 0]) / gamma < 1e-15);

fprintf('Stage 4 LLGS tests passed; tangent error %.3e.\n', tangent_error);
end
