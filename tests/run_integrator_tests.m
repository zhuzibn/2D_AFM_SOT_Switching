function run_integrator_tests()
%RUN_INTEGRATOR_TESTS Deterministic RK4, norm, energy, and time gates.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);

spins0 = [1, 0, 0];
layer = 1;
dynamics.gamma_rad_per_s_T = 1;
dynamics.alpha = 0;
torque.sigma = [0, 1, 0];
torque.dl_T = [0, 0];
torque.fl_T = [0, 0];
B = [0, 0, 1];
rhs = @(~, spins) cri3_llgs_rhs(spins, B, layer, dynamics, torque);

coarse = cri3_integrate_rk4(spins0, 0.203, 0.01, 7, rhs);
fine = cri3_integrate_rk4(spins0, 0.203, 0.005, 7, rhs);
exact = [cos(0.203), sin(0.203), 0];
coarse_error = norm(coarse.spins(:, :, end) - exact);
fine_error = norm(fine.spins(:, :, end) - exact);
assert(fine_error < coarse_error / 8);
assert(abs(coarse.time_s(1)) < eps);
assert(abs(coarse.time_s(end) - 0.203) < eps);
assert(size(coarse.spins, 3) == numel(coarse.time_s));
assert(coarse.max_norm_error < 5e-15);
assert(all(isfinite(coarse.spins(:))));

dynamics.alpha = 0.2;
rhs_damped = @(~, spins) cri3_llgs_rhs(spins, B, layer, dynamics, torque);
damped = cri3_integrate_rk4(spins0, 1.0, 0.005, 1, rhs_damped);
zeeman_energy = -squeeze(damped.spins(1, 3, :));
assert(all(diff(zeeman_energy) <= 5e-14));
assert(zeeman_energy(end) < zeeman_energy(1));

dynamics.alpha = 0;
rhs_undamped = @(~, spins) cri3_llgs_rhs(spins, B, layer, dynamics, torque);
undamped = cri3_integrate_rk4([sqrt(0.75), 0, 0.5], 1.0, 0.005, 1, rhs_undamped);
undamped_energy = -squeeze(undamped.spins(1, 3, :));
assert(max(abs(undamped_energy - undamped_energy(1))) < 1e-13);

fprintf(['Stage 5 RK4 tests passed; final-time error %.3e -> %.3e ' ...
    'under timestep refinement.\n'], coarse_error, fine_error);
end
