function run_open_sot_tests()
%RUN_OPEN_SOT_TESTS Deterministic short-time open-flake and size check.
%   This verifies boundary-aware dynamics, not a physical switching claim.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
config = cri3_default_config();
config.boundary = 'open';
config.alpha = 0.001;
config.t_final_s = 2e-14;
config.dt_s = 1e-15;
config.initial_canting = 1e-3;
config.torque.dl_T = [0.01, -0.01];
config.torque.fl_T = [0.002, -0.002];

config.nx = 5;
config.ny = 5;
small = cri3_run(config);
config.nx = 7;
config.ny = 7;
medium = cri3_run(config);
config.nx = 9;
config.ny = 9;
large = cri3_run(config);

small_center = find(all(small.system.cell == [2, 2], 2));
medium_center = find(all(medium.system.cell == [3, 3], 2));
large_center = find(all(large.system.cell == [4, 4], 2));
small_medium_difference = max(abs(small.spins(small_center, :, :) - ...
    medium.spins(medium_center, :, :)), [], 'all');
medium_large_difference = max(abs(medium.spins(medium_center, :, :) - ...
    large.spins(large_center, :, :)), [], 'all');
assert(small_medium_difference < 2e-7);
assert(medium_large_difference < 5e-9);
assert(medium_large_difference < small_medium_difference / 20);
assert(all(small.bonds.periodic_image(:) == 0));
assert(all(medium.bonds.periodic_image(:) == 0));
assert(all(large.bonds.periodic_image(:) == 0));
assert(small.observables.max_norm_error < 5e-15);
assert(medium.observables.max_norm_error < 5e-15);
assert(large.observables.max_norm_error < 5e-15);

fprintf(['Stage 7 open-SOT test passed; central-history differences ' ...
    '5x5/7x7 %.3e and 7x7/9x9 %.3e.\n'], ...
    small_medium_difference, medium_large_difference);
end
