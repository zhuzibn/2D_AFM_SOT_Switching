function run_interface_tests()
%RUN_INTERFACE_TESTS End-to-end configuration, provenance, and output gate.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
config = cri3_default_config();
config.nx = 1;
config.ny = 1;
config.alpha = 0.001; % source stability value, explicit in this test
config.t_final_s = 5e-15;
config.initial_canting = 1e-3;
config.dt_s = 1e-15;
config.save_every = 2;
result = cri3_run(config);

assert(isequal(size(result.spins), [4, 3, 4]));
assert(isequal(result.time_s, [0; 2e-15; 4e-15; 5e-15]));
assert(result.observables.max_norm_error < 5e-15);
assert(numel(result.bond_sha256) == 64);
assert(result.bonds.count == 34);
assert(strcmp(result.system.source.uid, '2CrI3-1-2-1_0_0_1-0_-0.69'));
assert(numel(result.energy.total_J) == numel(result.time_s));
assert(all(isfinite(result.energy.total_J)));
assert(isequal(size(result.observables.layer_m), [4, 3, 2]));
assert(isequal(size(result.observables.sublattice_m), [4, 3, 2, 2]));

missing_alpha = cri3_default_config();
missing_alpha.t_final_s = 1e-15;
assert_error(@() cri3_run(missing_alpha), 'CrI3:MissingDamping');
missing_time = cri3_default_config();
missing_time.alpha = 0.001;
assert_error(@() cri3_run(missing_time), 'CrI3:MissingFinalTime');

fprintf('Stage 7 interface test passed; bond SHA-256 %s.\n', result.bond_sha256);
end

function assert_error(action, expected_identifier)
try
    action();
catch exception
    assert(strcmp(exception.identifier, expected_identifier));
    return
end
error('CrI3:ExpectedErrorMissing', 'Expected error %s was not raised.', ...
    expected_identifier);
end
