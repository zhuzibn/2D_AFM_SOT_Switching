function cri3_result = cri3_run(config)
%CRI3_RUN Run the approved fitted-scalar bilayer CrI3 LLGS model.

validate_run_config(config);
system = cri3_build_supercell(config.nx, config.ny, config.boundary);
bonds = cri3_exchange_table(system);

parameters.mu_J_per_T = config.mu_J_per_T;
parameters.K_J = config.K_J;
parameters.B_external_T = config.B_external_T;
dynamics.alpha = config.alpha;
dynamics.gamma_rad_per_s_T = config.gamma_rad_per_s_T;
model.system = system;
model.bonds = bonds;
model.parameters = parameters;
model.dynamics = dynamics;
model.torque = config.torque;
model.torque_function = config.torque_function;

if isempty(config.initial_spins)
    spins0 = cri3_initial_state(system, config.initial_mode, config.initial_canting);
else
    spins0 = config.initial_spins;
end
rhs = @(time_s, spins) cri3_model_rhs(time_s, spins, model);
trajectory = cri3_integrate_rk4(spins0, config.t_final_s, config.dt_s, ...
    config.save_every, rhs);

n_time = numel(trajectory.time_s);
energy.exchange_J = zeros(n_time, 1);
energy.anisotropy_J = zeros(n_time, 1);
energy.zeeman_J = zeros(n_time, 1);
energy.total_J = zeros(n_time, 1);
for t = 1:n_time
    snapshot = cri3_energy_fields(trajectory.spins(:, :, t), system, bonds, parameters);
    energy.exchange_J(t) = snapshot.exchange_J;
    energy.anisotropy_J(t) = snapshot.anisotropy_J;
    energy.zeeman_J(t) = snapshot.zeeman_J;
    energy.total_J(t) = snapshot.total_J;
end
energy.units = 'J';

cri3_result.schema_version = '1.0';
cri3_result.model_id = config.model_id;
cri3_result.config = config;
cri3_result.system = system;
cri3_result.bonds = bonds;
cri3_result.bond_sha256 = cri3_bond_hash(bonds);
cri3_result.time_s = trajectory.time_s;
cri3_result.spins = trajectory.spins;
cri3_result.energy = energy;
cri3_result.observables = cri3_observables(trajectory.spins, system);
cri3_result.solver = rmfield(trajectory, 'spins');
cri3_result.provenance.created = char(datetime('now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
cri3_result.provenance.source_uid = system.source.uid;
cri3_result.provenance.source_cif_sha256 = system.source.cif_sha256;
cri3_result.provenance.mapping = system.source.coordinate_convention;
cri3_result.provenance.code_version = 'unversioned workspace (no Git metadata available)';
cri3_result.provenance.scope = system.source.fidelity_label;

if ~isempty(config.output_file)
    if isfile(config.output_file) && ~config.overwrite
        error('CrI3:OutputExists', ...
            'Output exists; choose another path or set overwrite=true: %s', ...
            config.output_file);
    end
    parent = fileparts(config.output_file);
    if ~isempty(parent) && ~isfolder(parent)
        error('CrI3:MissingOutputDirectory', ...
            'Output directory does not exist: %s', parent);
    end
    save(config.output_file, 'cri3_result', '-v7.3');
end
end

function validate_run_config(config)
required = {'model_id', 'nx', 'ny', 'boundary', 'mu_J_per_T', 'K_J', ...
    'B_external_T', 'gamma_rad_per_s_T', 'alpha', 'torque', ...
    'torque_function', 'dt_s', 't_final_s', 'save_every', 'initial_mode', ...
    'initial_canting', 'initial_spins', 'output_file', 'overwrite'};
for k = 1:numel(required)
    if ~isfield(config, required{k})
        error('CrI3:MissingConfigField', 'Missing configuration field: %s', required{k});
    end
end
if isempty(config.alpha)
    error('CrI3:MissingDamping', ...
        ['Set config.alpha explicitly. The source alpha=0.001 applies to ' ...
        'Soenen''s T=0 stability run and is not assumed universally.']);
end
if isempty(config.t_final_s)
    error('CrI3:MissingFinalTime', 'Set config.t_final_s explicitly.');
end
if ~isempty(config.torque_function) && ~isa(config.torque_function, 'function_handle')
    error('CrI3:InvalidTorqueFunction', 'torque_function must be empty or a function handle.');
end
end
