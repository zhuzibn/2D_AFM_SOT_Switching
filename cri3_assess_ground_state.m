function assessment = cri3_assess_ground_state(config, options)
%CRI3_ASSESS_GROUND_STATE Reproducible multistart local-minimum assessment.
%   The relaxation alpha is a numerical minimization choice, not a claimed
%   material damping value. The returned result never claims globality.

if nargin < 2
    options = struct();
end
options = apply_defaults(options);
system = cri3_build_supercell(config.nx, config.ny, 'periodic');
bonds = cri3_exchange_table(system);
parameters.mu_J_per_T = config.mu_J_per_T;
parameters.K_J = config.K_J;
parameters.B_external_T = [0, 0, 0];

dynamics.alpha = options.relaxation_alpha;
dynamics.gamma_rad_per_s_T = config.gamma_rad_per_s_T;
torque.sigma = [0, 1, 0];
torque.dl_T = [0, 0];
torque.fl_T = [0, 0];
model.system = system;
model.bonds = bonds;
model.parameters = parameters;
model.dynamics = dynamics;
model.torque = torque;
model.torque_function = [];

starts = {
    'collinear-afm', cri3_initial_state(system, 'afm', 0);
    'collinear-fm', cri3_initial_state(system, 'fm', 0);
    'canted-afm', cri3_initial_state(system, 'afm', options.canting)};
for seed = options.random_seeds(:)'
    random_spins = random_unit_spins(system.n_sites, seed);
    starts(end + 1, :) = {sprintf('random-seed-%d', seed), random_spins}; %#ok<AGROW>
end

n_start = size(starts, 1);
candidates = repmat(struct(), n_start, 1);
for k = 1:n_start
    rhs = @(time_s, spins) cri3_model_rhs(time_s, spins, model);
    trajectory = cri3_integrate_rk4(starts{k, 2}, options.t_final_s, ...
        options.dt_s, options.save_every, rhs);
    n_time = numel(trajectory.time_s);
    energy_J = zeros(n_time, 1);
    for t = 1:n_time
        snapshot = cri3_energy_fields(trajectory.spins(:, :, t), ...
            system, bonds, parameters);
        energy_J(t) = snapshot.total_J;
    end
    final_spins = trajectory.spins(:, :, end);
    [~, final_fields] = cri3_energy_fields(final_spins, system, bonds, parameters);
    residual = vecnorm(cross(final_spins, final_fields.total_T, 2), 2, 2);
    order = cri3_observables(reshape(final_spins, system.n_sites, 3, 1), system);

    candidates(k).label = starts{k, 1};
    candidates(k).initial_energy_J = energy_J(1);
    candidates(k).final_energy_J = energy_J(end);
    candidates(k).energy_trace_J = energy_J;
    candidates(k).time_s = trajectory.time_s;
    candidates(k).maximum_energy_increase_J = max([0; diff(energy_J)]);
    candidates(k).max_residual_torque_T = max(residual);
    candidates(k).final_spins = final_spins;
    candidates(k).m_bottom = reshape(order.layer_m(1, :, 1), 1, 3);
    candidates(k).m_top = reshape(order.layer_m(1, :, 2), 1, 3);
    candidates(k).M = order.M(1, :);
    candidates(k).L = order.L(1, :);
    candidates(k).layer_dot = dot(candidates(k).m_bottom, candidates(k).m_top);
end

[~, lowest] = min([candidates.final_energy_J]);
hessian = cri3_transverse_hessian(candidates(lowest).final_spins, system, ...
    bonds, parameters, options.hessian_epsilon);

assessment.model_id = config.model_id;
assessment.system_size = [config.nx, config.ny];
assessment.boundary = 'periodic';
assessment.options = options;
assessment.candidates = candidates;
assessment.lowest_candidate_index = lowest;
assessment.lowest_candidate_label = candidates(lowest).label;
assessment.lowest_transverse_hessian = hessian;
assessment.warning = [
    'Multistart damped LLGS plus a local transverse Hessian identifies ' ...
    'candidate local minima only; it is not proof of the global ground state.'];
assessment.units.energy = 'J';
assessment.units.residual_torque = 'T';
end

function options = apply_defaults(options)
defaults.relaxation_alpha = 1;
defaults.dt_s = 1e-15;
defaults.t_final_s = 5e-13;
defaults.save_every = 10;
defaults.canting = 0.05;
defaults.random_seeds = [17, 29];
defaults.hessian_epsilon = 1e-4;
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(options, names{k})
        options.(names{k}) = defaults.(names{k});
    end
end
validateattributes(options.relaxation_alpha, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(options.random_seeds, {'numeric'}, {'vector', 'integer', 'nonnegative'});
end

function spins = random_unit_spins(n_sites, seed)
previous = rng;
rng(seed, 'twister');
spins = randn(n_sites, 3);
spins = spins ./ vecnorm(spins, 2, 2);
rng(previous);
end
