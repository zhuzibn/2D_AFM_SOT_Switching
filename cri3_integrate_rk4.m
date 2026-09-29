function trajectory = cri3_integrate_rk4(spins0, t_final_s, dt_s, save_every, rhs_function)
%CRI3_INTEGRATE_RK4 Fixed-step RK4 with exact final-time and save alignment.
%   RHS_FUNCTION has signature dspindt=rhs_function(time_s, spins). Spins
%   are normalized only after each complete RK4 step.

validateattributes(spins0, {'numeric'}, {'2d', 'real', 'finite'});
if size(spins0, 2) ~= 3
    error('CrI3:InvalidSpinShape', 'Initial spins must be N-by-3.');
end
validateattributes(t_final_s, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(dt_s, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(save_every, {'numeric'}, {'scalar', 'integer', 'positive'});
if ~isa(rhs_function, 'function_handle')
    error('CrI3:InvalidRhsFunction', 'rhs_function must be a function handle.');
end

initial_norm = vecnorm(spins0, 2, 2);
if any(initial_norm == 0)
    error('CrI3:ZeroSpin', 'Initial spins must have nonzero length.');
end
spins = spins0 ./ initial_norm;

if t_final_s == 0
    n_steps = 0;
else
    n_steps = ceil(t_final_s / dt_s);
end
save_steps = unique([0:save_every:n_steps, n_steps]);
n_save = numel(save_steps);
saved_spins = zeros(size(spins, 1), 3, n_save);
saved_time = zeros(n_save, 1);
saved_spins(:, :, 1) = spins;
save_cursor = 2;
time_s = 0;

for step = 1:n_steps
    h = min(dt_s, t_final_s - time_s);
    k1 = rhs_function(time_s, spins);
    k2 = rhs_function(time_s + h / 2, spins + h * k1 / 2);
    k3 = rhs_function(time_s + h / 2, spins + h * k2 / 2);
    k4 = rhs_function(time_s + h, spins + h * k3);
    candidate = spins + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    candidate_norm = vecnorm(candidate, 2, 2);
    if any(~isfinite(candidate(:))) || any(candidate_norm == 0)
        error('CrI3:InvalidIntegratorState', ...
            'RK4 produced a nonfinite or zero-length spin at step %d.', step);
    end
    spins = candidate ./ candidate_norm;
    time_s = min(step * dt_s, t_final_s);

    if save_cursor <= n_save && step == save_steps(save_cursor)
        saved_spins(:, :, save_cursor) = spins;
        saved_time(save_cursor) = time_s;
        save_cursor = save_cursor + 1;
    end
end

if save_cursor ~= n_save + 1
    error('CrI3:SavedTimeMismatch', 'Saved spin count does not match saved times.');
end

trajectory.time_s = saved_time;
trajectory.spins = saved_spins;
trajectory.step_count = n_steps;
trajectory.save_steps = save_steps(:);
trajectory.dt_max_s = dt_s;
trajectory.max_norm_error = max(abs(vecnorm(saved_spins, 2, 2) - 1), [], 'all');
trajectory.solver = 'fixed-step RK4; normalize after complete step';
end
