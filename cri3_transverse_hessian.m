function analysis = cri3_transverse_hessian(spins, system, bonds, parameters, epsilon)
%CRI3_TRANSVERSE_HESSIAN Numerical unit-sphere Hessian about a spin state.
%   Eigenvalues diagnose local transverse stability only. They do not prove
%   that the supplied state is the global ground state.

validateattributes(epsilon, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
spins = spins ./ vecnorm(spins, 2, 2);
n = system.n_sites;
e1 = zeros(n, 3);
e2 = zeros(n, 3);
for site = 1:n
    reference = [0, 0, 1];
    if abs(dot(spins(site, :), reference)) > 0.9
        reference = [1, 0, 0];
    end
    e1(site, :) = cross(spins(site, :), reference);
    e1(site, :) = e1(site, :) / norm(e1(site, :));
    e2(site, :) = cross(spins(site, :), e1(site, :));
end

n_coordinate = 2 * n;
hessian_J = zeros(n_coordinate);
origin = zeros(n_coordinate, 1);
energy0 = energy_at(origin);
for a = 1:n_coordinate
    qa = origin;
    qa(a) = epsilon;
    energy_plus = energy_at(qa);
    qa(a) = -epsilon;
    energy_minus = energy_at(qa);
    hessian_J(a, a) = (energy_plus - 2 * energy0 + energy_minus) / epsilon^2;
    for b = (a + 1):n_coordinate
        q = origin;
        q([a, b]) = epsilon;
        energy_pp = energy_at(q);
        q(b) = -epsilon;
        energy_pm = energy_at(q);
        q(a) = -epsilon;
        energy_mm = energy_at(q);
        q(b) = epsilon;
        energy_mp = energy_at(q);
        value = (energy_pp - energy_pm - energy_mp + energy_mm) / (4 * epsilon^2);
        hessian_J(a, b) = value;
        hessian_J(b, a) = value;
    end
end

analysis.hessian_J = hessian_J;
analysis.eigenvalues_J = sort(real(eig((hessian_J + hessian_J') / 2)));
analysis.eigenvalues_meV = analysis.eigenvalues_J / 1.602176634e-22;
analysis.minimum_eigenvalue_meV = min(analysis.eigenvalues_meV);
analysis.negative_mode_count = sum(analysis.eigenvalues_meV < -1e-6);
analysis.epsilon = epsilon;
analysis.interpretation = [
    'Local transverse finite-difference curvature; assess only together ' ...
    'with residual torque and epsilon convergence. Not a global proof.'];

    function total_J = energy_at(coordinates)
        trial = spins;
        for inner_site = 1:n
            trial(inner_site, :) = trial(inner_site, :) + ...
                coordinates(2 * inner_site - 1) * e1(inner_site, :) + ...
                coordinates(2 * inner_site) * e2(inner_site, :);
        end
        trial = trial ./ vecnorm(trial, 2, 2);
        trial_energy = cri3_energy_fields(trial, system, bonds, parameters);
        total_J = trial_energy.total_J;
    end
end
