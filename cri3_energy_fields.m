function [energy, fields] = cri3_energy_fields(spins, system, bonds, parameters)
%CRI3_ENERGY_FIELDS Evaluate scalar energy and conservative fields.
%   SPINS is N-by-3. Energy coefficients are joules and returned fields are
%   flux densities in tesla using B_eff=-(1/mu)*dH/ds.

validateattributes(spins, {'numeric'}, {'2d', 'real', 'finite'});
if size(spins, 1) ~= system.n_sites || size(spins, 2) ~= 3
    error('CrI3:InvalidSpinShape', 'Spins must have shape system.n_sites-by-3.');
end

n = system.n_sites;
mu = expand_site_scalar(parameters.mu_J_per_T, n, 'mu_J_per_T');
K = expand_site_scalar(parameters.K_J, n, 'K_J');
B_external = expand_site_vector(parameters.B_external_T, n, 'B_external_T');
if any(mu <= 0)
    error('CrI3:InvalidMoment', 'Every magnetic moment must be positive.');
end

exchange_gradient = zeros(n, 3);
exchange_energy = 0;
for b = 1:bonds.count
    i = bonds.i(b);
    j = bonds.j(b);
    J = bonds.J_J(b);
    exchange_energy = exchange_energy + J * dot(spins(i, :), spins(j, :));
    exchange_gradient(i, :) = exchange_gradient(i, :) + J * spins(j, :);
    exchange_gradient(j, :) = exchange_gradient(j, :) + J * spins(i, :);
end

anisotropy_energy = -sum(K .* spins(:, 3).^2);
zeeman_energy = -sum(mu .* sum(spins .* B_external, 2));

fields.exchange_T = -exchange_gradient ./ mu;
fields.anisotropy_T = zeros(n, 3);
fields.anisotropy_T(:, 3) = 2 * K .* spins(:, 3) ./ mu;
fields.external_T = B_external;
fields.total_T = fields.exchange_T + fields.anisotropy_T + fields.external_T;

energy.exchange_J = exchange_energy;
energy.anisotropy_J = anisotropy_energy;
energy.zeeman_J = zeeman_energy;
energy.total_J = exchange_energy + anisotropy_energy + zeeman_energy;
energy.units = 'J';
fields.units = 'T';
end

function output = expand_site_scalar(value, n, name)
validateattributes(value, {'numeric'}, {'real', 'finite', 'nonempty'});
if isscalar(value)
    output = repmat(double(value), n, 1);
elseif isvector(value) && numel(value) == n
    output = double(value(:));
else
    error('CrI3:InvalidParameterShape', '%s must be scalar or N-by-1.', name);
end
end

function output = expand_site_vector(value, n, name)
validateattributes(value, {'numeric'}, {'2d', 'real', 'finite', 'nonempty'});
if isequal(size(value), [1, 3])
    output = repmat(double(value), n, 1);
elseif isequal(size(value), [n, 3])
    output = double(value);
else
    error('CrI3:InvalidParameterShape', '%s must be 1-by-3 or N-by-3.', name);
end
end
