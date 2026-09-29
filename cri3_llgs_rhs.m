function [dspindt, terms] = cri3_llgs_rhs(spins, B_effective_T, layer, dynamics, torque)
%CRI3_LLGS_RHS Explicit fixed-length LLGS right-hand side.
%   Implements Eq. eq:ll-with-flt from spintronic_methods.tex. Fields and
%   signed DL/FL amplitudes are in tesla; gamma is in rad/(s T).

validateattributes(spins, {'numeric'}, {'2d', 'real', 'finite'});
if size(spins, 2) ~= 3 || ~isequal(size(B_effective_T), size(spins))
    error('CrI3:InvalidRhsShape', 'Spins and effective field must be N-by-3.');
end
n = size(spins, 1);
if ~isvector(layer) || numel(layer) ~= n || any(~ismember(layer(:), [1, 2]))
    error('CrI3:InvalidLayerMap', 'Layer must contain one value (1 or 2) per site.');
end

alpha = expand_site_scalar(dynamics.alpha, n, 'alpha');
gamma = expand_site_scalar(dynamics.gamma_rad_per_s_T, n, ...
    'gamma_rad_per_s_T');
if any(alpha < 0) || any(gamma <= 0)
    error('CrI3:InvalidDynamicsParameter', ...
        'Alpha must be nonnegative and gamma must be positive.');
end

sigma = torque.sigma;
validateattributes(sigma, {'numeric'}, {'size', [1, 3], 'real', 'finite'});
if abs(norm(sigma) - 1) > 1e-12
    error('CrI3:InvalidPolarization', 'Torque polarization sigma must be a unit vector.');
end
sigma = repmat(sigma, n, 1);
dl_T = expand_layer_value(torque.dl_T, layer(:), 'dl_T');
fl_T = expand_layer_value(torque.fl_T, layer(:), 'fl_T');

s_cross_B = cross(spins, B_effective_T, 2);
s_cross_s_cross_B = cross(spins, s_cross_B, 2);
s_cross_sigma = cross(spins, sigma, 2);
s_cross_s_cross_sigma = cross(spins, s_cross_sigma, 2);

gamma3 = gamma(:, ones(1, 3));
alpha3 = alpha(:, ones(1, 3));
dl3 = dl_T(:, ones(1, 3));
fl3 = fl_T(:, ones(1, 3));
denominator = (1 + alpha.^2);
denominator3 = denominator(:, ones(1, 3));

terms.precession = -gamma3 .* s_cross_B;
terms.damping = -alpha3 .* gamma3 .* s_cross_s_cross_B;
terms.dl_cross = alpha3 .* gamma3 .* dl3 .* s_cross_sigma;
terms.dl_direct = -gamma3 .* dl3 .* s_cross_s_cross_sigma;
terms.fl_cross = alpha3 .* gamma3 .* fl3 .* s_cross_s_cross_sigma;
terms.fl_direct = gamma3 .* fl3 .* s_cross_sigma;

dspindt = (terms.precession + terms.damping + terms.dl_cross + ...
    terms.dl_direct + terms.fl_cross + terms.fl_direct) ./ denominator3;
terms.units = '1/s';
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

function output = expand_layer_value(value, layer, name)
validateattributes(value, {'numeric'}, {'real', 'finite', 'nonempty'});
if isscalar(value)
    layer_value = [value, value];
elseif isvector(value) && numel(value) == 2
    layer_value = double(value(:))';
else
    error('CrI3:InvalidTorqueShape', '%s must be scalar or [bottom, top].', name);
end
output = layer_value(layer(:));
output = output(:);
end
