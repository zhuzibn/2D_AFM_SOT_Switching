function spins = cri3_initial_state(system, mode, canting)
%CRI3_INITIAL_STATE Construct deterministic layerwise FM initial states.
%   CANTING is an angle in radians toward +x. AFM means bottom +z, top -z.

mode = validatestring(mode, {'afm', 'fm'});
validateattributes(canting, {'numeric'}, {'scalar', 'real', 'finite'});
if abs(canting) >= pi / 2
    error('CrI3:InvalidCanting', 'Initial canting magnitude must be below pi/2.');
end

spins = zeros(system.n_sites, 3);
spins(:, 1) = sin(canting);
spins(:, 3) = cos(canting);
if strcmp(mode, 'afm')
    spins(system.layer == 2, 3) = -cos(canting);
end
spins = spins ./ vecnorm(spins, 2, 2);
end
