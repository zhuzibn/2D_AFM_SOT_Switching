function observables = cri3_observables(spin_history, system)
%CRI3_OBSERVABLES Layer, net, Neel, and honeycomb averages.

if size(spin_history, 1) ~= system.n_sites || size(spin_history, 2) ~= 3
    error('CrI3:InvalidSpinHistory', 'Spin history must be N-by-3-by-time.');
end
n_time = size(spin_history, 3);
observables.layer_m = zeros(n_time, 3, 2);
observables.sublattice_m = zeros(n_time, 3, 2, 2);
for t = 1:n_time
    spins = spin_history(:, :, t);
    for layer = 1:2
        layer_mask = system.layer == layer;
        observables.layer_m(t, :, layer) = mean(spins(layer_mask, :), 1);
        for sublattice = 1:2
            mask = layer_mask & system.sublattice == sublattice;
            observables.sublattice_m(t, :, layer, sublattice) = mean(spins(mask, :), 1);
        end
    end
end
observables.M = (observables.layer_m(:, :, 1) + observables.layer_m(:, :, 2)) / 2;
observables.L = (observables.layer_m(:, :, 1) - observables.layer_m(:, :, 2)) / 2;
observables.max_norm_error = max(abs(vecnorm(spin_history, 2, 2) - 1), [], 'all');
observables.definitions = ['M=(m_bottom+m_top)/2; ' ...
    'L=(m_bottom-m_top)/2; sublattice 1=A, 2=B'];
end
