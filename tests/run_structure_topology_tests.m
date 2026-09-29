function run_structure_topology_tests()
%RUN_STRUCTURE_TOPOLOGY_TESTS Deterministic gates for stages 1 and 2.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);

source = cri3_source_structure();
assert(numel(source.label) == 16);
assert(sum(strcmp(source.species, 'Cr')) == 4);
assert(sum(strcmp(source.species, 'I')) == 12);
assert(isequal(source.cr_label, {'Cr1'; 'Cr2'; 'Cr3'; 'Cr4'}));

z_A = source.cr_cartesian_A(:, 3);
plane_gap_A = mean(z_A(3:4)) - mean(z_A(1:2));
assert(abs(plane_gap_A - 6.694495) < 1e-5);

bottom_delta = source.cr_fractional(2, :) - source.cr_fractional(1, :);
candidate_distances = zeros(9, 1);
row = 0;
for dx = -1:1
    for dy = -1:1
        row = row + 1;
        candidate_distances(row) = norm((bottom_delta + [dx, dy, 0]) * ...
            source.lattice_A);
    end
end
candidate_distances(candidate_distances < 1e-12) = inf;
assert(abs(min(candidate_distances) - 4.046038) < 1e-5);

periodic = cri3_build_supercell(2, 2, 'periodic');
assert(periodic.n_sites == 16);
assert(all(periodic.layer(1:4) == [1; 1; 2; 2]));
assert(isequal(periodic.source_label(1:4), {'Cr1'; 'Cr2'; 'Cr3'; 'Cr4'}));

bonds = cri3_exchange_table(periodic);
assert(bonds.count == 34 * 4);
assert(sum(strcmp(bonds.kind, 'intralayer')) == 18 * 4);
assert(sum(strcmp(bonds.kind, 'interlayer')) == 16 * 4);
assert(all(bonds.i >= 1 & bonds.i <= periodic.n_sites));
assert(all(bonds.j >= 1 & bonds.j <= periodic.n_sites));

cell_zero = all(bonds.source_cell == 0, 2);
assert(sum(cell_zero & strcmp(bonds.shell, 'inter-nn')) == 4);
assert(sum(cell_zero & strcmp(bonds.shell, 'inter-nnn')) == 4);
assert(sum(cell_zero & strcmp(bonds.shell, 'inter-3nn')) == 8);

nn_distance = bonds.distance_A(cell_zero & strcmp(bonds.shell, 'inter-nn'));
nnn_distance = bonds.distance_A(cell_zero & strcmp(bonds.shell, 'inter-nnn'));
third_distance = bonds.distance_A(cell_zero & strcmp(bonds.shell, 'inter-3nn'));
assert(all(nn_distance >= 7.043 & nn_distance <= 7.117));
assert(all(nnn_distance >= 8.121 & nnn_distance <= 8.249));
assert(all(third_distance >= 9.017 & third_distance <= 9.189));

inter_cell = cell_zero & strcmp(bonds.kind, 'interlayer');
assert(abs(sum(bonds.J_converted_meV(inter_cell)) - (-0.09)) < 1e-12);
assert(abs(sum(bonds.J_fit_meV(inter_cell)) - 0.12) < 1e-12);
assert(abs(sum(bonds.J_meV(inter_cell)) - 0.03) < 1e-12);
assert(all(bonds.J_meV(strcmp(bonds.kind, 'intralayer')) < 0));

template_rows = find(cell_zero);
template_keys = cell(numel(template_rows), 1);
for k = 1:numel(template_rows)
    b = template_rows(k);
    from_basis = periodic.source_cr_index(bonds.i(b));
    to_basis = periodic.source_cr_index(bonds.j(b));
    template_keys{k} = canonical_bond_key(from_basis, to_basis, ...
        bonds.image_offset(b, :));
end
assert(numel(unique(template_keys)) == numel(template_keys));

adjacency = cri3_bond_adjacency(bonds, periodic.n_sites);
assert(numel(adjacency.from) == 2 * bonds.count);
assert(all(adjacency.count_by_site == 17));
assert(isequal(adjacency.from(1:bonds.count), bonds.i));
assert(isequal(adjacency.to((bonds.count + 1):end), bonds.i));
intralayer_mask = strcmp(bonds.kind, 'intralayer');
intralayer_bonds = subset_bonds(bonds, intralayer_mask);
intralayer_adjacency = cri3_bond_adjacency(intralayer_bonds, periodic.n_sites);
for layer = 1:2
    layer_sites = find(periodic.layer == layer);
    assert(is_connected(layer_sites, intralayer_adjacency));
end

open_system = cri3_build_supercell(5, 5, 'open');
open_bonds = cri3_exchange_table(open_system);
assert(all(open_bonds.periodic_image(:) == 0));
assert(all(open_bonds.target_cell(:) >= 0));
assert(all(open_bonds.target_cell(:, 1) < open_system.nx));
assert(all(open_bonds.target_cell(:, 2) < open_system.ny));
open_adjacency = cri3_bond_adjacency(open_bonds, open_system.n_sites);
central = find(open_system.cell(:, 1) == 2 & open_system.cell(:, 2) == 2);
assert(all(open_adjacency.count_by_site(central) == 17));
edge = find(open_system.cell(:, 1) == 0 | open_system.cell(:, 1) == 4 | ...
    open_system.cell(:, 2) == 0 | open_system.cell(:, 2) == 4);
assert(all(open_adjacency.count_by_site(edge) < 17));

fprintf('Stage 1-2 structure/topology tests passed (%d periodic bonds).\n', ...
    bonds.count);
end

function key = canonical_bond_key(from_basis, to_basis, offset)
forward = sprintf('%d>%d@%d,%d', from_basis, to_basis, offset(1), offset(2));
reverse = sprintf('%d>%d@%d,%d', to_basis, from_basis, -offset(1), -offset(2));
ordered = sort({forward, reverse});
key = ordered{1};
end

function selected = subset_bonds(bonds, mask)
selected.i = bonds.i(mask);
selected.j = bonds.j(mask);
selected.J_J = bonds.J_J(mask);
selected.count = sum(mask);
end

function connected = is_connected(sites, adjacency)
visited = false(max(adjacency.from), 1);
visited(sites(1)) = true;
changed = true;
while changed
    changed = false;
    for edge = 1:numel(adjacency.from)
        if visited(adjacency.from(edge)) && ~visited(adjacency.to(edge))
            visited(adjacency.to(edge)) = true;
            changed = true;
        end
    end
end
connected = all(visited(sites));
end
