function bonds = cri3_exchange_table(system)
%CRI3_EXCHANGE_TABLE Build the approved unique-bond scalar exchange table.
%   Each template is one physical bond per source cell. Periodic endpoints
%   wrap to the supercell; open-boundary bonds crossing the flake are omitted.

templates = exchange_templates();
n_cells = system.nx * system.ny;
max_bonds = n_cells * numel(templates.source_pair);

bonds.i = zeros(max_bonds, 1);
bonds.j = zeros(max_bonds, 1);
bonds.source_cell = zeros(max_bonds, 2);
bonds.target_cell_unwrapped = zeros(max_bonds, 2);
bonds.target_cell = zeros(max_bonds, 2);
bonds.image_offset = zeros(max_bonds, 2);
bonds.periodic_image = zeros(max_bonds, 2);
bonds.source_pair = cell(max_bonds, 1);
bonds.shell = cell(max_bonds, 1);
bonds.kind = cell(max_bonds, 1);
bonds.provenance = cell(max_bonds, 1);
bonds.J_source_meV = zeros(max_bonds, 1);
bonds.J_converted_meV = zeros(max_bonds, 1);
bonds.J_fit_meV = zeros(max_bonds, 1);
bonds.J_meV = zeros(max_bonds, 1);
bonds.J_J = zeros(max_bonds, 1);
bonds.distance_A = zeros(max_bonds, 1);

row = 0;
for cy = 0:(system.ny - 1)
    for cx = 0:(system.nx - 1)
        origin_cell = [cx, cy];
        for k = 1:numel(templates.source_pair)
            target_unwrapped = origin_cell + templates.offset(k, :);
            if strcmp(system.boundary, 'open') && ...
                    (target_unwrapped(1) < 0 || target_unwrapped(1) >= system.nx || ...
                     target_unwrapped(2) < 0 || target_unwrapped(2) >= system.ny)
                continue
            end

            if strcmp(system.boundary, 'periodic')
                target_cell = [mod(target_unwrapped(1), system.nx), ...
                    mod(target_unwrapped(2), system.ny)];
            else
                target_cell = target_unwrapped;
            end

            from_basis = label_to_basis(templates.from{k});
            to_basis = label_to_basis(templates.to{k});
            row = row + 1;
            bonds.i(row) = site_index(origin_cell, from_basis, system.nx);
            bonds.j(row) = site_index(target_cell, to_basis, system.nx);
            bonds.source_cell(row, :) = origin_cell;
            bonds.target_cell_unwrapped(row, :) = target_unwrapped;
            bonds.target_cell(row, :) = target_cell;
            bonds.image_offset(row, :) = templates.offset(k, :);
            bonds.periodic_image(row, :) = target_unwrapped - target_cell;
            bonds.source_pair{row} = templates.source_pair{k};
            bonds.shell{row} = templates.shell{k};
            bonds.kind{row} = templates.kind{k};
            bonds.provenance{row} = templates.provenance{k};
            bonds.J_source_meV(row) = templates.J_source_meV(k);
            bonds.J_converted_meV(row) = templates.J_converted_meV(k);
            bonds.J_fit_meV(row) = templates.J_fit_meV(k);
            bonds.J_meV(row) = templates.J_meV(k);
            bonds.J_J(row) = templates.J_meV(k) * 1.602176634e-22;

            source = system.source;
            delta_fractional = source.cr_fractional(to_basis, :) + ...
                [templates.offset(k, :), 0] - source.cr_fractional(from_basis, :);
            bonds.distance_A(row) = norm(delta_fractional * source.lattice_A);
        end
    end
end

fields = fieldnames(bonds);
for k = 1:numel(fields)
    value = bonds.(fields{k});
    bonds.(fields{k}) = value(1:row, :);
end
bonds.count = row;
bonds.boundary = system.boundary;
bonds.mapping = system.source.coordinate_convention;
bonds.pair_counting = 'unique physical bonds; no Hamiltonian factor 1/2';
bonds.units.J = 'joule per unique unit-spin bond';
bonds.units.distance = 'angstrom';
end

function index = site_index(cell_xy, basis, nx)
index = (cell_xy(2) * nx + cell_xy(1)) * 4 + basis;
end

function basis = label_to_basis(label)
basis = sscanf(label, 'Cr%d');
if isempty(basis) || basis < 1 || basis > 4
    error('CrI3:InvalidSourceLabel', 'Unknown source-site label: %s', label);
end
end

function templates = exchange_templates()
top = {
    '1-2',   'Cr4', 'Cr3',  0,  0, 'intra-nn',  -4.62;
    '2-3',   'Cr3', 'Cr4',  1,  0, 'intra-nn',  -4.20;
    '2-7',   'Cr3', 'Cr4',  0, -1, 'intra-nn',  -4.64;
    '1-3',   'Cr4', 'Cr4',  1,  0, 'intra-nnn', -0.68;
    '4-8',   'Cr3', 'Cr3', -1, -1, 'intra-nnn', -0.68;
    '1-7',   'Cr4', 'Cr4',  0, -1, 'intra-nnn', -0.62;
    '2-4',   'Cr3', 'Cr3',  1,  0, 'intra-nnn', -0.62;
    '3-7',   'Cr4', 'Cr4', -1, -1, 'intra-nnn', -0.62;
    '2-8',   'Cr3', 'Cr3',  0, -1, 'intra-nnn', -0.62};

mirror = top;
for k = 1:size(mirror, 1)
    mirror{k, 1} = ['mirror ' top{k, 1}];
    mirror{k, 2} = mirror_label(top{k, 2});
    mirror{k, 3} = mirror_label(top{k, 3});
end

inter = {
    '7-20',  'Cr4', 'Cr1',  0,  1, 'inter-nn',  -0.30;
    '1-19',  'Cr4', 'Cr2',  0,  0, 'inter-nn',  -0.29;
    '2-19',  'Cr3', 'Cr2',  0,  0, 'inter-nn',  -0.28;
    '2-20',  'Cr3', 'Cr1',  0,  0, 'inter-nn',  -0.29;
    '9-21',  'Cr4', 'Cr2',  0,  1, 'inter-nnn',  0.13;
    '9-26',  'Cr4', 'Cr1', -1,  0, 'inter-nnn',  0.18;
    '10-22', 'Cr3', 'Cr1',  0,  1, 'inter-nnn',  0.13;
    '10-29', 'Cr3', 'Cr2',  1,  0, 'inter-nnn',  0.19;
    '9-20',  'Cr4', 'Cr1', -1,  1, 'inter-3nn',  0.12;
    '9-23',  'Cr4', 'Cr2',  1,  1, 'inter-3nn',  0.06;
    '9-25',  'Cr4', 'Cr2', -1,  0, 'inter-3nn',  0.06;
    '9-28',  'Cr4', 'Cr1',  0,  0, 'inter-3nn',  0.01;
    '10-23', 'Cr3', 'Cr2',  1,  1, 'inter-3nn',  0.12;
    '10-24', 'Cr3', 'Cr1',  1,  1, 'inter-3nn',  0.06;
    '10-26', 'Cr3', 'Cr1', -1,  0, 'inter-3nn',  0.05;
    '10-33', 'Cr3', 'Cr2',  0, -1, 'inter-3nn',  0.01};

all_rows = [top; mirror; inter];
templates.source_pair = all_rows(:, 1);
templates.from = all_rows(:, 2);
templates.to = all_rows(:, 3);
templates.offset = cell2mat(all_rows(:, 4:5));
templates.shell = all_rows(:, 6);
templates.J_source_meV = cell2mat(all_rows(:, 7));
templates.kind = [repmat({'intralayer'}, 18, 1); repmat({'interlayer'}, 16, 1)];
templates.J_converted_meV = 2.25 * templates.J_source_meV;
templates.J_fit_meV = [zeros(18, 1); repmat(0.0075, 16, 1)];
templates.J_meV = templates.J_converted_meV + templates.J_fit_meV;
templates.provenance = [
    repmat({'Soenen Table S.V, S^2 conversion'}, 9, 1);
    repmat({'layer-symmetry mirror of Table S.V, S^2 conversion'}, 9, 1);
    repmat({'Soenen Table S.IV Jzz, S^2 conversion + approved fit'}, 16, 1)];
end

function output = mirror_label(input)
switch input
    case 'Cr4'
        output = 'Cr2';
    case 'Cr3'
        output = 'Cr1';
    otherwise
        error('CrI3:InvalidMirrorLabel', 'Cannot mirror source label %s.', input);
end
end
