function system = cri3_build_supercell(nx, ny, boundary)
%CRI3_BUILD_SUPERCELL Replicate the four Cr spin sites of the BiDB cell.
%   BOUNDARY is 'periodic' or 'open'. Coordinates are never wrapped; the
%   boundary mode controls only which interaction bonds are retained.

validateattributes(nx, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(ny, {'numeric'}, {'scalar', 'integer', 'positive'});
boundary = validatestring(boundary, {'periodic', 'open'});

source = cri3_source_structure();
n_basis = numel(source.cr_source_index);
n_sites = nx * ny * n_basis;

system.nx = double(nx);
system.ny = double(ny);
system.boundary = boundary;
system.source = source;
system.n_sites = n_sites;
system.cell = zeros(n_sites, 2);
system.source_cr_index = zeros(n_sites, 1);
system.source_atom_index = zeros(n_sites, 1);
system.source_label = cell(n_sites, 1);
system.layer = zeros(n_sites, 1);
system.sublattice = zeros(n_sites, 1);
system.fractional = zeros(n_sites, 3);
system.cartesian_A = zeros(n_sites, 3);

row = 0;
for cy = 0:(ny - 1)
    for cx = 0:(nx - 1)
        for basis = 1:n_basis
            row = row + 1;
            global_fractional = source.cr_fractional(basis, :) + [cx, cy, 0];
            system.cell(row, :) = [cx, cy];
            system.source_cr_index(row) = basis;
            system.source_atom_index(row) = source.cr_source_index(basis);
            system.source_label{row} = source.cr_label{basis};
            system.layer(row) = source.cr_layer(basis);
            system.sublattice(row) = source.cr_sublattice(basis);
            system.fractional(row, :) = global_fractional;
            system.cartesian_A(row, :) = global_fractional * source.lattice_A;
        end
    end
end

system.site_order = '(cell y, cell x, source Cr1, Cr2, Cr3, Cr4)';
system.units.length = 'angstrom';
system.units.spin = 'unit vector';
end
