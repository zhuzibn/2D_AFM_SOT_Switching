function source = cri3_source_structure()
%CRI3_SOURCE_STRUCTURE Immutable BiDB source cell for AB-prime bilayer CrI3.
%   Lengths are in angstrom. Fractional coordinates and labels are copied
%   from the cited BiDB CIF without idealization or geometry relaxation.

source.uid = '2CrI3-1-2-1_0_0_1-0_-0.69';
source.url = ['https://bidb.fysik.dtu.dk/material/' source.uid];
source.access_date = '2026-09-27';
source.cif_sha256 = ...
    'd3386fbda276c27327d94820543824f29ed46571f5eca702a0b5226847e824cf';
source.fidelity_label = 'BiDB geometry + fitted Soenen scalar reduction';
source.lattice_A = [
     7.007944101379477,   0,                  0;
    -3.5039720506897396,  6.069057620095938,  0;
     0,                   0,                 27.197312819559507];

source.species = {
    'Cr'; 'Cr'; 'I'; 'I'; 'I'; 'I'; 'I'; 'I';
    'Cr'; 'Cr'; 'I'; 'I'; 'I'; 'I'; 'I'; 'I'};
source.label = {
    'Cr1'; 'Cr2'; 'I1'; 'I2'; 'I3'; 'I4'; 'I5'; 'I6';
    'Cr3'; 'Cr4'; 'I7'; 'I8'; 'I9'; 'I10'; 'I11'; 'I12'};
source.fractional = [
    0.6666666668402268, 0.3333333338772982, 0.3320831285767547;
    0.3333333330534472, 0.6666666661068942, 0.33208314549019036;
    0.6400805131903461, 0.9999977871200474, 0.27447123690217945;
    0.3599172743756347, 0.35991948614346325, 0.27447123690217945;
    2.2124182117120484e-06, 0.6400827267048738, 0.27447123690217945;
    0.3599174387432292, 0.9999980853541579, 0.38969498201225816;
    0.6400806465474529, 0.6400825602869448, 0.38969498201225816;
    1.914693510525997e-06, 0.3599193543272819, 0.38969498201225816;
    0.6663373607355287, 0.6455589525179144, 0.578228585828896;
    0.333004026948749, 0.9788922847475102, 0.5782286027423317;
    0.6397512069951293, 0.3122234057764714, 0.5206166937866373;
    0.35958796684398464, 0.6721451047840795, 0.5206166937866373;
    0.9996729061166689, 0.9523083453454899, 0.5206166937866373;
    0.3595881325480123, 0.31222370401058175, 0.635840438896716;
    0.6397513390158028, 0.9523081789275609, 0.635840438896716;
    0.9996726083919676, 0.672144972967898, 0.635840438896716];
source.cartesian_A = source.fractional * source.lattice_A;
source.occupancy = ones(16, 1);
source.cr_source_index = find(strcmp(source.species, 'Cr'));
source.cr_label = source.label(source.cr_source_index);
source.cr_fractional = source.fractional(source.cr_source_index, :);
source.cr_cartesian_A = source.cartesian_A(source.cr_source_index, :);
source.cr_layer = [1; 1; 2; 2]; % 1: bottom, 2: top
source.cr_sublattice = [2; 1; 2; 1]; % 1: A, 2: B
source.coordinate_convention = [
    'a_Soenen=-a2_BiDB, b_Soenen=a1_BiDB; ' ...
    'top s0=Cr4, top s1=Cr3, bottom s0=Cr2, bottom s1=Cr1'];
end
