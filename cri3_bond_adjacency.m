function adjacency = cri3_bond_adjacency(bonds, n_sites)
%CRI3_BOND_ADJACENCY Expand unique bonds into reciprocal directed entries.

adjacency.from = [bonds.i; bonds.j];
adjacency.to = [bonds.j; bonds.i];
adjacency.J_J = [bonds.J_J; bonds.J_J];
adjacency.unique_bond = [(1:bonds.count)'; (1:bonds.count)'];
adjacency.count_by_site = accumarray(adjacency.from, 1, [n_sites, 1]);
end
