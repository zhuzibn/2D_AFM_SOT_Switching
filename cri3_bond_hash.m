function digest = cri3_bond_hash(bonds)
%CRI3_BOND_HASH Stable SHA-256 audit hash of the physical bond ledger.

parts = cell(bonds.count, 1);
for b = 1:bonds.count
    parts{b} = sprintf(['%d|%d|%d,%d|%d,%d|%d,%d|%s|%s|' ...
        '%.17g|%.17g|%.17g|%.17g\n'], ...
        bonds.i(b), bonds.j(b), bonds.source_cell(b, :), ...
        bonds.target_cell_unwrapped(b, :), bonds.image_offset(b, :), ...
        bonds.source_pair{b}, bonds.shell{b}, bonds.J_source_meV(b), ...
        bonds.J_converted_meV(b), bonds.J_fit_meV(b), bonds.J_J(b));
end
payload = strjoin(parts, '');
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(uint8(payload));
raw = typecast(engine.digest(), 'uint8');
digest = lower(reshape(dec2hex(raw, 2).', 1, []));
end
