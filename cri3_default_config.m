function config = cri3_default_config()
%CRI3_DEFAULT_CONFIG Documented inputs for the fitted scalar CrI3 model.
%   ALPHA and T_FINAL_S intentionally have no production defaults. Set them
%   explicitly for a dynamical run. Soenen et al. used alpha=0.001 for their
%   separate T=0 stability calculations; that value is not silently assumed
%   to be a universal SOT-switching damping constant here.

mu_B = 9.2740100783e-24; % [J/T]
hbar = 1.054571817e-34; % [J s]
meV_to_J = 1.602176634e-22;

config.model_id = 'bidb-abprime-cri3-fitted-soenen-scalar-v1';
config.nx = 2;
config.ny = 2;
config.boundary = 'periodic';
config.mu_J_per_T = 3 * mu_B;
config.K_J = 0.1575 * meV_to_J;
config.B_external_T = [0, 0, 0];
config.gamma_rad_per_s_T = 2 * mu_B / hbar;
config.alpha = [];
config.torque.sigma = [0, 1, 0];
config.torque.dl_T = [0, 0]; % [bottom, top]
config.torque.fl_T = [0, 0]; % [bottom, top]
config.torque_function = []; % optional function handle torque=time-dependent torque(time_s)
config.dt_s = 1e-15;
config.t_final_s = [];
config.save_every = 1;
config.initial_mode = 'afm';
config.initial_canting = 0;
config.initial_spins = [];
config.output_file = '';
config.overwrite = false;
config.source_alpha_note = [
    'Soenen supplement stability simulation: alpha=0.001; ' ...
    'production value must be selected explicitly for the target dynamics.'];
config.units.energy = 'J';
config.units.field_and_torque = 'T';
config.units.time = 's';
config.units.length = 'angstrom';
end
