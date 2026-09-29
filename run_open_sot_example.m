function result = run_open_sot_example()
%RUN_OPEN_SOT_EXAMPLE Run the documented candidate open-flake simulation.

root = fileparts(mfilename('fullpath'));
output_dir = fullfile(root, 'results');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

config = cri3_default_config();
config.nx = 8;
config.ny = 8;
config.boundary = 'open';
config.alpha = 0.001; % Explicit candidate choice, not a production default
config.t_final_s = 2e-12;
config.torque.dl_T = [0.01, -0.01]; % T: [bottom, top]
config.torque.fl_T = [0, 0];
config.output_file = fullfile(output_dir, 'candidate.mat');

if isfile(config.output_file)
    error('CrI3:OutputExists', ...
        'Output exists; choose another path or remove it explicitly: %s', ...
        config.output_file);
end

result = cri3_run(config);
end
