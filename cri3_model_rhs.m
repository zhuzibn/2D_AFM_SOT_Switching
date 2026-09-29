function [dspindt, energy, fields, torque] = cri3_model_rhs(time_s, spins, model)
%CRI3_MODEL_RHS Compose verified fields and nonconservative LLGS torques.

[energy, fields] = cri3_energy_fields(spins, model.system, model.bonds, ...
    model.parameters);
if isfield(model, 'torque_function') && ~isempty(model.torque_function)
    torque = model.torque_function(time_s);
else
    torque = model.torque;
end
dspindt = cri3_llgs_rhs(spins, fields.total_T, model.system.layer, ...
    model.dynamics, torque);
end
