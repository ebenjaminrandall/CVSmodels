function [outputs,rout,J] = model_sol(adjpars,data)

tspan = data.tspan;  
DPbar = data.DPbar; 

ODE_TOL = data.gpars.ODE_TOL;

%% Parameters 

V_sau = adjpars(35);
V_svu = adjpars(36); 
V_pau = adjpars(37); 
V_pvu = adjpars(38); 

%% Get initial conditions

init = initialconditions(adjpars,data); 

%% Solve model 

M = speye(length(init));
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 

opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL);
sol  = ode15s(@model,[tspan(1) tspan(end)],init,opts,adjpars,data);
sols = deval(sol,tspan);

%% Calculate pressures

o = zeros(37,length(tspan)); 
for i = 1:length(tspan) 
    [~,o(:,i)] = model(tspan(i),sols(:,i),adjpars,data);
end 

%% Outputs 

% Convert m to cm 
displacements.xm_lv  = sols(1,:) * 1e2; 
displacements.xm_sep = sols(2,:) * 1e2; 
displacements.xm_rv  = sols(3,:) * 1e2; 
displacements.ym     = sols(4,:) * 1e2; 

% Convert m to micrometers
lengths.Lsc_lv  = sols(5,:) * 1e6;
lengths.Lsc_sep = sols(6,:) * 1e6; 
lengths.Lsc_rv  = sols(7,:) * 1e6; 

% Convert to m^3 to mL
volumes.V_la = sols(8,:) * 1e6; 
volumes.V_lv = sols(9,:) * 1e6; 
volumes.V_sa = (sols(10,:) + V_sau) * 1e6; 
volumes.V_sv = (sols(11,:) + V_svu) * 1e6; 
volumes.V_ra = sols(12,:) * 1e6; 
volumes.V_rv = sols(13,:) * 1e6; 
volumes.V_pa = (sols(14,:) + V_pau) * 1e6; 
volumes.V_pv = (sols(15,:) + V_pvu) * 1e6; 

% Convert kPa to mmHg
pressures.P_la = o(1,:) * 7.5; 
pressures.P_lv = o(2,:) * 7.5; 
pressures.P_sa = o(3,:) * 7.5; 
pressures.P_sv = o(4,:) * 7.5; 
pressures.P_ra = o(5,:) * 7.5; 
pressures.P_rv = o(6,:) * 7.5; 
pressures.P_pa = o(7,:) * 7.5; 
pressures.P_pv = o(8,:) * 7.5; 

% Convert m^3 to cm^3
wallvolumes.Vm_lv  = o(9,:)  * 1e6; 
wallvolumes.Vm_sep = o(10,:)  * 1e6; 
wallvolumes.Vm_rv  = o(11,:) * 1e6; 

% Convert m^2 to cm^2
areas.Am_lv  = o(12,:) * 1e4; 
areas.Am_sep = o(13,:) * 1e4; 
areas.Am_rv  = o(14,:) * 1e4; 

% Convert m^(-1) to cm^(-1)
curvatures.Cm_lv  = o(15,:) * 1e-2;
curvatures.Cm_sep = o(16,:) * 1e-2;
curvatures.Cm_rv  = o(17,:) * 1e-2; 

strains.eps_lv  = o(18,:); 
strains.eps_sep = o(19,:); 
strains.eps_rv  = o(20,:); 

stresses.passive.sigma_pas_lv  = o(21,:);
stresses.passive.sigma_pas_sep = o(22,:);
stresses.passive.sigma_pas_rv  = o(23,:);

stresses.active.sigma_act_lv  = o(24,:);
stresses.active.sigma_act_sep = o(25,:);
stresses.active.sigma_act_rv  = o(26,:);

stresses.total.sigma_lv  = o(27,:);
stresses.total.sigma_sep = o(28,:);
stresses.total.sigma_rv  = o(29,:);

% Convert m^3 s^(-1) to L min^(-1)
flows.Q_m_valve = o(30,:) * 1e3 * 60; 
flows.Q_a_valve = o(31,:) * 1e3 * 60; 
flows.Q_t_valve = o(32,:) * 1e3 * 60; 
flows.Q_p_valve = o(33,:) * 1e3 * 60; 

flows.Q_sa = o(34,:) * 1e3 * 60; 
flows.Q_sv = o(35,:) * 1e3 * 60; 
flows.Q_pa = o(36,:) * 1e3 * 60; 
flows.Q_pv = o(37,:) * 1e3 * 60; 

%outputs.beats         = beats; 
outputs.volumes       = volumes; 
outputs.pressures     = pressures; 
outputs.displacements = displacements; 
outputs.areas         = areas;
outputs.wallvolumes   = wallvolumes; 
outputs.curvatures    = curvatures; 
outputs.strains       = strains; 
outputs.stresses      = stresses;
outputs.lengths       = lengths; 
outputs.flows         = flows; 

rout = (min(pressures.P_sa) - DPbar) / DPbar; 
J    = rout' * rout; 


