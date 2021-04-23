function [outputs,rout,J] = model_sol(adjpars,data)

tspan = data.tspan;  
DPbar = data.DPbar; 

ODE_TOL = data.gpars.ODE_TOL;

%% Parameters 

% Compliances (m^3 kPa^(-1))
C_ao = adjpars(3); 
C_sa = adjpars(4); 
C_sv = adjpars(5); 
C_pa = adjpars(6); 
C_pv = adjpars(7); 

% Wall volume of ventricular wall segment (m^3)
Vw_lv  = adjpars(14); 
Vw_sep = adjpars(15); 
Vw_rv  = adjpars(16); 

% Reference midwall surface area (m^2)
Amref_lv  = adjpars(17); 
Amref_sep = adjpars(18); 
Amref_rv  = adjpars(19); 

% Sarcomere length parameters (m)
Lsref   = adjpars(23);
Lsc0    = adjpars(24); 
Lse_iso = adjpars(25); 

% Force scaling factors (kPa)
k_act = adjpars(26); 
k_pas = adjpars(27); 

V_aou = adjpars(30); 
V_sau = adjpars(31);
V_svu = adjpars(32); 
V_pau = adjpars(33); 
V_pvu = adjpars(34); 

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

V_lv = sols(5,:); 
V_ao = sols(6,:); 
V_sa = sols(7,:);
V_sv = sols(8,:); 
V_rv = sols(9,:);
V_pa = sols(10,:); 
V_pv = sols(11,:); 

% Axial distance
xm_lv  = sols(1,:); 
xm_sep = sols(2,:); 
xm_rv  = sols(3,:); 

% Radial distance
ym = sols(4,:); 

% Contractile element length 
Lsc_lv  = sols(12,:);
Lsc_sep = sols(13,:); 
Lsc_rv  = sols(14,:); 

% Mechanical activation 
Ca_lv  = sols(15,:);
Ca_sep = sols(16,:); 
Ca_rv  = sols(17,:); 

%% Calculate pressures

Vm_lv  = - (pi / 6) .* xm_lv  .* (xm_lv.^2  + 3 * ym.^2); 
Vm_sep = (pi / 6) .* xm_sep .* (xm_sep.^2 + 3 * ym.^2); 
Vm_rv  = (pi / 6) .* xm_rv  .* (xm_rv.^2  + 3 * ym.^2);

Am_lv  = pi * (xm_lv.^2  + ym.^2);
Am_sep = pi * (xm_sep.^2 + ym.^2); 
Am_rv  = pi * (xm_rv.^2  + ym.^2);

Cm_lv  = - 2 * xm_lv  ./ (xm_lv.^2  + ym.^2);
Cm_sep = 2 * xm_sep ./ (xm_sep.^2 + ym.^2);
Cm_rv  = 2 * xm_rv  ./ (xm_rv.^2  + ym.^2);

z_lv   = 3 * Cm_lv  .* Vw_lv   ./ (2 .* Am_lv);
z_sep  = 3 * Cm_sep .* Vw_sep  ./ (2 .* Am_sep); 
z_rv   = 3 * Cm_rv  .* Vw_rv   ./ (2 .* Am_rv);

eps_lv  = (1 / 2) * log(Am_lv   ./ Amref_lv)  - (1 / 12) * z_lv.^2  - 0.019 * z_lv.^4;
eps_sep = (1 /2 ) * log( Am_sep ./ Amref_sep) - (1 / 12) * z_sep.^2 - 0.019 * z_sep.^4; 
eps_rv  = (1 / 2) * log(Am_rv   ./ Amref_rv)  - (1 / 12) * z_rv.^2  - 0.019 * z_rv.^4;

Ls_lv  = Lsref * exp(eps_lv); 
Ls_sep = Lsref * exp(eps_sep); 
Ls_rv  = Lsref * exp(eps_rv); 

sigma_act_lv  = Ca_lv  .* (Lsc_lv  - Lsc0) ./ 1e-6 .* (Ls_lv  - Lsc_lv)  ./ Lse_iso;  
sigma_act_sep = Ca_sep .* (Lsc_sep - Lsc0) ./ 1e-6 .* (Ls_sep - Lsc_sep) ./ Lse_iso;
sigma_act_rv  = Ca_rv  .* (Lsc_rv  - Lsc0) ./ 1e-6 .* (Ls_rv  - Lsc_rv)  ./ Lse_iso; 

sigma_pas_lv  = 36 * max(0,eps_lv  - 0.1).^2 + 0.1 .* (eps_lv  - 0.1) + 0.0025 * exp(30 .* eps_lv); 
sigma_pas_sep = 36 * max(0,eps_sep - 0.1).^2 + 0.1 .* (eps_sep - 0.1) + 0.0025 * exp(30 .* eps_sep); 
sigma_pas_rv  = 36 * max(0,eps_rv  - 0.1).^2 + 0.1 .* (eps_rv  - 0.1) + 0.0025 * exp(30 .* eps_rv); 

sigma_lv  = k_act .* sigma_act_lv  + k_pas .* sigma_pas_lv; 
sigma_sep = k_act .* sigma_act_sep + k_pas .* sigma_pas_sep; 
sigma_rv  = k_act .* sigma_act_rv  + k_pas .* sigma_pas_rv; 

Tm_lv  = (Vw_lv .* sigma_lv ./ (2 .* Am_lv)) .* (1 + (z_lv.^2)/3  + (z_lv.^4)/5); 
Tm_rv  = (Vw_rv .* sigma_rv ./ (2 .* Am_rv)) .* (1 + (z_rv.^2)/3  + (z_rv.^4)/5);

Tx_lv = - Tm_lv .* 2 .* xm_lv .* ym ./ (xm_lv.^2 + ym.^2);
Tx_rv = Tm_rv .* 2 .* xm_rv .* ym ./ (xm_rv.^2 + ym.^2);

ptrans_lv = 2 * Tx_lv ./ ym;
ptrans_rv = 2 * Tx_rv ./ ym;

P_lv = -ptrans_lv;
P_ao = (V_ao - V_aou) ./ C_ao; 
P_sa = (V_sa - V_sau) ./ C_sa; 
P_sv = (V_sv - V_svu) ./ C_sv; 
P_rv = ptrans_rv;
P_pa = (V_pa - V_pau) ./ C_pa; 
P_pv = (V_pv - V_pvu) ./ C_pv; 

%% Outputs 

% Convert to m^3 to mL
volumes.V_lv = V_lv * 1e6; 
volumes.V_ao = V_ao * 1e6; 
volumes.V_sa = V_sa * 1e6; 
volumes.V_sv = V_sv * 1e6; 
volumes.V_rv = V_rv * 1e6; 
volumes.V_pa = V_pa * 1e6; 
volumes.V_pv = V_pv * 1e6; 

% Convert kPa to mmHg
pressures.P_lv = P_lv * 7.5; 
pressures.P_ao = P_ao * 7.5; 
pressures.P_sa = P_sa * 7.5; 
pressures.P_sv = P_sv * 7.5; 
pressures.P_rv = P_rv * 7.5; 
pressures.P_pa = P_pa * 7.5; 
pressures.P_pv = P_pv * 7.5; 

deformations.xm_lv  = xm_lv; 
deformations.xm_sep = xm_sep; 
deformations.xm_rv  = xm_rv; 
deformations.ym     = ym; 

areas.Am_lv = Am_lv; 
areas.Am_sep = Am_sep; 
areas.Am_rv  = Am_rv; 

wallvolumes.Vw_lv  = Vm_lv; 
wallvolumes.Vw_sep = Vm_sep; 
wallvolumes.Vw_rv  = Vm_rv; 

curvatures.Cm_lv  = Cm_lv;
curvatures.Cm_sep = Cm_sep;
curvatures.Cm_rv  = Cm_rv; 

stresses.passive.sigma_pas_lv  = sigma_pas_lv;
stresses.passive.sigma_pas_sep = sigma_pas_sep;
stresses.passive.sigma_pas_rv  = sigma_pas_rv;

stresses.active.sigma_act_lv  = sigma_act_lv;
stresses.active.sigma_act_sep = sigma_act_sep;
stresses.active.sigma_act_rv  = sigma_act_rv;

stresses.total.sigma_lv  = sigma_lv;
stresses.total.sigma_sep = sigma_sep;
stresses.total.sigma_rv  = sigma_rv;

lengths.Lsc_lv  = Lsc_lv;
lengths.Lsc_sep = Lsc_sep; 
lengths.Lsc_rv  = Lsc_rv; 

activation.Ca_lv  = Ca_lv;
activation.Ca_sep = Ca_sep;
activation.Ca_rv  = Ca_rv; 


outputs.volumes      = volumes; 
outputs.pressures    = pressures; 
outputs.deformations = deformations; 
outputs.areas        = areas;
outputs.wallvolumes  = wallvolumes; 
outputs.curvatures   = curvatures; 
outputs.stresses     = stresses;
outputs.lengths      = lengths; 
outputs.activation   = activation; 

rout = (min(P_sa) - DPbar) / DPbar; 
J    = rout' * rout; 


