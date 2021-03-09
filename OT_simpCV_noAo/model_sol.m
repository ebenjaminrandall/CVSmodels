function [outputs,rout,J] = model_sol(pars,data)


H = 75; 

%% Extract parameters
pars = exp(pars); 

k_TS = pars(30); 
k_TR = pars(31); 
tau_v = pars(32); 

%% Unpack data structure 

tspan = data.tspan; 
SPbar = data.SPbar;
DPbar = data.DPbar; 

dt = data.dt; 

EDV_LV = data.EDV_LV;
EDV_RV = EDV_LV;
ESV_LV = data.ESV_LV;
ESV_RV = ESV_LV;

ODE_TOL = data.gpars.ODE_TOL;

%% Get initial conditions

init = initialconditions(pars,data); 

%% Solve model 

M = speye(length(init));
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 

opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL);
ndone = 0; 
k1 = 1; 
tc_vec = []; 
beats   = []; 
x_keep = []; 
y_v_keep = []; 
y_a_keep = []; 
while ndone == 0 
    %clearvars t_per TS TR Y sols
    
    t_per = 60 / H; 
    k2 = k1 + round(t_per / dt) - 1; 
    if k2 > length(tspan) 
        k2 = length(tspan);
    end 
    t = tspan(k1:k2); 
    tc = t(1); 
    tc_vec = [tc_vec tc*ones(1,length(t)-1)]; 
    beats = [beats k1]; 
   
    TS = k_TS*t_per; 
    TR = k_TR*t_per; 
    
    % Atrial activation 
    y_a = zeros(size(t)); 
    for i = 1:length(t) 
        tn_a = t(i) - tc;
        if tn_a >= 0 && tn_a < TS 
            y_a(i) = 0.5*(1 - cos(pi*tn_a/TS)); 
        elseif tn_a >= TS && tn_a < TR + TS 
            y_a(i) = 0.5*(1 + cos(pi*(tn_a - TS)/TR)); 
        else
            y_a(i) = 0; 
        end 
    end 
    
    % Ventricular activation 
    y_v = zeros(size(t)); 
    for i = 1:length(t) 
        tn_v = t(i) - tc - tau_v;
        if tn_v >= 0 && tn_v < TS
            y_v(i) = 0.5*(1 - cos(pi*tn_v/TS)); 
        elseif tn_v >= TS  && tn_v < TR + TS 
            y_v(i) = 0.5*(1 + cos(pi*(tn_v - TS)/TR)); 
        else
            y_v(i) = 0; 
        end 
    end 
    
    y_a_keep = [y_a_keep y_a(1:end-1)]; 
    y_v_keep = [y_v_keep y_v(1:end-1)]; 
    
    Y_a = griddedInterpolant(t,y_a); 
    Y_v = griddedInterpolant(t,y_v); 
    
    sol  = ode15s(@model,[t(1) t(end)],init,opts,pars,data,Y_a,Y_v);
    sols = deval(sol,t);
    
    x_keep = [x_keep sols(:,1:end-1)]; 
    init = sols(:,end);
    
    k1 = k2; 
    
    if k1 >= length(tspan)
        ndone = 1; 
    end 
    
end 
sols   = [x_keep sols(:,end)]; 
tc_vec = [tc_vec tc_vec(end)];
y_a_keep = [y_a_keep y_a(end)]; 
y_v_keep = [y_v_keep y_v(end)]; 

Y_a = griddedInterpolant(tspan,y_a_keep); 
Y_v = griddedInterpolant(tspan,y_v_keep); 

%% Calculate pressures

o = zeros(37,length(tspan)); 
for i = 1:length(tspan) 
    [~,o(:,i)] = model(tspan(i),sols(:,i),pars,data,Y_a,Y_v);
end 

%% Outputs 

[~,fixpars] = parameters(data); 

V_sau = fixpars(9); 
V_svu = fixpars(10); 
V_pau = fixpars(11); 
V_pvu = fixpars(12); 

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

outputs.beats         = beats; 
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

%% Cost 

EDV_LV = EDV_LV * 1e6; 
ESV_LV = ESV_LV * 1e6; 
EDV_RV = EDV_RV * 1e6; 
ESV_RV = ESV_RV * 1e6; 

V_lv = outputs.volumes.V_lv; 
V_rv = outputs.volumes.V_rv;

P_lv = outputs.pressures.P_lv;
P_rv = outputs.pressures.P_rv;
P_sa = outputs.pressures.P_sa; 

beat = beats(end-3):beats(end-1); 

P_lvM = 125; 
P_rvM = 30; 

rout = [(max(V_lv(beat)) - EDV_LV)/EDV_LV; 
    (min(V_lv(beat)) - ESV_LV)/ESV_LV; 
    (max(V_rv(beat)) - EDV_RV)/EDV_RV; 
    (min(V_rv(beat)) - ESV_RV)/ESV_RV;
    (max(P_lv(beat)) - P_lvM)/P_lvM; 
    (max(P_rv(beat)) - P_rvM)/P_rvM; 
    (max(P_sa(beat)) - SPbar)/SPbar; 
    (min(P_sa(beat)) - DPbar)/DPbar;
    ];
    
J    = rout' * rout; 


