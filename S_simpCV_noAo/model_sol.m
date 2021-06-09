function [outputs,rout,J] = model_sol(pars,data)


%% Extract parameters

pars = exp(pars); 

k_TS  = pars(21); 
k_TR  = pars(22); 
tau_v = pars(23); 

%% Unpack data structure 

dt = data.dt; 

tspan = data.tspan; 
SPbar = data.SPbar;
DPbar = data.DPbar; 
H     = data.HR.HR_rest; 

P_mmHg2kPa = data.units.P_mmHg2kPa; 
V_mL2m3    = data.units.V_mL2m3; 

EDV_LV = data.EDV_LV;
EDV_RV = EDV_LV;
ESV_LV = data.ESV_LV;
ESV_RV = ESV_LV;

ODE_TOL = data.gpars.ODE_TOL;

%% Get initial conditions

init = initialconditions(pars,data); 

%% Solve model  

M = speye(length(init)); 
M(9,9) = 0; 
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

o = zeros(18,length(tspan)); 
for i = 1:length(tspan) 
    [~,o(:,i)] = model(tspan(i),sols(:,i),pars,data,Y_a,Y_v);
end 

%% Outputs 

[~,fixpars] = parameters(data); 

V_sau = fixpars(9); 
V_svu = fixpars(10); 
V_pau = fixpars(11); 
V_pvu = fixpars(12); 

% Convert to m^3 to mL
volumes.V_la = sols(1,:) / V_mL2m3;  
volumes.V_lv = sols(2,:) / V_mL2m3; 
volumes.V_sa = (sols(3,:) + V_sau) / V_mL2m3; 
volumes.V_sv = (sols(4,:) + V_svu) / V_mL2m3; 
volumes.V_ra = sols(5,:) / V_mL2m3; 
volumes.V_rv = sols(6,:) / V_mL2m3; 
volumes.V_pa = (sols(7,:) + V_pau) / V_mL2m3;  
volumes.V_pv = (sols(8,:) + V_pvu) / V_mL2m3; 

% Convert kPa to mmHg
pressures.P_la = o(1,:) / P_mmHg2kPa; 
pressures.P_lv = o(2,:) / P_mmHg2kPa; 
pressures.P_sa = o(3,:) / P_mmHg2kPa; 
pressures.P_sv = o(4,:) / P_mmHg2kPa; 
pressures.P_ra = o(5,:) / P_mmHg2kPa; 
pressures.P_rv = o(6,:) / P_mmHg2kPa; 
pressures.P_pa = o(7,:) / P_mmHg2kPa; 
pressures.P_pv = o(8,:) / P_mmHg2kPa; 

% Convert m^3 s^(-1) to L min^(-1)
flows.Q_m_valve = o(9,:)  * 1e3 * 60; 
flows.Q_a_valve = o(10,:) * 1e3 * 60; 
flows.Q_sa      = o(11,:) * 1e3 * 60; 
flows.Q_sv      = o(12,:) * 1e3 * 60; 
flows.Q_t_valve = o(13,:) * 1e3 * 60; 
flows.Q_p_valve = o(14,:) * 1e3 * 60; 
flows.Q_pa      = o(15,:) * 1e3 * 60; 
flows.Q_pv      = o(16,:) * 1e3 * 60;

elastances.E_lv = o(17,:) / 1e6 * 7.5; 
elastances.E_rv = o(18,:) / 1e6 * 7.5; 

outputs.beats      = beats; 
outputs.volumes    = volumes; 
outputs.pressures  = pressures; 
outputs.flows      = flows; 
outputs.elastances = elastances; 

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

beat = beats(end-1):beats(end); 

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


