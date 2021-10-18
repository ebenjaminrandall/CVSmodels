function [outputs,rout,J] = model_sol(pars,theta,data)

pars = exp(pars); 

%% Unpack data structure 

tspan    = data.time; 
stim_per = data.stim_per;

%% Get initial conditions 

init = initialconditions(pars); 

%% run the simulation 
try 
M = speye(19);
M(1,1) = 0; 
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0;
options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',1e-10);

[t,y] = ode15s(@model,[0 30*stim_per],init,options,pars,stim_per);
init_new = y(end,:); 

% figure(100)
% plot(t,y(:,9),t,y(:,14))

sol = ode15s(@model,[tspan(1) tspan(end)],init_new,options,pars,stim_per);
sols = deval(tspan,sol); 
o = zeros(16,length(tspan)); 
for i = 1:length(tspan) 
    [~,o(:,i)] = model(tspan(i),sols(:,i),pars,stim_per);
end 

%% Outputs

xm_LV  = sols(1,:); % LV heart geometry variable, cm
xm_SEP = sols(2,:); % septum heart geometry variable, cmsizs
xm_RV  = sols(3,:); % RV heart geometry variable, cm
ym     = sols(4,:); % Heart geometry variable, cm
SL_LV  = sols(5,:);
SL_SEP = sols(6,:);
SL_RV  = sols(7,:);

V_LA   = sols(8,:);  % volume of LA
V_LV   = sols(9,:);  % volume LV, mL
V_Ao   = sols(10,:); % volume of aorta
V_SA   = sols(11,:); % volume of systemic arterys
V_SV   = sols(12,:); % volume of systemic veins
V_RA   = sols(13,:); % volume of RA
V_RV   = sols(14,:); % volume RV, mL
V_PA   = sols(15,:); % volume of pulmonary arterys
V_PV   = sols(16,:); % volume of pulmonary veins

V_T = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao + V_RA + V_LA;

P_LA = o(1,:);
P_LV = o(2,:);
P_Ao = o(3,:);
P_SA = o(4,:);
P_SV = o(5,:);
P_RA = o(6,:);
P_RV = o(7,:);
P_PA = o(8,:);
P_PV = o(9,:);

Q_m = o(10,:);
Q_a = o(11,:); 
Q_t = o(12,:);
Q_p = o(13,:); 

sigmapas_LV  = o(14,:);
sigmapas_SEP = o(15,:);
sigmapas_RV  = o(16,:);

displacements.xm_LV  = xm_LV; 
displacements.xm_SEP = xm_SEP; 
displacements.xm_RV  = xm_RV; 
displacements.ym     = ym; 

lengths.SL_LV  = SL_LV; 
lengths.SL_SEP = SL_SEP; 
lengths.SL_RV  = SL_RV;

volumes.V_LA = V_LA; 
volumes.V_LV = V_LV; 
volumes.V_Ao = V_Ao; 
volumes.V_SA = V_SA; 
volumes.V_SV = V_SV; 
volumes.V_RA = V_RA; 
volumes.V_RV = V_RV; 
volumes.V_PA = V_PA;
volumes.V_PV = V_PV;
volumes.V_T  = V_T; 

pressures.P_LA = P_LA; 
pressures.P_LV = P_LV; 
pressures.P_Ao = P_Ao; 
pressures.P_SA = P_SA; 
pressures.P_SV = P_SV; 
pressures.P_RA = P_RA; 
pressures.P_RV = P_RV; 
pressures.P_PA = P_PA;
pressures.P_PV = P_PV;

flows.Q_m = Q_m; 
flows.Q_a = Q_a; 
flows.Q_t = Q_t; 
flows.Q_p = Q_p; 

stresses.sigmapas_LV  = sigmapas_LV; 
stresses.sigmapas_SEP = sigmapas_SEP; 
stresses.sigmapas_RV  = sigmapas_RV; 

outputs.time = tspan; 
outputs.displacements = displacements;
outputs.lengths = lengths;
outputs.volumes = volumes; 
outputs.pressures = pressures; 
outputs.flows = flows; 
outputs.stresses = stresses; 

%% Residual 

[rout,J] = calculatecost(outputs,theta,data); 
catch 
    outputs = []; 
    rout = 100; 
    J = rout'*rout; 
end 

