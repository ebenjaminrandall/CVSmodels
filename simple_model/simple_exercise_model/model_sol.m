function outputs = model_sol(pars,data)


stim_period = data.stim_period; 

%% Get initial conditions 

init = initialconditions(pars,data); 

%% Solve model 

M = speye(19);
M(1,1) = 0; 
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0;
options = odeset('Mass',M,'RelTol',1e-6,'AbsTol',1e-6);%,'MaxStep',stim_period/50);

% Run model to steady state 
[t,y] = ode15s(@model,[0 20*stim_period],init,options,pars,data);
init = y(end,:);

% Run model for 2 heart periods 
[t,y] = ode15s(@model,[0 2*stim_period],init,options,pars,data);
o = zeros(31,length(t)); 
for i = 1:length(t) 
    [~,o(:,i)] = model(t(i),y(i,:),pars,data);
end 

xm_LV  = y(:,1); % LV heart geometry variable, cm
xm_SEP = y(:,2); % septum heart geometry variable, cmsizs
xm_RV  = y(:,3); % RV heart geometry variable, cm
ym     = y(:,4); % Heart geometry variable, cm

SL_LV  = y(:,5);
SL_SEP = y(:,6);
SL_RV  = y(:,7);

V_LA   = y(:,8); 
V_LV   = y(:,9); 
V_Ao   = y(:,10); 
V_SA   = y(:,11);
V_SV   = y(:,12);
V_RA   = y(:,13); 
V_RV   = y(:,14);
V_PA   = y(:,15); 
V_PV   = y(:,16);

C_LV = y(:,17); 
C_SEP = y(:,18); 
C_RV = y(:,19); 

P_LA = o(1,:);
P_LV = o(2,:);
P_Ao = o(3,:);
P_SA = o(4,:);
P_SV = o(5,:);
P_RA = o(6,:);
P_RV = o(7,:);
P_PA = o(8,:);
P_PV = o(9,:);

Q_m  = o(10,:);
Q_a  = o(11,:);
Q_t  = o(12,:);
Q_p  = o(13,:); 

sigmapas_LV  = o(14,:);
sigmapas_SEP = o(15,:);
sigmapas_RV  = o(16,:);

sigmaact_LV  = o(17,:); 
sigmaact_SEP = o(18,:); 
sigmaact_RV  = o(19,:); 

SLo_LV  = o(20,:); 
SLo_SEP = o(21,:); 
SLo_RV  = o(22,:); 

z_LV  = o(23,:); 
z_SEP = o(24,:); 
z_RV  = o(25,:); 

Vm_LV  = o(26,:); 
Vm_SEP = o(27,:); 
Vm_RV  = o(28,:); 

Am_LV  = o(29,:); 
Am_SEP = o(30,:); 
Am_RV  = o(31,:); 

%% Outputs 

displacements.xm_LV  = xm_LV; 
displacements.xm_SEP = xm_SEP;
displacements.xm_RV  = xm_RV; 
displacements.ym     = ym; 

lengths.SL_LV   = SL_LV; 
lengths.SL_SEP  = SL_SEP;
lengths.SL_RV   = SL_RV; 
lengths.SLo_LV  = SLo_LV; 
lengths.SLo_SEP = SLo_SEP;
lengths.SLo_RV  = SLo_RV; 

volumes.V_LA = V_LA; 
volumes.V_LV = V_LV; 
volumes.V_Ao = V_Ao; 
volumes.V_SA = V_SA;
volumes.V_SV = V_SV; 
volumes.V_RA = V_RA; 
volumes.V_RV = V_RV; 
volumes.V_PA = V_PA; 
volumes.V_PV = V_PV; 

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

activation.C_LV  = C_LV; 
activation.C_SEP = C_SEP;
activation.C_RV  = C_RV; 

stresses.sigmapas_LV  = sigmapas_LV; 
stresses.sigmapas_SEP = sigmapas_SEP;
stresses.sigmapas_RV  = sigmapas_RV; 
stresses.sigmaact_LV  = sigmaact_LV; 
stresses.sigmaact_SEP = sigmaact_SEP; 
stresses.sigmaact_RV  = sigmaact_RV; 

Triseg.Vm_LV  = Vm_LV; 
Triseg.Vm_SEP = Vm_SEP;
Triseg.Vm_RV  = Vm_RV;
Triseg.Am_LV  = Am_LV;
Triseg.Am_SEP = Am_SEP;
Triseg.Am_RV  = Am_RV; 
Triseg.z_LV   = z_LV; 
Triseg.z_SEP  = z_SEP; 
Triseg.z_RV   = z_RV; 

outputs.time          = t; 
outputs.displacements = displacements;
outputs.lengths       = lengths; 
outputs.volumes       = volumes; 
outputs.pressures     = pressures; 
outputs.flows         = flows; 
outputs.activation    = activation; 
outputs.stresses      = stresses; 
outputs.Triseg        = Triseg; 
