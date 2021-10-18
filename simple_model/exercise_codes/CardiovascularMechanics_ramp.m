% This code is the main driver for the cardiovascular mechanics model 
clear; 

theta = 0.0; % exercise level, between 0 and 1

HR = 64*(1 + 1.9*theta);  % 1/sec

Vw = [80 38 28]; % Heart wall volumes (mL)
Amref_LV  = 0.975*80; % LV midwall reference surface area, cm^2
Amref_SEP = 0.975*45; % SEP midwall reference surface area, cm^2
Amref_RV  = 1.12*100; % RV midwall reference surface area, cm^2

V_LV  = 150;  % initial V_LV (mL)
V_RV  = 150; % initial V_LV (mL)

freq = HR/60; %Hz
stim_period = 1/freq;
para = [ stim_period Vw(1) Vw(2) Vw(3) Amref_LV Amref_SEP Amref_RV theta];

xm_LV  = -5.0;
xm_SEP = +2.5;
xm_RV  = +8.0;
ym     = +5.0;

SL_LV   = 2.2;
SL_SEP  = 2.2;
SL_RV   = 2.2;

x0 = [xm_LV ,xm_SEP ,xm_RV ,ym];

Vw_LV  = Vw(1); % LV wall volume, mL 
Vw_SEP = Vw(2); % Septal wall volume, mL 
Vw_RV  = Vw(3); % RV wall volume, mL 

% opts = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
opts = optimset('MaxFunEvals',10000,'MaxIter',1000);
TrisegEquations(x0,Vw_LV,Vw_SEP,Vw_RV, V_LV, V_RV);
x = fsolve(@TrisegEquations,x0,opts,Vw_LV,Vw_SEP,Vw_RV, V_LV, V_RV);


%% set up IC's and Parameters
% 
vfactor = 1.0; % control 20-yo
V_SA = vfactor*200;
V_SV = vfactor*1175;
V_PA = vfactor*35;
V_PV = vfactor*65; 
V_Ao = vfactor*60;
V_RA = 100; % mL
V_LA = 100; % mL

init = [x, SL_LV, SL_SEP, SL_RV, V_LV, V_RV, ...
       V_SV, V_PV ,V_SA ,V_PA, V_Ao V_RA V_LA, ...
       0, 0, 0]';


%% run the (resting) simulation 

% load init

M = speye(19);
M(1,1) = 0; 
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0;
options = odeset('Mass',M,'RelTol',1e-6,'AbsTol',1e-6,'MaxStep',stim_period/50);

[t,y] = ode15s(@dXdT_cardiovascular_mechanics_exercise,[0 20*stim_period],init,options,para);
init = y(end,:);
% save init init
[t,y] = ode15s(@dXdT_cardiovascular_mechanics_exercise,[0 2*stim_period],init,options,para);
o = zeros(14,length(t)); 
for i = 1:length(t) 
    [~,o(:,i)] = dXdT_cardiovascular_mechanics_exercise(t(i),y(i,:),para);
end 

% xm_LV  = y(:,1); % LV heart geometry variable, cm
% xm_SEP = y(:,2); % septum heart geometry variable, cmsizs
% xm_RV  = y(:,3); % RV heart geometry variable, cm
% ym     = y(:,4); % Heart geometry variable, cm
% SL_LV  = y(:,5);
% SL_SEP = y(:,6);
% SL_RV  = y(:,7);
% V_LV   = y(:,8); % volume LV, mL
% V_RV   = y(:,9); % volume RV, mL
% 
% V_SV   = y(:,10); % volume of systemic veins
% V_PV   = y(:,11); % volume of pulmonary veins
% V_SA   = y(:,12); % volume of systemic arterys
% V_PA   = y(:,13); % volume of pulmonary arterys
% V_Ao   = y(:,14); % volume of aorta
% V_RA   = y(:,15); % volume of RA
% V_LA   = y(:,16); % volume of LA
% V_T = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao + V_RA + V_LA;


%% Exercise ramp
theta = 0:0.1:1;
for i = 1:length(theta)
  theta(i)
  HR = 64*(1 + 1.9*theta(i));  % 1/sec
  freq = HR/60; %Hz
  stim_period = 1/freq;
  para = [ stim_period Vw(1) Vw(2) Vw(3) Amref_LV Amref_SEP Amref_RV theta(i)];
% run to steady state at theta(i)
  init = y(end,:);
  [t,y] = ode15s(@dXdT_cardiovascular_mechanics_exercise,[0 20*stim_period],init,options,para);
  init = y(end,:);
  % save init init
  [t,y] = ode15s(@dXdT_cardiovascular_mechanics_exercise,[0 2*stim_period],init,options,para);
  o = zeros(14,length(t)); 
  for j = 1:length(t) 
    [~,o(:,j)] = dXdT_cardiovascular_mechanics_exercise(t(j),y(j,:),para);
  end

  V_LV   = y(:,8); % volume LV, mL
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
  Q_t  = o(11,:);
  sigmapas_LV = o(12,:);
  sigmapas_SEP = o(13,:);
  sigmapas_RV = o(14,:);

  SV(i) = max(V_LV) - min(V_LV);
  EF(i) = SV(i)/max(V_LV);
  CO(i) = SV(i)*HR;
  SP(i) = max(P_Ao);
  DP(i) = min(P_Ao);
  
  Pwedge(i) = mean(interp1(t,P_PV,(0:0.1:1).*stim_period));
  
end

HR = 64*(1 + 1.9*theta);

figure(1); plot(HR,SP,'k--',HR,DP,'k--'); ylabel('systemic pressures'); xlabel('HR');
figure(2); plot(HR,CO*1e-3,'k--'); ylabel('CO'); xlabel('HR');
figure(3); plot(HR,Pwedge,'k--'); ylabel('P_{PV}'); xlabel('HR');

