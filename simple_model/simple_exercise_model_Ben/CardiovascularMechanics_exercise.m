% This code is the main driver for the cardiovascular mechanics model 
clear; 

%% Inputs

theta = 1; % exercise level, between 0 and 1
HR = 64*(1 + 1.9*theta)  % 1/min
freq = HR/60; %Hz
stim_period = 1/freq;

% Total blood volume scaling factor 
vfactor = 1.0; % control 20-yo

%% Get parameters and initialconditions

pars = parameters; 
init = initialconditions(pars,vfactor); 

%% run the simulation 

% load init

M = speye(16);
M(1,1) = 0; 
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0;
options = odeset('Mass',M,'RelTol',1e-6,'AbsTol',1e-6,'MaxStep',stim_period/50);

[t,y] = ode15s(@dXdT_cardiovascular_mechanics_exercise,[0 20*stim_period],init,options,pars,stim_period,theta);
init = y(end,:);
% save init init
[t,y] = ode15s(@dXdT_cardiovascular_mechanics_exercise,[0 2*stim_period],init,options,pars,stim_period,theta);
o = zeros(14,length(t)); 
for i = 1:length(t) 
    [~,o(:,i)] = dXdT_cardiovascular_mechanics_exercise(t(i),y(i,:),pars,stim_period,theta);
end 

xm_LV  = y(:,1); % LV heart geometry variable, cm
xm_SEP = y(:,2); % septum heart geometry variable, cmsizs
xm_RV  = y(:,3); % RV heart geometry variable, cm
ym     = y(:,4); % Heart geometry variable, cm

SL_LV  = y(:,5);
SL_SEP = y(:,6);
SL_RV  = y(:,7);

V_LA = y(:,8); % volume of LA
V_LV = y(:,9); % volume LV, mL
V_Ao = y(:,10); % volume of aorta
V_SA = y(:,11); % volume of systemic arterys
V_SV = y(:,12); % volume of systemic veins
V_RA = y(:,13); % volume of RA
V_RV = y(:,14); % volume RV, mL
V_PA = y(:,15); % volume of pulmonary arterys
V_PV = y(:,16); % volume of pulmonary veins
V_T  = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao + V_RA + V_LA;

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
sigmapas_LV  = o(12,:);
sigmapas_SEP = o(13,:);
sigmapas_RV  = o(14,:);

%% Plotting Baseline

figureson = 0; 

figure(1)
clf
hold on 
h1 = plot(t,V_LV,'b','linewidth',2);
h2 = plot(t,V_RV,'r','linewidth',2);
legend([h1 h2],'LV','RV')
xlabel('Time (s)')
ylabel('Volume (mL)') 
set(gca,'FontSize',20)

if figureson == 1
    print -dpng volumes.png 
end 

figure(3)
clf
hold on 
plot(t,P_RV,t,P_PA,t,P_SV,t,P_RA,'linewidth',2)
legend('P_{RV}','P_{PA}','P_{SV}','P_{RA}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)') 
set(gca,'FontSize',20)

if figureson == 1
    print -dpng pressures.png 
end 

figure(4)
clf
hold on 
plot(t,y(:,5:7),'linewidth',2)
title('SL')
legend('LV','SEP','RV')

figure(5)
clf
hold on 
plot(t,P_LV,t,P_Ao,t,P_SA,t,P_PV,t,P_LA,'linewidth',2)
legend('P_{LV}','P_{Ao}','P_{SA}','P_{PV}','P_{LA}')

figure(6)
clf
hold on 
h1 = plot(V_LV,P_LV,'b','linewidth',2);
h2 = plot(V_RV,P_RV,'r','linewidth',2);
set(gca,'Xlim',[0 200])
legend([h1 h2],'LV','RV')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)') 
set(gca,'FontSize',20)

if figureson == 1
    print -dpng pvloops.png 
end 

figure(7)
clf
hold on 
plot(t,Q_t,'b',t,Q_m,'r','linewidth',2)
legend('Tricuspid','Mitral')
xlabel('Flow (mL s^{-1})')
ylabel('Time (s)') 
set(gca,'FontSize',20)

figure(8)
clf
hold on 
plot(t,V_RA,t,V_LA,'linewidth',2)
title('Atria volumes')

SV = max(V_LV) - min(V_LV)
EF = SV/max(V_LV)
CO = SV*HR
SP = max(P_Ao)
DP = min(P_Ao)
max(V_LV)
