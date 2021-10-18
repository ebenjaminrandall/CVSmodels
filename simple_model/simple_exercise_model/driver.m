% This code is the main driver for the cardiovascular mechanics model 
clear all

%% Load pseudodata 

Vtot = 5000; 

%% Volume loading 

vfactor = 1;

%% Exercise 

theta = 0; % exercise level, between 0 and 1

HR = 64*(1 + 1.9*theta);  % 1/sec

freq = HR/60; %Hz
stim_period = 1/freq;

% exercise factors
a = 2.50; % inotropy factor
b = 2.00; % arterial vasodilation factor
c = 2.00; % venous vasoconstriction factor
d = 3.00; % arterial vasoconstriction factor
e = 5; % pulmonary vasodilation factor
f = 5; %1 + 5*theta; % calcium factor

exercise = [a; b; c; d; e; f]; 

%% Make data input structure 

data.vfactor = vfactor; 
data.Vtot = Vtot; 
data.theta = theta; 
data.HR = HR; 
data.stim_period = stim_period; 

%% Parameters 

[pars,data] = parameters(exercise,data); 

%% run the simulation 

outputs = model_sol(pars,data); 

t = outputs.time; 

V_LA = outputs.volumes.V_LA; 
V_LV = outputs.volumes.V_LV; 
V_RA = outputs.volumes.V_RA;
V_RV = outputs.volumes.V_RV;

P_LA = outputs.pressures.P_LA; 
P_LV = outputs.pressures.P_LV; 
P_Ao = outputs.pressures.P_Ao; 
P_SA = outputs.pressures.P_SA; 
P_SV = outputs.pressures.P_SV; 
P_RA = outputs.pressures.P_RA; 
P_RV = outputs.pressures.P_RV; 
P_PA = outputs.pressures.P_PA; 
P_PV = outputs.pressures.P_PV; 

SL_LV  = outputs.lengths.SL_LV; 
SL_SEP = outputs.lengths.SL_SEP; 
SL_RV  = outputs.lengths.SL_RV; 

Q_m = outputs.flows.Q_m; 
Q_a = outputs.flows.Q_a; 
Q_t = outputs.flows.Q_t; 
Q_p = outputs.flows.Q_p; 

%% Calculation of biomarkers 
 
SV  = max(V_LV) - min(V_LV)
EF  = SV/max(V_LV)
CO  = SV*HR
SBP = max(P_Ao)
DBP = min(P_Ao)
EDV = max(V_LV)

%save CVworkspace.mat 

%% Plots 

figure(1)
clf
plot(t,V_RV,'b',t,V_LV,'r')
legend('LV','RV')
xlabel('Time (s)')
ylabel('Volume (mL)')
set(gca,'FontSize',20)
title('Ventricular volumes')

figure(2)
clf
plot(t,P_LA,t,P_LV,t,P_Ao,t,P_SA)
legend('P_{LA}','P_{LV}','P_{Ao}','P_{SA}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)

figure(3)
clf
plot(t,P_RA,t,P_RV,t,P_PA,t,P_PV,t,P_SV)
legend('P_{RA}','P_{RV}','P_{PA}','P_{PV}','P_{SV}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)

figure(4)
clf
plot(t,SL_LV,t,SL_SEP,t,SL_RV) 
legend('LV','SEP','RV')
xlabel('Time (s)')
ylabel('Sarcomere length (\mu m)')
set(gca,'FontSize',20)

figure(5) 
clf
hold on 
plot(V_LV,P_LV,'b')
plot(V_RV,P_RV,'r')
set(gca,'Xlim',[0 200])
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)

figure(6)
clf
plot(t,Q_m,'b',t,Q_a,'b--',t,Q_t,'r',t,Q_p,'r--') 
legend('Q_m','Q_a','Q_t','Q_p')
xlabel('Time (s)')
ylabel('Flow (mL s^{-1})')
set(gca,'FontSize',20)

figure(7)
clf
plot(t,V_LA,'b',t,V_RA,'r')
legend('LA','RA')
title('Atrial volumes')
xlabel('Time (s)')
ylabel('Volume (mL)')
set(gca','FontSize',20)


