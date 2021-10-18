% This code is the main driver for the cardiovascular mechanics model 
clear all 
%close all

printon = 0; 

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


%% run the simulation 

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
o = zeros(22,length(t)); 
for i = 1:length(t) 
    [~,o(:,i)] = dXdT_cardiovascular_mechanics_exercise(t(i),y(i,:),para);
end 

xm_LV  = y(:,1); % LV heart geometry variable, cm
xm_SEP = y(:,2); % septum heart geometry variable, cmsizs
xm_RV  = y(:,3); % RV heart geometry variable, cm
ym     = y(:,4); % Heart geometry variable, cm
SL_LV  = y(:,5);
SL_SEP = y(:,6);
SL_RV  = y(:,7);
V_LV   = y(:,8); % volume LV, mL
V_RV   = y(:,9); % volume RV, mL

V_SV   = y(:,10); % volume of systemic veins
V_PV   = y(:,11); % volume of pulmonary veins
V_SA   = y(:,12); % volume of systemic arterys
V_PA   = y(:,13); % volume of pulmonary arterys
V_Ao   = y(:,14); % volume of aorta
V_RA   = y(:,15); % volume of RA
V_LA   = y(:,16); % volume of LA
V_T = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao + V_RA + V_LA;

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
sigmapas_LV = o(14,:);
sigmapas_SEP = o(15,:);
sigmapas_RV = o(16,:);

SLo_LV = o(17,:); 
SLo_SEP = o(18,:); 
SLo_RV = o(19,:); 

sigmaact_LV = o(20,:); 
sigmaact_SEP = o(21,:);
sigmaact_RV = o(22,:); 

SV = max(V_LV) - min(V_LV)
EF = SV/max(V_LV)
CO = SV*HR
SP = max(P_Ao)
DP = min(P_Ao)
EDV = max(V_LV)

save CVworkspace.mat 


%% Plotting Baseline

figure(1)
clf
plot(t,V_RV,'b',t,V_LV,'r')
legend('LV','RV')
xlabel('Time (s)')
ylabel('Volume (mL)')
set(gca,'FontSize',20)
title('Ventricular volumes')

if printon == 1
    print -dpng volumes.png
end 

figure(2)
clf
plot(t,P_LA,t,P_LV,t,P_Ao,t,P_SA)
legend('P_{LA}','P_{LV}','P_{Ao}','P_{SA}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)

if printon == 1
    print -dpng leftP.png 
end 

figure(3)
clf
plot(t,P_RA,t,P_RV,t,P_PA,t,P_PV,t,P_SV)
legend('P_{RA}','P_{RV}','P_{PA}','P_{PV}','P_{SV}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)

if printon == 1
    print -dpng rightP.png
end 

figure(4)
clf
plot(t,SL_LV,t,SL_SEP,t,SL_RV) 
legend('LV','SEP','RV')
xlabel('Time (s)')
ylabel('Sarcomere length (\mu m)')
set(gca,'FontSize',20)

if printon == 1
    print -dpng SL.png 
end 

figure(5) 
%clf
hold on 
plot(V_LV,P_LV,'b')
plot(V_RV,P_RV,'r')
set(gca,'Xlim',[0 200])
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)

if printon == 1
    print -dpng PVloop.png
end 

figure(6)
clf
plot(t,Q_m,'b',t,Q_a,'b--',t,Q_t,'r',t,Q_p,'r--') 
legend('Q_m','Q_a','Q_t','Q_p')
xlabel('Time (s)')
ylabel('Flow (mL s^{-1})')
set(gca,'FontSize',20)

if printon == 1
    print -dpng flows.png 
end 

figure(7)
clf
plot(t,V_LA,'b',t,V_RA,'r')
legend('LA','RA')
title('Atrial volumes')
xlabel('Time (s)')
ylabel('Volume (mL)')
set(gca','FontSize',20)

if printon == 1
    print -dpng atrialvolumes.png
end 


