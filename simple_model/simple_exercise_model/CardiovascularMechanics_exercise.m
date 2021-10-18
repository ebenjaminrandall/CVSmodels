% This code is the main driver for the cardiovascular mechanics model 
clear;

theta = 0; % exercise level, between 0 and 1


HR = 64*(1 + 1.9*theta);  % 1/sec


Vtot = 5000; 

Vw_LV = 80;
Vw_SEP = 38; 
Vw_RV = 28; % Heart wall volumes (mL)
Amref_LV  = 0.975*80; % LV midwall reference surface area, cm^2
Amref_SEP = 0.975*45; % SEP midwall reference surface area, cm^2
Amref_RV  = 1.12*100; % RV midwall reference surface area, cm^2

V_LV  = .03 * Vtot; % initial V_LV (mL)  150;  % 
V_RV  = .03 * Vtot; % initial V_LV (mL) 150; %

freq = HR/60; %Hz
stim_period = 1/freq;


para = [ stim_period Vw_LV Vw_SEP Vw_RV Amref_LV Amref_SEP Amref_RV theta];

xm_LV  = -5.0;
xm_SEP = +2.5;
xm_RV  = +8.0;
ym     = +5.0;

SL_LV   = 2.2;
SL_SEP  = 2.2;
SL_RV   = 2.2;

x0 = [xm_LV ,xm_SEP ,xm_RV ,ym];

% opts = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
opts = optimset('MaxFunEvals',10000,'MaxIter',1000);
%TrisegEquations(x0,Vw_LV,Vw_SEP,Vw_RV, V_LV, V_RV);

V = [Vw_LV Vw_SEP Vw_RV];
x = fsolve(@TrisegEquations,x0,opts,V, V_LV, V_RV);


%% set up IC's and Parameters

vfactor = 1; 
% 
% V_SA =  vfactor*200; %
% V_SV =  vfactor*1175; %
% V_PA = vfactor*35; %
% V_PV = vfactor*65; %
% V_Ao = vfactor*60; %
% V_RA =  35; %100; % mL 
% V_LA =  35; %100; % mL 


V_SA = .3*vfactor * .07 * Vtot; % vfactor*200; %
V_SV = .1*vfactor * .7 * Vtot; % vfactor*1175; %
V_PA = .4*vfactor * .03 * Vtot; %vfactor*35; %
V_PV = .1*vfactor * .10 * Vtot; %vfactor*65; %
V_Ao = .3*vfactor * .03 * Vtot; % vfactor*60; %
V_RA = .005 * Vtot; % 35; %100; % mL 
V_LA = .005 * Vtot; %  35; %100; % mL 

V_T = V_LV + V_RV + V_SA + V_SV + V_PA + V_PV + V_Ao + V_RA + V_LA

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

figure(200)
plot(t,y(:,8),'b',t,y(:,9),'r')

[t,y] = ode15s(@dXdT_cardiovascular_mechanics_exercise,[0 2*stim_period],init,options,para);
o = zeros(26,length(t)); 
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

sum(y(end,8:16))

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
Q_t  = o(11,:);
sigmapas_LV = o(12,:);
sigmapas_SEP = o(13,:);
sigmapas_RV = o(14,:);
sigmaact_LV = o(15,:); 
sigmaact_SEP = o(16,:); 
sigmaact_RV = o(17,:); 
SLo_LV = o(18,:); 
SLo_SEP = o(19,:); 
SLo_RV = o(20,:); 
z_LV = o(21,:); 
z_SEP = o(22,:); 
z_RV = o(23,:); 
Vm_LV = o(24,:); 
Vm_SEP = o(25,:); 
Vm_RV = o(26,:); 

%save CVworkspace.mat 


%% Plotting Baseline
% 
figure(1); clf; plot(t,V_RV,t,V_LV); title('ventricular volumes')
figure(3); clf; plot(t,P_RV,t,P_PA,t,P_SV,t,P_RA); title('pulmonary pressures'); legend('P_{RV}','P_{PA}','P_{SV}','P_{RA}');  %set(gca,'Ylim',[0 35]);
figure(4); clf; plot(t,y(:,5:7)); title('SL'); legend('LV','SEP','RV'); % SL's
figure(5); clf; plot(t,P_LV,t,P_Ao,t,P_SA,t,P_PV,t,P_LA); legend('P_{LV}','P_{Ao}','P_{SA}','P_{PV}','P_{LA}'); %set(gca,'Ylim',[0 125]);

figure(6); 
clf
hold on 
plot(V_LV,P_LV); 
plot(V_RV,P_RV); 
title('PV Loop');
set(gca,'Xlim',[0 200]); 

% figure(7); clf; plot(t,Q_t,t,Q_m); title('Mitral and Tric. valve flows');
% figure(8); clf; plot(t,V_RA,t,V_LA); title('Atria volumes');
% % 
% figure(100)
% clf
% plot3(V_LV,P_LV,sigmaact_LV)
% 
SV = max(V_LV) - min(V_LV)
EF = SV/max(V_LV)
CO = SV*HR
SP = max(P_Ao)
DP = min(P_Ao)
max(V_LV)
