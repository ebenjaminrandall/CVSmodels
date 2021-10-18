

load CVworkspace.mat 

%% Original Triseg passive curve 

%Amref = 0.975*80; 

k_pas = 22 ; 
L_0 = 1.6 ; 
L_C = 1.8 ;
g_2 = 3; 
g_1 = 22 ; 

k_act = 7.5*96; 
SLrest = 1.51; 
LSEiso = 0.04; 


SLo = [1.8:.05:2.2] ; 

xm = -1.96 * SLo - 0.26; 

%xm   = -4.6:.01:-3.8; 
xm_s = -0.5 * xm + 0.15; % Plotted xm_LV vs xm_SEP and got approximately this line 
y   = -0.8 * xm - 0.05; % Plotted xm_LV vs ym and got approximately this line

figure(110)
clf
hold on 
plot(xm_LV,xm_SEP,'b')
plot(xm,xm_s,'m','linewidth',2)
plot(xm_LV,ym,'b')
plot(xm,y,'m','linewidth',2)
title('xm_{SEP} and ym')

Vm = (pi/6) .* xm .* (xm.^2 + 3 * y.^2);
Vm_s = (pi/6) .* xm_s .* (xm_s.^2 + 3 * y.^2); 
Am = pi * (xm.^2 + y.^2); 
Cm = 2 * xm ./ (xm.^2 + y.^2); 
z  =  3 .* Cm .* Vw_l ./ (2 .* Am); 

%epsf = (1/2)*log(Am/Amref) - (1/12).*z.^2 - 0.019.*z.^4; 
epsf = (1/2)*log(Am) - (1/12).*z.^2 - 0.019.*z.^4; 

% SLo1 = Lsref * exp(epsf); 
% 
% figure(130) 
% plot(SLo,SLo1)
% return 

sigmapas1 = k_pas * (SLo  - L_0)  + g_1 * max(0,SLo - L_C).^g_2;

V_D = -1/2 * Vw_l - 1/2 * Vw_s + Vm_s - Vm; 


CC = -0.64 * (SLo - 1.75).^0.28 + 0.45 * (SLo - 1) ; %Used curve fitting tool to fit diastolic curve of activation function  

figure(111)
clf
hold on 
plot(SLo_LV,C_LV,'b')
plot(SLo,CC,'m','linewidth',2)
title('SL vs activation function')

SL = SLo - 0.04; % Plotted SLo_LV vs SL_LV and got this line 

figure(112)
clf
hold on 
plot(SLo_LV,SL_LV,'b')
plot(SLo,SL,'m','linewidth',2)
title('SL follower')

sigmaact = k_act .* CC .* (SL - SLrest) .* (SLo - SL) / LSEiso; 
sigma = sigmapas1 + sigmaact; 
Tm = (Vw_l .* sigma ./ (2 .* Am)) .* (1 + (1/3) * z.^2 + (1/5) * z.^4); 
Tx = Tm .* 2 .* xm .* y ./ (xm.^2  + y.^2); 
P_D = - 2 * Tx ./ y; 

figure(113)
clf
hold on 
plot(SLo_LV,sigmapas_LV,'b')
plot(SLo,sigmapas1,'m','linewidth',2) 
title('\sigma_{pas}')

figure(114) 
clf
hold on 
plot(SLo_LV,V_LV,'b')
plot(SLo,V_D,'m','linewidth',2)
title('SL vs V')

figure(115)
clf
hold on 
plot(SLo_LV,P_LV,'b')
plot(SLo,P_D,'m','linewidth',2)
title('SL vs P')

figure(116)
clf
hold on
plot(V_LV,P_LV,'b')
plot(V_D,P_D,'m','linewidth',2) 
title('PV loop')

figure(117)
clf
hold on 
plot(SLo_LV,sigmaact_LV,'b')
plot(SLo,sigmaact,'m','linewidth',2)
title('\sigma_{act}')



% return 
%  
% 
% %% Scaled passive curve parameters to Amref to predict the same passive curve
% 
% epsf_hat = (1/2)*log(Am) - (1/12).*z.^2 - 0.019.*z.^4; 
% 
% Lsref_hat = Lsref /sqrt(Amref); 
% 
% SLo_hat = Lsref_hat * exp(epsf_hat); 
% 
% % k_pas_hat = k_pas / sqrt(Amref); 
% % L_0_hat = L_0 * sqrt(Amref); 
% % g_1_hat = g_1 / sqrt(Amref)^g_2;
% % L_C_hat = L_C * sqrt(Amref); 
% 
% %sigmapas2 = k_pas_hat * (SLo_hat  - L_0_hat)  + g_1_hat * max(0,SLo_hat - L_C_hat).^g_2;
% sigmapas2 = k_pas * (SLo_hat  - L_0)  + g_1 * max(0,SLo_hat - L_C).^g_2;
% 
% 
% figure(113)
% hold on  
% plot(SLo_hat,sigmapas2,'r')








