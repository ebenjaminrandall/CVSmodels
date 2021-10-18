
load CVworkspace.mat 


% Normal EDPVR relationship 
EDP = 5; 
EDV = 130; 

%% Find EDPVR curve 

V_D = 1 * EDV;

An = 28; 
Bn = 3; 

V_0  = EDV * (0.6 - 0.006 * EDP); 
V_30 = V_0 + (EDV - V_0) / ((EDP / An)^(1/Bn));
beta = log(EDP/30) / log(EDV / V_30); 
alpha = 30 / V_30^beta; 

P_D = alpha * V_D^beta;

disp('Assigned V_D = ')
disp(V_D)
disp('Assigned P_D = ')
disp(P_D) 

%% Parameters 

% In diastole, SL is greatest 
SLo_D = 2.2; 

% Parameters 
Vw_l  = 80; 
Vw_s  = 38; 

Lsref = 1.9; 
k_act = 7.5*96; 
SLrest = 1.51; 
LSEiso = 0.04; 

k_pas = 22 ; 
L_0 = 1.6 ; 
L_C = 1.8 ;
g_2 = 3; 
g_1 = 22/k_pas; 
C_rest = 0.02; 

%% New Amref 

xm_l = -0.009 * V_D - 3.32; %-4.5; %
y    = -0.8 * xm_l - 0.05; %3.5; %

Vm_l = (pi/6) .* xm_l .* (xm_l.^2 + 3 * y.^2);
Am_l = pi * (xm_l.^2 + y.^2); 
Cm_l = 2 * xm_l ./ (xm_l.^2 + y.^2); 
z_l  = 3 .* Cm_l .* Vw_l ./ (2 .* Am_l); 

epsf_l = (1/2) * log(Am_l) - (1/12) .* z_l.^2 - 0.019 .* z_l.^4; 

Amref = (Lsref / SLo_D * exp(epsf_l)).^2

%% New kpas 

sinalpha_l = 2 .* xm_l .* y ./ (xm_l.^2 + y.^2); 
Gamma_l    = (Vw_l ./ (2 .* Am_l)) .* (1 + (1/3) .* z_l.^2 + (1/5) .* z_l.^4); 
sigmaact_l = k_act * C_rest * (SLo_D - SLrest) * 1;  
sigmapas_l = (SLo_D - L_0) + g_1 * (SLo_D - L_C).^g_2; 
k_pas = (P_D .* y ./ (-2 .* sinalpha_l .* Gamma_l) - sigmaact_l) ./ sigmapas_l 

%% Triseg model with new Amref 

%Amref = 0.975*80; 


%kpas = 22; 

%SLo = [1.8:.05:2.2] ; 

%xm = -1.96 * SLo - 0.26; 

xm   = [-.5:.01:.5] + xm_l; 
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

epsf = (1/2)*log(Am/Amref) - (1/12).*z.^2 - 0.019.*z.^4; 
%epsf = (1/2)*log(Am) - (1/12).*z.^2 - 0.019.*z.^4; 

SLo = Lsref * exp(epsf); 

xx = find(SLo >= 2.2,1,'last'); 
if isempty(xx)
    [~,xx] = max(SLo); 
end 
SLo22 = SLo(xx); 

% 
% figure(130) 
% plot(SLo,SLo1)
% return 

sigmapas1 = k_pas * ((SLo  - L_0)  + g_1 * max(0,SLo - L_C).^g_2);

V_D = -1/2 * Vw_l - 1/2 * Vw_s + Vm_s - Vm; 


CC = -0.64 * (SLo - 1.75).^0.28 + 0.45 * (SLo - 1) ; %Used curve fitting tool to fit diastolic curve of activation function  

figure(111)
clf
hold on 
plot(SLo22*ones(2,1),[min(C_LV) max(C_LV)],'k:')
plot(SLo_LV,C_LV,'b')
plot(SLo,CC,'m','linewidth',2)
title('SL vs activation function')

SL = SLo - 0.04; % Plotted SLo_LV vs SL_LV and got this line 

figure(112)
clf
hold on 
plot(SLo22*ones(2,1),[min(SL) max(SL)],'k:')
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
plot(SLo22*ones(2,1),[min(sigmapas1) max(sigmapas1)],'k:')
plot(SLo_LV,sigmapas_LV,'b')
plot(SLo,sigmapas1,'m','linewidth',2) 
title('\sigma_{pas}')

figure(114) 
clf
hold on 
plot(SLo22*ones(2,1),[min(V_D) max(V_D)],'k:')
plot(SLo_LV,V_LV,'b')
plot(SLo,V_D,'m','linewidth',2)
title('SL vs V')

figure(115)
clf
hold on 
plot(SLo22*ones(2,1),[min(P_LV) max(P_LV)],'k:')
plot(SLo_LV,P_LV,'b')
plot(SLo,P_D,'m','linewidth',2)
title('SL vs P')

figure(116)
clf
hold on
plot(V_D(xx)*ones(2,1),[min(P_LV) max(P_LV)],'k:')
plot(V_LV,P_LV,'b')
plot(V_D,P_D,'m','linewidth',2) 
title('PV loop')

figure(117)
clf
hold on 
plot(SLo22*ones(2,1),[min(sigmaact_LV) max(sigmaact_LV)],'k:')
plot(SLo_LV,sigmaact_LV,'b')
plot(SLo,sigmaact,'m','linewidth',2)
title('\sigma_{act}')

figure(118) 
%clf
hold on 
plot([V_0:V_30],alpha*[V_0:V_30].^beta,'b')
plot(V_D(xx),P_D(xx),'r*','MarkerSize',20)


disp('Model V_D = ')
disp(V_D(xx))
disp('Model P_D = ')
disp(P_D(xx)) 
