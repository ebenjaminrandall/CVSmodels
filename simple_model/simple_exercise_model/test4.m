
load ../exercise_codes/CVworkspace.mat 


% Normal EDPVR relationship 
EDP = 5; 
EDV = 140; 

%% Find EDPVR curve 

V_D_l = 1 * EDV;
V_D_r = V_D_l; 

An = 28; 
Bn = 3; 

V_0  = EDV * (0.6 - 0.006 * EDP); 
V_30 = V_0 + (EDV - V_0) / ((EDP / An)^(1/Bn));
beta = log(EDP/30) / log(EDV / V_30); 
alpha = 30 / V_30^beta; 

P_D_l = alpha * V_D_l^beta
P_D_r = P_D_l / 6

disp('Assigned V_D = ')
disp(V_D_l)
disp('Assigned P_D = ')
disp(P_D_l) 

%% Parameters 

% In diastole, SL is greatest 
SLo_D_l = .0046 * V_D_l + 1.55; %2.2; 

% Parameters 
Vw_l  = 80; 
Vw_s  = 38; 
Vw_r  = 28; 

Lsref = 1.9; 
k_act = 7.5*96; 
SLrest = 1.51; 
LSEiso = 0.04; 

k_pas_l = 22 ; 
L_0 = 1.6 ; 
L_C = 1.8 ;
g_2 = 3; 
g_1 = 22/k_pas_l; 
C_rest = 0.02; 

%% New Amref 

% Left 
xm_l = -0.009 * V_D_l - 3.32; %-4.5; %
y    = -0.8 * xm_l - 0.05; %3.5; %

Vm_l = (pi/6) .* xm_l .* (xm_l.^2 + 3 * y.^2);
Am_l = pi * (xm_l.^2 + y.^2); 
Cm_l = 2 * xm_l ./ (xm_l.^2 + y.^2); 
z_l  = 3 .* Cm_l .* Vw_l ./ (2 .* Am_l); 

epsf_l = (1/2) * log(Am_l) - (1/12) .* z_l.^2 - 0.019 .* z_l.^4; 

Amref_l = (Lsref / SLo_D_l * exp(epsf_l)).^2

% Right 
xm_r = 0.015 * V_D_r + 3.87; %-4.5; %

Vm_r = (pi/6) .* xm_r .* (xm_r.^2 + 3 * y.^2);
Am_r = pi * (xm_r.^2 + y.^2); 
Cm_r = 2 * xm_r ./ (xm_r.^2 + y.^2); 
z_r  = 3 .* Cm_r .* Vw_r ./ (2 .* Am_r); 

epsf_r = (1/2) * log(Am_r) - (1/12) .* z_r.^2 - 0.019 .* z_r.^4; 

Amref_r = (Lsref / SLo_D_l * exp(epsf_r)).^2

%% New kpas 

% Left 
sinalpha_l = 2 .* xm_l .* y ./ (xm_l.^2 + y.^2); 
Gamma_l    = (Vw_l ./ (2 .* Am_l)) .* (1 + (1/3) .* z_l.^2 + (1/5) .* z_l.^4); 
sigmaact_l = k_act * C_rest * (SLo_D_l - SLrest) * 1;  
sigmapas_l = (SLo_D_l - L_0) + g_1 * (SLo_D_l - L_C).^g_2; 
k_pas_l    = (P_D_l .* y ./ (-2 .* sinalpha_l .* Gamma_l) - sigmaact_l) ./ sigmapas_l 

% Right
sinalpha_r = 2 .* xm_r .* y ./ (xm_r.^2 + y.^2); 
Gamma_r    = (Vw_r ./ (2 .* Am_r)) .* (1 + (1/3) .* z_r.^2 + (1/5) .* z_r.^4); 
sigmaact_r = k_act * C_rest * (SLo_D_l - SLrest) * 1;  
sigmapas_r = (SLo_D_l - L_0) + g_1 * (SLo_D_l - L_C).^g_2; 
k_pas_r    = (P_D_r .* y ./ (2 .* sinalpha_r .* Gamma_r) - sigmaact_r) ./ sigmapas_r 

%% Triseg model with new Amref 

%Amref = 0.975*80; 


%kpas = 22; 

%SLo = [1.8:.05:2.2] ; 

%xm = -1.96 * SLo - 0.26; 

xm   = [-.5:.001:.75] + xm_l; 
xm_s = -0.5 * xm + 0.15; % Plotted xm_LV vs xm_SEP and got approximately this line 
xm_r = -1.57 * xm - 1.4; 
y   = -0.8 * xm - 0.05; % Plotted xm_LV vs ym and got approximately this line

figure(110)
clf
hold on 
plot(xm_LV,xm_SEP,'b','linewidth',2)
h1 = plot(xm,xm_s,'m','linewidth',2);
plot(xm_LV,ym,'b','linewidth',2)
h2 = plot(xm,y,'g','linewidth',2);
xlabel('xm_{LV}')
ylabel('Displacement (\mu m)')
legend([h1 h2],'xm_{SEP}','ym')
text(-4.8,4,'ym = -0.8 xm_{LV} - 0.05','FontSize',20)
text(-4.8,2.75,'xm_{SEP} = -0.5 xm_{LV} + 0.15','FontSize',20)
set(gca,'FontSize',20)

print -dpng displacements.png 

%%

Vm = (pi/6) .* xm .* (xm.^2 + 3 * y.^2);
Vm_s = (pi/6) .* xm_s .* (xm_s.^2 + 3 * y.^2); 
Vm_r = (pi/6) .* xm_r .* (xm_r.^2 + 3 * y.^2); 
Am = pi * (xm.^2 + y.^2); 
Cm = 2 * xm ./ (xm.^2 + y.^2); 
z  =  3 .* Cm .* Vw_l ./ (2 .* Am); 

epsf = (1/2)*log(Am/Amref_l) - (1/12).*z.^2 - 0.019.*z.^4; 
%epsf = (1/2)*log(Am) - (1/12).*z.^2 - 0.019.*z.^4; 

SLo = Lsref * exp(epsf); 



% 
% figure(130) 
% plot(SLo,SLo1)
% return 

sigmapas1 = k_pas_l * ((SLo  - L_0)  + g_1 * max(0,SLo - L_C).^g_2);

V_D = -1/2 * Vw_l - 1/2 * Vw_s + Vm_s - Vm; 
V_D_RV = -1/2 * Vw_r - 1/2 * Vw_s - Vm_s + Vm_r; 


xx = find(V_D >= max(V_LV),1,'last');  
SLo_at_EDV = SLo(xx); 


CC = -0.64 * (SLo - 1.75).^0.28 + 0.45 * (SLo - 1) ; %Used curve fitting tool to fit diastolic curve of activation function  

figure(111)
clf
hold on 
plot(SLo_at_EDV*ones(2,1),[min(C_LV) max(C_LV)],'k:')
h1 = plot(SLo_LV,C_LV,'b','linewidth',2);
h2 = plot(SLo,CC,'m','linewidth',2);
legend([h1 h2],'C_{LV} from Triseg','Approximation')
xlabel('Sarcomere length (\mu m)') 
ylabel('Activation function')
set(gca,'FontSize',20)

print -dpng activationfunction.png 

SL = SLo - 0.04; % Plotted SLo_LV vs SL_LV and got this line 

figure(112)
clf
hold on 
plot(SLo_at_EDV*ones(2,1),[min(SL) max(SL)],'k:')
h1 = plot(SLo_LV,SL_LV,'b','linewidth',2);
h2 = plot(SLo,SL,'m','linewidth',2);
text(2,1.9,'SL = SLo - 0.04','FontSize',20)
xlabel('Sarcomere length (\mu m)')
ylabel('Sarcomere length (\mu m)') 
legend([h1 h2],'SL follower','Approximation')
set(gca,'FontSize',20)

print -dpng SLfollower.png 

%%

sigmaact = k_act .* CC .* (SL - SLrest) .* (SLo - SL) / LSEiso; 
sigma = sigmapas1 + sigmaact; 
Tm = (Vw_l .* sigma ./ (2 .* Am)) .* (1 + (1/3) * z.^2 + (1/5) * z.^4); 
Tx = Tm .* 2 .* xm .* y ./ (xm.^2  + y.^2); 
P_D = - 2 * Tx ./ y; 

figure(113)
clf
hold on 
plot(SLo_at_EDV*ones(2,1),[min(sigmapas1) max(sigmapas1)],'k:')
h1 = plot(SLo_LV,sigmapas_LV,'b','linewidth',2);
h2 = plot(SLo,sigmapas1,'m','linewidth',2);
legend([h1 h2],'\sigma_{pas,LV} from Triseg','Approximation')
xlabel('Sarcomere length (\mum)')
ylabel('Stress (kPa)')
set(gca,'FontSize',20)

print -dpng sigmapas_LV.png 

%% 

V_Dl = 50:150; 

m = (2.2 - 1.85)/(142.1 - 65.5); 
b = 2.2 - m * 142.1; 

figure(114) 
clf
hold on 
%plot([min(V_D) max(V_D)],SLo22*ones(2,1),'k:')
h1 = plot(V_LV,SLo_LV,'b','linewidth',2);
%h2 = plot(V_D,SLo,'m','linewidth',2);
h2 = plot(V_Dl,m * V_Dl + b,'m','linewidth',2);
%h3 = plot(V_RV,SLo_RV,'r','linewidth',2);
h3 = plot(V_Dl,.0046 * V_Dl + 1.55,'g','linewidth',2);  

xlabel('Volume (mL)')
ylabel('Sarcomere length (\mum)')
set(gca,'FontSize',20)
legend([h1 h2],'SL_{LV} from Triseg','Approximation','location','northwest')

%print -dpng SLvsV.png 

figure(115) 
clf
hold on 
%plot([min(V_D) max(V_D)],SLo22*ones(2,1),'k:')
h1 = plot(V_RV,SLo_RV,'b','linewidth',2);
%h2 = plot(V_D,SLo,'m','linewidth',2);
h2 = plot(V_Dl,m * V_Dl + b,'m','linewidth',2);
%h3 = plot(V_RV,SLo_RV,'r','linewidth',2);
h3 = plot(V_Dl,.0046 * V_Dl + 1.55,'g','linewidth',2);  

xlabel('Volume (mL)')
ylabel('Sarcomere length (\mum)')
set(gca,'FontSize',20)
legend([h1 h2],'SL_{LV} from Triseg','Approximation','location','northwest')

%print -dpng SLvsV.png 

%% 

figure(116)
clf
hold on 
plot(SLo_at_EDV*ones(2,1),[min(P_LV) max(P_LV)],'k:')
h1 = plot(SLo_LV,P_LV,'b','linewidth',2);
h2 = plot(SLo,P_D,'m','linewidth',2);
xlabel('Sarcomere length (\mum)')
ylabel('Pressure (mmHg)')
legend([h1 h2],'P_{LV}','Approximation during diastole')
set(gca,'FontSize',20)

print -dpng PvsSL.png 

%%

figure(117)
clf
hold on
plot(V_D(xx)*ones(2,1),[min(P_LV) max(P_LV)],'k:')
plot(V_LV,P_LV,'b','linewidth',2)
plot(V_D,P_D,'m','linewidth',2) 
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend([h1 h2],'Model','Approximation during diastole')
set(gca,'FontSize',20)

print -dpng PVloop.png 

%% 

figure(118)
clf
hold on 
plot(SLo_at_EDV*ones(2,1),[min(sigmaact_LV) max(sigmaact_LV)],'k:')
h1 = plot(SLo_LV,sigmaact_LV,'b','linewidth',2); 
h2 = plot(SLo,sigmaact,'m','linewidth',2); 
legend([h1 h2],'\sigma_{act,LV} from Triseg','Approximationg during diastole')
xlabel('Sarcomere length (\mum)')
ylabel('Stress (kPa)')
set(gca,'FontSize',20)

print -dpng sigmaact_LV.png 


%%

figure(119) 
%clf
hold on 
plot([V_0:V_30],alpha*[V_0:V_30].^beta,'b')
plot(V_D(xx),P_D(xx),'r*','MarkerSize',20)


[mm,ii] = max(V_LV); 

disp('Model V_LV = ')
disp(mm)
disp('Model P_LV = ')
disp(P_LV(ii)) 


disp('Nominal V_D = ')
disp(V_D(xx))
disp('Nominal P_D = ')
disp(P_D(xx)) 

