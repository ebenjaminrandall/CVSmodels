%% This code is the main driver for the cardiovascular mechanics model 
clear all

figure(10)
clf

%% Load pseudodata 

Vtot = 5000; 

%% Volume loading 

vfactor = [.8:.05:1.2];

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

data.Vtot = Vtot; 
data.theta = theta; 
data.HR = HR; 
data.stim_period = stim_period; 

%% Parameters 

EDV_LV = zeros(size(vfactor)); 
EDV_RV = zeros(size(vfactor)); 
EDP_LV = zeros(size(vfactor)); 
EDP_RV = zeros(size(vfactor)); 
SV     = zeros(size(vfactor)); 
CO     = zeros(size(vfactor)); 
for i = 1:length(vfactor) 
    data.vfactor = vfactor(i); 
    %data.Vtot = Vtot * vfactor(i); 

    
    [pars,data] = parameters(exercise,data); 
    
    outputs = model_sol(pars,data); 

    [~,msgid] = lastwarn; 
    
    if ~strcmp(msgid,'MATLAB:ode15s:IntegrationTolNotMet')
        
        t = outputs.time; 

        V_LV = outputs.volumes.V_LV; 
        V_RV = outputs.volumes.V_RV;

        P_LV = outputs.pressures.P_LV; 
        P_RV = outputs.pressures.P_RV; 

        [Vlv,i_Vlv] = max(V_LV); 
        [Vrv,i_Vrv] = max(V_RV);
        Plv = P_LV(i_Vlv); 
        Prv = P_RV(i_Vrv); 

        EDV_LV(i) = Vlv; 
        EDP_LV(i) = Plv; 
        EDV_RV(i) = Vrv; 
        EDP_RV(i) = Prv; 
        
        
        SV(i) = max(V_LV) - min(V_LV); 
        CO(i) = SV(i) * HR * 1e-3;    

        figure(10) 
        hold on 
        plot(V_LV,P_LV,'b')
        plot(V_RV,P_RV,'r')
        set(gca,'Xlim',[0 200])
        xlabel('Volume (mL)')
        ylabel('Pressure (mmHg)')
        set(gca,'FontSize',20)
    end 
    
    lastwarn('')

end 

print -dpng PVloopprogression.png 

%% Plot EDPVR curve 

ii = find(vfactor == 1); 

% LV
EDP_lv_normal = EDP_LV(ii); 
EDV_lv_normal = EDV_LV(ii); 

An = 28; 
Bn = 3; 

V_0  = EDV_lv_normal * (0.6 - 0.006 * EDP_lv_normal); 
V_30 = V_0 + (EDV_lv_normal - V_0) / ((EDP_lv_normal / An)^(1/Bn));
beta_lv = log(EDP_lv_normal/30) / log(EDV_lv_normal / V_30); 
alpha_lv = 30 / V_30^beta_lv; 

% RV
EDP_rv_normal = EDP_RV(ii); 
EDV_rv_normal = EDV_RV(ii); 

V_0  = EDV_rv_normal * (0.6 - 0.006 * EDP_rv_normal); 
V_30 = V_0 + (EDV_rv_normal - V_0) / ((EDP_rv_normal / An)^(1/Bn));
beta_rv = log(EDP_rv_normal/30) / log(EDV_rv_normal / V_30); 
alpha_rv = 30 / V_30^beta_rv; 

V_EDPVR = [80:150]; 
P_LV_EDPVR = alpha_lv * V_EDPVR.^beta_lv;
P_RV_EDPVR = alpha_rv * V_EDPVR.^beta_rv;
Plims = [min([P_LV_EDPVR P_RV_EDPVR]) max([P_LV_EDPVR P_RV_EDPVR])]; 
Vlims = [min(V_EDPVR) max(V_EDPVR)]; 

figure(11)
clf
hold on 
plot(EDV_LV(ii)*ones(2,1),Plims,'b:')
plot(EDV_RV(ii)*ones(2,1),Plims,'r:')
plot(V_EDPVR,P_LV_EDPVR,'b','linewidth',2)
plot(V_EDPVR,P_RV_EDPVR,'r','linewidth',2)
h1 = plot(EDV_LV,EDP_LV,'bo','MarkerSize',10);
h2 = plot(EDV_RV,EDP_RV,'ro','MarkerSize',10);
legend([h1 h2],'LV','RV')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg')
set(gca,'FontSize',20)
xlim(Vlims)
ylim(Plims)

print -dpng EDPVR.png 


xx = SV ~= 0; 

EDP_FS = [0:20];

mu = 120; 
nu = -log(1 - SV(ii)/mu)/EDP_LV(ii); 

SV_FS  = mu * (1 - exp(-nu * EDP_FS)); 

figure(12)
clf
hold on 
plot(EDP_LV(xx),SV(xx),'k*','MarkerSize',20)
plot(EDP_LV(xx),SV(xx),'k')
plot(EDP_FS,SV_FS,'b')
xlabel('LV EDP (mmHg)')
ylabel('Stroke Volume (mL)')
set(gca,'FontSize',20)

print -dpng SV.png 

figure(13)
clf
hold on 
plot(EDP_LV(xx),CO(xx),'k*','MarkerSize',20)
plot(EDP_LV(xx),CO(xx),'k')
xlabel('LV EDP (mmHg)')
ylabel('Cardiac Output (mL s^{-1})')
set(gca,'FontSize',20)

print -dpng CO.png 



return 

