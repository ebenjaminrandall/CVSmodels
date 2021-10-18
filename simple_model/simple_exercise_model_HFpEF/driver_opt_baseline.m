% This code is the main driver for the cardiovascular mechanics model 
clear all 

%% Exercise input 

theta = 0; % exercise level, between 0 and 1

HR   = 64*(1 + 1.9*theta);  % 1/sec
freq = HR/60; %Hz
stim_per = 1/freq;

t = 0:.01:3*stim_per; 

data.time     = t; 
data.HR       = HR; 
data.theta    = theta; 
data.stim_per = stim_per; 

% Nominal values 
a = 2.50; % inotropy factor
b = 2.00; % arterial vasodilation factor
c = 2.00; % venous vasoconstriction factor
d = 3.00; % arterial vasoconstriction factor
e = 0.50; % pulmonary vasodilation factor
f = 5;    % calcium factor

exercise0 = [a; b; c; d; e; f]; 

%% HF input 

z = 1/3;    % C_Ao and C_SA scalar
y = 1/8;    % C_PA scalar
x = 1.5;   % R_SA scalar 
w = 2;    % R_vlv scalar
v = 3;   % k_pas_lv scalar
u = 2;   % k_pas_rv scalar

HF = [z; y; x; w; v; u]; 

%% Global variables 

[pars,lo,hi] = parameters(exercise0,HF,theta);

INDMAP  = [4 5 6 8 10 14 16 19 27];
ALLPARS = pars; 

gvars.INDMAP = INDMAP; 
gvars.ALLPARS = ALLPARS; 

data.gvars = gvars; 

%% Optimization 

p0 = pars(INDMAP); 
lb = lo(INDMAP); 
ub = hi(INDMAP); 

options = optimoptions('fmincon','Display','iter','MaxIterations',40); 
fun = @(x) model_wrap_baseline(x,data); 
popt = fmincon(fun,p0,[],[],[],[],lb,ub,[],options); 

%% Solve optimized model 

pars_opt = pars; 
pars_opt(INDMAP) = popt; 
outputs  = model_sol(pars_opt,theta,data); 

V_LV = outputs.volumes.V_LV; 
V_LA = outputs.volumes.V_LA;
V_RV = outputs.volumes.V_RV;
V_RA = outputs.volumes.V_RA; 

P_LA = outputs.pressures.P_LA; 
P_LV = outputs.pressures.P_LV; 
P_Ao = outputs.pressures.P_Ao; 
P_SA = outputs.pressures.P_SA; 
P_SV = outputs.pressures.P_SV; 
P_RA = outputs.pressures.P_RA; 
P_RV = outputs.pressures.P_RV; 
P_PA = outputs.pressures.P_PA;
P_PV = outputs.pressures.P_PV;

Q_m = outputs.flows.Q_m;
Q_a = outputs.flows.Q_a; 
Q_t = outputs.flows.Q_t; 
Q_p = outputs.flows.Q_p; 

SL_LV  = outputs.lengths.SL_LV; 
SL_SEP = outputs.lengths.SL_SEP; 
SL_RV  = outputs.lengths.SL_RV; 

%% Calculations 

% Volumes
V_LV_max = max(V_LV)
V_LV_min = min(V_LV)
V_LA_max = max(V_LA)
V_LA_min = min(V_LA)

    
%Pressures 
P_SA_s = max(P_SA)
P_SA_d = min(P_SA)
P_PA_s = max(P_PA) 
P_PA_d = min(P_PA)
P_PV_m = mean(P_PV) 

%Flows 
a = findpeaks(Q_m); 
E2A = a(1)/a(2)

SV = V_LV_max - V_LV_min
EF = SV/max(V_LV)
CO = SV*HR * 1e-3 

save opt_baseline.mat


%% Plotting Baseline

figure(1)
clf
hold on 
h1 = plot(t,V_LV,'b','linewidth',2);
h2 = plot(t,V_RV,'r','linewidth',2);
legend([h1 h2],'LV','RV')
xlabel('Time (s)')
ylabel('Volume (mL)')
set(gca,'FontSize',20)

figure(2)
clf
hold on 
h1 = plot(t,P_LV,'b','linewidth',2); 
h2 = plot(t,P_Ao,'c','linewidth',2);
h3 = plot(t,P_SA,'m','linewidth',2);
h4 = plot(t,P_LA,'k','linewidth',2);
legend([h1 h2 h3 h4],'P_{LV}','P_{Ao}','P_{SA}','P_{LV}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)

figure(3)
clf
plot(t,SL_LV,t,SL_SEP,t,SL_RV)
legend('LV','SEP','RV')
xlabel('Time (s)')
ylabel('Sarcomere length (\mu m)')
set(gca,'FontSize',20)

figure(4)
clf
hold on 
h1 = plot(t,P_SV,'k','linewidth',2); 
h2 = plot(t,P_RV,'r','linewidth',2);
h3 = plot(t,P_PA,'g','linewidth',2);
h4 = plot(t,P_PV,'m','linewidth',2);
h5 = plot(t,P_RA,'b','linewidth',2);
legend([h1 h2 h3 h4 h5],'P_{SV}','P_{RV}','P_{PA}','P_{PV}','P_{RA}')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)

figure(5)
clf
hold on 
h1 = plot(V_LV,P_LV,'b','linewidth',2);
h2 = plot(V_RV,P_RV,'r','linewidth',2);
set(gca,'Xlim',[0 200])
legend([h1 h2],'LV','RV')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
set(gca,'FontSize',20)

figure(6)
clf
hold on 
h1 = plot(t,Q_m,'b','linewidth',2); 
h2 = plot(t,Q_a,'b:','linewidth',2);
h3 = plot(t,Q_t,'r','linewidth',2);
h4 = plot(t,Q_p,'r:','linewidth',2);
legend([h1 h2 h3 h4],'Q_m','Q_a','Q_t','Q_p')
xlabel('Time (s)')
ylabel('Flow (mL s^{-1})')
set(gca,'FontSize',20)

figure(7)
clf
hold on 
h1 = plot(t,V_LA,'b','linewidth',2);
h2 = plot(t,V_RA,'r','linewidth',2);
legend([h1 h2],'LA','RA')
xlabel('Time (s)')
ylabel('Volume (mL)')
set(gca,'FontSize',20)
