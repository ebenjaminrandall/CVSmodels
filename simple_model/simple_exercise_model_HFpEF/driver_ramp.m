% Ramps up HR from baseline to 100% exercise for normal 

load opt_exercise.mat 

theta = 0:0.1:1; 
HR = 64*(1 + 1.9*theta);  % 1/sec

a = 10; 
b = 2.5; 
c = 3.5; 
d = 5; 
e = 1; 
f = 2.5; 

exercise = [a; b; c; d; e; f]; 

% Pre-allocate vectors 
SV  = zeros(size(theta)); 
EF  = zeros(size(theta)); 
CO  = zeros(size(theta)); 
SBP = zeros(size(theta)); 
DBP = zeros(size(theta)); 
PCW = zeros(size(theta)); 

for i = 1:length(theta)
    
    freq  = HR(i)/60; %Hz
    stim_per = 1/freq;

    t = 0:.01:5*stim_per; 

    data.time     = t; 
    data.HR       = HR(i); 
    data.stim_per = stim_per; 
    
    pars = parameters_opt(exercise,pars_opt,theta(i));
    outputs = model_sol(pars,theta(i),data); 
    
    V_LV = outputs.volumes.V_LV; 
    P_SA = outputs.pressures.P_SA; 
    P_PV = outputs.pressures.P_PV; 
    
    V_LV_dias = max(V_LV); 
    V_LV_syst = min(V_LV); 
    
    SV(i) = V_LV_dias - V_LV_syst; 
    EF(i) = SV(i) / V_LV_dias; 
    CO(i) = SV(i) * HR(i); 
    SBP(i) = max(P_SA); 
    DBP(i) = min(P_SA); 
    PCW(i) = mean(P_PV); 
    
end 

%% Plots 

figure(10)
clf
hold on 
plot(HR(1),120,'k*','MarkerSize',10)
plot(HR(1),80,'k*','MarkerSize',10)
plot(HR(8),194,'k*','MarkerSize',10)
plot(HR(8),81,'k*','MarkerSize',10) 
h1 = plot(HR,SBP,'b--','linewidth',2);
h2 = plot(HR,DBP,'r--','linewidth',2);
xlabel('HR (bpm)')
ylabel('Pressure (mmHg)')
legend([h1 h2],'SBP','DBP')
set(gca,'FontSize',20)

print -dpng SBP_DBP_exercise.png 

figure(11)
clf
hold on 
plot(HR(1),5.76,'k*','MarkerSize',10)
plot(HR(8),16.2,'k*','MarkerSize',10)
plot(HR(11),20.6,'k*','MarkerSize',10)
plot(HR,CO * 1e-3,'k--','linewidth',2)
xlabel('HR (bpm)')
ylabel('Cardiac output (mL s^{-1})')
set(gca,'FontSize',20)

print -dpng CO_exercise.png

figure(12)
clf
hold on
plot(HR(1),8,'k*','MarkerSize',10)
plot(HR(8),19.5,'k*','MarkerSize',10)
h1 = plot(HR,PCW,'k--','linewidth',2);
xlabel('HR (bpm)')
ylabel('Pressure (mmHg)')
legend('PCW')
set(gca,'FontSize',20)

print -dpng PCW_exercise.png



