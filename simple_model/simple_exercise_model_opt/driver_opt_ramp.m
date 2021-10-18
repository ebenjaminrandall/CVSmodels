% Ramps up HR from baseline to 100% exercise for normal 

%% Set number of steady state heart periods

cycles = 5; 

%% Heart rate ramp 

theta = 0:0.1:1; 
HR = 64*(1 + 1.9*theta);  % 1/sec

%% Exercise input 

a = 10; 
b = 2.5; 
c = 3.5; 
d = 5; 
e = 1; 
f = 2.5; 

exercise = [a; b; c; d; e; f]; 

%% Optimized parameters

pars_opt = [4.38202663467388;3.63758615972639;3.33220451017520;3.96452244030183;3.90731493124816;4.85015526551576;-0.430782916092454;-0.101966371056364;5.85793315448346;1.31942850715666;3.21887582486820;-6.21460809842219;-3.72970144863419;-4.67518438663743;-6.21460809842219;0.101174582643829;-2.99573227355399;-3.72970144863419;-3.25796810133788;-2.59026716544583;-1.89711998488588;0.641853886172395;-3.21887582486820;0.412109650826833;1.94591014905531;3.09104245335832;6.37386430009458;-3.91202300542815]; 

%% Compute ramp 

% Pre-allocate vectors 
SV  = zeros(size(theta)); 
EF  = zeros(size(theta)); 
CO  = zeros(size(theta)); 
SBP = zeros(size(theta)); 
MBP = zeros(size(theta)); 
DBP = zeros(size(theta)); 
PCW = zeros(size(theta)); 
CPO = zeros(size(theta)); 

for i = 1:length(theta)
    
    freq     = HR(i)/60; %Hz
    stim_per = 1/freq;

    t = 0:.01:cycles*stim_per; 

    data.time     = t; 
    data.HR       = HR(i); 
    data.stim_per = stim_per; 
    
    pars = parameters_opt(exercise,pars_opt,theta(i));
    outputs = model_sol(pars,theta(i),data); 
    
    V_LV = outputs.volumes.V_LV; 
    P_LV = outputs.pressures.P_LV; 
    P_SA = outputs.pressures.P_SA; 
    P_PV = outputs.pressures.P_PV; 
    
    V_LV_dias = max(V_LV); 
    V_LV_syst = min(V_LV); 
    
    SV(i) = V_LV_dias - V_LV_syst; 
    EF(i) = SV(i) / V_LV_dias; 
    CO(i) = SV(i) * HR(i); 
    SBP(i) = max(P_SA); 
    MBP(i) = mean(P_SA); 
    DBP(i) = min(P_SA); 
    PCW(i) = mean(P_PV); 
    CPO(i) = trapz(P_LV,V_LV) / 7.5 * 1e-3 * HR(i)/60 / cycles;
    
end 

figure(1); clf; axes('position',[0.15 0.15 0.75 0.75]); hold on; box on;
h1 = plot(HR,SBP,'k--','linewidth',1.5);
h2 = plot(HR,MBP,'k--','linewidth',1.5);
h3 = plot(HR,DBP,'k--','linewidth',1.5);
plot(HR(1),120,'ko','MarkerSize',10,'linewidth',1.5,'markerfacecolor',[1 1 1])
plot(HR(1),80,'ko','MarkerSize',10,'linewidth',1.5,'markerfacecolor',[1 1 1])
plot(HR(8),194,'ko','MarkerSize',10,'linewidth',1.5,'markerfacecolor',[1 1 1])
plot(HR(8),81,'ko','MarkerSize',10,'linewidth',1.5,'markerfacecolor',[1 1 1])
xlabel('HR (bpm)','interpreter','latex','fontsize',16)
ylabel('Pressure (mmHg)','interpreter','latex','fontsize',16)
legend([h1 h2 h3],'SBP','MBP','DBP','location','northwest')
set(gca,'FontSize',16)
 
figure(2); clf; axes('position',[0.15 0.15 0.75 0.75]); hold on; box on;
plot(HR,CO * 1e-3,'k--','linewidth',1.5)
plot(HR(1),5.76,'ko','MarkerSize',10,'linewidth',1.5,'markerfacecolor',[1 1 1])
plot(HR(8),16.2,'ko','MarkerSize',10,'linewidth',1.5,'markerfacecolor',[1 1 1])
plot(HR(11),20.6,'ko','MarkerSize',10,'linewidth',1.5,'markerfacecolor',[1 1 1])
xlabel('HR (bpm)','interpreter','latex','fontsize',16)
ylabel('Cardiac output (L$\cdot$min$^{-1}$)','interpreter','latex','fontsize',16)
set(gca,'FontSize',16)
axis([50 200 0 25]);
 
figure(3); clf; axes('position',[0.15 0.15 0.75 0.75]); hold on; box on;
h1 = plot(HR,PCW,'k--','linewidth',1.5);
plot(HR(1),8,'ko','MarkerSize',10,'linewidth',1.5,'markerfacecolor',[1 1 1])
plot(HR(8),19.5,'ko','MarkerSize',10,'linewidth',1.5,'markerfacecolor',[1 1 1])
xlabel('HR (bpm)','interpreter','latex','fontsize',16)
ylabel('Pulmonary Venous Pressure (mmHg)','interpreter','latex','fontsize',16)
set(gca,'FontSize',16)
 
figure(4); clf; axes('position',[0.15 0.15 0.75 0.75]); hold on; box on;
plot(HR,CPO,'k--','linewidth',1.5);
xlabel('HR (bpm)','interpreter','latex','fontsize',16)
ylabel('Cardiac power (W)','interpreter','latex','fontsize',16)
set(gca,'FontSize',16)

return 

%% Plots 

figure(10)
clf
hold on 
plot(HR(1),120,'k*','MarkerSize',10)
plot(HR(1),80,'k*','MarkerSize',10)
plot(HR(8),194,'k*','MarkerSize',10)
plot(HR(8),81,'k*','MarkerSize',10) 
h1 = plot(HR,SBP,'b--','linewidth',2);
h2 = plot(HR,MBP,'k--','linewidth',2);
h3 = plot(HR,DBP,'r--','linewidth',2);
xlabel('HR (bpm)')
ylabel('Pressure (mmHg)')
legend([h1 h2 h3],'SBP','MBP','DBP')
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
legend(h1,'PCW')
set(gca,'FontSize',20)

print -dpng PCW_exercise.png

figure(13)
clf
hold on
plot(HR,CPO,'k--','linewidth',2);
xlabel('HR (bpm)')
ylabel('Cardiac power (W)')
set(gca,'FontSize',20)

print -dpng CPO_exercise.png






