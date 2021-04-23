% driver

clear all

%% Load pseudodata

% Mean values 
SPbar = 120; 
DPbar = 80; 
Pbar  = 100; %SPbar/3 + DPbar * 2/3; 

SPbar = SPbar / 7.5; % convert from mmHg to kPa 
DPbar = DPbar / 7.5; % convert from mmHg to kPa 
Pbar  = Pbar  / 7.5; % convert from mmHg to kPa 

HR    = 60  / 60;  % convert from beats per min to beats per sec 

% Total blodo volume % convert from mL to m^3
Vtot = 4500 * 1e-6; 

% Cardiac Output 
CO = Vtot / 60;     % m^3 s^(-1)

% Individual wall volumes % convert from mL to m^3
Vw_LV_and_SEP = 115 * 1e-6; 
Vw_RV         = 30  * 1e-6;  

% Individual reference surface area % convert from cm^2 to m^2
Am_LV_and_SEP = 125  * 1e-4; 
Am_RV         = 100 * 1e-4; 

% Evaluation time 
dt = 0.001; 
tspan = 0:dt:20;

% Deformation of heart wall (convert cm to m)
xm_lv0  = 4.5 * 1e-2; 
xm_sep0 = 2   * 1e-2; 
xm_rv0  = 5.5 * 1e-2; 
ym0     = 3.3 * 1e-2; 

deformation.xm_lv0  = xm_lv0; 
deformation.xm_sep0 = xm_sep0; 
deformation.xm_rv0  = xm_rv0; 
deformation.ym0     = ym0; 

% Make data/input structure
data.Pbar  = Pbar; 
data.SPbar = SPbar; 
data.DPbar = DPbar; 
data.HR    = HR; 
data.Vtot  = Vtot; 
data.CO    = CO; 
data.tspan = tspan; 

data.Vw_LV_and_SEP = Vw_LV_and_SEP; 
data.Vw_RV         = Vw_RV; 
data.Am_LV_and_SEP = Am_LV_and_SEP; 
data.Am_RV         = Am_RV; 
data.deformation   = deformation; 

% Global parameters
ODE_TOL = 1e-8; 
data.gpars.ODE_TOL = ODE_TOL; 

%% Get parameters 

[adjpars,fixpars] = parameters(data); 
data.fixpars      = fixpars; 

%% Solve model 

outputs = model_sol(adjpars,data); 

V_lv = outputs.volumes.V_lv; 
V_sa = outputs.volumes.V_sa; 
V_rv = outputs.volumes.V_rv; 

P_lv = outputs.pressures.P_lv; 
P_sa = outputs.pressures.P_sa; 
P_rv = outputs.pressures.P_rv; 

save nom.mat

%% Plot

beat = find(tspan >= tspan(end) - HR); 

figure(1)
clf
plot(tspan,V_lv,'b',tspan,V_rv,'r')
xlabel('Time (s)')
ylabel('Volume (mL)')
legend('V_{lv}','V_{rv}')
set(gca,'FontSize',20)

figure(2)
clf
plot(tspan,P_lv,'b',tspan,P_rv,'r')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
legend('P_{lv}','P_{rv}')
set(gca,'FontSize',20)

figure(3)
clf
plot(tspan,P_lv,'b',tspan,P_sa, 'm')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
legend('P_{lv}','P_{sa}')
set(gca,'FontSize',20)

figure(4)
clf
plot(V_lv(beat), P_lv(beat), 'b')
hold on 
plot(V_rv(beat), P_rv(beat), 'r')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend('V_{lv}','V_{rv}')
set(gca,'FontSize',20)
xlim([0 140])

figure(5) 
clf
hold on 
plot([tspan(1) tspan(end)],(SPbar * 7.5) * ones(2,1),'k:','linewidth',0.5)
plot([tspan(1) tspan(end)],(DPbar * 7.5) * ones(2,1),'k:','linewidth',0.5)
h1 = plot(tspan,P_sa, 'm'); 
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
legend(h1,'P_{sa}')
set(gca,'FontSize',20)
ylim([60 140])



