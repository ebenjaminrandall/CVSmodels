% driver

clear all

printon = 0; 

%% Load pseudodata

% Mean values 
SPbar = 120; 
DPbar = 80; 
Pbar  = 100; %SPbar/3 + DPbar * 2/3; 

SPbar = SPbar / 7.5; % convert from mmHg to kPa 
DPbar = DPbar / 7.5; % convert from mmHg to kPa 
Pbar  = Pbar  / 7.5; % convert from mmHg to kPa 

HR    = 75;  % convert from beats per min to beats per sec 

% Total blodo volume % convert from mL to m^3
Vtot = 4500 * 1e-6; 

% Cardiac Output 
CO = Vtot / 60;     % m^3 s^(-1)

% Individual wall volumes % convert from mL to m^3
Vw_LV_and_SEP = 165 * 1e-6; %115 * 1e-6;
Vw_RV         = 30 * 1e-6;  %30  * 1e-6;  

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

V_la = outputs.volumes.V_la; 
V_lv = outputs.volumes.V_lv; 
V_sa = outputs.volumes.V_sa; 
V_sv = outputs.volumes.V_sv; 
V_ra = outputs.volumes.V_ra; 
V_rv = outputs.volumes.V_rv; 
V_pa = outputs.volumes.V_pa; 
V_pv = outputs.volumes.V_pv; 

P_la = outputs.pressures.P_la; 
P_lv = outputs.pressures.P_lv; 
P_sa = outputs.pressures.P_sa; 
P_sv = outputs.pressures.P_sv; 
P_ra = outputs.pressures.P_ra; 
P_rv = outputs.pressures.P_rv; 
P_pa = outputs.pressures.P_pa; 
P_pv = outputs.pressures.P_pv; 

Q_m_valve = outputs.flows.Q_m_valve; 
Q_a_valve = outputs.flows.Q_a_valve; 
Q_sa      = outputs.flows.Q_sa; 
Q_sv      = outputs.flows.Q_sv; 
Q_t_valve = outputs.flows.Q_t_valve; 
Q_p_valve = outputs.flows.Q_p_valve; 
Q_pa      = outputs.flows.Q_pa; 
Q_pv      = outputs.flows.Q_pv; 

b = mod(tspan,round(60/HR,3)); 
beats = find(round(b,3) ==  0); 
beat = beats(end-3):beats(end-1); 

t_beat = tspan(beat) - tspan(beat(1)); 

SV = max(V_lv(beat)) - min(V_lv(beat)) % mL
EF = SV / max(V_lv(beat)) % dimensionless
CO = trapz(t_beat/60,Q_a_valve(beat))/(t_beat(end)/60 - t_beat(1)/60) %SV * HR_end * 1e-3 % L min^(-1)
CP = trapz(P_lv(beat),V_lv(beat)) / 7.5 * 1e-3 * HR/60; %mean(P_sa(beat)) / 7.5 * 1e3 * SV * 1e-6 * HR_end/60 % W 
CP = CP / 2 %average over 2 beats 

save nom.mat


%% Plot

vlims = [0 175]; 
plims = [0 145]; 

hfig1 = figure(1); 
clf
hold on 
h1 = plot(V_lv(beat), P_lv(beat), 'b','linewidth',2);
h2 = plot(V_rv(beat), P_rv(beat), 'r','linewidth',2);
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend([h1 h2],'LV','RV')
set(gca,'FontSize',25)
xlim(vlims)

hfig2 = figure(2);  
clf
hold on 
plot([tspan(1) tspan(end)],(SPbar * 7.5) * ones(2,1),'k:','linewidth',0.5)
plot([tspan(1) tspan(end)],(DPbar * 7.5) * ones(2,1),'k:','linewidth',0.5)
h1 = plot(tspan,P_sa, 'm','linewidth',2); 
xlabel('Time (s)')
ylabel('Pressure (mmHg)','linewidth',2)
legend([h1],'P_{sa}')
set(gca,'FontSize',25)
ylim([50 140])

hfig3 = figure(3);
clf 
hold on
plot(t_beat,P_lv(beat),'b','linewidth',2)
plot(t_beat,P_la(beat),'r','linewidth',2)
plot(t_beat,P_sa(beat),'m','linewidth',2)
set(gca,'FontSize',25)
legend('P_{lv}','P_{la}','P_{sa}','orientation','horizontal')%
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
xlim([t_beat(1) t_beat(end)])

hfig4 = figure(4);
clf 
hold on
plot(t_beat,P_sv(beat),'k','linewidth',2)
plot(t_beat,P_rv(beat),'b','linewidth',2)
plot(t_beat,P_ra(beat),'r','linewidth',2)
plot(t_beat,P_pa(beat),'c','linewidth',2)
plot(t_beat,P_pv(beat),'m','linewidth',2)
set(gca,'FontSize',25)
legend('P_{sv}','P_{rv}','P_{ra}','P_{pa}','P_{pv}','orientation','horizontal')%
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
xlim([t_beat(1) t_beat(end)])
ylim([0 30])

hfig5 = figure(5); 
clf
hold on 
plot(t_beat,V_lv(beat),'b','linewidth',2)
plot(t_beat,V_rv(beat),'r','linewidth',2)
set(gca,'FontSize',25)
legend('V_{lv}','V_{rv}','orientation','horizontal')
xlabel('Time (s)')
ylabel('Volume (mmHg)')
xlim([t_beat(1) t_beat(end)])

hfig6 = figure(6); 
clf
hold on 
plot(t_beat,Q_t_valve(beat),'r','linewidth',2)
plot(t_beat,Q_a_valve(beat),'b--','linewidth',2)
plot(t_beat,Q_p_valve(beat),'r--','linewidth',2)
plot(t_beat,Q_m_valve(beat),'b','linewidth',2)
set(gca,'FontSize',25)
legend('Q_{t}','Q_{a}','Q_{p}','Q_{m}')
xlabel('Time (s)')
ylabel('Flow (L min^{-1})')

if printon == 1
    print(hfig1,'-dpng','F1_PVloops.png')
    print(hfig2,'-dpng','F2_Psa.png')
    print(hfig3,'-dpng','F3_highP.png')
    print(hfig4,'-dpng','F4_lowP.png')
    print(hfig5,'-dpng','F5_ventrV.png')
    print(hfig6,'-dpng','F6_valveF.png')
end 



