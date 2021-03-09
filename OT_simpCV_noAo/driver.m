% driver

clear all

printon = 0; 
exerciseon = 0; 

% Evaluation time 
dt    = 0.001; 
tstart = 0; 
tend  = 20; 

tspan = tstart:dt:tend;

%% Load pseudodata

% Mean pressure values (mmHg)
SPbar = 120; 
Pbar  = 100;
DPbar = 80; 

% Heart rate (bpm)
HR_rest = 75; 

if exerciseon == 1
    HR_ex = 180; 
else 
    HR_ex   = 75;  
end 

N_beats_restHR = 10;
if tend <= 20 
    N_beats_ramp = 10; 
else 
    N_beats_ramp = 75; 
end  

% Total blood volume (mL)
Vtot = 4.5 * 1e3; 

% Stroke volume (mL)
SV_nom = .015 * Vtot; 

% Cardiac output (mL min^(-1))
CO_nom = SV_nom * HR_rest;  

% Ejection fraction (dimensionless) 
EF_nom = .6; 

% End-diastolic and end systolic left ventricular volume (mL)
EDV_LV = SV_nom / EF_nom;  
ESV_LV = EDV_LV - SV_nom; 

% Total systemic vascular resistance (mmHg min ml^(-1))
SVR_nom = Pbar / CO_nom; 

%% Conversions

% Convert from mmHg to kPa
SPbar = SPbar / 7.5;   
DPbar = DPbar / 7.5; 
Pbar  = Pbar  / 7.5; 

% Convert from mL to m^3
Vtot   = Vtot * 1e-6; 
SV_nom = SV_nom * 1e-6; 
EDV_LV = EDV_LV * 1e-6; 
ESV_LV = ESV_LV * 1e-6; 

% Convert from mL min^(-1) to m^3 s^(-1)
CO_nom = CO_nom * 1e-6 / 60;     

% Convert from mmHg min ml^(-1) to kPa s m^(-3)
SVR_nom = SVR_nom / 7.5 / 1e-6 * 60; 

%% Model specifications 

% Make data/input structure
data.Pbar   = Pbar; 
data.SPbar  = SPbar; 
data.DPbar  = DPbar; 
data.EF     = EF_nom; 
data.Vtot   = Vtot; 
data.SV     = SV_nom; 
data.EDV_LV = EDV_LV; 
data.ESV_LV = ESV_LV;
data.CO     = CO_nom; 
data.tspan  = tspan; 
data.dt     = dt; 

data.HR.HR_rest = HR_rest; 
data.HR.HR_ex   = HR_ex; 
data.HR.N_beats_restHR = N_beats_restHR; 
data.HR.N_beats_ramp   = N_beats_ramp; 

% Global specifications
ODE_TOL = 1e-8; 

data.gpars.ODE_TOL = ODE_TOL; 
data.gpars.exerciseon = exerciseon; 

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

beats = outputs.beats; 
beat = beats(end-3):beats(end-1); 
HR_end = HR_rest; 

t_beat = tspan(beat) - tspan(beat(1)); 

SV = max(V_lv(beat)) - min(V_lv(beat)) % mL
EF = SV / max(V_lv(beat)) % dimensionless
CO = trapz(t_beat/60,Q_a_valve(beat))/(t_beat(end)/60 - t_beat(1)/60) %SV * HR_end * 1e-3 % L min^(-1)
CP = trapz(P_lv(beat),V_lv(beat)) / 7.5 * 1e-3 * HR_end/60; %mean(P_sa(beat)) / 7.5 * 1e3 * SV * 1e-6 * HR_end/60 % W 
CP = CP / 2 %average over 2 beats 

if exerciseon == 1 
    save exercise.mat
else
    save nom.mat
end 

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
legend('P_{lv}','P_{la}','P_{sa}','orientation','horizontal')
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
legend('P_{sv}','P_{rv}','P_{ra}','P_{pa}','P_{pv}','orientation','horizontal')
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
