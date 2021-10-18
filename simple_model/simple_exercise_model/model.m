% function DXDT = dXdT_cardiovascular_mechanics(t,x,...)
%
%
% Output parameters:
%   DXDT  time derivatives of the model
%
% Mandatory input parameters:
%   t     time
%   x     state variables at time t%   
%
% State Variables:
% []
%
% Parameters:
% []

function [dXdT,outputs] = model(t,x,pars,data)

stim_period = data.stim_period; 

%% Parameters

Vw_LV  = pars(1); 
Vw_SEP = pars(2); 
Vw_RV  = pars(3); 

Amref_LV  = pars(4);
Amref_SEP = pars(5); 
Amref_RV  = pars(6); 

C_Ao = pars(7); 
C_SA = pars(8); 
C_SV = pars(9); 
C_PA = pars(10); 
C_PV = pars(11); 

R_vlv = pars(12);
R_Ao  = pars(13);
R_tAo = pars(14);
R_SA  = pars(15);
R_tSA = pars(16); 
R_SV  = pars(17);
R_PA  = pars(18);
R_PV  = pars(19); 

Lsref  = pars(20);
L_0    = pars(21);
L_C    = pars(22); 
LSEiso = pars(23); 
SLrest = pars(24);

k_pas_LV = pars(25); 
k_pas_RV = pars(26); 
k_act    = pars(27); 

gamma_1 = pars(28); 
gamma_2 = pars(29); 
v_max  = pars(30); 
C_rest = pars(31); 

%% Variables

xm_LV  = x(1); % LV heart geometry variable, cm
xm_SEP = x(2); % septum heart geometry variable, cm
xm_RV  = x(3); % RV heart geometry variable, cm
ym     = x(4); % Heart geometry variable, cm
SL_LV  = x(5); % sarcomere length, LV, micron
SL_SEP = x(6); % sarcomere length, septum, micron
SL_RV  = x(7); % sarcomere length, RV, micron

V_LA   = x(8);
V_LV   = x(9); 
V_Ao   = x(10); 
V_SA   = x(11); 
V_SV   = x(12);
V_RA   = x(13); 
V_RV   = x(14); 
V_PA   = x(15); 
V_PV   = x(16); 

C_LV   = x(17); % LV activtion function
C_SEP  = x(18); % SEP activtion function
C_RV   = x(19); % RV activtion function

%% Equations 

% ventricular mechanics
Vm_LV  = (pi/6) * xm_LV  * (xm_LV^2  + 3*ym^2);
Vm_SEP = (pi/6) * xm_SEP * (xm_SEP^2 + 3*ym^2);
Vm_RV  = (pi/6) * xm_RV  * (xm_RV^2  + 3*ym^2);

Am_LV  = pi * (xm_LV^2  + ym^2);
Am_SEP = pi * (xm_SEP^2 + ym^2);
Am_RV  = pi * (xm_RV^2  + ym^2);

Cm_LV  = 2 * xm_LV  / (xm_LV^2  + ym^2);
Cm_SEP = 2 * xm_SEP / (xm_SEP^2 + ym^2);
Cm_RV  = 2 * xm_RV  / (xm_RV^2  + ym^2);

z_LV  = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV);
z_SEP = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP);
z_RV  = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV);

epsf_LV  = (1/2) * log(Am_LV  / Amref_LV)  - (1/12) * z_LV^2  - 0.019 * z_LV^4;
epsf_SEP = (1/2) * log(Am_SEP / Amref_SEP) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4;
epsf_RV  = (1/2) * log(Am_RV  / Amref_RV)  - (1/12) * z_RV^2  - 0.019 * z_RV^4;

SLo_LV  = Lsref * exp(epsf_LV); 
SLo_SEP = Lsref * exp(epsf_SEP); 
SLo_RV  = Lsref * exp(epsf_RV);

% Passive forces 
sigmapas_LV  = (SLo_LV  - L_0) + gamma_1 * max(0,SLo_LV  - L_C)^gamma_2; %(SLo_LV - L_0) + (exp(gamma * (SLo_LV - L_C)) - 1); %
sigmapas_SEP = (SLo_SEP - L_0) + gamma_1 * max(0,SLo_SEP - L_C)^gamma_2; %(SLo_LV - L_0) + (exp(gamma * (SLo_LV - L_C)) - 1); %
sigmapas_RV  = (SLo_RV  - L_0) + gamma_1 * max(0,SLo_RV  - L_C)^gamma_2; %(SLo_LV - L_0) + (exp(gamma * (SLo_LV - L_C)) - 1); %

% Active forces
sigmaact_LV  = C_LV  * (SL_LV  - SLrest) * (SLo_LV  - SL_LV)  / LSEiso;
sigmaact_SEP = C_SEP * (SL_SEP - SLrest) * (SLo_SEP - SL_SEP) / LSEiso;
sigmaact_RV  = C_RV  * (SL_RV  - SLrest) * (SLo_RV  - SL_RV)  / LSEiso;

% Total forces
sigmaM_LV  = k_act * sigmaact_LV  + k_pas_LV * sigmapas_LV;
sigmaM_SEP = k_act * sigmaact_SEP + k_pas_LV * sigmapas_SEP;
sigmaM_RV  = k_act * sigmaact_RV  + k_pas_RV * sigmapas_RV;

% equilibrium of forces at junction circle
Tm_LV  = (Vw_LV  * sigmaM_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5);
Tm_SEP = (Vw_SEP * sigmaM_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5);
Tm_RV  = (Vw_RV  * sigmaM_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);

Tx_LV  = Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2);
Tx_SEP = Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2);
Tx_RV  = Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2);

Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2);
Ty_SEP = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2);
Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);

% ventricular pressure
ptrans_LV = 2 * Tx_LV / ym ;
ptrans_RV = 2 * Tx_RV / ym ;
P_LV = -ptrans_LV ;
P_RV = +ptrans_RV ;

%% Cardiac activation 

% calcium activation factors for TriSeg
tauD  = 0.032*stim_period; % sec
tauR  = 0.048*stim_period; % sec
tauSC = 0.425*stim_period; % sec


% heart timing 
phi = mod(t-0.1,stim_period)/stim_period;

x  = min(8,max(0,phi*stim_period/tauR));
Fr = 0.02*(x^3)*((8-x)^2)*exp(-x);

CL_LV  = tanh(4*(SL_LV - SLrest)^2);
CL_SEP = tanh(4*(SL_SEP - SLrest)^2);
CL_RV  = tanh(4*(SL_RV - SLrest)^2);

T_LV  = tauSC*(0.29 + 0.3*SL_LV);
T_SEP = tauSC*(0.29 + 0.3*SL_SEP);
T_RV  = tauSC*(0.29 + 0.3*SL_RV);

%% Atria 

Tact    = +0.015;
Emin    = 1.5*0.050 ;
Emax    = 0.150 ;
sigma_a = 0.100;

phi_atria = phi-Tact - 1*((phi-Tact)>(0.5)) + 1*((phi-Tact)<(-0.5));
act  = exp( -(phi_atria/sigma_a)^2 );
E_LA = 3 * (Emin + Emax * act); 
E_RA = Emin + Emax * act; 

P_LA  = E_LA*V_LA;
P_RA  = E_RA*V_RA;

%% Lumped circulatory model
P_SV  = V_SV/C_SV;
P_PV  = V_PV/C_PV;
P_PA  = V_PA/C_PA;

% Ao valves closed equations
Q_a = 0;
P_Ao = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
P_SA = (C_Ao*R_Ao*R_SA*V_SA + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
Q_Ao = -(C_Ao*R_SA*V_SA - C_SA*R_SA*V_Ao - C_SA*R_tSA*V_Ao + C_Ao*C_SA*P_SV*R_tSA)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
if (P_Ao < P_LV)*(V_LV>0) 
  % Ao valve open equations 
  P_SA    = (C_Ao*R_Ao*R_SA*R_tAo*V_SA + C_Ao*R_Ao*R_SA*R_vlv*V_SA + C_SA*R_SA*R_tSA*R_vlv*V_Ao + C_Ao*R_SA*R_tAo*R_vlv*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_vlv + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_vlv)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_vlv + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_vlv + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_vlv + R_SA*R_tAo*R_vlv + R_tSA*R_tAo*R_vlv));
  Q_a = -(C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA - C_Ao*C_SA*P_LV*R_Ao*R_SA - C_Ao*C_SA*P_LV*R_Ao*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_vlv + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_vlv + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_vlv + R_SA*R_tAo*R_vlv + R_tSA*R_tAo*R_vlv));
  Q_Ao    = -(C_Ao*R_SA*R_tAo*V_SA + C_Ao*R_SA*R_vlv*V_SA - C_SA*R_SA*R_vlv*V_Ao - C_SA*R_tSA*R_vlv*V_Ao - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_vlv)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_vlv + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_vlv + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_vlv + R_SA*R_tAo*R_vlv + R_tSA*R_tAo*R_vlv));
end
Q_m = max((P_LA - P_LV) / R_vlv,0)*(V_LA>0);
Q_t = max((P_RA - P_RV) / R_vlv,0)*(V_RA>0);
Q_p = max((P_RV - P_PA) / R_vlv,0)*(V_RV>0);

Q_SA = (P_SA - P_SV) / R_SA; 
Q_PA = (P_PA - P_PV) / R_PA; 

Q_SV = (P_SV - P_RA) / R_PV * (P_SV>P_RA) + (P_SV - P_RA) / R_PV * (P_SV<P_RA)*(V_RA>0);
Q_PV = (P_PV - P_LA) / R_SV * (P_PV>P_LA) + (P_PV - P_LA) / R_SV * (P_PV<P_LA)*(V_LA>0);


%% Differential equations 

% TriSeg DAEs
dxm_LVdt  = (-V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV) / V_LV; 
dxm_SEPdt = (Tx_LV + Tx_SEP + Tx_RV); 
dxm_RVdt  = (+V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV) / V_LV; 
dymdt     = (Ty_LV + Ty_SEP + Ty_RV); 

% Sliding velocities -- Eq. (B2) Lumens et al.
dSL_LVdt  = ((SLo_LV  - SL_LV)  / LSEiso - 1) * v_max;
dSL_SEPdt = ((SLo_SEP - SL_SEP) / LSEiso - 1) * v_max;
dSL_RVdt  = ((SLo_RV  - SL_RV)  / LSEiso - 1) * v_max;

% Circulatory model ODEs
dV_LAdt = Q_PV - Q_m; 
dV_LVdt = Q_m  - Q_a; 
dV_Aodt = Q_a  - Q_Ao; 
dV_SAdt = Q_Ao - Q_SA; 
dV_SVdt = Q_SA - Q_SV;
dV_RAdt = Q_SV - Q_t;
dV_RVdt = Q_t  - Q_p;
dV_PAdt = Q_p  - Q_PA;
dV_PVdt = Q_PA - Q_PV;

% Activation ODEs
dCLVdt  = CL_LV  * Fr / tauR + (C_rest - C_LV)  / (1 + exp((T_LV  - phi * stim_period) / tauD)) / tauD;
dCSEPdt = CL_SEP * Fr / tauR + (C_rest - C_SEP) / (1 + exp((T_SEP - phi * stim_period) / tauD)) / tauD;
dCRVdt  = CL_RV  * Fr / tauR + (C_rest - C_RV)  / (1 + exp((T_RV  - phi * stim_period) / tauD)) / tauD;

dXdT = [dxm_LVdt; dxm_SEPdt; dxm_RVdt; dymdt; 
    dSL_LVdt; dSL_SEPdt; dSL_RVdt; 
    dV_LAdt; dV_LVdt; dV_Aodt; dV_SAdt; dV_SVdt; 
    dV_RAdt; dV_RVdt; dV_PAdt; dV_PVdt; 
    dCLVdt; dCSEPdt; dCRVdt; 
    ]; 

outputs = [P_LA; P_LV; P_Ao; P_SA; P_SV;
    P_RA; P_RV; P_PA; P_PV;
    Q_m; Q_a; Q_t; Q_p; 
    sigmapas_LV; sigmapas_SEP; sigmapas_RV;
    sigmaact_LV; sigmaact_SEP; sigmaact_RV;
    SLo_LV; SLo_SEP; SLo_RV; 
    z_LV; z_SEP; z_RV; 
    Vm_LV; Vm_SEP; Vm_RV; 
    Am_LV; Am_SEP; Am_RV; 
    ];


