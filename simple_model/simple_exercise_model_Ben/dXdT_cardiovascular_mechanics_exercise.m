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

function [dXdT,outputs] = dXdT_cardiovascular_mechanics_exercise(t,x,pars,stim_period,theta)

%% Parameters
Vw_LV     = pars(1); % LV wall volume, mL 
Vw_SEP    = pars(2); % Septal wall volume, mL 
Vw_RV     = pars(3); % RV wall volume, mL 
Amref_LV  = pars(4) ; % LV midwall reference surface area, cm^2
Amref_SEP = pars(5) ; % SEP midwall reference surface area, cm^2
Amref_RV  = pars(6) ; % RV midwall reference surface area, cm^2

tau_v = 0; 
tau_a =  min(- 0.2 * (stim_period),.035); 

% exercise factors
a = 2.50; % inotropy factor
b = 2.00; % arterial vasodilation factor
c = 2.00; % venous vasoconstriction factor
d = 3.00; % arterial vasoconstriction factor
e = 0.50; % pulmonary vasodilation factor
f = 1 + 5*theta; % calcium factor

% Triseg parameters
Lsref     = 1.9; % Resting SL, micron
vmax      = 7; % micron/sec
LSEiso    = 0.04; % micron
sigma_act = 7.5*96*(1 + a*theta); % mmHg 
SLrest    = 1.51; % microns

% Lumped circulatory parameters
C_Ao = 0.65;                    % Proximal aortic compliance, mL/mmHg
C_SA = 1.65/(1 + d*theta);      % Systemic arterial compliance, mL/mmHg
C_SV = 1.4*250/(1 + c*theta);   % Systemic venous compliance, mL/mmHg 
C_PV = 25;                      % Pulmonary venous compliance, mL/mmHg
C_PA = 5.4;                     % Pulmonary arterial compliance, mL/mmHg

R_Ao  = 0.01;                  % resistance of aorta , mmHg*sec/mL
R_SA  = 0.965/(1 + b*theta);   % Systemic vasculature resistance, mmHg*sec/mL
R_PA  = 0.05/(1 + e*theta); %0.05*(1 - e*theta);    % Pulmonary vasculature resistance, mmHg*sec/mL 
R_vlv = 0.002;                 % valve resistance, mmHg*sec/mL
R_tAo = 0.0020;
R_tSA = 0.05;
R_RA  = 0.024;
R_LA  = 0.024;

%% Variables

xm_LV  = x(1); % LV heart geometry variable, cm
xm_SEP = x(2); % septum heart geometry variable, cm
xm_RV  = x(3); % RV heart geometry variable, cm
ym     = x(4); % Heart geometry variable, cm

SL_LV  = x(5); % sarcomere length, LV, micron
SL_SEP = x(6); % sarcomere length, septum, micron
SL_RV  = x(7); % sarcomere length, RV, micron

V_LA   = x(8); % volume of right atrium
V_LV   = x(9); % volume LV, mL
V_Ao   = x(10); % volume of proximal aorta
V_SA   = x(11); % volume of systemic arterys
V_SV   = x(12); % volume of systemic veins
V_RA   = x(13); % volume of left atrium
V_RV   = x(14); % volume RV, mL
V_PA   = x(15); % volume of pulmonary arterys
V_PV   = x(16); % volume of pulmonary veins

%% Ventricular activation 

TS_v = max(.4 * stim_period,.2); 
TR_v = max(.3 * stim_period,.15); 

tn_v = mod(t,stim_period)/stim_period; 
if tn_v >= 0 && tn_v < TS_v 
    y_v = 0.5 * (1 - cos(pi * tn_v / TS_v)); 
elseif tn_v >= TS_v && tn_v < TR_v + TS_v 
    y_v = 0.5 * (1 + cos(pi * (tn_v - TS_v) / TR_v)); 
else
    y_v = 0; 
end 

EM_v = 5; 
Em_v = .02; 
E_v  = (EM_v - Em_v)/2 * y_v + Em_v; 

%% Atrial activation

TS_a = .4 * stim_period; 
TR_a = .3 * stim_period; 

tn_a = mod(t - tau_a,stim_period)/stim_period;     
if tn_a >= 0 && tn_a < TS_a
    y_a = 0.5 * (1 - cos(pi * tn_a / TS_a)); 
elseif tn_a >= TS_a && tn_a < TR_a + TS_a 
    y_a = 0.5 * (1 + cos(pi * (tn_a - TS_a) / TR_a)); 
else
    y_a = 0; 
end 

EM_a = 0.375; 
Em_a = 0.075; 
E_a  = (EM_a - Em_a)/2 * y_a + Em_a; 

E_LA = 3 * E_a; 
E_RA = E_a; 

%% Heart model 

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

z_LV   = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV);
z_SEP  = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP);
z_RV   = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV);

epsf_LV  = (1/2) * log(Am_LV  / Amref_LV)  - (1/12) * z_LV^2  - 0.019 * z_LV^4;
epsf_SEP = (1/2) * log(Am_SEP / Amref_SEP) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4;
epsf_RV  = (1/2) * log(Am_RV  / Amref_RV)  - (1/12) * z_RV^2  - 0.019 * z_RV^4;

SLo_LV  = Lsref * exp(epsf_LV); 
SLo_SEP = Lsref * exp(epsf_SEP); 
SLo_RV  = Lsref * exp(epsf_RV);

% sigmapas_LV  = 22 * (SLo_LV  - 1.6)  + 22 * max(0,SLo_LV  - 1.8)^3 ;
% sigmapas_SEP = 22 * (SLo_SEP  - 1.6) + 22 * max(0,SLo_SEP - 1.8)^3 ;
% sigmapas_RV  = 22 * (SLo_RV  - 1.6)  + 22 * max(0,SLo_RV  - 1.8)^3 ;


P_C = 26.4; 
P_expC = 2.48; 
SL_C = 2.09; 

k_P = 5; 
L_0 = 1.6; 

sigmaC_LV  = P_C * max(0, SLo_LV  - SL_C)^P_expC; 
sigmaC_SEP = P_C * max(0, SLo_SEP - SL_C)^P_expC; 
sigmaC_RV  = P_C * max(0, SLo_RV  - SL_C)^P_expC; 

sigmapas_LV  = k_P * (SLo_LV  - L_0) + sigmaC_LV; 
sigmapas_SEP = k_P * (SLo_SEP - L_0) + sigmaC_SEP; 
sigmapas_RV  = k_P * (SLo_RV  - L_0) + sigmaC_RV; 


% Active forces
sigmaact_LV  = sigma_act * E_v * (SL_LV  - SLrest) * (SLo_LV  - SL_LV)  / LSEiso;
sigmaact_SEP = sigma_act * E_v * (SL_SEP - SLrest) * (SLo_SEP - SL_SEP) / LSEiso;
sigmaact_RV  = sigma_act * E_v * (SL_RV  - SLrest) * (SLo_RV  - SL_RV)  / LSEiso;

% Total forces
sigmaM_LV  = sigmaact_LV  + sigmapas_LV;
sigmaM_SEP = sigmaact_SEP + sigmapas_SEP;
sigmaM_RV  = sigmaact_RV  + sigmapas_RV;

% equilibrium of forces at junction circle
Tm_LV  = (Vw_LV  * sigmaM_LV  / (2 * Am_LV))  * (1 + (1/3) * z_LV^2  + (1/5) * z_LV^4);
Tm_SEP = (Vw_SEP * sigmaM_SEP / (2 * Am_SEP)) * (1 + (1/3) * z_SEP^2 + (1/5) * z_SEP^4);
Tm_RV  = (Vw_RV  * sigmaM_RV  / (2 * Am_RV))  * (1 + (1/3) * z_RV^2  + (1/5) * z_RV^4);

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

%% Lumped circulatory model

% V0_col = 40; 
% s      = 10; 
% P_LA_col = exp(s * (V_LA / V0_col - 1));
% P_RA_col = exp(s * (V_RA / V0_col - 1));

P_LA = E_LA * V_LA;% + P_LA_col; 
P_RA = E_RA * V_RA;% + P_RA_col; 
P_SV = V_SV / C_SV;
P_PV = V_PV / C_PV;
P_PA = V_PA / C_PA ;

% Ao valves closed equations
Q_avlv = 0;
P_Ao = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
P_SA = (C_Ao*R_Ao*R_SA*V_SA + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
Q_Ao = -(C_Ao*R_SA*V_SA - C_SA*R_SA*V_Ao - C_SA*R_tSA*V_Ao + C_Ao*C_SA*P_SV*R_tSA)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
if (P_Ao < P_LV) 
  % Ao valve open equations 
  P_SA    = (C_Ao*R_Ao*R_SA*R_tAo*V_SA + C_Ao*R_Ao*R_SA*R_vlv*V_SA + C_SA*R_SA*R_tSA*R_vlv*V_Ao + C_Ao*R_SA*R_tAo*R_vlv*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_vlv + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_vlv)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_vlv + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_vlv + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_vlv + R_SA*R_tAo*R_vlv + R_tSA*R_tAo*R_vlv));
  Q_avlv = -(C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA - C_Ao*C_SA*P_LV*R_Ao*R_SA - C_Ao*C_SA*P_LV*R_Ao*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_vlv + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_vlv + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_vlv + R_SA*R_tAo*R_vlv + R_tSA*R_tAo*R_vlv));
  Q_Ao    = -(C_Ao*R_SA*R_tAo*V_SA + C_Ao*R_SA*R_vlv*V_SA - C_SA*R_SA*R_vlv*V_Ao - C_SA*R_tSA*R_vlv*V_Ao - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_vlv)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_vlv + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_vlv + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_vlv + R_SA*R_tAo*R_vlv + R_tSA*R_tAo*R_vlv));
end
Q_SA = (P_SA - P_SV)/R_SA; 
Q_SV = (P_SV - P_RA)/R_RA;
Q_PA = (P_PA - P_PV)/R_PA; 
Q_PV = (P_PV - P_LA)/R_LA;
Q_mvlv = max((P_LA - P_LV)/R_vlv,0);
Q_tvlv = max((P_RA - P_RV)/R_vlv,0);
Q_pvlv = max((P_RV - P_PA)/R_vlv,0);

%% Differential and algebraic equations 

% TriSeg
dxm_LVdt  = (-V_LV - 0.5*Vw_LV - 0.5*Vw_SEP + Vm_SEP - Vm_LV)/V_LV; %
dxm_SEPdt = (Tx_LV + Tx_SEP + Tx_RV); % 
dxm_RVdt  = (+V_RV + 0.5*Vw_RV + 0.5*Vw_SEP + Vm_SEP - Vm_RV)/V_LV; % 
dymdt     = (Ty_LV + Ty_SEP + Ty_RV);  %

% Sliding velocities -- Eq. (B2) Lumens et al.
dSL_LVdt  = ((SLo_LV  - SL_LV) /LSEiso - 1)*vmax;
dSL_SEPdt = ((SLo_SEP - SL_SEP)/LSEiso - 1)*vmax;
dSL_RVdt  = ((SLo_RV  - SL_RV) /LSEiso - 1)*vmax;

% Circulatory model
dV_LAdt = Q_PV   - Q_mvlv; 
dV_LVdt = Q_mvlv - Q_avlv; 
dV_Aodt = Q_avlv - Q_Ao; 
dV_SAdt = Q_Ao   - Q_SA;  
dV_SVdt = Q_SA   - Q_SV;  
dV_RAdt = Q_SV   - Q_tvlv; 
dV_RVdt = Q_tvlv - Q_pvlv; 
dV_PAdt = Q_pvlv - Q_PA;  
dV_PVdt = Q_PA   - Q_PV; 

dXdT = [dxm_LVdt; dxm_SEPdt; dxm_RVdt; dymdt; 
    dSL_LVdt; dSL_SEPdt; dSL_RVdt; 
    dV_LAdt; dV_LVdt; dV_Aodt; dV_SAdt; dV_SVdt; 
    dV_RAdt; dV_RVdt; dV_PAdt; dV_PVdt;      
    ]; 

outputs = [P_LA; P_LV; P_Ao; P_SA; P_SV;
    P_RA; P_RV; P_PA; P_PV;
    Q_mvlv; Q_tvlv;
    sigmapas_LV; sigmapas_SEP; sigmapas_RV;
    ];


