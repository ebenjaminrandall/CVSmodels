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

function [dXdT,outputs] = dXdT_cardiovascular_mechanics(t,x,para)

%% Parameters
stim_period = para(1); % 
Vw_LV     = para(2); % LV wall volume, mL 
Vw_SEP    = para(3); % Septal wall volume, mL 
Vw_RV     = para(4); % RV wall volume, mL 
Amref_LV  = para(5) ; % LV midwall reference surface area, cm^2
Amref_SEP = para(6) ; % SEP midwall reference surface area, cm^2
Amref_RV  = para(7) ; % RV midwall reference surface area, cm^2
theta     = para(8); % exercise level

% exercise factors
a = 2.50; % inotropy factor
b = 2.00; % arterial vasodilation factor
c = 2.00; % venous vasoconstriction factor
d = 3.00; % arterial vasoconstriction factor
e = 10*0.50; % pulmonary vasodilation factor
f = 1 + 5*theta; % calcium factor

% Triseg parameters
Lsref   = 1.9; % Resting SL, micron
vmax    = 7; % micron/sec
LSEiso  = 0.04; % micron
sigma_act = 7.5*96*(1 + a*theta); % mmHg 
SLrest  = 1.51; % microns

% % Active force parameters
% TS = 0.20*stim_period; % time from the end of diastole to maximal systole
% TR = 0.30*stim_period; % time for relaxation from systole

% Lumped circulatory parameters
C_Ao = 0.65;  % Proximal aortic compliance, mL/mmHg
C_SA = 1.65/(1 + d*theta); % Systemic arterial compliance, mL/mmHg
C_SV = 1.4*250/(1 + c*theta); % Systemic venous compliance, mL/mmHg 
C_PV =25; % Pulmonary venous compliance, mL/mmHg
C_PA =5.4; % Pulmonary arterial compliance, mL/mmHg
R_Ao   = 0.01; % resistance of aorta , mmHg*sec/mL
R_SA   =  0.965/(1 + b*theta);% mmHg*sec/mL; % Systemic vasculature resistance, mmHg*sec/mL
R_PA   = 0.05*(1 - e*theta); % Pulmonary vasculature resistance, mmHg*sec/mL 
R_vlv  = 0.002; %  valve resistance, mmHg*sec/mL
R_tAo  = 0.0020;
R_tSA  = 0.05;
R_PV   = 0.024;
R_SV   = 0.024;

C_Ao = 0.36; 
C_SA = 0.875 / (1 + d * theta); 
C_SV = 87.5 / (1 + c * theta);
C_PV = 12.5; 
C_PA = 6; 

R_SA = 1.4 / (1 + b * theta); 
R_PA = 0.25 / (1 + e * theta); %
R_PV = 0.036; 
R_SV = 0.036; 

%% Variables

xm_LV  = x(1); % LV heart geometry variable, cm
xm_SEP = x(2); % septum heart geometry variable, cm
xm_RV  = x(3); % RV heart geometry variable, cm
ym     = x(4); % Heart geometry variable, cm
SL_LV  = x(5); % sarcomere length, LV, micron
SL_SEP = x(6); % sarcomere length, septum, micron
SL_RV  = x(7); % sarcomere length, RV, micron
V_LV   = x(8); % volume LV, mL
V_RV   = x(9); % volume RV, mL

V_SV   = x(10); % volume of systemic veins
V_PV   = x(11); % volume of pulmonary veins
V_SA   = x(12); % volume of systemic arterys
V_PA   = x(13); % volume of pulmonary arterys
V_Ao   = x(14); % volume of proximal aorta
V_RA   = x(15); % volume of left atrium
V_LA   = x(16); % volume of right atrium

C_LV   = x(17); % LV activtion function
C_SEP  = x(18); % SEP activtion function
C_RV   = x(19); % RV activtion function

% ventricular mechanics
Vm_LV  = (pi/6)*xm_LV*(xm_LV^2 + 3*ym^2);
Vm_SEP = (pi/6)*xm_SEP*(xm_SEP^2 + 3*ym^2);
Vm_RV  = (pi/6)*xm_RV*(xm_RV^2 + 3*ym^2);
Am_LV  = pi*(xm_LV^2 + ym^2);
Am_SEP = pi*(xm_SEP^2 + ym^2);
Am_RV  = pi*(xm_RV^2 + ym^2);
Cm_LV  = 2*xm_LV/(xm_LV^2 + ym^2);
Cm_SEP = 2*xm_SEP/(xm_SEP^2 + ym^2);
Cm_RV  = 2*xm_RV/(xm_RV^2 + ym^2);
z_LV   = 3*Cm_LV*Vw_LV/(2*Am_LV);
z_SEP  = 3*Cm_SEP*Vw_SEP/(2*Am_SEP);
z_RV   = 3*Cm_RV*Vw_RV/(2*Am_RV);

epsf_LV = (1/2)*log(Am_LV/Amref_LV) - (1/12)*z_LV^2 - 0.019*z_LV^4;
epsf_SEP = (1/2)*log(Am_SEP/Amref_SEP) - (1/12)*z_SEP^2 - 0.019*z_SEP^4;
epsf_RV = (1/2)*log(Am_RV/Amref_RV) - (1/12)*z_RV^2 - 0.019*z_RV^4;
SLo_LV = Lsref*exp(epsf_LV); SLo_SEP = Lsref*exp(epsf_SEP); SLo_RV = Lsref*exp(epsf_RV);

% heart timing 
phi = mod(t-0.1,stim_period)/stim_period;

% calcium activation factors for TriSeg
tauD  = 0.032*stim_period; % sec
tauR  = 0.048*stim_period; % sec
tauSC = 0.425*stim_period; % sec
Crest = 0.02/f;
x = min(8,max(0,phi*stim_period/tauR));
Fr = 0.02*(x^3)*((8-x)^2)*exp(-x);
CL_LV = tanh(4*(SL_LV - SLrest)^2);
T_LV = tauSC*(0.29 + 0.3*SL_LV);
dCLVdt = CL_LV*Fr/tauR + (Crest - C_LV)/(1 + exp((T_LV-phi*stim_period)/tauD))/tauD;

CL_SEP = tanh(4*(SL_SEP - SLrest)^2);
T_SEP = tauSC*(0.29 + 0.3*SL_SEP);
dCSEPdt = CL_SEP*Fr/tauR + (Crest - C_SEP)/(1 + exp((T_SEP-phi*stim_period)/tauD))/tauD;

CL_RV = tanh(4*(SL_RV - SLrest)^2);
T_RV = tauSC*(0.29 + 0.3*SL_RV);
dCRVdt = CL_RV*Fr/tauR + (Crest - C_RV)/(1 + exp((T_RV-phi*stim_period)/tauD))/tauD;

% Passive forces
% sigmapas_LV  = (22) * (SLo_LV  - 1.8)  + 10*exp(25*(SLo_LV-2.2)) ;
% sigmapas_SEP = (22) * (SLo_SEP  - 1.8) + 10*exp(25*(SLo_SEP-2.2)) ;
% sigmapas_RV  = (22) * (SLo_RV  - 1.8)  + 10*exp(25*(SLo_RV-2.2)) ;

sigmapas_LV  = (SLo_LV  - 1.6)  + max(0,SLo_LV-1.8)^3 ;
sigmapas_SEP = (SLo_SEP  - 1.6) + max(0,SLo_SEP-1.8)^3 ;
sigmapas_RV  = (SLo_RV  - 1.6) + max(0,SLo_RV-1.8)^3 ;

% Active forces
sigmaact_LV  = sigma_act*C_LV*(SL_LV-SLrest)*(SLo_LV - SL_LV)/LSEiso;
sigmaact_SEP = sigma_act*C_SEP*(SL_SEP-SLrest)*(SLo_SEP - SL_SEP)/LSEiso;
sigmaact_RV  = sigma_act*C_RV*(SL_RV-SLrest)*(SLo_RV - SL_RV)/LSEiso;

% sigmaact_LV  = 500 * y_v  * (SL_LV-SLrest)  ; 
% sigmaact_SEP = 500 * y_v  * (SL_SEP-SLrest) ;
% sigmaact_RV  = 500 * y_v  * (SL_RV-SLrest)  ;

% Total forces
sigmaM_LV  = sigmaact_LV  + 22 * sigmapas_LV;
sigmaM_SEP = sigmaact_SEP + 22 * sigmapas_SEP;
sigmaM_RV  = sigmaact_RV  + 22 * sigmapas_RV;

% equilibrium of forces at junction circle
Tm_LV  = (Vw_LV*sigmaM_LV/(2*Am_LV))*(1 + (z_LV^2)/3 + (z_LV^4)/5);
Tm_SEP = (Vw_SEP*sigmaM_SEP/(2*Am_SEP))*(1 + (z_SEP^2)/3 + (z_SEP^4)/5);
Tm_RV  = (Vw_RV*sigmaM_RV/(2*Am_RV))*(1 + (z_RV^2)/3 + (z_RV^4)/5);

Tx_LV  = Tm_LV*2*xm_LV*ym/(xm_LV^2 + ym^2);
Tx_SEP = Tm_SEP*2*xm_SEP*ym/(xm_SEP^2 + ym^2);
Tx_RV  = Tm_RV*2*xm_RV*ym/(xm_RV^2 + ym^2);

Ty_LV  = Tm_LV*(-xm_LV^2 + ym^2)/(xm_LV^2 + ym^2);
Ty_SEP = Tm_SEP*(-xm_SEP^2 + ym^2)/(xm_SEP^2 + ym^2);
Ty_RV  = Tm_RV*(-xm_RV^2 + ym^2)/(xm_RV^2 + ym^2);

% ventricular pressure
ptrans_LV = 2*Tx_LV/ym ;
ptrans_RV = 2*Tx_RV/ym ;
P_LV = -ptrans_LV ;
P_RV = +ptrans_RV ;

%% Lumped circulatory model
P_SV  = V_SV/C_SV;
P_PV  = V_PV/C_PV;
P_PA  = V_PA/C_PA ;

% atria
Tact = +0.015;
phi_atria = phi-Tact - 1*((phi-Tact)>(0.5)) + 1*((phi-Tact)<(-0.5));
% Emin = 0.05 + 0.20;
% Emax = 0.15 + 0.20; 
Emin = 1.5*0.050 ;
Emax = 0.150 ;
sigma_a = 0.100;
act = exp( -(phi_atria/sigma_a)^2 );

P_LA  = 3.0*(Emin + Emax*act)*V_LA ;
P_RA  = (Emin + Emax*act)*V_RA ;

% Ao valves closed equations
QOUT_LV = 0;
P_Ao = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
P_SA = (C_Ao*R_Ao*R_SA*V_SA + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
Q_Ao = -(C_Ao*R_SA*V_SA - C_SA*R_SA*V_Ao - C_SA*R_tSA*V_Ao + C_Ao*C_SA*P_SV*R_tSA)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
if (P_Ao < P_LV)*(V_LV>0) 
  % Ao valve open equations 
  P_SA    = (C_Ao*R_Ao*R_SA*R_tAo*V_SA + C_Ao*R_Ao*R_SA*R_vlv*V_SA + C_SA*R_SA*R_tSA*R_vlv*V_Ao + C_Ao*R_SA*R_tAo*R_vlv*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_vlv + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_vlv)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_vlv + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_vlv + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_vlv + R_SA*R_tAo*R_vlv + R_tSA*R_tAo*R_vlv));
  QOUT_LV = -(C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA - C_Ao*C_SA*P_LV*R_Ao*R_SA - C_Ao*C_SA*P_LV*R_Ao*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_vlv + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_vlv + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_vlv + R_SA*R_tAo*R_vlv + R_tSA*R_tAo*R_vlv));
  Q_Ao    = -(C_Ao*R_SA*R_tAo*V_SA + C_Ao*R_SA*R_vlv*V_SA - C_SA*R_SA*R_vlv*V_Ao - C_SA*R_tSA*R_vlv*V_Ao - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_vlv)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_vlv + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_vlv + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_vlv + R_SA*R_tAo*R_vlv + R_tSA*R_tAo*R_vlv));
end
QIN_LV  = max((P_LA - P_LV)/R_vlv,0)*(V_LA>0);
QIN_RV  = max((P_RA - P_RV)/R_vlv,0)*(V_RA>0);
QOUT_RV = max((P_RV - P_PA)/R_vlv,0)*(V_RV>0);

QIN_RA = (P_SV - P_RA)/R_PV*(P_SV>P_RA) + (P_SV - P_RA)/R_PV*(P_SV<P_RA)*(V_RA>0);
QIN_LA = (P_PV - P_LA)/R_SV*(P_PV>P_LA) + (P_PV - P_LA)/R_SV*(P_PV<P_LA)*(V_LA>0);
% QIN_RA = (P_SV - P_RA)/R_RA;
% QIN_LA = (P_PV - P_LA)/R_LA;

% TriSeg
dXdT(1) = (-V_LV - 0.5*Vw_LV - 0.5*Vw_SEP + Vm_SEP - Vm_LV)/V_LV; %
dXdT(2) = (Tx_LV + Tx_SEP + Tx_RV); % 
dXdT(3) = (+V_RV + 0.5*Vw_RV + 0.5*Vw_SEP + Vm_SEP - Vm_RV)/V_LV; % 
dXdT(4) = (Ty_LV + Ty_SEP + Ty_RV);  %

% Sliding velocities -- Eq. (B2) Lumens et al.

dSL_LV  = ((SLo_LV - SL_LV)/LSEiso   - 1)*vmax;
dSL_SEP = ((SLo_SEP - SL_SEP)/LSEiso - 1)*vmax;
dSL_RV  = ((SLo_RV - SL_RV)/LSEiso   - 1)*vmax;
dXdT(5) = dSL_LV;
dXdT(6) = dSL_SEP;
dXdT(7) = dSL_RV;

% Circulatory model
dXdT(8)  = QIN_LV - QOUT_LV; % V_LV
dXdT(9)  = QIN_RV - QOUT_RV; % V_RV
dXdT(10) = (P_SA - P_SV)/R_SA - QIN_RA;  % V_SV
dXdT(11) = (P_PA - P_PV)/R_PA - QIN_LA;  % V_PV
dXdT(12) = Q_Ao - (P_SA - P_SV)/R_SA; % V_SA 
dXdT(13) = QOUT_RV - (P_PA - P_PV)/R_PA; % V_PA 
dXdT(14) = QOUT_LV - Q_Ao; % V_Ao
dXdT(15) = QIN_RA - QIN_RV; % V_RA
dXdT(16) = QIN_LA - QIN_LV; % V_LA

dXdT(17) = dCLVdt;
dXdT(18) = dCSEPdt;
dXdT(19) = dCRVdt;

dXdT = dXdT(:);

outputs = [P_LA; P_LV; P_Ao; P_SA; P_SV;
           P_RA; P_RV; P_PA; P_PV;
           QIN_LV; QIN_RV;
           sigmapas_LV; sigmapas_SEP; sigmapas_RV;
           sigmaact_LV; sigmaact_SEP; sigmaact_RV;
           SLo_LV; SLo_SEP; SLo_RV; 
           z_LV; z_SEP; z_RV; 
           Vm_LV; Vm_SEP; Vm_RV; 
            ];


