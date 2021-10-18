function f = TrisegEquations(x,V_LV,V_RV,pars)


%% Unpack parameters

Vw_LV  = pars(1); % LV wall volume, mL 
Vw_SEP = pars(2); % Septal wall volume, mL 
Vw_RV  = pars(3); % RV wall volume, mL 

Amref_LV  = pars(4); 
Amref_SEP = pars(5); 
Amref_RV  = pars(6); 

Lsref = pars(22); 
k_pas = pars(26); 

%% Variables 

xm_LV  = x(1); % LV heart geometry variable, cm
xm_SEP = x(2); % septum heart geometry variable, cm
xm_RV  = x(3); % RV heart geometry variable, cm
ym     = x(4); % Heart geometry variable, cm

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

sigmapas_LV  = k_pas * ((SLo_LV  - 1.6) + max(0,SLo_LV  - 1.8)^3);
sigmapas_SEP = k_pas * ((SLo_SEP - 1.6) + max(0,SLo_SEP - 1.8)^3);
sigmapas_RV  = k_pas * ((SLo_RV  - 1.6) + max(0,SLo_RV  - 1.8)^3);

% Active forces
sigmaact_LV  = 0; %k_act * C_LV  * (SL_LV  - SLrest) * (SLo_LV  - SL_LV)  / LSEiso;
sigmaact_SEP = 0; %k_act * C_SEP * (SL_SEP - SLrest) * (SLo_SEP - SL_SEP) / LSEiso;
sigmaact_RV  = 0; %k_act * C_RV  * (SL_RV  - SLrest) * (SLo_RV  - SL_RV)  / LSEiso;

% Total forces
sigmaM_LV  = sigmaact_LV  + sigmapas_LV;
sigmaM_SEP = sigmaact_SEP + sigmapas_SEP;
sigmaM_RV  = sigmaact_RV  + sigmapas_RV;

% equilibrium of forces at junction circle
Tm_LV  = (Vw_LV  * sigmaM_LV  / (2 * Am_LV))  * (1 + (z_LV^2)  / 3 + (z_LV^4)  / 5);
Tm_SEP = (Vw_SEP * sigmaM_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2) / 3 + (z_SEP^4) / 5);
Tm_RV  = (Vw_RV  * sigmaM_RV  / (2 * Am_RV))  * (1 + (z_RV^2)  / 3 + (z_RV^4)  / 5);

Tx_LV  = Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2);
Tx_SEP = Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2);
Tx_RV  = Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2);

Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2);
Ty_SEP = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2);
Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);

f(1) = (-V_LV - 0.5*Vw_LV - 0.5*Vw_SEP + Vm_SEP - Vm_LV)/V_LV; %xm_LV
f(2) = (Tx_LV + Tx_SEP + Tx_RV); % xm_SEP
f(3) = (V_RV + 0.5*Vw_RV + 0.5*Vw_SEP + Vm_SEP - Vm_RV)/V_LV; %xm_RV
f(4) = (Ty_LV + Ty_SEP + Ty_RV);  %ym

