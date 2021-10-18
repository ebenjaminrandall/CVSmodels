function [pars,data] = parameters(exercise,data)

Vtot  = data.Vtot; 
theta = data.theta;
vfactor = data.vfactor; 

CO = Vtot / 60; 

%% Cardiac shape parameters

Vw_LV  = 80;
Vw_SEP = 38; 
Vw_RV  = 28; 

Amref_LV  = 0.975*80; % LV midwall reference surface area, cm^2
Amref_SEP = 0.975*45; % SEP midwall reference surface area, cm^2
Amref_RV  = 1.12*100; % RV midwall reference surface area, cm^2

%% Cardiovascular parameters 

% Volume fractions 
d_LA = .005; 
d_LV = .03; 
d_Ao = .03; 
d_SA = .07; 
d_SV = .70; 
d_RA = .005; 
d_RV = .03; 
d_PA = .03; 
d_PV = .10; 

% Chamber volumes 
V_LA = d_LA * Vtot; 
V_LV = d_LV * Vtot; 
V_Ao = d_Ao * Vtot; 
V_SA = d_SA * Vtot; 
V_SV = d_SV * Vtot; 
V_RA = d_RA * Vtot; 
V_RV = d_RV * Vtot; 
V_PA = d_PA * Vtot; 
V_PV = d_PV * Vtot; 
V_T = V_LA + V_LV + V_Ao + V_SA + V_SV + V_RA + V_RV + V_PA + V_PV;

% Stressed volumes 
V_Ao_s = 0.3 * V_Ao; 
V_SA_s = 0.3 * V_SA; 
V_SV_s = 0.1 * V_SV; 
V_PA_s = 0.4 * V_PA; 
V_PV_s = 0.1 * V_PV; 
V_Ts = V_LA + V_LV + V_Ao_s + V_SA_s + V_SV_s + V_RA + V_RV + V_PA_s + V_PV_s; 

V0 = [V_LA; V_LV; V_Ao_s; V_SA_s; V_SV_s; V_RA; V_RV; V_PA_s; V_PV_s]; 
data.V0 = V0; 

% Normal pressures 
P_LA_M = 10;    P_LV_M = 130;   P_Ao_M = 125;   P_SA_M = 120; 
P_LA_m = 1;     P_LV_m = 3;     P_Ao_m = 85;    P_SA_m = 80;
 
P_SV_M = 8;     P_RA_M = 6;     P_RV_M = 25;    P_PA_M = 25;    P_PV_M = 8; 
P_SV_m = 4;     P_RA_m = 1;     P_RM_m = 5;     P_PA_m = 15;    P_PV_m = 4; 

% Compliances 
C_Ao = V_Ao_s / P_Ao_M;
C_SA = V_SA_s / P_SA_M;
C_SV = V_SV_s / (P_SV_M - P_SV_m);
C_PA = V_PA_s / (P_PA_M - P_PA_m);
C_PV = V_PV_s / (P_PV_M - P_PV_m); 

% Resistances 
R_vlv = 0.002; 
R_Ao  = 0.01; 
R_tAo = 0.002; 
R_SA  = (P_SA_M - P_SV_m) / CO; 
R_tSA = 0.05; 
R_SV  = (P_SV_m - P_RA_m) / CO;
R_PA  = (P_PA_M - P_PV_m) / CO;
R_PV  = (P_PV_m - P_LA_m) / CO;

%% Triseg parameters 

Lsref    = 1.9; 
L_0      = 1.6; 
L_C      = 1.8; 
LSEiso   = 0.04; 
SLrest   = 1.51; 
gamma_1  = 1;
gamma_2  = 3 ; %.5; 
v_max    = 7; 
k_pas_LV = 22;
k_pas_RV = 22; 
k_act    = 7.5 * 96 ; 
C_rest   = 0.02; %* (2 * vfactor - 1)

%% New Amref 

EDP = 5; 
EDV = 140; 

V_D = vfactor * EDV; 

An = 28; 
Bn = 3; 

V_0  = EDV * (0.6 - 0.006 * EDP); 
V_30 = V_0 + (EDV - V_0) / ((EDP / An)^(1/Bn));
beta = log(EDP/30) / log(EDV / V_30); 
alpha = 30 / V_30^beta; 

P_D_l = alpha * V_D^beta;
P_D_r = P_D_l / 6;

SLo_D = 2.2; %0.0046 * V_D + 1.55; %   

xm_l = -0.009 * V_D - 3.32; %4-4.5; %
xm_s = -0.5 * xm_l + 0.15;
xm_r = 0.015 * V_D + 3.4; %3.87; %-4.5; %
y    = -0.8 * xm_l - 0.05; %3.5; %

% LV
Am_l = pi * (xm_l.^2 + y.^2); 
Cm_l = 2 * xm_l ./ (xm_l.^2 + y.^2); 
z_l  = 3 .* Cm_l .* Vw_LV ./ (2 .* Am_l); 

epsf_l = (1/2) * log(Am_l) - (1/12) .* z_l.^2 - 0.019 .* z_l.^4; 

Amref_LV = (Lsref / SLo_D * exp(epsf_l)).^2;

% SEP 
Am_s = pi * (xm_s.^2 + y.^2); 
Cm_s = 2 * xm_s ./ (xm_s.^2 + y.^2); 
z_s  = 3 .* Cm_s .* Vw_SEP ./ (2 .* Am_s); 

epsf_s = (1/2) * log(Am_s) - (1/12) .* z_s.^2 - 0.019 .* z_s.^4;

Amref_SEP = (Lsref / SLo_D * exp(epsf_s)).^2;

% RV 
Am_r = pi * (xm_r.^2 + y.^2); 
Cm_r = 2 * xm_r ./ (xm_r.^2 + y.^2); 
z_r  = 3 .* Cm_r .* Vw_RV ./ (2 .* Am_r); 

epsf_r = (1/2) * log(Am_r) - (1/12) .* z_r.^2 - 0.019 .* z_r.^4; 

Amref_RV = (Lsref / SLo_D * exp(epsf_r)).^2;

%% New k_pas


% LV and SP
sinalpha_l = 2 .* xm_l .* y ./ (xm_l.^2 + y.^2); 
Gamma_l    = (Vw_LV ./ (2 .* Am_l)) .* (1 + (1/3) .* z_l.^2 + (1/5) .* z_l.^4); 
sigmaact_l = k_act * C_rest * (SLo_D - SLrest) * 1;  
sigmapas_l = (SLo_D - L_0) + gamma_1 * (SLo_D - L_C).^gamma_2; %(SLo_D - L_0) + (exp(gamma * (SLo_D - L_C)) - 1); %
k_pas_LV   = (P_D_l .* y ./ (-2 .* sinalpha_l .* Gamma_l) - sigmaact_l) ./ sigmapas_l;

% RV
sinalpha_r = 2 .* xm_r .* y ./ (xm_r.^2 + y.^2); 
Gamma_r    = (Vw_RV ./ (2 .* Am_r)) .* (1 + (1/3) .* z_r.^2 + (1/5) .* z_r.^4); 
sigmaact_r = k_act * C_rest * (SLo_D - SLrest) * 1;  
sigmapas_r = (SLo_D - L_0) + gamma_1 * (SLo_D - L_C).^gamma_2; 
k_pas_RV   = (P_D_r .* y ./ (2 .* sinalpha_r .* Gamma_r) - sigmaact_r) ./ sigmapas_r;

%% Exercise 

a = exercise(1); 
b = exercise(2); 
c = exercise(3); 
d = exercise(4); 
e = exercise(5); 
f = exercise(6); 

k_act  = k_act  * (1 + a * theta); 
R_SA   = R_SA   / (1 + b * theta); 
C_SV   = C_SV   / (1 + c * theta); 
C_SA   = C_SA   / (1 + d * theta); 
R_PA   = R_PA   / (1 + e * theta); 
C_rest = C_rest / (1 + f * theta); 


pars = [Vw_LV; Vw_SEP; Vw_RV; 
    Amref_LV; Amref_SEP; Amref_RV; 
    C_Ao; C_SA; C_SV; C_PA; C_PV; 
    R_vlv; R_Ao; R_tAo; R_SA; R_tSA; R_SV; R_PA; R_PV; 
    Lsref; L_0; L_C; LSEiso; SLrest; 
    k_pas_LV; k_pas_RV; k_act; 
    gamma_1; gamma_2; v_max; C_rest]; 



