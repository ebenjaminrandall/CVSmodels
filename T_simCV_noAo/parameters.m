function [adjpars,fixpars] = parameters(data)

Pbar  = data.Pbar;
DPbar = data.DPbar; 
SPbar = data.SPbar; 
Vtot  = data.Vtot; 
CO    = data.CO;

Vw_LV_and_SEP = data.Vw_LV_and_SEP; 
Vw_RV         = data.Vw_RV;

Am_LV_and_SEP = data.Am_LV_and_SEP; 
Am_RV         = data.Am_RV; 

xm_lv0 = data.deformation.xm_lv0; 
ym0    = data.deformation.ym0; 

%% Stressed Volumes (m^3)

d_lv = 0.025;
d_ao = 0.04;
d_sa = .1 - d_ao; 
d_sv = 0.75; 
d_rv = 0.025; 
d_pa = 0.04;
d_pv = .1 - d_pa; %0.075; 

V_lv0 = d_lv*Vtot; 
V_ao0 = d_ao*Vtot; 
V_sa0 = d_sa*Vtot;
V_sv0 = d_sv*Vtot;
V_rv0 = d_rv*Vtot; 
V_pa0 = d_pa*Vtot; 
V_pv0 = d_pv*Vtot;

%% Unstressed Volumes (m^3)

vaperc = .72; %Vu for arteries is ~72% of Vtot - Beneken
V_aou = V_ao0*vaperc;
V_sau = V_sa0*vaperc;
V_pau = V_pa0*vaperc; 

vvperc = .92; %Vu for veins is ~92% of Vtot - Beneken
V_svu = V_sv0*vvperc; 
V_pvu = V_pv0*vvperc; 

%% Pressures (kPa)
%Ratios from Boron book (mmHg converted to kPa)

SBP = 120 / 7.5;
DBP = 80 / 7.5; 
MBP = 100 / 7.5;

% LV
k_lvM = (125 / 7.5) / SBP;
k_lvm = (3   / 7.5) / DBP;

P_lvM   = k_lvM * SPbar;
P_lvm   = k_lvm * DPbar;

% AO
k_aoM = (125 / 7.5) / SBP; 
k_ao  = (105 / 7.5) / MBP; 
k_aom = (85  / 7.5) / DBP;

P_aoM   = k_aoM * SPbar; 
P_aobar = k_ao  * Pbar; 
P_aom   = k_aom * DPbar; 

% SA
k_saM = (120 / 7.5) / SBP; 
k_sa  = (100 / 7.5) / MBP;
k_sam = (80  / 7.5) / DBP; 

P_saM   = k_saM * SPbar; 
P_sabar = k_sa  * Pbar; 
P_sam   = k_sam * DPbar; 

% SV
k_svM = (15 / 7.5) / SBP;
k_sv  = (10 / 7.5) / MBP;
k_svm = (5 / 7.5) / DBP;

P_svM   = k_svM * SPbar;
P_svbar = k_sv  * Pbar; 
P_svm   = k_svm * DPbar;

% RV
k_rvM = (25 / 7.5) / SBP; 
k_rvm = (3  / 7.5) / DBP; 

P_rvM   = k_rvM * SPbar;
P_rvm   = k_rvm * DPbar;

% PA
k_paM = (25 / 7.5) / SBP; 
k_pa  = (15 / 7.5) / MBP;
k_pam = (5 / 7.5) / DBP; 

P_paM   = k_paM * SPbar; 
P_pabar = k_pa  * Pbar; 
P_pam   = k_pam * DPbar; 

% PV
k_pvM = (5  / 7.5) / SBP; 
k_pv  = (4  / 7.5) / MBP;
k_pvm = (3  / 7.5) / DBP; 

P_pvM   = k_pvM * SPbar; 
P_pvbar = k_pv  * Pbar; 
P_pvm   = k_pvm * DPbar; 

%% Compliances (mL mmHg^(-1) converted to m^3 kPa^(-1))

E_lvm = P_lvm / V_lv0; % Minimal elastance of the left ventricle 
E_rvm = P_rvm / V_rv0; 

C_ao = (V_ao0 - V_aou)/P_aobar; 
C_sa = (V_sa0 - V_sau)/P_sabar; 
C_sv = (V_sv0 - V_svu)/P_svbar; 
C_pa = (V_pa0 - V_pau)/P_pabar; 
C_pv = (V_pv0 - V_pvu)/P_pvbar; 

%% Resistances (mmHg s mL^(-1) converted to kPa s m^(-3))

R_ao = (P_aobar - P_sabar)/CO; 
R_sa = (P_sabar - P_svbar)/CO;
R_sv = (P_svbar - P_rvm)/CO;
R_pa = (P_pabar - P_pvbar)/CO; 
R_pv = (P_pvbar - P_lvm)/CO;
R_v  = 1e-4 / 7.5 / 1e-6;

%% Heart model 

% Wall volume of ventricular wall segment (m^3)
Vw_lv  = Vw_LV_and_SEP * 2/3; 
Vw_sep = Vw_LV_and_SEP * 1/3; 
Vw_rv  = Vw_RV; 

% Reference midwall surface area (m^2)
Amref_lv  = Am_LV_and_SEP * 2/3;
Amref_sep = Am_LV_and_SEP * 1/3; 
Amref_rv  = Am_RV; 

% Time scale (convet from ms to s) 
tauR  = 48 * 1e-3; 
tauD  = 32 * 1e-3; 
tausc = 425 * 1e-3; 

% Sarcomere length parameters (convert from �m to m)
Lsref   = 2 * 1e-6; 
Lsc0    = 1.51 * 1e-6; 
Lse_iso = 0.04 * 1e-6; 

v_max   = 7 * 1e-6;    % convert from um s^(-1) to m s^(-1) % sarcomere length shortening velocity
Ca_rest = 0.02;        % dimensionless %diastolic resting level of activation

%% Calculate force parameters 

% Diastole
xm_lv_d = xm_lv0; 
ym_d    = ym0; 

Am_lv_d  = pi * (xm_lv_d^2  + ym_d^2);
Cm_lv_d  = - 2 * xm_lv_d  / (xm_lv_d^2  + ym_d^2);
z_lv_d   = 3 * Cm_lv_d  * Vw_lv  / (2 * Am_lv_d); 
eps_lv_d = 0.5 * log( Am_lv_d  / Amref_lv  ) - (1/12) * z_lv_d^2  - 0.019 * z_lv_d^4; 

sigma_pas_lv_d = 36 * max(0,eps_lv_d - 0.1)^2  + 0.1 * (eps_lv_d  - 0.1) + 0.0025 * exp(30 * eps_lv_d); 
Gamma_lv_d     = - (2 / 3) * z_lv_d * (1 + (1 / 3) * z_lv_d^2 + (1 / 5) * z_lv_d^4);

% Force scaling factors (kPa)
k_pas = P_lvm / (Gamma_lv_d * sigma_pas_lv_d);
k_act = 120; 

%% Outputs

adjpars = [E_lvm; E_rvm;
    C_ao; C_sa; C_sv; C_pa; C_pv; 
    R_ao; R_sa; R_sv; R_pa; R_pv; R_v; 
    Vw_lv; Vw_sep; Vw_rv; 
    Amref_lv; Amref_sep; Amref_rv; 
    tauR; tauD; tausc; 
    Lsref; Lsc0; Lse_iso; 
    k_act; k_pas; 
    v_max; 
    Ca_rest; 
    V_aou; V_sau; V_svu; V_pau; V_pvu; 
    ]; 

fixpars = [k_lvm; k_ao; k_sa; k_sv; k_rvm; k_pa; k_pv; 
    ]; 
    

end 