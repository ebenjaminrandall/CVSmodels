function [adjpars,fixpars,low,hi] = parameters(data)

Pbar  = data.Pbar;
DPbar = data.DPbar; 
SPbar = data.SPbar; 
Vtot  = data.Vtot; 
CO    = data.CO;

%% Stressed Volumes (m^3)

d_la = 0.005; 
d_lv = 0.02;
d_ao = 0.04;
d_sa = 0.15;
d_sv = 0.55;
d_ra = 0.005; 
d_rv = 0.02; 
d_pa = 0.06;
d_pv = 0.15; 

V_la0 = d_la * Vtot; 
V_lv0 = d_lv * Vtot; 
V_ao0 = d_ao * Vtot; 
V_sa0 = d_sa * Vtot;
V_sv0 = d_sv * Vtot;
V_ra0 = d_ra * Vtot;
V_rv0 = d_rv * Vtot; 
V_pa0 = d_pa * Vtot; 
V_pv0 = d_pv * Vtot;

%% Unstressed Volumes (m^3)

vaperc = .7; %Vu for arteries is ~72% of stressed volume - Beneken
V_aou = V_ao0 * vaperc;
V_sau = V_sa0 * vaperc;
V_pau = V_pa0 * vaperc; 

vvperc = .9; %Vu for veins is ~92% of stressed volume - Beneken
V_svu = V_sv0 * vvperc; 
V_pvu = V_pv0 * vvperc; 

% Stressed volumes
V_aos = V_ao0 - V_aou; 
V_sas = V_sa0 - V_sau; 
V_svs = V_sv0 - V_svu; 
V_pas = V_pa0 - V_pau; 
V_pvs = V_pv0 - V_pvu; 

%% Pressures (kPa)
%Ratios from Boron book (mmHg converted to kPa)

SBP = 120 / 7.5;
MBP = 100 / 7.5;
DBP = 80  / 7.5; 

% LA
k_laM = (10 / 7.5) / SBP; 
k_lam = (2 / 7.5) / DBP; 

P_laM = k_laM * SPbar; 
P_lam = k_lam * DPbar; 

% LV
k_lvM = (125 / 7.5) / SBP;
k_lvm = (5 / 7.5) / DBP;

P_lvM = k_lvM * SPbar;
P_lvm = k_lvm * DPbar;

% AO
k_aoM = (121 / 7.5) / SBP; 
k_ao  = (101 / 7.5) / MBP; 
k_aom = (81  / 7.5) / DBP;

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
k_svM = (4.5 / 7.5) / SBP;
k_sv  = (4 / 7.5) / MBP;
k_svm = (3.5 / 7.5) / DBP;

P_svM   = k_svM * SPbar;
P_svbar = k_sv  * Pbar; 
P_svm   = k_svm * DPbar;

% RA 
k_raM = (10 / 7.5) / SBP; 
k_ram = (2.5 / 7.5) / DBP; 

P_raM = k_raM * SPbar; 
P_ram = k_ram * DPbar; 

% RV
k_rvM = (25 / 7.5) / SBP; 
k_rvm = (3  / 7.5) / DBP; 

P_rvM = k_rvM * SPbar;
P_rvm = k_rvm * DPbar;

% PA
k_paM = (20 / 7.5) / SBP; 
k_pa  = (15 / 7.5) / MBP;
k_pam = (10 / 7.5) / DBP; 

P_paM   = k_paM * SPbar; 
P_pabar = k_pa  * Pbar; 
P_pam   = k_pam * DPbar; 

% PV
k_pvM = (4.5  / 7.5) / SBP; 
k_pv  = (4  / 7.5) / MBP;
k_pvm = (3.5  / 7.5) / DBP; 

P_pvM   = k_pvM * SPbar; 
P_pvbar = k_pv  * Pbar; 
P_pvm   = k_pvm * DPbar; 

%% Compliances (mL mmHg^(-1) converted to m^3 kPa^(-1))

E_laM = P_laM / V_la0; 
E_lam = P_lam / V_la0; 
E_raM = P_raM / V_ra0;
E_ram = P_ram / V_ra0; 

E_lvm = P_lvm / V_lv0; % Minimal elastance of the left ventricle 
E_rvm = P_rvm / V_rv0; 

C_ao = V_aos/P_aoM; 
C_sa = V_sas/P_saM; 
C_sv = V_svs/P_svM; 
C_pa = V_pas/P_paM; 
C_pv = V_pvs/P_pvM; 

%% Resistances (mmHg s mL^(-1) converted to kPa s m^(-3))

R_v   = 1e-3 / 7.5e-6; 
R_tao = .0001 / C_ao; 
R_ao  = (P_aom - P_sam) / CO; 
R_tsa = .1 / C_sa; 
R_sa  = (P_sam - P_svm) / CO; 
R_sv  = (P_svm - P_ram) / CO; 
R_pa  = (P_pam - P_pvm) / CO; 
R_pv  = (P_pvm - P_lam) / CO; 

R_m_valve = 5e-3 / 7.5e-6; 
R_a_valve = 5e-3 / 7.5e-6; 
R_t_valve = 1e-3 / 7.5e-6; 
R_p_valve = 1e-4 / 7.5e-6; 

%% Heart model 

% Individual wall volumes % convert from mL to m^3
Vw_LV_SEP = 80 * 1e-6; 
Vw_RV     = 10  * 1e-6;  

% Sarcomere length parameters (convert from µm to m)
Lsref   = 2 * 1e-6; 
Lsc0    = 1.51 * 1e-6; 
Lse_iso = 0.04 * 1e-6; 

v_max   = 7 * 1e-6;    % convert from um s^(-1) to m s^(-1) % sarcomere length shortening velocity

%% Calculate force parameters from diastolic measurements

Vw_lv    = Vw_LV_SEP * 2/3; 
Vw_rv    = Vw_RV; 

xm_lv_d = 4.5 * 1e-2; 
xm_rv_d = 5 * 1e-2; 
ym_d    = 3.3 * 1e-2; 

Am_lv_d  = pi * (xm_lv_d^2  + ym_d^2);
Am_rv_d  = pi * (xm_rv_d^2  + ym_d^2);

Cm_lv_d  = - 2 * xm_lv_d  / (xm_lv_d^2  + ym_d^2);
Cm_rv_d  = - 2 * xm_rv_d  / (xm_rv_d^2  + ym_d^2);

z_lv_d   = 3 * Cm_lv_d  * Vw_lv  / (2 * Am_lv_d); 
z_rv_d   = 3 * Cm_rv_d  * Vw_rv  / (2 * Am_rv_d); 

A1 = Am_lv_d * exp( -2 * (1/12 * z_lv_d^2 + 0.019 * z_lv_d^4));
A2 = Am_rv_d * exp( -2 * z_rv_d^2 * (1/12 + 0.019 * z_rv_d^2));

% Individual reference surface area % convert from cm^2 to m^2
Amref_lv  = .9 * A1; 
Am_LV_SEP = Amref_lv * 3/2;
Am_RV     = .9 * A2; 

eps_lv_d = 0.5 * log( Am_lv_d  / Amref_lv  ) - (1/12) * z_lv_d^2  - 0.019 * z_lv_d^4; 
sigma_pas_lv_d = eps_lv_d + exp(40 * eps_lv_d); 
Gamma_lv_d     = - (2 / 3) * z_lv_d * (1 + (1 / 3) * z_lv_d^2 + (1 / 5) * z_lv_d^4);

% Force scaling factors (kPa)
k_pas = P_lvm / (Gamma_lv_d * sigma_pas_lv_d);

k_act_lv = 120 * 8;  
k_act_rv = 120 * 8;  

%% Activation function parameters 

k_TS = .1; % percentage of heart period from the end of diastole to maximal systole
k_TR = .3; % percentage of heart period due to relaxation from systole 

tau_v = .2; % ventricles delayed by 0.2 s after atrial contraction  

%% Outputs

adjpars = [E_laM; E_lam; E_raM; E_ram; E_lvm; E_rvm; 
    C_ao; C_sa; C_sv; C_pa; C_pv; 
    R_ao; R_sa; R_sv; R_pa; R_pv; 
    R_m_valve; R_a_valve; R_t_valve; R_p_valve; 
    Vw_LV_SEP; Vw_RV; 
    Am_LV_SEP; Am_RV; 
    Lsref; Lsc0; Lse_iso; 
    k_pas; k_act_lv; k_act_rv;  
    v_max; 
    k_TS; k_TR; tau_v; 
    ]; 

low = adjpars * (1 - .5); 
hi  = adjpars * (1 + .5); 

adjpars = log(adjpars); 
low     = log(low);
hi      = log(hi);

fixpars = [k_lam; k_lvm; k_aom; k_sam; k_svm; k_ram; k_rvm; k_pam; k_pvm; 
    V_aou; V_sau; V_svu; V_pau; V_pvu; 
    ]; 
    

end 