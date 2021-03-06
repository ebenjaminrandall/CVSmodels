function [adjpars,fixpars] = parameters(data)


Vtot  = data.Vtot; 
CO    = data.CO;

Vw_LV_and_SEP = data.Vw_LV_and_SEP; 
Vw_RV         = data.Vw_RV;

Am_LV_and_SEP = data.Am_LV_and_SEP; 
Am_RV         = data.Am_RV; 

xm_lv0 = data.deformation.xm_lv0; 
ym0    = data.deformation.ym0; 

%% Stressed Volumes (m^3)

d_la = 0.005; 
d_lv = 0.02;
d_sa = .12; 
d_sv = 0.62;
d_ra = 0.005;
d_rv = 0.02; 
d_pa = 0.06;
d_pv = .15; 

V_la0 = d_la*Vtot; 
V_lv0 = d_lv*Vtot; 
V_sa0 = d_sa*Vtot;
V_sv0 = d_sv*Vtot;
V_ra0 = d_ra*Vtot; 
V_rv0 = d_rv*Vtot; 
V_pa0 = d_pa*Vtot; 
V_pv0 = d_pv*Vtot;

%% Unstressed Volumes (m^3)

vaperc = .7; %Vu for arteries is ~72% of Vtot - Beneken
V_sau = V_sa0*vaperc;
V_pau = V_pa0*vaperc; 

vvperc = .9; %Vu for veins is ~92% of Vtot - Beneken
V_svu = V_sv0*vvperc; 
V_pvu = V_pv0*vvperc; 

% Stressed volumes
V_sas = V_sa0 - V_sau; 
V_svs = V_sv0 - V_svu; 
V_pas = V_pa0 - V_pau; 
V_pvs = V_pv0 - V_pvu; 

%% Pressures (kPa)
%Ratios from Boron book (mmHg converted to kPa)

% LA
P_laM = 10 / 7.5; 
P_lam = 2 / 7.5; 

% LV
P_lvM   = 125 / 7.5;
P_lvm   = 5   / 7.5;

% SA
P_saM   = 120 / 7.5; 
P_sabar = 100 / 7.5; 
P_sam   = 80  / 7.5; 

% SV
P_svM   = 6 / 7.5;
P_svbar = 4 / 7.5; 
P_svm   = 2 / 7.5;

% RA 
P_raM = 10 / 7.5; 
P_ram = 2 / 7.5; 

% RV
P_rvM   = 25 / 7.5;
P_rvm   = 5  / 7.5;

% PA
P_paM   = 20 / 7.5; 
P_pabar = 15 / 7.5; 
P_pam   = 10 / 7.5; 

% PV
P_pvM   = 6  / 7.5; 
P_pvbar = 4  / 7.5; 
P_pvm   = 2  / 7.5; 

%% Compliances (mL mmHg^(-1) converted to m^3 kPa^(-1))

E_laM = P_laM / V_la0;  
E_lam = P_lam / V_la0; 
E_raM = P_raM / V_ra0;
E_ram = P_ram / V_ra0; 

E_lvm = P_lvm / V_lv0; % Minimal elastance of the left ventricle 
E_rvm = P_rvm / V_rv0; 

C_sa = V_sas/P_saM;
C_sv = V_svs/P_svM; 
C_pa = V_pas/P_paM; 
C_pv = V_pvs/P_pvM; 

%% Resistances (mmHg s mL^(-1) converted to kPa s m^(-3))

R_sa = (P_sabar - P_svbar)/CO;
R_sv = (P_svbar - P_ram)/CO;
R_pa = (P_pabar - P_pvbar)/CO; 
R_pv = (P_pvbar - P_lam)/CO;

% % Parameter values in Wood units 
% R_sa * 7.5e-6 * 1000/60
% R_sv * 7.5e-6 * 1000/60
%R_pa * 7.5e-6 * 1000/60
% R_pv * 7.5e-6 * 1000/60

R_m_valve = 1e-3 / 7.5e-6; %1e-4 / 7.5 / 1e-6;
R_a_valve = 5e-3 / 7.5e-6; 
R_t_valve = 1e-3 / 7.5e-6; 
R_p_valve = 1e-3 / 7.5e-6; 

%% Heart model 

% Wall volume of ventricular wall segment (m^3)
Vw_lv  = Vw_LV_and_SEP * 2/3; 
Vw_sep = Vw_LV_and_SEP * 1/3; 
Vw_rv  = Vw_RV; 

% Reference midwall surface area (m^2)
Amref_lv  = Am_LV_and_SEP * 2/3;
Amref_sep = Am_LV_and_SEP * 1/3; 
Amref_rv  = Am_RV; 

% Sarcomere length parameters (convert from ?m to m)
Lsref   = 2 * 1e-6; 
Lsc0    = 1.51 * 1e-6; 
Lse_iso = 0.04 * 1e-6; 

v_max   = 7 * 1e-6;    % convert from um s^(-1) to m s^(-1) % sarcomere length shortening velocity

% Collagen parameters
alpha_c = 100; %8*19; 
beta_c  = 10; %3*2.93; 

%% Calculate force parameters 

% Diastole
xm_lv_d = xm_lv0; 
ym_d    = ym0; 

Am_lv_d  = pi * (xm_lv_d^2  + ym_d^2);
Cm_lv_d  = - 2 * xm_lv_d  / (xm_lv_d^2  + ym_d^2);
z_lv_d   = 3 * Cm_lv_d  * Vw_lv  / (2 * Am_lv_d); 
eps_lv_d = 0.5 * log( Am_lv_d  / Amref_lv  ) - (1/12) * z_lv_d^2  - 0.019 * z_lv_d^4; 

Ls_lv_d  = Lsref * exp(eps_lv_d); 

sigma_coll_lv  =  alpha_c * ((Ls_lv_d - Lsc0)/1e-6)^beta_c; 
sigma_pas_lv_d =  (Ls_lv_d  - Lsc0)/1e-6 + sigma_coll_lv; 

Gamma_lv_d     = - (2 / 3) * z_lv_d * (1 + (1 / 3) * z_lv_d^2 + (1 / 5) * z_lv_d^4);

% Force scaling factors (kPa)
k_pas    = P_lvm / (Gamma_lv_d * sigma_pas_lv_d);
k_act_lv = 120; 
k_act_rv = 120; 

%% Atrial offset 

tau_a = 0.2; 

k_TS = 0.1; 
k_TR = 0.3; 

%% Outputs

adjpars = [E_laM; E_lam; E_lvm; E_raM; E_ram; E_rvm; 
    C_sa; C_sv; C_pa; C_pv; 
    R_sa; R_sv; R_pa; R_pv; 
    R_m_valve; R_a_valve; R_t_valve; R_p_valve;  
    Vw_lv; Vw_sep; Vw_rv; 
    Amref_lv; Amref_sep; Amref_rv; 
    Lsref; Lsc0; Lse_iso; 
    k_act_lv; k_act_rv; k_pas; 
    v_max; 
    V_sau; V_svu; V_pau; V_pvu; 
    tau_a; 
    k_TS; k_TR; 
    alpha_c; beta_c; 
    ]; 

fixpars = [P_lam; P_lvm; P_sam; P_svm; P_ram; P_rvm; P_pam; P_pvm; 
    ]; 
    

end 