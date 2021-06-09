function [adjpars,fixpars,low,hi] = parameters(data)

P_mmHg2kPa = data.units.P_mmHg2kPa; 

Pbar  = data.Pbar;
DPbar = data.DPbar; 
SPbar = data.SPbar; 
Vtot  = data.Vtot; 
CO    = data.CO;

%% Stressed Volumes (m^3)

d_la = 0.005; 
d_lv = 0.02;
d_sa = 0.15;
d_sv = 0.59;
d_ra = 0.005; 
d_rv = 0.02; 
d_pa = 0.06;
d_pv = 0.15; 

V_la0 = d_la * Vtot; 
V_lv0 = d_lv * Vtot; 
V_sa0 = d_sa * Vtot;
V_sv0 = d_sv * Vtot;
V_ra0 = d_ra * Vtot;
V_rv0 = d_rv * Vtot; 
V_pa0 = d_pa * Vtot; 
V_pv0 = d_pv * Vtot;

%% Unstressed Volumes (m^3)

vaperc = .7; %Vu for arteries is ~72% of stressed volume - Beneken
V_sau = V_sa0 * vaperc;
V_pau = V_pa0 * vaperc; 

vvperc = .9; %Vu for veins is ~92% of stressed volume - Beneken
V_svu = V_sv0 * vvperc; 
V_pvu = V_pv0 * vvperc; 

%%  Stressed volumes

V_sas = V_sa0 - V_sau; 
V_svs = V_sv0 - V_svu; 
V_pas = V_pa0 - V_pau; 
V_pvs = V_pv0 - V_pvu; 

%% Pressures (kPa)
%Ratios from Boron book (mmHg converted to kPa)

SBP = 120 * P_mmHg2kPa;
MBP = 100 * P_mmHg2kPa;
DBP = 80  * P_mmHg2kPa; 

% LA
k_laM = (10 * P_mmHg2kPa) / SBP; 
k_lam = (2 * P_mmHg2kPa) / DBP; 

P_laM = k_laM * SPbar; 
P_lam = k_lam * DPbar; 

% LV
k_lvM = (125 * P_mmHg2kPa) / SBP;
k_lvm = (5 * P_mmHg2kPa) / DBP;

P_lvM = k_lvM * SPbar;
P_lvm = k_lvm * DPbar;

% SA
k_saM = (120 * P_mmHg2kPa) / SBP; 
k_sa  = (100 * P_mmHg2kPa) / MBP;
k_sam = (80  * P_mmHg2kPa) / DBP; 

P_saM   = k_saM * SPbar; 
P_sabar = k_sa  * Pbar; 
P_sam   = k_sam * DPbar; 

% SV
k_svM = (4.5 * P_mmHg2kPa) / SBP;
k_sv  = (4   * P_mmHg2kPa) / MBP;
k_svm = (3.5 * P_mmHg2kPa) / DBP;

P_svM   = k_svM * SPbar;
P_svbar = k_sv  * Pbar; 
P_svm   = k_svm * DPbar;

% RA 
k_raM = (10 * P_mmHg2kPa) / SBP; 
k_ram = (2.5  * P_mmHg2kPa) / DBP; 

P_raM = k_raM * SPbar; 
P_ram = k_ram * DPbar; 

% RV
k_rvM = (25 * P_mmHg2kPa) / SBP; 
k_rvm = (3  * P_mmHg2kPa) / DBP; 

P_rvM = k_rvM * SPbar;
P_rvm = k_rvm * DPbar;

% PA
k_paM = (20 * P_mmHg2kPa) / SBP; 
k_pa  = (15 * P_mmHg2kPa) / MBP;
k_pam = (10 * P_mmHg2kPa) / DBP; 

P_paM   = k_paM * SPbar; 
P_pabar = k_pa  * Pbar; 
P_pam   = k_pam * DPbar; 

% PV
k_pvM = (4.5 * P_mmHg2kPa) / SBP; 
k_pv  = (4   * P_mmHg2kPa) / MBP;
k_pvm = (3.5 * P_mmHg2kPa) / DBP; 

P_pvM   = k_pvM * SPbar; 
P_pvbar = k_pv  * Pbar; 
P_pvm   = k_pvm * DPbar; 

%% Compliances (mL mmHg^(-1) converted to m^3 kPa^(-1))

E_laM = P_laM / V_la0; 
E_lam = P_lam / V_la0; 
E_lvM = P_lvM / (.25 * V_lv0); % Assume maximal elastance occurs at 25% stressed volume 
E_lvm = P_lvm / V_lv0; 
E_raM = P_raM / V_ra0;
E_ram = P_ram / V_ra0; 
E_rvM = P_rvM / (.25 * V_rv0); % Assume maximal elastance occurs at 25% stressed volume 
E_rvm = P_rvm / V_rv0;  

C_sa = V_sas/P_saM / 1.25; 
C_sv = V_svs/P_svM; 
C_pa = V_pas/P_paM; 
C_pv = V_pvs/P_pvM; 

%% Resistances (mmHg s mL^(-1) converted to kPa s m^(-3))

R_sa  = (P_sam - P_svm) / CO * 1.25; 
R_sv  = (P_svm - P_ram) / CO; 
R_pa  = (P_pam - P_pvm) / CO; 
R_pv  = (P_pvm - P_lam) / CO; 

R_m_valve = 1e-3 / 7.5e-6; 
R_a_valve = 5e-3 / 7.5e-6;  
R_t_valve = 1e-3 / 7.5e-6;  
R_p_valve = 1e-3 / 7.5e-6; 

%% Elastance parameters  

k_TS = .1; % percentage of heart period from the end of diastole to maximal systole
k_TR = .3; % percentage of heart period due to relaxation from systole 

tau_v = .2; 

%% Outputs

adjpars = [E_laM; E_lam; E_lvM; E_lvm; E_raM; E_ram; E_rvM; E_rvm; 
    C_sa; C_sv; C_pa; C_pv; 
    R_sa; R_sv; R_pa; R_pv; 
    R_m_valve; R_a_valve; R_t_valve; R_p_valve; 
    k_TS; k_TR; tau_v; 
    ]; 

low = adjpars * (1 - .5); 
hi  = adjpars * (1 + .5); 

adjpars = log(adjpars); 
low     = log(low);
hi      = log(hi);

fixpars = [k_lam; k_lvm; k_sam; k_svm; k_ram; k_rvm; k_pam; k_pvm; 
    V_sau; V_svu; V_pau; V_pvu; 
    ]; 
    

end 