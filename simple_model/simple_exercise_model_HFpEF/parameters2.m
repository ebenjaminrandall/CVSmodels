function [pars,lb,ub] = parameters2(exercise,HF,theta,HFpEFdata,data) 

TotBV = 5e3; 

%% Exercise factors 

a = exercise(1); 
b = exercise(2); 
c = exercise(3); 
d = exercise(4); 
e = exercise(5); 
f = exercise(6); 

%% HF 

z = HF(1); 
y = HF(2); 
x = HF(3); 
w = HF(4); 

%% HFpEF data 

P_SA_syst_data = HFpEFdata.P_SAs; 
P_SA_dias_data = HFpEFdata.P_SAd; 
P_RA_syst_data = HFpEFdata.P_RAs; 
P_RA_dias_data = HFpEFdata.P_RAd; 
P_RV_syst_data = HFpEFdata.P_RVs; 
P_RV_dias_data = HFpEFdata.P_RVd; 
P_PA_syst_data = HFpEFdata.P_PAs;
P_PA_dias_data = HFpEFdata.P_PAd; 
P_PV_mean_data = HFpEFdata.PCWP; 

if strcmp(data.procedure,'RHC') == 1 
    CO = HFpEFdata.CO_RHC * 1e3 / 60; 
else 
    CO = HFpEFdata.LVCO_MRI * 1e3 / 60; 
end 

%% Recalculate stressed volumes 

vfactor = 0.4; 

VsB_sa = 160;           VuB_sa = 425;           VtB_sa = 585;
VsB_sv = 219;           VuB_sv = 2697;          VtB_sv = 2916;
VsB_pa = 69;            VuB_pa = 50;            VtB_pa = 119;
VsB_pv = 54;            VuB_pv = 460;           VtB_pv = 514;
VsB_lv = 125;           VuB_lv = 0;             VtB_lv = 125;
VsB_la = 50;            VuB_la = 30;            VtB_la = 80;
VsB_rv = 125;           VuB_rv = 0;             VtB_rv = 125;
VsB_ra = 50;            VuB_ra = 30;            VtB_ra = 80;
VsB_tot = 852;          VuB_tot = 3692;         VtB_tot = 4544;    

VsBNew_lv = 0.70 * VtB_lv;
VsBNew_la = 0.50 * VtB_la;
VsBNew_rv = 0.70 * VtB_rv;
VsBNew_ra = 0.50 * VtB_ra;

VsB_nh    = VsB_tot - VsB_lv - VsB_la - VsB_rv - VsB_ra;
VsBP_tot  = vfactor * VtB_tot;
VsBP_nh   = VsBP_tot - VsBNew_lv - VsBNew_la - VsBNew_rv - VsBNew_ra;
VsRec_tot = VsBP_nh - VsB_nh;

VsBRec_sa = VsRec_tot * (VuB_sa / VuB_tot);
VsBRec_sv = VsRec_tot * (VuB_sv / VuB_tot);
VsBRec_pa = VsRec_tot * (VuB_pa / VuB_tot);
VsBRec_pv = VsRec_tot * (VuB_pv / VuB_tot);

VsBNew_sa = VsB_sa + VsBRec_sa;
VsBNew_sv = VsB_sv + VsBRec_sv;
VsBNew_pa = VsB_pa + VsBRec_pa;
VsBNew_pv = VsB_pv + VsBRec_pv;

Frac_sa = VsBNew_sa / VtB_sa;
Frac_sv = VsBNew_sv / VtB_sv;
Frac_pa = VsBNew_pa / VtB_pa;
Frac_pv = VsBNew_pv / VtB_pv;

CircBV = vfactor * TotBV;

Vsa_tot = CircBV * (VtB_sa / VtB_tot);
Vsv_tot = CircBV * (VtB_sv / VtB_tot);
Vpa_tot = CircBV * (VtB_pa / VtB_tot);
Vpv_tot = CircBV * (VtB_pv / VtB_tot);

Vs_SA = Frac_sa * Vsa_tot;
Vs_SV = Frac_sv * Vsv_tot;
Vs_PA = Frac_pa * Vpa_tot;
Vs_PV = Frac_pv * Vpv_tot;

%% Heart parameters 

Vw_LV  = 80;
Vw_SEP = 38; 
Vw_RV  = 28; % Heart wall volumes (mL)

Amref_LV  = 52; %0.975*80; % LV midwall reference surface area, cm^2
Amref_SEP = 50; %0.975*45; % SEP midwall reference surface area, cm^2
Amref_RV  = 128; %1.12*100; % RV midwall reference surface area, cm^2

%% Lumped circulatory parameters

P_RA_pp = P_RA_syst_data - P_RA_dias_data;
P_SA_pp = P_SA_syst_data - P_SA_dias_data; 
P_SV_pp = 0.05 * P_SA_pp; 
P_PA_pp = P_PA_syst_data - P_PA_dias_data; 
P_PV_pp = 3/7 * P_PV_mean_data; 


P_RA_m = 1/3 * P_RA_syst_data + 2/3 * P_RA_dias_data; 
P_SA_m = 1/3 * P_SA_syst_data + 2/3 * P_SA_dias_data; 
P_SV_m = 1.025 * P_RA_m; 
P_PA_m = 1/3 * P_PA_syst_data + 2/3 * P_PA_dias_data; 
P_PV_m = P_PV_mean_data; 
P_LA_m = .975 * P_PV_m; 

C_Ao = 0.65; 
C_SA = Vs_SA / P_SA_pp; 
C_SV = Vs_SV / P_SV_pp; 
C_PA = Vs_PA / P_PA_pp; 
C_PV = Vs_PV / P_PV_pp; 

R_Ao = 0.01; 
R_SA = (P_SA_m - P_SV_m) / CO; 
R_RA = (P_SV_m - P_RA_m) / CO; 
R_PA = (P_PA_m - P_PV_m) / CO; 
R_LA = (P_PV_m - P_LA_m) / CO; 


% % Lumped circulatory parameters
% C_Ao = 0.65;  % Proximal aortic compliance, mL/mmHg
% C_SA = 1 / (1 + d * theta); %1.65/(1 + d*theta)/1.2; % Systemic arterial compliance, mL/mmHg
% C_SV = 1.4*250/(1 + c*theta); % Systemic venous compliance, mL/mmHg 
% C_PA = 4; %5.4; % Pulmonary arterial compliance, mL/mmHg
% C_PV = 25; % Pulmonary venous compliance, mL/mmHg

%R_Ao  = 0.01; % resistance of aorta , mmHg*sec/mL
%R_SA  = (1.1 / (1 + b * theta); %1.1*0.965/(1 + b*theta);% mmHg*sec/mL; % Systemic vasculature resistance, mmHg*sec/mL
%R_PA  = 0.05 / (1 + e * theta); % Pulmonary vasculature resistance, mmHg*sec/mL 
R_vlv = 0.002; %  valve resistance, mmHg*sec/mL
R_tAo = 0.0020;
R_tSA = 0.05;
%R_RA  = 0.024;
%R_LA  = 0.024;

% Atrial elastances 
Emin = P_RA_pp / 35; %1.5*0.050 ;
Emax = P_RA_pp / 10; %0.150 ;


% Triseg parameters
Lsref    = 1.9; % Resting SL, micron
vmax     = 7; % micron/sec
LSEiso   = 0.04; % micron
k_pas_lv = 22; 
k_pas_rv = 22; 
k_act    = 600; %7.5*96*(1 + a*theta); % mmHg 
SLrest   = 1.51; % microns
Crest    = 0.02;

%% HF 

C_SA = C_SA * z; 
R_SA = R_SA * y; 
k_pas_lv = k_pas_lv * x; 
k_pas_rv = k_pas_rv * w; 

%% Exercise 

C_SA  = C_SA  / (1 + d * theta); 
C_SV  = C_SV  / (1 + c * theta); 
R_SA  = R_SA  / (1 + b * theta); 
R_PA  = R_PA  / (1 + e * theta); 
k_act = k_act * (1 + a * theta); 
Crest = Crest / (1 + f * theta); 

%% Outputs 

pars = [Vw_LV; Vw_SEP; Vw_RV;                           %1-3
    Amref_LV; Amref_SEP; Amref_RV;                      %4-6
    C_Ao; C_SA; C_SV; C_PA; C_PV;                       %7-11
    R_vlv; R_LA; R_Ao; R_tAo; R_SA; R_tSA; R_RA; R_PA;  %12-19
    Emin; Emax;                                         %20-21
    Lsref; LSEiso; SLrest; vmax; k_pas_lv; k_pas_rv; k_act; Crest;   %22-29
    1; 
    ];

pars = log(pars); 
lb   = pars - log(10); 
ub   = pars + log(10); 

end 