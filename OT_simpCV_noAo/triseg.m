function f = triseg(x,pars,data,init,P_lv,P_rv)

HR = data.HR.HR_rest; 

%% Parameters 

% Wall volume of ventricular wall segment (m^3)
Vw_LV_SEP = pars(19); 
Vw_RV     = pars(20); 

% Reference midwall surface area (m^2)
Am_LV_SEP = pars(21); 
Am_RV     = pars(22); 

% Sarcomere length parameters (m)
Lsref   = pars(23);
Lsc0    = pars(24); 
Lse_iso = pars(25); 

% Force scaling factors (kPa) 
k_pas    = pars(26); 
k_act_lv = pars(27); 
k_act_rv = pars(28); 

v_max = pars(29); 

%% Variables 

xm_lv  = x(1); 
xm_sep = x(2); 
xm_rv  = x(3); 
ym     = x(4); 

% Contractile element length 
Lsc_lv  = x(5); 
Lsc_sep = x(6);  
Lsc_rv  = x(7);

% Volumes 
V_lv = init(8); 
V_rv = init(12);

%% Activation 

y = 0; 

%% Heart and sarcomere model 

% Wall volume of ventricular wall segment (m^3)
Vw_lv  = Vw_LV_SEP * 2/3; 
Vw_sep = Vw_LV_SEP * 1/3; 
Vw_rv  = Vw_RV; 

% Reference midwall surface area (m^2)
Amref_lv  = Am_LV_SEP * 2/3;
Amref_sep = Am_LV_SEP * 1/3; 
Amref_rv  = Am_RV; 

% Volume of spherical cap formed by midwall surface (m^3)
Vm_lv  = - (pi / 6) * xm_lv  * (xm_lv^2  + 3 * ym^2); 
Vm_sep = (pi / 6) * xm_sep * (xm_sep^2 + 3 * ym^2); 
Vm_rv  = (pi / 6) * xm_rv  * (xm_rv^2  + 3 * ym^2); 

% Surface area of midwall surface (m^2) 
Am_lv  = pi * (xm_lv^2  + ym^2);
Am_sep = pi * (xm_sep^2 + ym^2); 
Am_rv  = pi * (xm_rv^2  + ym^2); 

% Curvature of midwall surface (m^(-1))
Cm_lv  = - 2 * xm_lv  / (xm_lv^2  + ym^2);
Cm_sep = 2 * xm_sep / (xm_sep^2 + ym^2);
Cm_rv  = 2 * xm_rv  / (xm_rv^2  + ym^2);

% Ratio of wall thickness to midwall radius of curvature (dimensionless)
z_lv  = 3 * Cm_lv  * Vw_lv  / (2 * Am_lv); 
z_sep = 3 * Cm_sep * Vw_sep / (2 * Am_sep); 
z_rv  = 3 * Cm_rv  * Vw_rv  / (2 * Am_rv);

% Myofiber strain (dimensionless)
eps_lv  = 0.5 * log( Am_lv  / Amref_lv  ) - (1/12) * z_lv^2  - 0.019 * z_lv^4; 
eps_sep = 0.5 * log( Am_sep / Amref_sep ) - (1/12) * z_sep^2 - 0.019 * z_sep^4; 
eps_rv  = 0.5 * log( Am_rv  / Amref_rv  ) - (1/12) * z_rv^2  - 0.019 * z_rv^4; 

% Sarcomere length (m)
Ls_lv  = Lsref * exp(eps_lv); 
Ls_sep = Lsref * exp(eps_sep); 
Ls_rv  = Lsref * exp(eps_rv); 

% Active stress (kPa)
sigma_act_lv  = k_act_lv * y * (Lsc_lv  - Lsc0) / 1e-6 * (Ls_lv  - Lsc_lv)  / Lse_iso;  
sigma_act_sep = k_act_lv * y * (Lsc_sep - Lsc0) / 1e-6 * (Ls_sep - Lsc_sep) / Lse_iso;
sigma_act_rv  = k_act_rv * y * (Lsc_rv  - Lsc0) / 1e-6 * (Ls_rv  - Lsc_rv)  / Lse_iso; 

% Passive stress (kPa)
sigma_pas_lv  = k_pas * (eps_lv  + exp(40 * eps_lv)); 
sigma_pas_sep = k_pas * (eps_sep + exp(40 * eps_sep)); 
sigma_pas_rv  = k_pas * (eps_rv  + exp(40 * eps_rv)); 

% Total stress (kPa)
sigma_lv  = sigma_act_lv  + sigma_pas_lv; 
sigma_sep = sigma_act_sep + sigma_pas_sep; 
sigma_rv  = sigma_act_rv  + sigma_pas_rv; 

% Representative midwall tension (mmHg cm^(-2))
Tm_lv  = (Vw_lv  * sigma_lv  / (2 * Am_lv))  * (1 + (z_lv^2)/3  + (z_lv^4)/5); 
Tm_sep = (Vw_sep * sigma_sep / (2 * Am_sep)) * (1 + (z_sep^2)/3 + (z_sep^4)/5); 
Tm_rv  = (Vw_rv  * sigma_rv  / (2 * Am_rv))  * (1 + (z_rv^2)/3  + (z_rv^4)/5);

% Axial midwall tension component 
Tx_lv  = - Tm_lv  * 2 * xm_lv  * ym / (xm_lv^2  + ym^2); 
Tx_sep = Tm_sep * 2 * xm_sep * ym / (xm_sep^2 + ym^2); 
Tx_rv  = Tm_rv  * 2 * xm_rv  * ym / (xm_rv^2  + ym^2); 

% Radial midwall tension component 
Ty_lv  = Tm_lv  * (-xm_lv^2  + ym^2) / (xm_lv^2  + ym^2); 
Ty_sep = Tm_sep * (-xm_sep^2 + ym^2) / (xm_sep^2 + ym^2); 
Ty_rv  = Tm_rv  * (-xm_rv^2  + ym^2) / (xm_rv^2  + ym^2);

% Ventricular pressure (kPa)
ptrans1 = 2 * Tx_lv / ym; 
ptrans3 = 2 * Tx_rv / ym;  

f1 = -V_lv - 0.5 * Vw_lv - 0.5 * Vw_sep + Vm_sep - Vm_lv; 
f2 = Tx_lv + Tx_sep + Tx_rv;
f3 = V_rv + 0.5 * Vw_rv + 0.5 * Vw_sep + Vm_sep - Vm_rv;
f4 = Ty_lv + Ty_sep + Ty_rv; 

f5 = ((Ls_lv  - Lsc_lv)  /Lse_iso - 1) * v_max;
f6 = ((Ls_sep - Lsc_sep) /Lse_iso - 1) * v_max;
f7 = ((Ls_rv  - Lsc_rv)  /Lse_iso - 1) * v_max;

f8 = P_rv - ptrans3;  
f9 = P_lv + ptrans1; 

f = [f1; f2; f3; f4; f5; f6; f7];% f8; f9]; 
end 

