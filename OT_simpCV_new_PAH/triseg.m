function f = triseg(x,pars,data,init)

HR = data.HR; 

%% Parameters 

% Wall volume of ventricular wall segment (m^3)
Vw_lv  = pars(19); 
Vw_sep = pars(20); 
Vw_rv  = pars(21); 

% Reference midwall surface area (m^2)
Amref_lv  = pars(22); 
Amref_sep = pars(23); 
Amref_rv  = pars(24); 

% Sarcomere length parameters (?m)
Lsref   = pars(25);
Lsc0    = pars(26); 
Lse_iso = pars(27); 

% Force scaling factors (mmHg)
k_act = pars(28); 
k_pas = pars(29);

k_TS = pars(36); 
k_TR = pars(37); 

%% Variables 

xm_lv  = x(1); 
xm_sep = x(2); 
xm_rv  = x(3); 
ym     = x(4); 

% Contractile element length 
Lsc_lv  = init(5); 
Lsc_sep = init(6); 
Lsc_rv  = init(7); 

% Volumes 
V_lv = init(9); 
V_rv = init(13);

%% Activation function 

T = 60/HR; 
TS = k_TS * T; 
TR = k_TR * T; 

t = 0; 
tc_v = mod(t,T);
if tc_v >= 0 && tc_v < TS 
    y_v = 0.5*(1 - cos(pi*tc_v/TS)); 
elseif tc_v >= TS && tc_v < TR + TS 
    y_v = 0.5*(1 + cos(pi*(tc_v - TS)/TR)); 
else
    y_v = 0; 
end 

%% Heart and sarcomere model 

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
sigma_act_lv  = y_v * (Lsc_lv  - Lsc0) / 1e-6 * (Ls_lv  - Lsc_lv)  / Lse_iso; 
sigma_act_sep = y_v * (Lsc_sep - Lsc0) / 1e-6 * (Ls_sep - Lsc_sep) / Lse_iso;
sigma_act_rv  = y_v * (Lsc_rv  - Lsc0) / 1e-6 * (Ls_rv  - Lsc_rv)  / Lse_iso;

% Passive stress (kPa)
% sigma_pas_lv  = 36 * max(0,eps_lv  - 0.1)^2 + 0.1 * (eps_lv  - 0.1) + 0.0025 * exp(30 * eps_lv); 
% sigma_pas_sep = 36 * max(0,eps_sep - 0.1)^2 + 0.1 * (eps_sep - 0.1) + 0.0025 * exp(30 * eps_sep); 
% sigma_pas_rv  = 36 * max(0,eps_rv  - 0.1)^2 + 0.1 * (eps_rv  - 0.1) + 0.0025 * exp(30 * eps_rv);  

PConcollagen = 19; 
PExpcollagen = 2.93; 
Ls_0 = 1.8e-6; %1.94 * 1e-6; 

sigma_coll_lv  =  ((Ls_lv - Ls_0)/1e-6)^PExpcollagen; 
sigma_coll_sep =  ((Ls_sep - Ls_0)/1e-6)^PExpcollagen; 
sigma_coll_rv  =  ((Ls_rv - Ls_0)/1e-6)^PExpcollagen; 


sigma_pas_lv  = (Ls_lv  - Ls_0)/1e-6 + PConcollagen * sigma_coll_lv; 
sigma_pas_sep = (Ls_sep - Ls_0)/1e-6 + PConcollagen * sigma_coll_sep; 
sigma_pas_rv  = (Ls_rv  - Ls_0)/1e-6 + PConcollagen * sigma_coll_rv; 


% Total stress (kPa)
sigma_lv  = k_act * sigma_act_lv  + k_pas * sigma_pas_lv; 
sigma_sep = k_act * sigma_act_sep + k_pas * sigma_pas_sep; 
sigma_rv  = k_act * sigma_act_rv  + k_pas * sigma_pas_rv; 

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

f1 = -V_lv - 0.5 * Vw_lv - 0.5 * Vw_sep + Vm_sep - Vm_lv; 
f2 = Tx_lv + Tx_sep + Tx_rv;
f3 = V_rv + 0.5 * Vw_rv + 0.5 * Vw_sep + Vm_sep - Vm_rv;
f4 = Ty_lv + Ty_sep + Ty_rv; 

f = [f1; f2; f3; f4]; 
end 

