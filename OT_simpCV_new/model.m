function [dxdt, outputs] = model(t,x,pars,data) 

HR = data.HR; 

%% Parameters

E_laM = pars(1);
E_lam = pars(2);

E_raM = pars(4);
E_ram = pars(5); 

% Compliance (mL mmHg^(-1))
C_sa = pars(7); 
C_sv = pars(8); 
C_pa = pars(9); 
C_pv = pars(10); 

% Resistance (mmHg s mL^(-1))
R_sa = pars(11); 
R_sv = pars(12); 
R_pa = pars(13); 
R_pv = pars(14); 

R_m_valve = pars(15); 
R_a_valve = pars(16); 
R_t_valve = pars(17); 
R_p_valve = pars(18); 

% Wall volume of ventricular wall segment (m^3)
Vw_lv  = pars(19); 
Vw_sep = pars(20); 
Vw_rv  = pars(21); 

% Reference midwall surface area (m^2)
Amref_lv  = pars(22); 
Amref_sep = pars(23); 
Amref_rv  = pars(24); 

% Sarcomere length parameters (m)
Lsref   = pars(25);
Lsc0    = pars(26); 
Lse_iso = pars(27); 

% Force scaling factors (kPa) 
k_act_lv = pars(28); 
k_act_rv = pars(29); 
k_pas    = pars(30); 

v_max   = pars(31); % m s^(-1) sarcomere length shortening velocity

tau_a = pars(36); 

k_TS = pars(37); 
k_TR = pars(38); 

alpha_c = pars(39); 
beta_c  = pars(40); 

%% Variables 

% Axial distance of midwall junction 
xm_lv  = x(1); 
xm_sep = x(2); 
xm_rv  = x(3);

% Radial distance of midwall junction 
ym = x(4); 

% Contractile element length 
Lsc_lv  = x(5); 
Lsc_sep = x(6); 
Lsc_rv  = x(7); 

% Volumes 
V_la = x(8); 
V_lv = x(9); 
V_sa = x(10); 
V_sv = x(11);
V_ra = x(12); 
V_rv = x(13);
V_pa = x(14); 
V_pv = x(15); 

%% Activation function

T = 60/HR; 
TS = k_TS * T; 
TR = k_TR * T; 

tc_v = mod(t,T);
if tc_v >= 0 && tc_v < TS 
    y_v = 0.5*(1 - cos(pi*tc_v/TS)); 
elseif tc_v >= TS && tc_v < TR + TS 
    y_v = 0.5*(1 + cos(pi*(tc_v - TS)/TR)); 
else
    y_v = 0; 
end 

tc_a = mod(t + tau_a,T);
if tc_a >= 0 && tc_a < TS 
    y_a = 0.5*(1 - cos(pi*tc_a/TS)); 
elseif tc_a >= TS && tc_a < TR + TS 
    y_a = 0.5*(1 + cos(pi*(tc_a - TS)/TR)); 
else
    y_a = 0; 
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
sigma_act_lv  = y_v  * (Lsc_lv  - Lsc0) / 1e-6 * (Ls_lv  - Lsc_lv)  / Lse_iso; 
sigma_act_sep = y_v  * (Lsc_sep - Lsc0) / 1e-6 * (Ls_sep - Lsc_sep) / Lse_iso;
sigma_act_rv  = y_v  * (Lsc_rv  - Lsc0) / 1e-6 * (Ls_rv  - Lsc_rv)  / Lse_iso;

% Stress due to collagen (kPa) 
sigma_coll_lv  =  alpha_c * ((Ls_lv - Lsc0)/1e-6)^beta_c; 
sigma_coll_sep =  alpha_c * ((Ls_sep - Lsc0)/1e-6)^beta_c; 
sigma_coll_rv  =  alpha_c * ((Ls_rv - Lsc0)/1e-6)^beta_c; 

% Passive stress (kPa)
sigma_pas_lv  = (Ls_lv  - Lsc0)/1e-6 + sigma_coll_lv; 
sigma_pas_sep = (Ls_sep - Lsc0)/1e-6 + sigma_coll_sep; 
sigma_pas_rv  = (Ls_rv  - Lsc0)/1e-6 + sigma_coll_rv; 

% Total stress (kPa)
sigma_lv  = k_act_lv * sigma_act_lv  + k_pas * sigma_pas_lv; 
sigma_sep = k_act_lv * sigma_act_sep + k_pas * sigma_pas_sep; 
sigma_rv  = k_act_rv * sigma_act_rv  + k_pas * sigma_pas_rv; 

% Representative midwall tension (kPa m)
Tm_lv  = (Vw_lv  * sigma_lv  / (2 * Am_lv))  * (1 + (z_lv^2)/3  + (z_lv^4)/5); 
Tm_sep = (Vw_sep * sigma_sep / (2 * Am_sep)) * (1 + (z_sep^2)/3 + (z_sep^4)/5); 
Tm_rv  = (Vw_rv  * sigma_rv  / (2 * Am_rv))  * (1 + (z_rv^2)/3  + (z_rv^4)/5);

% Axial midwall tension component (kPa m)
Tx_lv  = - Tm_lv  * 2 * xm_lv  * ym / (xm_lv^2  + ym^2); 
Tx_sep = Tm_sep * 2 * xm_sep * ym / (xm_sep^2 + ym^2); 
Tx_rv  = Tm_rv  * 2 * xm_rv  * ym / (xm_rv^2  + ym^2); 

% Radial midwall tension component (kPa m)
Ty_lv  = Tm_lv  * (-xm_lv^2  + ym^2) / (xm_lv^2  + ym^2); 
Ty_sep = Tm_sep * (-xm_sep^2 + ym^2) / (xm_sep^2 + ym^2); 
Ty_rv  = Tm_rv  * (-xm_rv^2  + ym^2) / (xm_rv^2  + ym^2);

% Ventricular pressure (kPa)
ptrans_lv = 2 * Tx_lv / ym; 
ptrans_rv = 2 * Tx_rv / ym; 
P_lv = -ptrans_lv; 
P_rv = ptrans_rv; 

%% Lumped circulatory model 

% Elastance
E_la = (E_laM - E_lam) * y_a + E_lam; 
E_ra = (E_raM - E_ram) * y_a + E_ram; 

% Pressure (kPa)
P_la = V_la * E_la; 
P_sa = V_sa / C_sa; 
P_sv = V_sv / C_sv; 
P_ra = V_ra * E_ra; 
P_pa = V_pa / C_pa; 
P_pv = V_pv / C_pv; 

% Flow (m^3 s^1) 
Q_m_valve = max((P_la - P_lv) / R_m_valve, 0); 
Q_a_valve = max((P_lv - P_sa) / R_a_valve, 0); 
Q_sa      = (P_sa - P_sv) / R_sa;
Q_sv      = (P_sv - P_ra) / R_sv; 
Q_t_valve = max((P_ra - P_rv) / R_t_valve, 0); 
Q_p_valve = max((P_rv - P_pa) / R_p_valve, 0); 
Q_pa      = (P_pa - P_pv) / R_pa; 
Q_pv      = (P_pv - P_la) / R_pv; 

%% ODEs

% 1 - 4
dxm_lv  = -V_lv - 0.5 * Vw_lv - 0.5 * Vw_sep + Vm_sep - Vm_lv; 
dxm_sep = Tx_lv + Tx_sep + Tx_rv;
dxm_rv  = V_rv + 0.5 * Vw_rv + 0.5 * Vw_sep + Vm_sep - Vm_rv;
dym     = Ty_lv + Ty_sep + Ty_rv; 

% 5 - 7
dLsc_lv  = ((Ls_lv  - Lsc_lv)  /Lse_iso - 1) * v_max;
dLsc_sep = ((Ls_sep - Lsc_sep) /Lse_iso - 1) * v_max;
dLsc_rv  = ((Ls_rv  - Lsc_rv)  /Lse_iso - 1) * v_max;

% 8 - 14
dV_la = Q_pv      - Q_m_valve; 
dV_lv = Q_m_valve - Q_a_valve; 
dV_sa = Q_a_valve - Q_sa; 
dV_sv = Q_sa      - Q_sv; 
dV_ra = Q_sv      - Q_t_valve;
dV_rv = Q_t_valve - Q_p_valve; 
dV_pa = Q_p_valve - Q_pa; 
dV_pv = Q_pa      - Q_pv; 

dxdt = [dxm_lv; dxm_sep; dxm_rv; dym;
    dLsc_lv; dLsc_sep; dLsc_rv; 
    dV_la; dV_lv; dV_sa; dV_sv; dV_ra; dV_rv; dV_pa; dV_pv; 
    ]; 

outputs = [P_la; P_lv; P_sa; P_sv; P_ra; P_rv; P_pa; P_pv; 
    Vm_lv; Vm_sep; Vm_rv; 
    Am_lv; Am_sep; Am_rv; 
    Cm_lv; Cm_sep; Cm_rv; 
    eps_lv; eps_sep; eps_rv; 
    sigma_pas_lv; sigma_pas_sep; sigma_pas_rv; 
    sigma_act_lv; sigma_act_sep; sigma_act_rv; 
    sigma_lv; sigma_sep; sigma_rv; 
    Q_m_valve; Q_a_valve; Q_t_valve; Q_p_valve; 
    Q_sa; Q_sv; Q_pa; Q_pv;   
    ];

end 