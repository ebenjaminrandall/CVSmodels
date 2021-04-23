function [dxdt, outputs] = model(t,x,pars,data) 

HR = data.HR; 

%% Parameters

% Compliance (mL mmHg^(-1))
C_sa = pars(3); 
C_sv = pars(4); 
C_pa = pars(5); 
C_pv = pars(6); 

% Resistance (mmHg s mL^(-1))
R_sa = pars(7); 
R_sv = pars(8); 
R_pa = pars(9); 
R_pv = pars(10); 
R_v  = pars(11); 

% Wall volume of ventricular wall segment (m^3)
Vw_lv  = pars(12); 
Vw_sep = pars(13); 
Vw_rv  = pars(14); 

% Reference midwall surface area (m^2)
Amref_lv  = pars(15); 
Amref_sep = pars(16); 
Amref_rv  = pars(17); 

% Time-scale (s)
tauR  = pars(18); 
tauD  = pars(19); 
tausc = pars(20); 

% Sarcomere length parameters (m)
Lsref   = pars(21);
Lsc0    = pars(22); 
Lse_iso = pars(23); 

% Force scaling factors (kPa) 
k_act = pars(24); 
k_pas = pars(25); 

v_max   = pars(26); % m s^(-1) sarcomere length shortening velocity
Ca_rest = pars(27); % dimensionless Diastolic resting level of activation

V_sau = pars(28);
V_svu = pars(29); 
V_pau = pars(30); 
V_pvu = pars(31); 

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
V_lv = x(8); 
V_sa = x(9); 
V_sv = x(10); 
V_rv = x(11);
V_pa = x(12); 
V_pv = x(13); 

% Mechanical activation 
Ca_lv  = x(14); 
Ca_sep = x(15);
Ca_rv  = x(16); 

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
sigma_act_lv  = Ca_lv  * (Lsc_lv  - Lsc0) / 1e-6 * (Ls_lv  - Lsc_lv)  / Lse_iso; 
sigma_act_sep = Ca_sep * (Lsc_sep - Lsc0) / 1e-6 * (Ls_sep - Lsc_sep) / Lse_iso;
sigma_act_rv  = Ca_rv  * (Lsc_rv  - Lsc0) / 1e-6 * (Ls_rv  - Lsc_rv)  / Lse_iso;

% Passive stress (kPa)
sigma_pas_lv  = 36 * max(0,eps_lv - 0.1)^2  + 0.1 * (eps_lv  - 0.1) + 0.0025 * exp(30 * eps_lv); 
sigma_pas_sep = 36 * max(0,eps_sep - 0.1)^2 + 0.1 * (eps_sep - 0.1) + 0.0025 * exp(30 * eps_sep); 
sigma_pas_rv  = 36 * max(0,eps_rv - 0.1)^2  + 0.1 * (eps_rv  - 0.1) + 0.0025 * exp(30 * eps_rv); 

% Total stress (kPa)
sigma_lv  = k_act * sigma_act_lv  + k_pas * sigma_pas_lv; 
sigma_sep = k_act * sigma_act_sep + k_pas * sigma_pas_sep; 
sigma_rv  = k_act * sigma_act_rv  + k_pas * sigma_pas_rv; 

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
ptrans1 = 2 * Tx_lv / ym; 
ptrans3 = 2 * Tx_rv / ym; 
P_lv = -ptrans1; 
P_rv = ptrans3; 

%% Calcium handling 

% Increase of activation with sarcomere length (dimensionless)
CaL_lv  = tanh(4 * ((Lsc_lv  - Lsc0) / 1e-6)^2); 
CaL_sep = tanh(4 * ((Lsc_sep - Lsc0) / 1e-6)^2); 
CaL_rv  = tanh(4 * ((Lsc_rv  - Lsc0) / 1e-6)^2); 

% Rise of mechanical activation 
T     = 1 / HR; 
tc    = mod(t, T); 
x     = min(8, max(0, tc / tauR)); 
Frise = 0.02 * x^3 * (8 - x)^2 * exp(-x);

% decrease of activation duration with decrease of sarcomere length (s)
T_lv  = tausc * (0.29 + 0.3 * (Lsc_lv  / 1e-6)); 
T_sep = tausc * (0.29 + 0.3 * (Lsc_sep / 1e-6));
T_rv  = tausc * (0.29 + 0.3 * (Lsc_rv  / 1e-6));

%% Lumped circulatory model 

% Pressure (kPa)
P_sa = (V_sa - V_sau) / C_sa; 
P_sv = (V_sv - V_svu) / C_sv; 
P_pa = (V_pa - V_pau) / C_pa; 
P_pv = (V_pv - V_pvu) / C_pv; 

% Flow (m^3 s^1) 
Q_m_valve = max((P_pv - P_lv) / R_pv, 0); 
Q_a_valve = max((P_lv - P_sa) / R_v, 0); 
Q_sa      = (P_sa - P_sv) / R_sa; 
Q_t_valve = max((P_sv - P_rv) / R_sv, 0); 
Q_p_valve = max((P_rv - P_pa) / R_v, 0); 
Q_pa      = (P_pa - P_pv) / R_pa; 

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
dV_lv = Q_m_valve - Q_a_valve; 
dV_sa = Q_a_valve - Q_sa; 
dV_sv = Q_sa      - Q_t_valve; 
dV_rv = Q_t_valve - Q_p_valve; 
dV_pa = Q_p_valve - Q_pa; 
dV_pv = Q_pa      - Q_m_valve; 

% 15 - 17
dCa_lv  = 1/tauR * CaL_lv  * Frise + 1/tauD * (Ca_rest - Ca_lv)  / (1 + exp((T_lv  - tc) / tauD)); 
dCa_sep = 1/tauR * CaL_sep * Frise + 1/tauD * (Ca_rest - Ca_sep) / (1 + exp((T_sep - tc) / tauD)); 
dCa_rv  = 1/tauR * CaL_rv  * Frise + 1/tauD * (Ca_rest - Ca_rv)  / (1 + exp((T_rv  - tc) / tauD)); 

dxdt = [dxm_lv; dxm_sep; dxm_rv; dym;
    dLsc_lv; dLsc_sep; dLsc_rv; 
    dV_lv; dV_sa; dV_sv; dV_rv; dV_pa; dV_pv; 
    dCa_lv; dCa_sep; dCa_rv; 
    ]; 

outputs = [P_lv; P_sa; P_sv; P_rv; P_pa; P_pv; 
    Vm_lv; Vm_sep; Vm_rv; 
    Am_lv; Am_sep; Am_rv; 
    Cm_lv; Cm_sep; Cm_rv; 
    eps_lv; eps_sep; eps_rv; 
    sigma_pas_lv; sigma_pas_sep; sigma_pas_rv; 
    sigma_act_lv; sigma_act_sep; sigma_act_rv; 
    sigma_lv; sigma_sep; sigma_rv; 
    Q_m_valve; Q_a_valve; Q_t_valve; Q_p_valve; 
    Q_sa; Q_pa;  
    ];

end 