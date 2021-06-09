function [dxdt, outputs] = model(t,x,pars,data,Y_a,Y_v) 

%% Parameters

% Elastance (kPa m^(-3))
E_laM = pars(1); 
E_lam = pars(2); 
E_lvM = pars(3); 
E_lvm = pars(4); 
E_raM = pars(5); 
E_ram = pars(6);
E_rvM = pars(7); 
E_rvm = pars(8); 

% Compliance (m^3 kPa^(-1))
C_sa = pars(9); 
C_sv = pars(10); 
C_pa = pars(11); 
C_pv = pars(12); 

% Resistance (kPa s m^(-3))
R_sa      = pars(13); 
R_sv      = pars(14); 
R_pa      = pars(15); 
R_pv      = pars(16); 
R_m_valve = pars(17); 
R_a_valve = pars(18); 
R_t_valve = pars(19); 
R_p_valve = pars(20);

%% Variables 

% Volumes (m^3)
V_la = x(1); 
V_lv = x(2); 
V_sa = x(3); 
V_sv = x(4);
V_ra = x(5); 
V_rv = x(6);
V_pa = x(7); 
V_pv = x(8); 
V_sep = x(9); 

%% Lumped circulatory model 

E_la = (E_laM - E_lam) * Y_a(t) + E_lam; 
E_ra = (E_raM - E_ram) * Y_a(t) + E_ram; 

E_lv = (E_lvM - E_lvm) * Y_v(t) + E_lvm; 
E_rv = (E_rvM - E_rvm) * Y_v(t) + E_lvm; 

% Left ventricular pressure 
e = exp(-80*(t - 0.27)^2); 

P_es_lv  = E_es_lv * V_lv; 
P_es_rv  = E_es_rv * V_rv; 
P_es_sep = E_es_sep * V_sep; 

P_ed_lv  = exp(lambda_lv  * V_lv)  - 1; 
P_ed_rv  = exp(lambda_rv  * V_rv)  - 1; 
P_ed_sep = exp(lambda_sep * V_sep) - 1; 

% Pressure (kPa)
P_la = V_la * E_la;
P_lv = V_lv * E_lv; 
P_sa = V_sa / C_sa; 
P_sv = V_sv / C_sv; 
P_ra = V_ra * E_ra; 
P_rv = V_rv * E_rv; 
P_pa = V_pa / C_pa; 
P_pv = V_pv / C_pv; 

% Flow (m^3 s^(-1)) 
if P_la > P_lv 
    Q_m_valve = (P_la - P_lv) / R_m_valve; 
else 
    Q_m_valve = 0; 
end 
if P_lv > P_sa 
    Q_a_valve = (P_lv - P_sa) / R_a_valve; 
else 
    Q_a_valve = 0; 
end 
Q_sa = (P_sa - P_sv) / R_sa; 
if P_sv > P_ra 
    Q_sv = max((P_sv - P_ra) / R_sv, 0); 
else
    Q_sv = 0; 
end 
if P_ra > P_rv 
    Q_t_valve = (P_ra - P_rv) / R_t_valve; 
else 
    Q_t_valve = 0; 
end 
if P_rv > P_pa 
    Q_p_valve = (P_rv - P_pa) / R_p_valve; 
else 
    Q_p_valve = 0; 
end 
Q_pa = (P_pa - P_pv) / R_pa; 
Q_pv = (P_pv - P_la) / R_pv; 

%% ODEs

% 1 - 8
dV_la = Q_pv      - Q_m_valve; 
dV_lv = Q_m_valve - Q_a_valve; 
dV_sa = Q_a_valve - Q_sa; 
dV_sv = Q_sa      - Q_sv; 
dV_ra = Q_sv      - Q_t_valve; 
dV_rv = Q_t_valve - Q_p_valve; 
dV_pa = Q_p_valve - Q_pa; 
dV_pv = Q_pa      - Q_pv; 
dV_sep = P_sep - P_lv + P_rv; 

dxdt = [dV_la; dV_lv; dV_sa; dV_sv; dV_ra; dV_rv; dV_pa; dV_pv; dV_sep]; 

outputs = [P_la; P_lv; P_sa; P_sv; P_ra; P_rv; P_pa; P_pv; 
    Q_m_valve; Q_a_valve; Q_sa; Q_sv; Q_t_valve; Q_p_valve; Q_pa; Q_pv; 
    E_lv; E_rv; 
    ];



end 