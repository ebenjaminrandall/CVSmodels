function f = smith(x,pars,data)

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

E_es_lv = ; 



%% Variables 

V_lv  = x(2); 
V_rv  = x(6);
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

P_lv  = e * P_es_lv  + (1 - e) * P_ed_lv; 
P_rv  = e * P_es_rv  + (1 - e) * P_ed_rv; 
P_sep = e * P_es_sep + (1 - e) * P_ed_sep; 

%% Outputs

f = P_sep - P_lv + P_rv; 

end 