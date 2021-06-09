function init = initialconditions(pars,data) 

DPbar = data.DPbar; 

%% Fixed parameters 

fixpars = data.fixpars; 

k_lam = fixpars(1); 
k_lvm = fixpars(2); 
k_sam = fixpars(3); 
k_svm = fixpars(4);
k_ram = fixpars(5); 
k_rvm = fixpars(6); 
k_pam = fixpars(7); 
k_pvm = fixpars(8); 

%% Adjustable parameters

% Elastance (kPa m^(-3))
E_lam = pars(2); 
E_lvm = pars(4); 
E_ram = pars(6); 
E_rvm = pars(8); 

% Compliance (m^3 kPa^(-1))
C_sa = pars(9); 
C_sv = pars(10); 
C_pa = pars(11); 
C_pv = pars(12); 

%% Pressures (kPa)

P_la = k_lam * DPbar; 
P_lv = k_lvm * DPbar;
P_sa = k_sam * DPbar; 
P_sv = k_svm * DPbar; 
P_ra = k_ram * DPbar;
P_rv = k_rvm * DPbar; 
P_pa = k_pam * DPbar; 
P_pv = k_pvm * DPbar; 

%% Initial conditions 

% 1 - 8 Volumes (m^3)
V_la0 = P_la / E_lam; 
V_lv0 = P_lv / E_lvm;
V_sa0 = C_sa * P_sa;  
V_sv0 = C_sv * P_sv;  
V_ra0 = P_ra / E_ram; 
V_rv0 = P_rv / E_rvm; 
V_pa0 = C_pa * P_pa;  
V_pv0 = C_pv * P_pv;  
V_sep0 = ; 

init = [V_la0; V_lv0; V_sa0; V_sv0; V_ra0; V_rv0; V_pa0; V_pv0; V_sep0];

x0 = [V_lv0; V_rv0; V_sep0]; 
opts = optimoptions('fsolve','Display','iter'); 
xopt = fsolve(@(x) smith(x,pars,data,init),x0,opts); 


init(2) = x0(1); 
init(6) = x0(2); 
init(9) = x0(3); 


end 