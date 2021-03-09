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
E_ram = pars(4); 
E_lvm = pars(5); 
E_rvm = pars(6); 

% Compliance (m^3 kPa^(-1))
C_sa = pars(7); 
C_sv = pars(8); 
C_pa = pars(9); 
C_pv = pars(10); 

%% Pressures (kPa)

P_la = k_lam * DPbar; 
P_lv = k_lvm * DPbar;
P_sa = k_sam * DPbar; 
P_sv = k_svm * DPbar; 
P_ra = k_ram * DPbar;
P_rv = k_rvm * DPbar; 
P_pa = k_pam * DPbar; 
P_pv = k_pvm * DPbar; 

P_lvm = k_lvm * DPbar; 
P_rvm = k_rvm * DPbar; 

%% Initial conditions 

% 1 - 4 Displacements (convert cm to m)
xm_lv0  = 4.5 * 1e-2; 
xm_sep0 = 2   * 1e-2; 
xm_rv0  = 5   * 1e-2; 
ym0     = 3.3 * 1e-2; 

% 5 - 7 Sarcomere lengths (convert um to m)
Lsc_lv0  = 2 * 1e-6;
Lsc_sep0 = 2 * 1e-6;
Lsc_rv0  = 2 * 1e-6;

% 8 - 14 Volumes (m^3)
V_la0 = P_la / E_lam; 
V_lv0 = P_lv / E_lvm;
V_sa0 = C_sa * P_sa;  
V_sv0 = C_sv * P_sv;  
V_ra0 = P_ra / E_ram; 
V_rv0 = P_rv / E_rvm; 
V_pa0 = C_pa * P_pa;  
V_pv0 = C_pv * P_pv;  

init = [xm_lv0; xm_sep0; xm_rv0; ym0;
    Lsc_lv0; Lsc_sep0; Lsc_rv0; 
    V_la0; V_lv0; V_sa0; V_sv0; V_ra0; V_rv0; V_pa0; V_pv0;  
    ]; 

% Solve the Triseg model DAE for consistent ICS 
x0   = init(1:7);  
opts = optimoptions('fsolve','Display','off'); 
xopt = fsolve(@(x) triseg(x,pars,data,init,P_lvm, P_rvm),x0,opts); 

init(1:7) = xopt(1:7); 

end 