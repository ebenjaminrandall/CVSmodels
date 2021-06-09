function init = initialconditions(pars,data) 

Pbar  = data.Pbar;
DPbar = data.DPbar; 

%% Fixed parameters 

fixpars = data.fixpars; 

k_lam = fixpars(1); 
k_lvm = fixpars(2); 
k_sa  = fixpars(3); 
k_sv  = fixpars(4);
k_ram = fixpars(5); 
k_rvm = fixpars(6); 
k_pa  = fixpars(7); 
k_pv  = fixpars(8); 

%% Adjustable parameters

E_lam = pars(2); 
E_lvm = pars(3);

E_ram = pars(5);
E_rvm = pars(6); 

C_sa = pars(7); 
C_sv = pars(8); 
C_pa = pars(9); 
C_pv = pars(10); 

Ca_rest = pars(34); 

%% Pressures 

P_sa = k_sa * Pbar; 
P_sv = k_sv * Pbar; 
P_pa = k_pa * Pbar; 
P_pv = k_pv * Pbar; 

P_lam = k_lam * DPbar;
P_lvm = k_lvm * DPbar;
P_ram = k_ram * DPbar; 
P_rvm = k_rvm * DPbar; 

%% Initial conditions 

% 1 - 4 (convert cm to m)
xm_lv0  = data.deformation.xm_lv0;
xm_sep0 = data.deformation.xm_sep0;
xm_rv0  = data.deformation.xm_rv0;
ym0     = data.deformation.ym0; 

% 5 - 7 (convert um to m)
Lsc_lv0  = 2 * 1e-6;
Lsc_sep0 = 2 * 1e-6;
Lsc_rv0  = 2 * 1e-6;

% 8 - 15 (m^3)
V_la0 = P_lam / E_lam; 
V_lv0 = P_lvm / E_lvm;
V_sa0 = C_sa * P_sa;  
V_sv0 = C_sv * P_sv;  
V_ra0 = P_ram / E_ram; 
V_rv0 = P_rvm / E_rvm; 
V_pa0 = C_pa * P_pa;  
V_pv0 = C_pv * P_pv;   

% 16 - 18
Ca_lv0  = Ca_rest;
Ca_sep0 = Ca_rest;
Ca_rv0  = Ca_rest;

init = [xm_lv0; xm_sep0; xm_rv0; ym0;
    Lsc_lv0; Lsc_sep0; Lsc_rv0; 
    V_la0; V_lv0; V_sa0; V_sv0; V_ra0; V_rv0; V_pa0; V_pv0; 
    Ca_lv0; Ca_sep0; Ca_rv0; 
    ]; 

x0   = init(1:4);  
opts = optimoptions('fsolve','Display','iter'); 
xopt = fsolve(@(x) triseg(x,pars,data,init),x0,opts); 
init(1:4) = xopt(1:4); 

end 