function init = initialconditions(pars,data) 

Pbar  = data.Pbar;
DPbar = data.DPbar; 

%% Fixed parameters 

fixpars = data.fixpars; 

k_lvm = fixpars(1); 
k_sa  = fixpars(2); 
k_sv  = fixpars(3); 
k_rvm = fixpars(4); 
k_pa  = fixpars(5); 
k_pv  = fixpars(6); 

%% Adjustable parameters

E_lvm = pars(1);
E_rvm = pars(2); 

C_sa = pars(3); 
C_sv = pars(4); 
C_pa = pars(5); 
C_pv = pars(6); 

Ca_rest = pars(27); 

V_sau = pars(28);
V_svu = pars(29); 
V_pau = pars(30); 
V_pvu = pars(31); 

%% Pressures 

P_sa = k_sa * Pbar; 
P_sv = k_sv * Pbar; 
P_pa = k_pa * Pbar; 
P_pv = k_pv * Pbar; 

P_lvm = k_lvm * DPbar; 
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

% 8 - 14 (m^3)
V_lv0 = P_lvm / E_lvm;
V_sa0 = C_sa * P_sa + V_sau;  
V_sv0 = C_sv * P_sv + V_svu;  
V_rv0 = P_rvm / E_rvm; 
V_pa0 = C_pa * P_pa + V_pau;  
V_pv0 = C_pv * P_pv + V_pvu;  

% 15 - 17
Ca_lv0  = Ca_rest;
Ca_sep0 = Ca_rest;
Ca_rv0  = Ca_rest;

init = [xm_lv0; xm_sep0; xm_rv0; ym0;
    Lsc_lv0; Lsc_sep0; Lsc_rv0; 
    V_lv0; V_sa0; V_sv0; V_rv0; V_pa0; V_pv0; 
    Ca_lv0; Ca_sep0; Ca_rv0; 
    ]; 

x0   = init(1:4);  
opts = optimoptions('fsolve','Display','iter'); 
xopt = fsolve(@(x) triseg(x,pars,data,init),x0,opts); 
init(1:4) = xopt(1:4); 

end 