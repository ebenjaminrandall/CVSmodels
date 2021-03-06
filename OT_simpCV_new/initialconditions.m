function init = initialconditions(pars,data) 

%% Fixed parameters 

fixpars = data.fixpars; 

P_lam = fixpars(1); 
P_lvm = fixpars(2); 
P_sam = fixpars(3); 
P_svm = fixpars(4);
P_ram = fixpars(5); 
P_rvm = fixpars(6); 
P_pam = fixpars(7); 
P_pvm = fixpars(8); 

%% Adjustable parameters

E_lam = pars(2); 
E_lvm = pars(3);

E_ram = pars(5);
E_rvm = pars(6); 

C_sa = pars(7); 
C_sv = pars(8); 
C_pa = pars(9); 
C_pv = pars(10); 

%% Initial conditions 

% 1 - 4 (convert cm to m)
xm_lv0  = data.deformation.xm_lv0;
xm_sep0 = data.deformation.xm_sep0;
xm_rv0  = data.deformation.xm_rv0;
ym0     = data.deformation.ym0; 

% 5 - 7 (convert um to m)
Lsc_lv0  = 1.6 * 1e-6;
Lsc_sep0 = 1.6 * 1e-6;
Lsc_rv0  = 1.6 * 1e-6;

% 8 - 15 (m^3)
V_la0 = P_lam / E_lam; 
V_lv0 = P_lvm / E_lvm;
V_sa0 = C_sa * P_sam;  
V_sv0 = C_sv * P_svm;  
V_ra0 = P_ram / E_ram; 
V_rv0 = P_rvm / E_rvm; 
V_pa0 = C_pa * P_pam;  
V_pv0 = C_pv * P_pvm; 

init = [xm_lv0; xm_sep0; xm_rv0; ym0;
    Lsc_lv0; Lsc_sep0; Lsc_rv0; 
    V_la0; V_lv0; V_sa0; V_sv0; V_ra0; V_rv0; V_pa0; V_pv0; 
    ]; 

x0   = log(init(1:7)); 
opts = optimoptions('fsolve','Display','iter',...
    'MaxFunctionEvaluations',1e3); 
xopt = fsolve(@(x) triseg(x,pars,data,init),x0,opts); 
init(1:7) = exp(xopt(1:7)); 

end 