function init = initialconditions(pars)


xm_LV  = -5.0;
xm_SEP = +2.5;
xm_RV  = +8.0;
ym     = +5.0;

SL_LV   = 2.2;
SL_SEP  = 2.2;
SL_RV   = 2.2;

vfactor = 1.0; % control 20-yo

V_LA = 35; %100; % mL
V_LV = 150;  % initial V_LV (mL)
V_Ao = vfactor*60;
V_SA = vfactor*200;
V_SV = vfactor*1175;
V_RA = 35; %100; % mL
V_RV = 150; % initial V_LV (mL)
V_PA = vfactor*35;
V_PV = vfactor*65; 

%V_tot = V_LA + V_LV + V_Ao + V_SA + V_RA + V_RV + V_PA + V_PV

x0 = [xm_LV ,xm_SEP ,xm_RV ,ym];

opts = optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');
x = fsolve(@TrisegEquations,x0,opts,V_LV,V_RV,pars);

init = [x, SL_LV, SL_SEP, SL_RV, ... 
    V_LA, V_LV, V_Ao, V_SA ,V_SV ,V_RA, V_RV, V_PA, V_PV, ...
    0, 0, 0]';
