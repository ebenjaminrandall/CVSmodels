function init = initialconditions(pars,data)

vfactor = data.vfactor;
V0      = data.V0; 

C_rest = pars(end); 

%% Volumes 

V_LA0 = V0(1);% * vfactor;  
V_LV0 = V0(2);% * vfactor; 
V_Ao0 = V0(3) * vfactor; 
V_SA0 = V0(4) * vfactor; 
V_SV0 = V0(5) * vfactor; 
V_RA0 = V0(6);% * vfactor; 
V_RV0 = V0(7);% * vfactor;
V_PA0 = V0(8) * vfactor; 
V_PV0 = V0(9) * vfactor; 

xm_LV0  = -5.0;
xm_SEP0 = +2.5;
xm_RV0  = +8.0;
ym0     = +5.0;

SL_LV0   = 2.2;
SL_SEP0  = 2.2;
SL_RV0   = 2.2;

C_LV0  = C_rest; 
C_SEP0 = C_rest; 
C_RV0  = C_rest; 

init = [xm_LV0; xm_SEP0; xm_RV0; ym0; 
    SL_LV0; SL_SEP0; SL_RV0; 
    V_LA0; V_LV0; V_Ao0; V_SA0; V_SV0; V_RA0; V_RV0; V_PA0; V_PV0;
    C_LV0; C_SEP0; C_RV0; 
    ];

xopt0 = [xm_LV0 ,xm_SEP0 ,xm_RV0 ,ym0];

% opts = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
%opts = optimset('MaxFunEvals',100000,'MaxIter',10000);

opts = optimset('Display','none','MaxFunEvals',10000,'MaxIter',1000);
xopt = fsolve(@TrisegEquations, xopt0, opts, init, pars);



init(1) = xopt(1);
init(2) = xopt(2); 
init(3) = xopt(3); 
init(4) = xopt(4); 


