function init = initialconditions2(pars,HFpEFdata)



%% Parameters

C_Ao = pars(7); 
C_SA = pars(8); 
C_SV = pars(9); 
C_PA = pars(10); 
C_PV = pars(11); 

%% Data 

P_Ao = HFpEFdata.P_SAd; 
P_SA = HFpEFdata.P_SAd; 
P_SV = 5;%HFpEFdata.P_SA_dias_data; 

P_PA = HFpEFdata.P_PAd; 
P_PV = HFpEFdata.PCWP; 


%% Initial conditions 

xm_LV  = -5.0;
xm_SEP = +2.5;
xm_RV  = +8.0;
ym     = +5.0;

SL_LV   = 2.2;
SL_SEP  = 2.2;
SL_RV   = 2.2;

%vfactor = 2 * 1.0; % control 20-yo

V_LA = 35; %100; % mL
V_LV = 150;  % initial V_LV (mL)
V_Ao = C_Ao * P_Ao; %vfactor*60;
V_SA = C_SA * P_SA; %vfactor*200;
V_SV = C_SV * P_SV; %vfactor*1175;
V_RA = 35; %100; % mL
V_RV = 150; % initial V_LV (mL)
V_PA = C_PA * P_PA; %vfactor*35;
V_PV = C_PV * P_PV; %vfactor*65; 

% V_tot = V_LA + V_LV + V_Ao + V_SA + V_RA + V_RV + V_PA + V_PV
% V_tot / 5e3

%% Solve for DAE consistency 

x0 = [xm_LV ,xm_SEP ,xm_RV ,ym];

opts = optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');
x = fsolve(@TrisegEquations,x0,opts,V_LV,V_RV,pars);

%% Outputs 

init = [x, SL_LV, SL_SEP, SL_RV, ... 
    V_LA, V_LV, V_Ao, V_SA ,V_SV ,V_RA, V_RV, V_PA, V_PV, ...
    0, 0, 0]';
