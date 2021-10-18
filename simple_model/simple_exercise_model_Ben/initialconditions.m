function ICs = initialconditions(pars, vfactor) 


%% Unpack parameters

Vw_LV  = pars(1); 
Vw_SEP = pars(2); 
Vw_RV  = pars(3); 

%% Heart submodel 

% Deflections (mu m)
xm_LV  = -5.0;
xm_SEP = +2.5;
xm_RV  = +8.0;
ym     = +5.0;

% Sarcomere lengths (mu m)
SL_LV   = 2.2;
SL_SEP  = 2.2;
SL_RV   = 2.2;

% Ventricular volumes
V_LV  = 150;  % initial V_LV (mL)
V_RV  = 150; % initial V_LV (mL)

x0 = [xm_LV ,xm_SEP ,xm_RV ,ym];

% opts = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
opts = optimset('MaxFunEvals',10000,'MaxIter',1000);
TrisegEquations(x0,Vw_LV,Vw_SEP,Vw_RV, V_LV, V_RV);
x = fsolve(@TrisegEquations,x0,opts,Vw_LV,Vw_SEP,Vw_RV, V_LV, V_RV);

%% Circulation submodel 

V_LA = 100; % mL
V_Ao = vfactor*60;
V_SA = vfactor*200;
V_SV = vfactor*1175;
V_RA = 100; % mL
V_PA = vfactor*35;
V_PV = vfactor*65; 

%% Outputs

ICs = [x'; SL_LV; SL_SEP; SL_RV; 
    V_LA; V_LV; V_Ao; V_SA; V_SV; V_RA; V_RV; V_PA; V_PV; 
    ]';
end 