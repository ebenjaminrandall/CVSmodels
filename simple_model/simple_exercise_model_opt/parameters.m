function [pars,lb,ub] = parameters(exercise,theta) 


%% Exercise factors 

a = exercise(1); 
b = exercise(2); 
c = exercise(3); 
d = exercise(4); 
e = exercise(5); 
f = exercise(6); 

%% Parameters 

Vw_LV  = 80;
Vw_SEP = 38; 
Vw_RV  = 28; % Heart wall volumes (mL)

Amref_LV  = 0.975*80; % LV midwall reference surface area, cm^2
Amref_SEP = 0.975*45; % SEP midwall reference surface area, cm^2
Amref_RV  = 1.12*100; % RV midwall reference surface area, cm^2

% Lumped circulatory parameters
C_Ao = 0.65;  % Proximal aortic compliance, mL/mmHg
C_SA = 1.65/(1 + d*theta)/1.2; % Systemic arterial compliance, mL/mmHg
C_SV = 1.4*250/(1 + c*theta); % Systemic venous compliance, mL/mmHg 
C_PA = 5.4; % Pulmonary arterial compliance, mL/mmHg
C_PV = 25; % Pulmonary venous compliance, mL/mmHg

R_Ao  = 0.01; % resistance of aorta , mmHg*sec/mL
R_SA  = 1.1*0.965/(1 + b*theta);% mmHg*sec/mL; % Systemic vasculature resistance, mmHg*sec/mL
R_PA  = 0.05/(1 + e*theta); % Pulmonary vasculature resistance, mmHg*sec/mL 
R_vlv = 0.002; %  valve resistance, mmHg*sec/mL
R_tAo = 0.0020;
R_tSA = 0.05;
R_RA  = 0.024;
R_LA  = 0.024;

% Atrial elastances 
Emin = 1.5*0.050 ;
Emax = 0.150 ;

% Triseg parameters
Lsref  = 1.9; % Resting SL, micron
vmax   = 7; % micron/sec
LSEiso = 0.04; % micron
k_pas  = 22; 
k_act  = 7.5*96*(1 + a*theta); % mmHg 
SLrest = 1.51; % microns
Crest  = 0.02/(1 + f*theta) ;

%% Outputs 

pars = [Vw_LV; Vw_SEP; Vw_RV;                           %1-3
    Amref_LV; Amref_SEP; Amref_RV;                      %4-6
    C_Ao; C_SA; C_SV; C_PA; C_PV;                       %7-11
    R_vlv; R_LA; R_Ao; R_tAo; R_SA; R_tSA; R_RA; R_PA;  %12-19
    Emin; Emax;                                         %20-21
    Lsref; LSEiso; SLrest; vmax; k_pas; k_act; Crest;   %22-28
    ];

pars = log(pars); 
lb   = pars - log(10); 
ub   = pars + log(10); 

end 