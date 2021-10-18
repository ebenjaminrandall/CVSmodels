function pars = parameters


Vw_LV  = 80; % LV wall volume, mL 
Vw_SEP = 38; % Septal wall volume, mL 
Vw_RV  = 28; % RV wall volume, mL 

Amref_LV  = 0.975*80; % LV midwall reference surface area, cm^2
Amref_SEP = 0.975*45; % SEP midwall reference surface area, cm^2
Amref_RV  = 1.12*100; % RV midwall reference surface area, cm^2


pars = [Vw_LV; Vw_SEP; Vw_RV;
    Amref_LV; Amref_SEP; Amref_RV;
    ];
