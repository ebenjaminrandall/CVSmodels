function f = TrisegEquations(x,Vw_LV,Vw_SEP,Vw_RV, V_LV, V_RV)

xm_LV  = x(1); % LV heart geometry variable, cm
xm_SEP = x(2); % septum heart geometry variable, cm
xm_RV  = x(3); % RV heart geometry variable, cm
ym     = x(4); % Heart geometry variable, cm

%% ventricular mechanics
Vm_LV  = (pi/6)*xm_LV*(xm_LV^2 + 3*ym^2);
Vm_SEP = (pi/6)*xm_SEP*(xm_SEP^2 + 3*ym^2);
Vm_RV  = (pi/6)*xm_RV*(xm_RV^2 + 3*ym^2);
Am_LV  = pi*(xm_LV^2 + ym^2);
Am_SEP = pi*(xm_SEP^2 + ym^2);
Am_RV  = pi*(xm_RV^2 + ym^2);
Cm_LV  = 2*xm_LV/(xm_LV^2 + ym^2);
Cm_SEP = 2*xm_SEP/(xm_SEP^2 + ym^2);
Cm_RV  = 2*xm_RV/(xm_RV^2 + ym^2);
z_LV   = 3*Cm_LV*Vw_LV/(2*Am_LV);
z_SEP  = 3*Cm_SEP*Vw_SEP/(2*Am_SEP);
z_RV   = 3*Cm_RV*Vw_RV/(2*Am_RV);

sigmaf_LV = 25;
sigmaf_SEP = 25;
sigmaf_RV = 25;

% equilibrium of forces at junction circle
Tm_LV = (Vw_LV*sigmaf_LV/(2*Am_LV))*(1 + (z_LV^2)/3 + (z_LV^4)/5);
Tm_SEP = (Vw_SEP*sigmaf_SEP/(2*Am_SEP))*(1 + (z_SEP^2)/3 + (z_SEP^4)/5);
Tm_RV = (Vw_RV*sigmaf_RV/(2*Am_RV))*(1 + (z_RV^2)/3 + (z_RV^4)/5);

Tx_LV = Tm_LV*2*xm_LV*ym/(xm_LV^2 + ym^2);
Tx_SEP = Tm_SEP*2*xm_SEP*ym/(xm_SEP^2 + ym^2);
Tx_RV = Tm_RV*2*xm_RV*ym/(xm_RV^2 + ym^2);

Ty_LV = Tm_LV*(-xm_LV^2 + ym^2)/(xm_LV^2 + ym^2);
Ty_SEP = Tm_SEP*(-xm_SEP^2 + ym^2)/(xm_SEP^2 + ym^2);
Ty_RV = Tm_RV*(-xm_RV^2 + ym^2)/(xm_RV^2 + ym^2);

f(1) = (-V_LV - 0.5*Vw_LV - 0.5*Vw_SEP + Vm_SEP - Vm_LV)/V_LV; %xm_LV
f(2) = (Tx_LV + Tx_SEP + Tx_RV); % xm_SEP
f(3) = (V_RV + 0.5*Vw_RV + 0.5*Vw_SEP + Vm_SEP - Vm_RV)/V_LV; %xm_RV
f(4) = (Ty_LV + Ty_SEP + Ty_RV);  %ym

