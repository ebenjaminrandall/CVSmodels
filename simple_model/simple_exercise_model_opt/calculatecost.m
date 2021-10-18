function [rout,J] = calculatecost(outputs,theta,data) 

HR       = data.HR;

if theta == 0 

    % Baseline target values 
    P_SA_syst_data = 120; 
    P_SA_dias_data = 80; 
    P_PA_dias_data = 8.8; 
    P_PA_syst_data = 20.5; 
    P_PV_mean_data = 8; 

    CO_data        = 5.76; 
    V_LV_dias_data = 150; 
    V_LV_syst_data = 60; 
    V_RV_dias_data = 150; 
    V_RV_syst_data = 60; 
    V_LA_dias_data = 41; 
    V_LA_syst_data = 12; 

    E_A_ratio_data = 1.7; 
%     ejection_time_data = 0.292; 
    
    % Residuals 
    P_SA = outputs.pressures.P_SA; 
    P_PA = outputs.pressures.P_PA; 
    P_PV = outputs.pressures.P_PV; 

    P_SA_syst = max(P_SA); 
    P_SA_dias = min(P_SA); 
    P_PA_syst = max(P_PA); 
    P_PA_dias = min(P_PA); 
    P_PV_mean = mean(P_PV); 

    rout_1 = (P_SA_syst - P_SA_syst_data)/P_SA_syst_data; 
    rout_2 = (P_SA_dias - P_SA_dias_data)/P_SA_dias_data; 
    rout_3 = (P_PA_syst - P_PA_syst_data)/P_PA_syst_data; 
    rout_4 = (P_PA_dias - P_PA_dias_data)/P_PA_dias_data;
    rout_5 = (P_PV_mean - P_PV_mean_data)/P_PV_mean_data; 

    V_LV = outputs.volumes.V_LV; 
    V_RV = outputs.volumes.V_RV; 
    V_LA = outputs.volumes.V_LA; 

    V_LV_dias = max(V_LV); 
    V_LV_syst = min(V_LV); 
    V_RV_dias = max(V_RV); 
    V_RV_syst = min(V_RV); 
    V_LA_dias = max(V_LA); 
    V_LA_syst = min(V_LA); 

    rout_6  = (V_LV_dias - V_LV_dias_data)/V_LV_dias_data; 
    rout_7  = (V_LV_syst - V_LV_syst_data)/V_LV_syst_data; 
    rout_8  = (V_RV_dias - V_RV_dias_data)/V_RV_dias_data; 
    rout_9  = (V_RV_syst - V_RV_syst_data)/V_RV_syst_data; 
    rout_10 = (V_LA_dias - V_LA_dias_data)/V_LA_dias_data; 
    rout_11 = (V_LA_syst - V_LA_syst_data)/V_LA_syst_data; 

    SV = V_LV_dias - V_LV_syst; 
    CO = SV * HR * 1e-3; 
    rout_12 = (CO - CO_data)/CO_data; 

    Q_m = outputs.flows.Q_m; 
    a = findpeaks(Q_m); 
    E_A_ratio = a(1)/a(2); 
    rout_13 = (E_A_ratio - E_A_ratio_data)/E_A_ratio_data; 
    
%     % Average ejection time across multiple heart periods 
%     Q_a = outputs.flows.Q_a; 
%     b = find(Q_a > 0); 
%     c = mod(t(b),stim_per); 
%     d = find(diff(c) < 0); 
%     avlv_open = t([b(1) b(d + 1)]); 
%     avlv_closed = t([b(d) b(end)]); 
%     ejections = avlv_closed - avlv_open; 
%     ejection_time = mean(ejections); 
% 
%     rout_12 = (ejection_time - ejection_time_data)/ejection_time_data; 

    %Output
    rout = [rout_1; rout_2; rout_3; rout_4; rout_5; 
        rout_6; rout_7; rout_8; rout_9; rout_10; rout_11; 
        rout_12; rout_13;
        ]; 
    
elseif theta == 0.7
    % 70% exercise target values 
    P_SA_syst_data = 194; 
    P_SA_dias_data = 81;  
    P_PA_syst_data = 34; 
    P_PV_mean_data = 19.5; 
    
    V_LV_dias_data = 200; 
    V_RV_dias_data = 200; 

    CO_data = 16.2;
    EF_data = 0.78; 
%     ejection_time_data = 196; 
    
    % Residuals 
    P_SA = outputs.pressures.P_SA; 
    P_PA = outputs.pressures.P_PA; 
    P_PV = outputs.pressures.P_PV; 

    P_SA_syst = max(P_SA); 
    P_SA_dias = min(P_SA); 
    P_PA_syst = max(P_PA); 
    P_PV_mean = mean(P_PV); 

    rout_1 = (P_SA_syst - P_SA_syst_data)/P_SA_syst_data; 
    rout_2 = (P_SA_dias - P_SA_dias_data)/P_SA_dias_data; 
    rout_3 = (P_PA_syst - P_PA_syst_data)/P_PA_syst_data; 
    rout_4 = (P_PV_mean - P_PV_mean_data)/P_PV_mean_data; 

    V_LV = outputs.volumes.V_LV; 
    V_RV = outputs.volumes.V_RV; 
    
    V_LV_dias = max(V_LV);
    V_RV_dias = max(V_RV);
    
%     if V_LV_dias > V_LV_dias_data || V_RV_dias > V_RV_dias_data 
%         rout_5 = (V_LV_dias - V_LV_dias_data)/V_LV_dias_data; 
%         rout_6 = (V_RV_dias - V_RV_dias_data)/V_RV_dias_data; 
%     else 
        rout_5 = 0; 
        rout_6 = 0; 
%     end 
    
    V_LV_dias = max(V_LV); 
    V_LV_syst = min(V_LV); 
    SV = V_LV_dias - V_LV_syst; 
    CO = SV * HR * 1e-3; 
    EF = SV / V_LV_dias; 
    
    rout_7 = (CO - CO_data)/CO_data; 
    rout_8 = (EF - EF_data)/EF_data; 
%     
%     % Average ejection time across multiple heart periods 
%     Q_a = outputs.flows.Q_a; 
%     b = find(Q_a > 0); 
%     c = mod(t(b),stim_per); 
%     d = find(diff(c) < 0); 
%     avlv_open = t([b(1) b(d + 1)]); 
%     avlv_closed = t([b(d) b(end)]); 
%     ejections = avlv_closed - avlv_open; 
%     ejection_time = mean(ejections); 

%     rout_7 = (ejection_time - ejection_time_data)/ejection_time_data; 
    
    %Output 
    rout = [rout_1; rout_2; rout_3; rout_4; 
        rout_5; rout_6; rout_7; rout_8]; 
    
else 
    % 100% exercise target values 
    CO_data = 20.6; 
    
    % Residuals
    V_LV = outputs.volumes.V_LV; 
    V_LV_dias = max(V_LV); 
    V_LV_syst = min(V_LV); 
    SV = V_LV_dias - V_LV_syst; 
    CO = SV * HR * 1e-3; 
    
    rout = (CO - CO_data)/CO_data; 
end 

J = rout'*rout; 

