function [rout,J] = calculatecost(outputs,theta,HFpEFdata,data) 


%% Calculate residuals 

HR        = data.HR;
procedure = data.procedure; 

if theta == 0 

    if strcmp(procedure,'RHC') == 1 
        % RHC data 
        P_SA_syst_data = HFpEFdata.P_SAs; 
        P_SA_dias_data = HFpEFdata.P_SAd; 
        P_RA_syst_data = HFpEFdata.P_RAs; 
        P_RA_dias_data = HFpEFdata.P_RAd; 
        P_RV_syst_data = HFpEFdata.P_RVs; 
        P_RV_dias_data = HFpEFdata.P_RVd; 
        P_PA_syst_data = HFpEFdata.P_PAs;
        P_PA_dias_data = HFpEFdata.P_PAd; 
        P_PV_mean_data = HFpEFdata.PCWP; 

        CO_data = HFpEFdata.CO_RHC;  

        % Residuals 
        P_SA = outputs.pressures.P_SA; 
        P_RA = outputs.pressures.P_RA; 
        P_RV = outputs.pressures.P_RV; 
        P_PA = outputs.pressures.P_PA; 
        P_PV = outputs.pressures.P_PV; 

        P_SA_syst = max(P_SA); 
        P_SA_dias = min(P_SA); 
        P_RA_syst = max(P_RA); 
        P_RA_dias = min(P_RA); 
        P_RV_syst = max(P_RV); 
        P_RV_dias = min(P_RV); 
        P_PA_syst = max(P_PA); 
        P_PA_dias = min(P_PA); 
        P_PV_mean = mean(P_PV); 

        rout_1 = (P_SA_syst - P_SA_syst_data)/P_SA_syst_data; 
        rout_2 = (P_SA_dias - P_SA_dias_data)/P_SA_dias_data; 
        rout_3 = (P_RA_syst - P_RA_syst_data)/P_RA_syst_data; 
        rout_4 = (P_RA_dias - P_RA_dias_data)/P_RA_dias_data;
        rout_5 = (P_RV_syst - P_RV_syst_data)/P_RV_syst_data; 
        rout_6 = (P_RV_dias - P_RV_dias_data)/P_RV_dias_data; 
        rout_7 = (P_PA_syst - P_PA_syst_data)/P_PA_syst_data; 
        rout_8 = (P_PA_dias - P_PA_dias_data)/P_PA_dias_data;
        rout_9 = (P_PV_mean - P_PV_mean_data)/P_PV_mean_data; 

        V_LV = outputs.volumes.V_LV; 
        
        V_LV_dias = max(V_LV); 
        V_LV_syst = min(V_LV); 

        SV = V_LV_dias - V_LV_syst; 
        CO = SV * HR * 1e-3; 
        rout_10 = (CO - CO_data)/CO_data; 

        % Output
        rout = [rout_1; rout_2; rout_3; rout_4; rout_5; 
            rout_6; rout_7; rout_8; rout_9; rout_10; 
            ]; 
    else
        % MRI data 
        V_LV_dias_data = HFpEFdata.LVEDV; 
        V_LV_syst_data = HFpEFdata.LVESV; 
        V_RV_dias_data = HFpEFdata.RVEDV; 
        V_RV_syst_data = HFpEFdata.RVESV; 
        
        CO_LV_data = HFpEFdata.LVCO_MRI; 
        CO_RV_data = HFpEFdata.RVCO_MRI; 
        
        % Residuals 
        V_LV = outputs.volumes.V_LV;
        V_RV = outputs.volumes.V_RV; 
        
        V_LV_dias = max(V_LV); 
        V_LV_syst = min(V_LV); 
        V_RV_dias = max(V_RV); 
        V_RV_syst = min(V_RV);
        
        rout_1 = (V_LV_dias - V_LV_dias_data)/V_LV_dias_data; 
        rout_2 = (V_LV_syst - V_LV_syst_data)/V_LV_syst_data; 
        rout_3 = (V_RV_dias - V_LV_dias_data)/V_RV_dias_data; 
        rout_4 = (V_RV_syst - V_LV_syst_data)/V_RV_syst_data; 
        
        SV_LV = V_LV_dias - V_LV_syst; 
        CO_LV = SV_LV * HR * 1e-3; 
        rout_5 = (CO_LV - CO_LV_data)/CO_LV_data; 
        
        SV_RV = V_RV_dias - V_RV_syst; 
        CO_RV = SV_RV * HR * 1e-3; 
        rout_6 = (CO_RV - CO_RV_data)/CO_RV_data; 
        
        % Output
        rout = [rout_1; rout_2; rout_3; rout_4; rout_5; rout_6]; 
    end 
        
    
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
    % 100% exercise 
    P_SA_syst_data = 173; %120; 
    P_SA_dias_data = 113; %80; 
    P_PA_dias_data = 10; %8.8; 
    P_PA_syst_data = 28; %20.5; 
    P_PV_mean_data = 13; %8; 

    CO_data        = 4; %5.76; 
    V_LV_dias_data = 118.5; %150; 
    V_LV_syst_data = 45.5; %60; 
    V_RV_dias_data = 108; %150; 
    V_RV_syst_data = 34; %60; 
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
    rout_10 = 0; %(V_LA_dias - V_LA_dias_data)/V_LA_dias_data; 
    rout_11 = 0; %(V_LA_syst - V_LA_syst_data)/V_LA_syst_data; 

    SV = V_LV_dias - V_LV_syst; 
    CO = SV * HR * 1e-3; 
    rout_12 = (CO - CO_data)/CO_data; 

    Q_m = outputs.flows.Q_m; 
    a = findpeaks(Q_m); 
    E_A_ratio = a(1)/a(2); 
    rout_13 = 0; %(E_A_ratio - E_A_ratio_data)/E_A_ratio_data; 
    
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
end 

J = rout'*rout; 

