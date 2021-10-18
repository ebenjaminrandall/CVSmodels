function [J,rout,outputs] = model_wrap_baseline_HFpEF(p,exercise,data,HFpEFdata)

theta = data.theta;

if theta == 0 

    % RHC    
    procedure = 'RHC'; 

    HR = (HFpEFdata.HR_RHC_max - HFpEFdata.HR_RHC) * theta + HFpEFdata.HR_RHC; %bpm
    freq = HR/60; %Hz
    stim_per = 1/freq; 
    t = 0:.01:5*stim_per; 

    data.time      = t; 
    data.HR        = HR; 
    data.stim_per  = stim_per; 
    data.procedure = procedure; 

    pars_RHC = parameters2(exercise,p,theta,HFpEFdata,data); 
    [o_RHC,rout_RHC,J_RHC] = model_sol(pars_RHC,theta,HFpEFdata,data); 

    % MRI 
    procedure = 'MRI'; 

    theta = 0; 
    HR = HFpEFdata.HR_MRI; 
    freq = HR/60; %Hz
    stim_per = 1/freq; 
    t = 0:.01:5*stim_per; 

    data.time     = t; 
    data.HR       = HR; 
    data.stim_per = stim_per; 
    data.procedure = procedure;

    pars_RHC = parameters2(exercise,p,theta,HFpEFdata,data); 
    [o_MRI,rout_MRI,J_MRI] = model_sol(pars_RHC,theta,HFpEFdata,data); 

    %% Outputs 

    outputs.o_RHC = o_RHC; 
    outputs.o_MRI = o_MRI; 

    rout = [rout_RHC; rout_MRI]; 

    J = J_RHC+ J_MRI; 
    
else 
    % RHC    
    procedure = 'RHC'; 

    HR = (HFpEFdata.HR_RHC_max - HFpEFdata.HR_RHC) * theta + HFpEFdata.HR_RHC; %bpm
    freq = HR/60; %Hz
    stim_per = 1/freq; 
    t = 0:.01:5*stim_per; 

    data.time     = t; 
    data.HR       = HR; 
    data.stim_per = stim_per; 
    data.procedure = procedure; 

    pars_RHC = parameters(exercise,p,theta); 
    [o_RHC,rout_RHC,J_RHC] = model_sol(pars_RHC,theta,data); 
    
    outputs = o_RHC; 
    rout    = rout_RHC; 
    J       = J_RHC; 
end 