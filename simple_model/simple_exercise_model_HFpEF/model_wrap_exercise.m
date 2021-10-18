function [J,rout,outputs] = model_wrap_exercise(p,pars_opt,data)


%% 70% exercise 

theta_70 = 0.7;
HR_70   = 64*(1 + 1.9*theta_70);  % 1/sec
freq_70 = HR_70/60; %Hz
stim_per_70 = 1/freq_70;

t_70 = 0:0.01:5*stim_per_70; 

data.t        = t_70; 
data.HR       = HR_70; 
data.stim_per = stim_per_70; 

pars_70 = parameters_opt(p,pars_opt,theta_70); 
[o_70,rout_70,J_70] = model_sol(pars_70,theta_70,data); 

%% 100% exercise

theta_100 = 1; 
HR_100    = 64*(1 + 1.9*theta_100);  % 1/sec
freq_100  = HR_100/60; %Hz
stim_per_100 = 1/freq_100;

t_100 = 0:0.01:5*stim_per_100; 

data.t        = t_100; 
data.HR       = HR_100; 
data.stim_per = stim_per_100; 

pars_100 = parameters_opt(p,pars_opt,theta_100); 
[o_100,rout_100,J_100] = model_sol(pars_100,theta_100,data); 

%% Outputs 

outputs.o_70  = o_70; 
outputs.o_100 = o_100; 

rout = [rout_70; rout_100]; 

J = J_70 + J_100; 

end 