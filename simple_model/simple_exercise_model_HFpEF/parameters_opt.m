function pars = parameters_opt(exercise,pars_opt,theta) 



%% Exercise factors 

a = exercise(1); 
b = exercise(2); 
c = exercise(3); 
d = exercise(4); 
e = exercise(5); 
f = exercise(6); 

%% Parameters 

pars = pars_opt; 

pars(8)  = pars_opt(8)  - log(1 + d * theta); %C_SA
pars(9)  = pars_opt(9)  - log(1 + c * theta); %C_SV
pars(16) = pars_opt(16) - log(1 + b * theta); %R_SA
pars(19) = pars_opt(19) - log(1 + e * theta); %R_PA
pars(27) = pars_opt(27) + log(1 + a * theta);  %k_act
pars(28) = pars_opt(28) - log(1 + f * theta); %Crest 


end 