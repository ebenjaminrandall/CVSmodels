function [J,rout,outputs] = model_wrap_baseline(p0,data)

INDMAP = data.gvars.INDMAP; 
ALLPARS = data.gvars.ALLPARS; 

tpars = ALLPARS; 
tpars(INDMAP) = p0; 

% Baseline 
theta = 0; 
[outputs,rout,J] = model_sol(tpars,theta,data); 

end 