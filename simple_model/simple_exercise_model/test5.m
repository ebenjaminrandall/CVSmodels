

% Normal EDPVR relationship 
EDP = 5; 
EDV = 140; 

ESP = 125; 
ESV = 50; 

%% Find EDPVR curve 

V_D_l = 1 * EDV;
V_D_r = V_D_l; 

An = 28; 
Bn = 3; 

V_0  = EDV * (0.6 - 0.006 * EDP); 
V_30 = V_0 + (EDV - V_0) / ((EDP / An)^(1/Bn));
beta = log(EDP/30) / log(EDV / V_30); 
alpha = 30 / V_30^beta; 

P_D_l = alpha * V_D_l^beta;
%P_D_r = P_D_l / 6

P = @(V) alpha .* V.^beta; 

x0 = [1; 1; 1; .05]; 
fun = @(x) passiveparameters(x,P,EDV,ESV,V_0,V_30);

options = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1e6,'MaxIterations',1e6); 
xopt = fsolve(fun,x0,options); 



function y = passiveparameters(x,P,EDV,ESV,V_0,V_30) 


a = x(1); 
b = x(2); 
c = x(3); 
gamma = x(4); 

L1 = 1.6; 
L2 = 1.51; 
%Crest = 0.02; 
Vw_LV = 80; 

f = @(SL) (SL - L1)./L1;
g = @(SL,gamma) exp(gamma .* (SL - L1)/L1) - 1; 
h = @(SL,Crest) Crest .* (SL - L2); 

xm = @(V) -0.009 .* V - 3.32; 
ym = @(V) .0072 .* V + 2.6; 

Am = @(V) pi .* (xm(V).^2 + ym(V).^2); 
Cm = @(V) 2 .* xm(V) / (xm(V).^2 + ym(V).^2); 
z  = @(V) 3 .* Cm(V) .* Vw_LV ./ (2 .* Am(V)); 

sinalpha = @(V) 2 .* xm(V) .* ym(V) ./ (xm(V)^2 + ym(V).^2); 
Gamma    = @(V) Vw_LV ./ (2 .* Am(V)) .* (1 + (1/3) .* z(V).^2 + (1/5) .* z(V).^4); 

point_0   = -P(V_0)  .* ym(V_0)  ./ (2 .* sinalpha(V_0)  .* Gamma(V_0));
point_5   = -P(EDV)  .* ym(EDV)  ./ (2 .* sinalpha(EDV)  .* Gamma(EDV)); 
point_10  = -P(P(10)) .* ym(P(10)) ./ (2 .* sinalpha(P(10)) .* Gamma(P(10))); 
point_125 = -P(ESV)  .* ym(ESV)  ./ (2 .* sinalpha(ESV)  .* Gamma(ESV));

y(1) = -point_0   + a * f(2.2) + b * g(2.2,gamma) + c * h(2.2,0.02); 
y(2) = -point_5   + a * f(2.2) + b * g(2.2,gamma) + c * h(2.2,0.02); 
y(3) = -point_10  + a * f(2.2) + b * g(2.2,gamma) + c * h(2.2,0.02); 
y(4) = -point_125 + a * f(1.8) + b * g(1.8,gamma) + c * h(1.8,2); 

y*y'

end 








