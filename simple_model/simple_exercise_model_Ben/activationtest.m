

clear all

load trisegC.mat 

%% Ventricles 

TS = .35 * stim_period; 
TR = .35 * stim_period; 


y = zeros(size(t)); 
for i = 1:length(t)
    
    tn_v  = mod(t(i) - 0.1,stim_period)/stim_period; 
    if tn_v >= 0 && tn_v < TS 
        y(i) = 0.5*(1 - cos(pi*tn_v/TS)); 
    elseif tn_v >= TS && tn_v < TR + TS 
        y(i) = 0.5*(1 + cos(pi*(tn_v - TS)/TR)); 
    else
        y(i) = 0; 
    end 
    
    
end 


E_M = 5; 
E_m = min(C_LV); 
E = (E_M - E_m)/2 * y + E_m; 

figure(10)
clf 
hold on 
h1 = plot(t,C_LV,'b'); 
h2 = plot(t,E,'r'); 
legend([h1 h2],'C_{LV}','E')

%% Atria 


Tact = +0.015;

for i = 1:length(t)
    
    tn_v  = mod(t(i) - 0.1,stim_period)/stim_period; 
    phi_atria(i) = tn_v - Tact - ((tn_v - Tact)>(0.5)) + ((tn_v - Tact)<(-0.5));
end 
% Emin = 0.05 + 0.20;
% Emax = 0.15 + 0.20; 
Emin    = 1.5 * 0.050 ;
Emax    = 0.375; %0.150 ;
sigma_a = 0.100;
act     = exp( -(phi_atria / sigma_a).^2 );

E_LA = 3.0*(Emin + (Emax-Emin)/2*act); 
E_RA = (Emin + (Emax-Emin)/2*act); 


figure(11) 
clf
hold on 
plot(t,E_RA)

TS_a = .25 * stim_period; 
TR_a = .2 * stim_period; 

for i = 1:length(t)
    tn_a = mod(t(i) + 0.1,stim_period)/stim_period;     
    
    if tn_a >= 0 && tn_a < TS_a
        y_atria(i) = 0.5 * (1 - cos(pi * tn_a / TS_a)); 
    elseif tn_a >= TS_a && tn_a < TR_a + TS_a 
        y_atria(i) = 0.5 * (1 + cos(pi * (tn_a - TS_a) / TR_a)); 
    else
        y_atria(i) = 0; 
    end 
end 

E_atria   = (Emax - Emin)/2 * y_atria + Emin; 


figure(11)
hold on
plot(t,E_atria)


