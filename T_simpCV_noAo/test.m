




t = 0:.01:2; 
HR = 60; 

T     = 60 / HR; 
tauR  = 48 * 1e-3; 


Frise = zeros(size(t)); 
for i = 1:length(t) 
    tc    = mod(t(i), T); 
    x     = min(8, max(0, tc / tauR)); 
    Frise(i) = 0.02 * x^3 * (8 - x)^2 * exp(-x);
end 

figure(10)
clf
plot(t,Frise,'b')


k_TS = .1; 
k_TR = .3; 
TS = k_TS * T; 
TR = k_TR * T; 

y_a = zeros(size(t)); 
for i = 1:length(t) 
    tc    = mod(t(i)+.2, T); 
    x     = min(8, max(0, tc / tauR)); 
    if tc >= 0 && tc < TS 
        y_a(i) = 0.5*(1 - cos(pi*tc/TS)); 
    elseif tc >= TS && tc < TR + TS 
        y_a(i) = 0.5*(1 + cos(pi*(tc - TS)/TR)); 
    else
        y_a(i) = 0; 
    end 
end 


figure(10)
hold on 
plot(t,y_a,'r')