% computation of global relativistic correction to Bremsstrahlung
% from: Fusion reactivity of the pB11 plasma revisited,  
% S.V. Putvinski et al 2019 Nucl. Fusion 59 076018
% https://doi.org/10.1088/1741-4326/ab1a60
% formula 7 divide bu non relativistic formula
function cor = brem_rel_cor(te,Zeff)

% physical constants
phys = cphys;
%
t = phys.e .* te ./ phys.me ./ phys.c .^2;

% ion part fully ionised + electron part
cor = 1 + 1.78 .* t .^ 1.34 + ...
      2.12 .* t .* (1 + 1.1 .* t + t .^ 2  - 1.25 .* t .^ 2.5) ./ max(1,Zeff);
  
% for comparaison only
return
  
% previous formula from P.E. Stott, Plasma. Physics and Controlled Fusion 47 (2005) 1305-1338 ref [15
xrel   = (1 + 2 .* t) .* (1 + (2 ./ Zeff) .* (1 - (1 + t) .^ (-1)));

figure;
subplot(2,2,1)
plot(te,cor,'r',te,xrel,'b');

subplot(2,2,2)
semilogx(te,cor,'r',te,xrel,'b');
  
subplot(2,2,3)
semilogy(te,cor,'r',te,xrel,'b');
    
subplot(2,2,4)
loglog(te,cor,'r',te,xrel,'b');

    

