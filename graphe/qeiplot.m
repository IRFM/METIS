% script : trace le graphe de Qei(ne,Te) 
x   = linspace(0,1,101);
fte = trapz(x, x.* (1-x.^2) .^2); 
fne = trapz(x, x.* ((1-x.^2) .^0.5) .^2); 
te = linspace(0.1,10)';
ne = [1e18,2e18,5e18,1e19,2e19,5e19,1e20];
tem = te * ones(size(ne));
nem = ones(size(te)) * ne;
me  =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
mp  =   1.6726485e-27;            % masse au repos du proton (kg)
ee  =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
vol =  25;                        % volume plasma
taue = 6.4e14 .* tem .^ (3/2) ./ nem;
tauex =  taue .* mp ./ me ./ 2;
pei   = (3/2) .* ee .* 1e3 .* nem ./ tauex .* vol .* fte;

h1 =findobj(0,'type','figure','tag','qeiplot');
if isempty(h1)
	  h1=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
  	          'times','defaultlinelinewidth',3,'tag','qeiplot','name','ECRH');
else
	  figure(h1);
end
clf


semilogy(te,pei./1e6);
xlabel('<Te> (keV)');
ylabel('Pei/(Te -Ti) (MW/keV)');
title('Valeur du terme d''equipartition'); 
axis([0 10 0.001 10])
grid on

hold on
try
   p = mean(data.gene.ploss(isfinite(data.gene.ploss)))/1e6;
   q = mean(data.gene.qei(isfinite(data.gene.qei)))/1e6;
   t = mean(data.gene.temoy(isfinite(data.gene.temoy)))/1e3;
   n = mean(data.gene.nbar(isfinite(data.gene.nbar)));
   tw = 6.4e14 .* t .^ (3/2) ./ n;
   tw =  tw .* mp ./ me ./ 2;
   w   = (3/2) .* ee .* 1e3 .* n ./ tw .* vol/1e6 .* fte;
   plot([0,10],[p,p],'g');
   plot(t,w,'ko');
   plot(t,q,'ro');
   legend('n_b_a_r = 1e18','n_b_a_r = 2e18','n_b_a_r = 5e18','n_b_a_r = 1e19', ...
       'n_b_a_r = 2e19','n_b_a_r = 5e19','n_b_a_r = 1e20', ...
       'Ploss',sprintf('# %d',fix(param.from.shot.num)),sprintf('Pei reel de # %d',fix(param.from.shot.num)));
   
catch   
   legend('n_b_a_r = 1e18','n_b_a_r = 2e18','n_b_a_r = 5e18','n_b_a_r = 1e19', ...
       'n_b_a_r = 2e19','n_b_a_r = 5e19','n_b_a_r = 1e20');
   
end
hold off
