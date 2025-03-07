% calcul du courant du au runaway pour le debut du plasma
function jrun = zjrun(epar,ne,te,zeff,ftrap,tauj,epr_in)

% constante physique (phys)
c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
ee           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
epsi0       =   1./c.^2./mu0;             % permitivite du vide (F/m)  (definition)
me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)

%epar  = max(eps,abs(epar));
kzeff = max(0.05,pchip([1 2 10],[0.43 0.32 0.14],zeff));

ld    = 14.9 - 0.5 .* log(ne ./ 1e20) + log(te./ 1e3);
vth   = sqrt(ee .* te ./ me);
ec    = (1 + ftrap) .* ne .* ee .^ 3 .* ld ./ 4 ./ pi ./ epsi0 .^ 2 ./me ./ vth .^ 2;
er    = (1 + ftrap) .* ne .* ee .^ 3 .* ld ./ 4 ./ pi ./ epsi0 .^ 2 ./me ./ c .^ 2;
vc    = min(c,sqrt(ne .*ee .^ 3 .* ld ./ 4 ./ pi ./ epsi0 .^ 2 ./ me ./ max(1e-6,epar) .* (2 + zeff)));
vc    = min(c,(vc ./c) .* sqrt(sqrt(vc .^ 4 + 4 .* c .^ 4 ) - vc .^ 2)  ./ sqrt(2));
tau    = 4 .* pi .* epsi0 .^ 2 .* me .^ 2 .* vth .^ 3 ./ ne ./ ee .^ 4 ./ ld;
% on prend un minimum pour eviter les instabilites numerique
%epr    = max(epar ./ ec, 0.01);
% epr(1) = epr_in;
%epr    = 0.03;
%  if epr_in < 0
%  	epr  = max(0.03,min(abs(epr_in ./ ec),abs(epar ./ec))); 
%  else
%  	epr  = max(0.03,min(epr_in,abs(epar ./ec))); 
%  end 
%epr  = max(0.03,epar ./ ec);
epr  = max(eps,epar ./ ec);
% au moins un point doit atteindre epr = 1 pour qu'il y ait des runaways
%  ind_low = find(all(epr < 1,2));
%  if ~isempty(ind_low)
%      ve = ones(1,size(epr,2));
%      epr(ind_low,:) = epr(ind_low,:) ./ (max(epr(ind_low,:),[],2) * ve);
%  end
ve  = ones(1,size(epr,2));
epr = epr  ./ (max(epr,[],2) * ve);

se     = epr .^ -(3 ./ 16 .* (1 + zeff)) .* exp( - 1 ./ epr ./ 4  - sqrt((1+zeff) ./ epr));
jrun   = max(eps,ee .* ne .* (1 - ftrap) .* vc .* kzeff .* se ./ tau .* (tauj * ones(1,size(ne,2))));
% can't be a Dirac
jrun(:,2) = max((jrun(:,1) + jrun(:,2)) ./ 2,jrun(:,2));
jrun(:,3) = max((jrun(:,2) + jrun(:,3)) ./ 2,jrun(:,3));
jrun(:,end) = 0; 

%  figure(21);clf;
%  subplot(2,1,1)
%  plot(epr);
%  subplot(2,1,2);
%  plot(jrun');
%  drawnow
%  keyboard
%  %  
%  ind = find((jrun(:,1).*eps) > max(jrun,[],2));
%  if ~isempty(ind)
%    keyboard
%  end
