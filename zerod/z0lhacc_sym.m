% diagramme LH
function [x,fpout,xlh,dlh,xmem,lc,hc,acc,landau] = z0lhacc(flh,npar0,width,agaz,zgaz,temps,x,nep,tep,qp,Raxe,rmx,spr,vpr,fdia,plh,xlhin,dlhin)


phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)


% extention des donnees cote LFS et HFS
x    = cat(2,-x(end:-1:2),x);
nep  = cat(2,nep(:,end:-1:2),nep);
tep  = cat(2,tep(:,end:-1:2),tep);
qp   = cat(2,qp(:,end:-1:2),qp);
Raxe = cat(2,Raxe(:,end:-1:2),Raxe);
fdia = cat(2,fdia(:,end:-1:2),fdia);
rmx  = cat(2,rmx(:,end:-1:2),rmx);
spr  = cat(2,spr(:,end:-1:2),spr);
vpr  = cat(2,vpr(:,end:-1:2),vpr);

% vecture utils
ve = ones(size(x));
vt = ones(size(tep,1),1);

% champ total
btor      = fdia ./ (Raxe + sign(vt * x) .* rmx);
bpol      = sign(vt * x) .* rmx .* btor ./ Raxe ./ qp;
btot      = sqrt( btor .^2  + bpol .^2);
qp        = max(1,min(qp(:,end) *ve,qp));

% frequence
w    = 2 .* pi .* flh;
wpe  = sqrt(nep .* phys.e .^ 2 ./ phys.me ./ phys.epsi0);
wpi  = sqrt(nep ./ (zgaz *ve) .* phys.e .^ 2 ./ phys.mp ./ (agaz * ve) ./ phys.epsi0);
wce  = phys.e ./ phys.me .* btot;
% SDP
S = 1 + wpe .^2 ./ wce .^2 - wpi .^ 2 ./ w .^ 2;
D = wpe .^ 2 ./ w ./ wce;
P = 1 - wpe .^2 ./ w .^2 - wpi .^2 ./ w .^ 2;
% courbes
lc     = npar0 ./ (1 + rmx ./ qp ./ Raxe .* sqrt(- P ./ S)); 
hc     = npar0 ./ (1 - rmx ./ qp ./ Raxe .* sqrt(- P ./ S));
acc    = wpe ./ wce + sqrt(1 + wpe .^ 2 ./ wce .^ 2 - wpi .^2 ./ w .^2);
landau = 6.5 ./ sqrt(tep ./1e3);
indbad  = find((hc <lc) | (hc > (6.5./sqrt(0.03))));
hc(indbad) = (6.5./sqrt(0.03));
hc(find(imag(hc))) = npar0;
hc(find(~isfinite(hc))) = npar0;

% largeur naturelle
dn0  = phys.c ./flh ./ width; % depend du grill

% boucle de recherche de l'acrochage
mul = 1;
fin = 0;
nb = 100;
indlh =  ones(size(temps));
landau_mat = landau;
vali  =  zeros(size(temps));
mulvec = NaN .* indlh;
while (fin == 0 ) & ( nb >0)

   landau_mul  = 6.5 ./ sqrt(max(30,tep) .* mul ./ 1e3);

   % l'onde ne peut pas se coupler sous acc mais se propage
   % recherche du n// de premiere absorbtion
   landau2 = landau_mul;
   landau2((landau_mul < acc) | (landau_mul < lc) | (landau_mul > hc)) = Inf;
   dd  = abs(landau2 - npar0);
   dd  = abs(landau2 - (vt * sign(x)) .* npar0 .* sqrt((qp(:,end)*ve) ./ qp));
   [ddmin,indlhs] = min(dd,[],2);
   % indice  ok
   indok = find(isfinite(ddmin) & (vali == 0));
   indlh(indok) = indlhs(indok);
   landau_mat(indok,:) = landau_mul(indok,:);
   vali(indok) = 1;
   mulvec(indok) = mul;


   nb = nb - 1;
   mul = mul .* 1.05;
   if all(vali ==1)
      fin =1;
   end
end
% nouvelle valeur en sortie
xlh  = x(indlh)';
landau = landau_mat;

nparabs = diag(landau(:,indlh));
indnonok = find(~isfinite(nparabs));
indok = find(isfinite(nparabs));
if ~isempty(indnonok)
   nparabs(indnonok) = mean(nparabs(indnonok));
end
nparabs = nparabs *ve;

% gaussienne d'absorbsion
% gaussienne ecart a l'optimum . la largeur est change par la propagation de l'onde (dispersion)
dlhabs  = dn0 .* (nparabs ./ npar0);
%dlhabs  = max(dn0,abs(hc - max(lc,acc)))./4;
pabs   = nep .* exp(-0.5 .* ((nparabs - landau) ./ dlhabs).^ 2);
pabs   = pabs ./ max(eps,(max(pabs,[],2) * ve));

% repliement
xlhmem = xlh;
xmem   = x;
x      = x(21:end);
ve     = ve(21:end);
xlh = abs(xlh);
pabsmem = pabs;
pabsr(:,2:21) = (pabs(:,20:-1:1) + pabs(:,22:end)) ./ 2;
pabsr(:,1) = pabs(:,21);
pabs = pabsr;

indpbad = find(sum(pabs,2) == 0);
fplh    = exp(-(vt*x -(xlhin * ve)).^ 2 ./ ((max(dlhin,0.05) * ve) .^ 2) ./ 2);
pabs(indpbad,:)  = fplh(indpbad,:);

% la valeur de dlh en sortie est differente de cette definition
%dlh      = max(0.05,sqrt(trapz(x,pabs .* (vt * x - xlh * ve) .^ 2,2) ./  max(eps,trapz(x,pabs,2))));
dlh      = max(0.05,sqrt(trapz(x,abs(vt*x) .* pabs .* (vt * x - xlh * ve) .^ 2,2) ./  max(eps,trapz(x,abs(vt*x) .* pabs,2))));
%dlh       = max(0.05, mean(sqrt(abs((vt * x - xlh * ve) .^ 2 ./ 2 ./ max(eps,-log(max(eps,pabs ./ (max(pabs,[],2)*ve)))))),2));
% anti oscillation ?
fpout    = pabs;
%dlh   = 0.5 .* (dlh + mean(dlh));



indp = find(plh >1e5);
if isempty(indp)
       return
end
mref =mean(pabs(indp,:),1);
mfp  = mean(fplh(indp,:),1);
figure(51);
clf
subplot(2,2,1)
plot(xmem,lc(indp,:),'b',xmem,hc(indp,:),'c', ...
     xmem,acc(indp,:),'r',xmem,landau_mat(indp,:), ...
     'g',xmem(indlh(indp)),nparabs(indp),'om',1,npar0,'*k');
axis([-1,1,0,30])
subplot(2,2,2)
plot(x,pabs(indp,:),'b',x,mref./max(mref),'r')
subplot(2,2,3)
plot(temps(indp),xlh(indp),'r',temps(indp),dlh(indp),'b',temps(indp),dlhabs(indp),'g');
set(gca,'ylim',[0,1])
subplot(2,2,4)
plot(temps(indp),mean(nparabs(indp,:),2),'b',temps(indp),mean(mulvec(indp,:),2),'r');
set(gca,'ylim',[0,10])
drawnow
%keyboard
