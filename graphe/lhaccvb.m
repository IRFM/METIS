function [x,fpout,effacc,xlh,dlh,lc,hc,acc,landau] = lhaccvb(flh,npar0,width,agaz,zgaz,temps,x,nep,tep,qp,equi,plh,xlhin,dlhin,phys);

% securite sur npar0
npar0    = max(1,min(10,npar0));

% vecture utils
ve = ones(size(x));
vt = ones(size(temps,1),1);
va = ones(size(temps,1),1);


% champ total
btot      = sqrt(equi.b2);
% rmx
rmx       = vt * (x .* equi.rhomax);
% Raxe
Raxe      = vt * equi.raxe;
% q
qp      = vt * qp;

% frequence
w    = (2 .* pi .* flh) * ve;
wpe  = va * (sqrt(nep .* phys.e .^ 2 ./ phys.me ./ phys.epsi0));
wpi  = va * (sqrt(nep ./ (zgaz *ve) .* phys.e .^ 2 ./ phys.mp ./ (agaz * ve) ./ phys.epsi0));
wce  = va * (phys.e ./ phys.me .* btot);
% SDP
S = 1 + wpe .^2 ./ wce .^2 - wpi .^ 2 ./ w .^ 2;
D = wpe .^ 2 ./ w ./ wce;
P = 1 - wpe .^2 ./ w .^2 - wpi .^2 ./ w .^ 2;
% courbes
lc     = (npar0*ve) ./ (1 + rmx ./ qp ./ Raxe .* sqrt(max(eps,- P ./ S)));
hc     = (npar0*ve) ./ (1 - rmx ./ qp ./ Raxe .* sqrt(max(eps,- P ./ S)));
acc    = wpe ./ wce + sqrt(1 + wpe .^ 2 ./ wce .^ 2 - wpi .^2 ./ w .^2);

effacc = min(1,exp((npar0 * ve  -  acc)));


landau = 6.5 ./ sqrt(max(30,tep) ./1e3);
indbad  = find((hc <lc) | (hc > (6.5./sqrt(0.03))));
hc(indbad) = (6.5./sqrt(0.03));
hc(find(imag(hc))) = npar0(find(imag(hc)));
hc(find(~isfinite(hc))) = npar0(find(~isfinite(hc)));

% largeur naturelle
dn0  = (phys.c ./flh ./ width) * ve; % depend du grill

% boucle de recherche de l'acrochage
mul = 1;
fin = 0;
nb = 100;
indlh =  ones(size(temps));
landau_mat = va * landau;
vali  =  zeros(size(temps));
mulvec = NaN .* indlh;
while (fin == 0 ) & ( nb >0)
   % condition de landau + effet de largeur de spectre
   landau_mul  = va * (6.5 ./ sqrt(max(30,tep) .* mul ./ 1e3));

   % l'onde ne peut pas se coupler sous acc mais se propage
   % recherche du n// de premiere absorbtion
   landau2 = landau_mul;
   landau2((landau_mul < acc) | (landau_mul < lc) | (landau_mul > hc)) = Inf;
   dd  = abs(landau2 - npar0*ve);
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
nparabs = nparabs * ve;

% gaussienne d'absorbsion
% gaussienne ecart a l'optimum . la largeur est change par la propagation de l'onde (dispersion)
%dlhabs  = max(0.1,dn0 .* max(1,nparabs ./ (npar0*ve)));
dlhabs  = max(0.1,dn0 .* max(medfilt1(pdederive(x,landau,0,2,2,1)./landau,3),max(1,nparabs ./ (npar0*ve))));

pabs   = (va*nep) .* exp(-0.5 .* ((nparabs - landau) ./ dlhabs).^ 2);
pabs   = pabs ./ max(eps,(max(pabs,[],2) * ve));


indpbad = find(sum(pabs,2) == 0);
fplh    = exp(-(vt*x -(xlhin .*(va * ve))).^ 2 ./ ((max(dlhin,0.05) .* (va * ve)) .^ 2) ./ 2);
pabs(indpbad,:)  = fplh(indpbad,:);

% la valeur de dlh en sortie est differente de cette definition
dlh      = max(0.05,sqrt(trapz(x,abs(vt*x) .* pabs .* (vt * x - xlh * ve) .^ 2,2) ./  max(eps,trapz(x,abs(vt*x) .* pabs,2))));
fpout    = pabs;


h = findobj(0,'type','figure','tag','z0acclh');
if isempty(h)
       h=figure('tag','z0acclh');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(2,2,1)
plot(x,lc,'b',x,hc,'b',x,acc,'r',x,landau_mat,'g',xlh,nparabs,'om',1,npar0,'*k');
axis([0,1,0,10])
subplot(2,2,2)
plot(x,pabs,'b')
subplot(2,2,3)
plot(temps,xlh,'or',temps,dlh,'ob',temps,dlhabs,'og');
set(gca,'ylim',[0,1])
subplot(2,2,4)
plot(temps,mean(nparabs,2),'ob',temps,mean(mulvec,2),'or');
set(gca,'ylim',[0,10])
drawnow
