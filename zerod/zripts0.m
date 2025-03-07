% perte ripple : fci = couche centrale, 10% H dans D par defaut
function [pfci,plh,ptherm,einjfci,einjlh] = zripts0(pfci,plh,nbar,a,r0,b0,betap,li,ip,tem,ate,nem,ane,nDm,cmin,mino,freq,tite,qa,qmin,meff,tebord,nebord,vp,R1H)


% constante physique (phys)
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)


nbar = nbar ./ 1e19;
nl   = nbar .* 2 .* a;
ip   = ip  ./ 1e6;
te0  = tem .* (1+ate) ./ 1e3;
ne0  = nem .* (1+ane) ./ 1e19;
if cmin == 0 
   nH   = (1e13  + nDm .* (1+ane)) ./ 1e19; 
else
   nH   = (1e13  + cmin .* nDm .* (1+ane)) ./ 1e19; 
end
% initialisation
Ilost = zeros(size(pfci));
Emean       = [];
beli = betap + li ./ 2;
ra   = a ./ r0;
d0   = a .^ 2 ./ r0 ./ 2 .* beli .* (1 - ra);
indfci = find(pfci > 4e5);
if ~isempty(indfci)
  geo.a = a(indfci);
  geo.r0 = r0(indfci);
  geo.b0 = b0(indfci);
  geo.d0 = d0(indfci);
  pos  = geo.r0 + geo.a + 0.02;
  % sert juste a determiner si le scenario est bien celui pour lequel le fit fonctionne
  [scenar,scenstr,Rres]=scenarzerodfci(geo,freq,pos,mino);
  if isfinite(Rres)

  % perte ripple de TS pour icrh
    Ilost       = 18.*5.72e-2 .* (pfci/1e6).^1.01 .* nl.^(-0.87) .* R1H .^ 9.1 .* ip.^0.2; 
    val         = pfci ./ 1e6 ./ nH .*  te0 .^ 1.5 ./ ne0;
    Emean       = min(8 .* val + 150,300);
    frlost      = Emean .* Ilost  ./ (1+pfci);
    frlost(~isfinite(frlost)) = 0;
    pfci = pfci .* (1  - min(1,(pfci > 350e3) .* frlost));
  end
else
   Ilost = zeros(size(pfci));
   Emean       = [];
end

% changement d'unite
einjfci  = Emean .* 1e3;
 
% flux ripple d'ions D/H 
flux_fci = Ilost./phys.e ./ 1000;

% perte pour lh
% Parametres du Ripple (modif 02/09/03)  (loi d'�helle pour le diagnostic DRIPPLE)
alpha = [0.6662639631967 -0.503759394462734 0.235509558581612;0.942854084818663 -0.212233746620332 -0.606046698438498;1.17449344302278 0.406956986489841 -1.55647263246365;1.26490295439731 0.807360991632855 -1.77323814523948;1.34391626203478 0.829563873473516 -1.67640737283367;1.27263426767301 0.445323609124571 -1.06744934260616;0.886779152061443 -0.4217528815393 -0.069881000409826;0.146970707366509 -0.758932144063924 0.917539880313246];
dalpha = [0.00155681373475998 0.00518733657582794 0.00765538869327886;0.00104073495105042 0.00346775105896577 0.00511765177752616;0.000831525873545265 0.00277066194964786 0.00408889877342242;0.000773591595767396 0.00257762370017632 0.00380401599961882;0.000934569549592637 0.00311400567647508 0.00459559997660044;0.00123139407027466 0.00410303131156422 0.00605518825539183;0.00152775646084226 0.00509051711925342 0.0075125040814332;0.00186667541971748 0.00621980232040562 0.0091790851937274];
c0 = [2.85482634661752;2.24535175106449;8.06969332838054;14.8139712278853;10.4570952588182;3.32534745291915;0.533723654679241;0.0812120305892144];
%
% premire version en attendant une loi d'�helle sur le pertes par calorim�rie
% loi d'�helle pour l'�ergie
escal         = [0.6 0.4];
%
% loi d'�helle pour l'��ation en temp�ature
%
cT            = 3.45;
aT            = [1.74 -1.92 1.06];
%  
% estimation nl voie 2 de TS
%
nl2            = nl * 0.72;
ve             = ones(size(c0'));
vt             = ones(size(nl2));
c0             = vt * c0';
profil         = c0 .* ((plh * ve)./ 1e6) .^ (vt *alpha(:,1)') .* ...
                 (ip*ve) .^  (vt * alpha(:,2)') .* ...
		 (nl2*ve) .^  (vt*alpha(:,3)');
Erip          = 130 .* (ip*ve) .^ 0.6 .* (nl2*ve) .^ 0.4;
ploss         = Erip .*  profil .* 18;
%
% facteur 1.4 tenant compte du fait que DRIPPLE ne voit que 70 % des pertes 
%
ptot          = sum(ploss,2)*1.4;
plh           = plh - ptot;

% flux d'electron
flux_lh       = ptot ./ sum(Erip,2) ./ phys.e;
einjlh        =  mean(Erip,2) .*1e3 ./ 1.4 ;

% calcul du ripple thermique de TS
x   = linspace(0,1,21); 
ux  = (1  - x .^ 2);
ve  = ones(size(x));
vt  = ones(size(tem));
nebord  = max(1e13,min(nem/2,nebord));
nep = ((nem .* (1 + ane)-nebord) * ve)  .* (vt * ux)  .^ (ane * ve) + nebord * ve;   
%nDp = (((nDm - nebord) .* (1 + ane)) * ve)  .* (vt * ux)  .^ (ane * ve) + nebord * ve;
nDp = nep; % une seule espece;
tebord  = max(13.6,min(300,tebord));
tep = ((tem .* (1 + ate)-tebord) * ve)  .* (vt * ux)  .^ (ate *ve) +tebord *ve;
tip = (tite * ve) .* tep;
meffp = meff * ve;
qa  = min(qa,30);
%qp  = max(1,min(qmin,qa-1)) * ve + (qa -qmin) * (x .^ 2);
qp  = z0qp(x,max(1,min(qmin,qa-1)),qa);
rp  = (r0 * ve) + a * x;
rm  = (r0 * ve) - a * x;
vpr = (2 .* vp) * x + eps;

% Er neoclassique + proportionnel a Ti
erp = pdederive(x,nep.*tip,0,2,2,1) ./ nep ;  % charge de l'ion  P. Maget these

%**************************************
% COLLISIONS entre les diverses especes -> restriction a une espece d'ion
%**************************************
% Logarithme coulombien (NRL)
Lnee = 24-0.5*log(nep./1e6)+log(tep);
Lnei = 24.3-0.5*log(nep./1e6)+log(tep);  % changement de la constante d'apres Wesson
Lnii = 23 - log( sqrt(nDp./1e6) + (nDp./1e6)) + 1.5 .* log(tip);

% Temps intercollision (Braginskii)
Tauee = (3.*phys.epsi0.^2.*phys.me.^2.*(2.*pi.*phys.e./phys.me.*tep).^1.5)./(phys.e.^4.*nep.*Lnee);	
Tauei = (3.*phys.epsi0.^2.*phys.me.^2.*(2.*pi.*phys.e./phys.me.*tep).^1.5)./(phys.e.^4.*nDp.*Lnei);	
Tauie = (3.*phys.epsi0.^2.*(phys.mp.*meffp).^2.*(2.*pi.*phys.e./(phys.mp.* meffp).* tip).^1.5)./(phys.e.^4 .*nep.*Lnei);	
Tauii = (3.*phys.epsi0.^2.*(phys.mp.*meffp).^2.*(2.*pi.*phys.e./(phys.mp.*meffp).* tip).^1.5)./(phys.e.^4.*nDp.*Lnii);
 
warning off  % a cause des divisions par 0 au centre du plasma
% calcul pour chaque espece des flux de particules sortants
Nbobine = 18; % pour TS et ITER
delta = (deltaTS(rp,0)+deltaTS(rm,0))/2; % delta "moyen" de la surface magnetique consideree
nhu_el = 1./Tauee + 1./Tauei;
nhu_s = 1./Tauie+ 1./ Tauii ;

% calcul des longueurs de gradient
gne = pdederive(x,nep,0,2,2,1) ./ (a * ve);
gte = pdederive(x,tep,0,2,2,1) ./ (a * ve);
gti = pdederive(x,tip,0,2,2,1) ./ (a * ve);
gnD = pdederive(x,nDp,0,2,2,1) ./ (a * ve);
%
lmax = mean(2.* r0(isfinite(r0)));
% minimum
lmin = mean(a(isfinite(a)) ./ 100);
% longueur de gradient
lne = lcalc(nep,gne,lmin,lmax);
lte = lcalc(tep,gte,lmin,lmax);
lti = lcalc(tip,gti,lmin,lmax);
lnD = lcalc(nDp,gnD,lmin,lmax);

ar     = 2 .* delta .* Nbobine .* qp ./ (a *ve) .* abs( erp)./(b0*ve)./nhu_el;
ar(1)  = ar(2); % pour enlever le NaN du au fait que equi.a est nul en 0

flux_el = nep .* 16 .*(delta./(2*pi)).^1.5 .* (tep./(b0*ve)./(r0*ve)).^2 .*J4(ar) ...
       .* (1./lne + (J5(ar)./J4(ar)-1.5)./lte + 1./tep .* erp) ...
       ./ nhu_el;      % Flux en particules /m^2 et /s;   facteur 16, car avant integration en theta

ar = 2.*delta .*Nbobine.* qp ./ (a*ve) .* abs(erp)./(b0*ve)./nhu_s;
ar(1) = ar(2); % pour enlever le NaN du au fait que equi.a est nul en 0

flux_ion  = nDp .* 16 .*(delta./(2*pi)).^1.5 .* (tip./(b0*ve)./(r0*ve)).^2 .*J4(ar) ...
       .* (1./lnD + (J5(ar)./J4(ar)-1.5)./lti - 1./tip.* erp) ./  nhu_s;   % Flux en particules /m^2 et /s;   facteur 16, car avant integration en theta


% facteur d'integration sur l'angle poloidal (int(sin(theta)^2))
un_sur_alpha = (r0 * ve).* Nbobine .* qp .* delta ./ (a*ve);
un_sur_alpha(1)=un_sur_alpha(2);
un_sur_alpha(find(abs(un_sur_alpha)>1))=1;  % si un_sur_alpha > 1, les particules sont piegees quelque soit theta
u = asin(un_sur_alpha);
integ = u/pi - sin(2*u)/2/pi;
integ(integ<0) = 0;

% on applique l'integration sur l'angle poloidal aux flux
flux_el  = flux_el   .* integ;
flux_ion = flux_ion  .* integ;


% attention : les sources dues au ripple ont un signe : pour el,ne,ion : negatif pour des pertes
%             la rotation suit la convention cronos
%             les donnees de nrip sont negatives pour des pertes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sources de particules equivalentes; on convertit le flux radial en source volumique avec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% l'equation de conservation dn/dt = -div(flux)
nerip = - pdederive(x,vpr.*flux_el,0,2,2,1) ./ (a*ve) ./ vpr; % dn/dt = -div(flux)
nerip(1)=0; % point central a zero 
nirip = - pdederive(x,vpr.*flux_ion,0,2,2,1) ./ (a*ve) ./ vpr; % dn/dt = -div(flux)
nirip(1)=0; % point central a zero
pel  = nerip .* tep .* phys.e;
pion = nirip .* tip .* phys.e;
%source.w = gene.signe.ip .* abs(prof.psid1) .* equi.grho .* phys.e .* (-flux_el + sum(squeeze(flux_ion)'.*(compo.z'*ones(1,gene.nbrho))));

% fin de la zone protegee
warning on

% securite
pel(~isfinite(pel)) = 0;
pion(~isfinite(pion)) = 0;


% calcul de la puissance de perte du plasma au bord
ptherm      =  - trapz(x,vpr .* (pel .* (pel <0)  + pion  .* (pion <0)),2);


% securite
pfci = real(pfci);
plh = real(plh);
ptherm = real(ptherm);
einjfci = real(einjfci);
einjlh = real(einjlh);

% fin du calcul
%disp('in zripts0')
%keyboard

% calcul de la longueur de gradient
function out = lcalc(v,gv,lmin,lmax)
ind = find( abs(gv) < eps);
if ~isempty(ind)
	 gv(ind) =eps;
end
out = abs(v ./ gv);
if ~isempty(ind)
	 out(ind) = lmax;
end
out = min(lmax,max(lmin,out));

function out = J4(a)
out = -21.333 ./(1+3.4.*abs(a)+55.5.*a.^2).^0.633;

function out = J5(a)
out = -106.667 ./(1+3.37.*abs(a)+107.75.*a.^2).^0.643;

function out = deltaTS(Rconst,Zconst)
      % calcule la profondeur du ripple delta pour TS en fonction de la position (R,Z) 
      % parametres ripple pour TS (cf these Arslanbekov IV-22)
      a1=2.2e-4;
      a2=5;
      a3=1.6;
      b1=0.32;
      b2=0.26;
      Rco=2.36; % Position en R du centre de la chambre de TS (m)
      petitr=sqrt((Rconst-Rco).^2+Zconst.^2);
      costeta=(Rconst-Rco)./petitr;
      Aeq=b2*b2;
      Beq=-(2*b1*b2+2*petitr.*costeta*b2+1);
      Ceq=petitr.*petitr+2*petitr.*costeta*b1+b1*b1;
      delteq=sqrt(Beq.*Beq-4*Aeq*Ceq);
      rprime=sqrt((-Beq-delteq)/2/Aeq);
      eli=find(imag(rprime)~=0);
      rprime(eli)=zeros(size(eli));
      out=a1*exp(a2*rprime+a3*rprime.*rprime);
