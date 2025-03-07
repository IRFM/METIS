% ZRIPPLE_THERM calcule le ripple thermique et les sources associees
%-------------------------------------------------------------------------------
% fichier zripple_therm.m ->  zripple_therm
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule la distribution des pertes ripples thermiques
% du plasma. Elle calcule aussi les sources associees.
% Elle correspond a la formule analytique donnee dans la these de C. Bourdelle 
% (rapport EUR-CEA-FC-1706, octobre 2000), revisitee par P. Maget (note informelle du 14/11/2001)
% Cette formule suppose une amplitude de ripple "delta" independante de theta,
% ce qui n'est pas tout a fait realiste. Neanmoins, on l'utilise ici telle qu'elle,
% en prenant pour delta(r) la moyenne de delta en theta = 0 et pi (Z = 0, moyenne cote bas champ et fort champ)
% de la surface magnetique de rayon r. 
%  
% syntaxe  :
%  
%     [source_rip,nrip] = zripple_therm(par,source,geo,equi,prof,neo, ...
%                                  impur,bord,phys,compo,gene,premiertemps);
%
% entrees :
%
%     par        = struture des consignes de calcul de la fonction
%     source     = portotype a 0 de la structure source remplie par la fonction
%     geo        = geometrie de la derniere surface magnetique (datak.geo)
%     equi       = structure equilibre (datak.equi)
%     prof       = structure des profils (datak.prof)
%     neo        = structure des grandeurs neoclassiques (datak.neo)
%     impur      = structure des impuretes et du rayonnement (datak.impur)
%     bord       = structure decrivant le bord et le mur (datak.bord)
%     phys       = constante physiques (param.phys)
%     compo      = composition du plasma (param.compo)
%     gene       = parametre generiques (param.gene)
%                  pas utiliser dans cette fonction, reserve pour d'autres modules)
%     premiertemps = si 1, force le calcul pour avoir des donnees realiste
%
% sorties :
% 
%     source _rip    = structure source  remplie, pour le ripple
%     nrip           = profil de perte ripple par espece d'ions
% 
% parametres :
% 
%    cette fonction a des parametres 
%
% fonction ecrite par F.Imbeaux , poste 63-26
% version 2.2, du 03/03/2004.
% 
% 
% liste des modifications : 
%
%   * 03/03/2004 -> correction bug sur indice dans impur.impur
%
%--------------------------------------------------------------
%
function [source,nrip] = zripple_therm(par,source,geo,equi,prof,neo, ...
                                      impur,bord,phys,compo,gene,premiertemps)

% fonction auto declarante                             
if nargin <=1 
	langue                  = getappdata(0,'langue_cronos');
	
	% parametres de la fonction
	valeur.modecalcul = 0; 
	type.modecalcul   = 'integer'; 
	borne.modecalcul  = {0};      
	defaut.modecalcul = 0;           
	info.modecalcul   = 'si 0: amplitude du ripple prise a sa valeur bas champ (max); 1: moyenne poloidale tient compte de la dependance du ripple en theta';
        if strcmp(langue,'anglais')
          info.modecalcul = 'if 0, ripple amplitude taken at low field side (max); 1 : poloidal integration takes into account the theta dependence of the ripple';	
        end

	interface.ts        = 'a faire ts';   
	interface.jet       = 'a faire jet';  
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;

	sortie.description = 'Calcul des flux ripple thermique + sources associes';   % description (une ligne) de la fonction
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	sortie.resume ='';                           % nom de la fonction qui ecrit le resume du parametrage
	
	source = sortie;
	return

elseif nargin < 12
   premiertemps = 0;
elseif isempty(premiertemps)
   premiertemps = 0;
end

if premiertemps == 1
end

nrip = NaN .* ones(1,gene.nbrho,gene.nbg);


%securite
impur.impur = max(1e11,impur.impur);

%**************************************
% COLLISIONS entre les diverses especes
%**************************************
% Logarithme coulombien (NRL)
Lnee=24-0.5*log(prof.ne./1e6)+log(prof.te);
Lnei=24.3-0.5*log(prof.ne./1e6)+log(prof.te);  % changement de la constante d'apres Wesson
Lnii = zeros(gene.nbg,gene.nbg,gene.nbrho);
for i = 1:gene.nbg
   for j = 1:gene.nbg
      Lnii(i,j,:) = 23-log(compo.z(i).*compo.z(j).*sqrt((squeeze(impur.impur(:,:,i))./1e6).*compo.z(i).^2+(squeeze(impur.impur(:,:,j))./1e6).*compo.z(j).^2))+1.5*log(prof.ti);
      % matrice symetrique
   end
end

% Temps intercollision (Braginskii)
Tauee=(3.*phys.epsi0.^2.*phys.me.^2.*(2.*pi.*phys.e./phys.me.*prof.te).^1.5)./(phys.e.^4.*prof.ne.*Lnee);	
Tauei = zeros(gene.nbg,gene.nbrho);
Tauie = zeros(gene.nbg,gene.nbrho);
Tauii = zeros(gene.nbg,gene.nbg,gene.nbrho);
for i = 1:5
   Tauei(i,:)=(3.*phys.epsi0.^2.*phys.me.^2.*(2.*pi.*phys.e./phys.me.*prof.te).^1.5)./(phys.e.^4.*compo.z(i).^2.*max(1e11,squeeze(impur.impur(:,:,i))).*Lnei);	
   Tauie(i,:)=(3.*phys.epsi0.^2.*(phys.mp.*compo.a(i)).^2.*(2.*pi.*phys.e./(phys.mp.*compo.a(i)).*prof.ti).^1.5)./(phys.e.^4.*compo.z(i).^2.*prof.ne.*Lnei);	
   for j = 1:5
      Tauii(i,j,:) = (3.*phys.epsi0.^2.*(phys.mp.*compo.a(i)).^2.*(2.*pi.*phys.e./(phys.mp.*compo.a(i)).*prof.ti).^1.5)./(phys.e.^4.*compo.z(i).^2.*compo.z(j).^2.*max(1e11,squeeze(impur.impur(:,:,j))).*Lnii(i,j));
      % matrice non symetrique
   end
end

%Taueimp=(3*phys.epsi0^2*me^2*(2*pi*qe/me*TE).^1.5)./(qe^4*Zimp^2*NImp.*LneimpNRL);		
%Tauimpe=(3*phys.epsi0^2*mimp^2*(2*pi*qe/mimp*TI).^1.5)./(qe^4*Zimp^2*NE.*LneimpNRL);
%Tauiimp=(3*phys.epsi0^2*mi^2*(2*pi*qe/mi*TI).^1.5)./(qe^4*Zi^2*Zimp^2*NImp.*LniimpNRL);
%Tauimpi=(3*phys.epsi0^2*mimp^2*(2*pi*qe/mimp*TI).^1.5)./(qe^4*Zi^2*Zimp^2*NI.*LniimpNRL);
%Tauimpimp=(3*phys.epsi0^2*mimp^2*(2*pi*qe/mimp*TI).^1.5)./(qe^4*Zimp^4*NImp.*LnimpimpNRL);
%Tauei=(3*phys.epsi0^2*me^2*(2*pi*qe/me*TE).^1.5)./(qe^4*Zi^2*NI.*LneiNRL);	
%Tauie=(3*phys.epsi0^2*mi^2*(2*pi*qe/mi*TI).^1.5)./(qe^4*Zi^2*NE.*LneiNRL);	
%Tauii=(3*phys.epsi0^2*mi^2*(2*pi*qe/mi*TI).^1.5)./(qe^4*Zi^4*NI.*LniiNRL);	
%#############################################################################
warning off  % a cause des divisions par 0 au centre du plasma
% calcul pour chaque espece des flux de particules sortants
Nbobine = 18; % pour TS et ITER
if ~par.modecalcul
   delta = (deltaTS(squeeze(double(equi.R(:,:,1))),0)+deltaTS(squeeze(double(equi.R(:,:,ceil(size(equi.R,3)/2)))),0))/2; % delta "moyen" de la surface magnetique consideree
   nhu_el = 1./Tauee + sum(1./Tauei);
   a = 2.*delta .*Nbobine.* prof.q ./ equi.a .* abs(neo.er)./sqrt(prof.bphi.^2+prof.bpol.^2)./nhu_el;
   a(1) = a(2); % pour enlever le NaN du au fait que equi.a est nul en 0

   flux_el = prof.ne .* 16 .*(delta./(2*pi)).^1.5 .* (prof.te./geo.b0./geo.r0).^2 .*J4(a) ...
          .* (1./prof.lne + (J5(a)./J4(a)-1.5)./prof.lte + 1./prof.te.*neo.er) ...
          ./ nhu_el;      % Flux en particules /m^2 et /s;   facteur 16, car avant integration en theta

   flux_ion = zeros(1,gene.nbrho,gene.nbg);
   % boucle sur les especes ioniques de Cronos
   for i = 1:gene.nbg
      % calcul de la longueur de gradient de densite pour l'espece consideree :
      % -----------------------------------------------------------------------
      % gradient
      gns = rpdederive(gene.x,impur.impur(1,:,i),0,2,2,1) ./ equi.rhomax .* equi.grho;
      % bornes pour les longueurs de gradient
      % maximum
      lmax = 2.* equi.rmoy(1);
      % minimum
      lmin = equi.rhomax ./ length(gene.x);
      % longueur de gradient
      lns = lcalc(impur.impur(i),gns,lmin,lmax);

      nhu_s = (1./squeeze(Tauie(i,:))+squeeze(sum(1./squeeze(Tauii(i,:,:)))));
      a = 2.*delta .*Nbobine.* prof.q ./ equi.a .* abs(neo.er)./sqrt(prof.bphi.^2+prof.bpol.^2)./nhu_s;
      a(1) = a(2); % pour enlever le NaN du au fait que equi.a est nul en 0

      flux_ion(:,:,i) = impur.impur(:,:,i) .* 16 .*(delta./(2*pi)).^1.5 .* (prof.ti./compo.z(i)./geo.b0./geo.r0).^2 .*J4(a) ...
          .* (1./lns + (J5(a)./J4(a)-1.5)./prof.lti - compo.z(i)./prof.ti.*neo.er) ...
          ./  nhu_s;   % Flux en particules /m^2 et /s;   facteur 16, car avant integration en theta
   end
else
   disp('L''option modecalcul = 1 n''est pas active pour le moment');
end

% facteur d'integration sur l'angle poloidal (int(sin(theta)^2))
un_sur_alpha = equi.rmoy .* Nbobine .* prof.q .* delta ./ equi.a;
un_sur_alpha(1)=un_sur_alpha(2);
un_sur_alpha(find(un_sur_alpha>1))=1;  % si un_sur_alpha > 1, les particules sont piegees quelque soit theta
u = asin(un_sur_alpha);
integ = u/pi - sin(2*u)/2/pi;

% on applique l'integration sur l'angle poloidal aux flux
flux_el = flux_el .* integ;
for i = 1:5
   flux_ion(:,:,i) = squeeze(flux_ion(:,:,i)) .* integ;
end

% attention : les sources dues au ripple ont un signe : pour el,ne,ion : negatif pour des pertes
%             la rotation suit la convention cronos
%             les donnees de nrip sont negatives pour des pertes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sources de particules equivalentes; on convertit le flux radial en source volumique avec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% l'equation de conservation dn/dt = -div(flux)
source.ne = - 1./equi.vpr .*rpdederive(gene.x.*equi.rhomax,equi.vpr.*flux_el.*equi.grho,0,2,2,1); % dn/dt = -div(flux)
source.ne(1)=0; % point central a zero 
for i=1:5
   nrip(1,:,i) = - 1./equi.vpr .*rpdederive(gene.x.*equi.rhomax,equi.vpr.*squeeze(flux_ion(1,:,i)).*equi.grho,0,2,2,1); % dn/dt = -div(flux)
end
nrip(1,1,:)=0; % point central a zero
source.el = source.ne .* prof.te .* phys.e;
source.ion = source.ne .* prof.ti .* phys.e;   % pourquoi nrip n'intervient pas ?
source.w = gene.signe.ip .* abs(prof.psid1) .* equi.grho .* phys.e .* (-flux_el + sum(squeeze(flux_ion)'.*(compo.z'*ones(1,gene.nbrho))));

warning on

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
