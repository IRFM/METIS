% ZDDSCRASH calcul du crash des DDS
%--------------------------------------------------------------
% fichier zddscrash.m ->  zddscrash
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule la reorganisation du plasma lors
% du crash d'une dent de scie selon le modele proposer dans la reference :
% ref :  F. porcelli et all, Plasma Phys. Control Fusion, vol 38, p 2163-2186, 1996
%
% 
% syntaxe  :
%  
%       [cr,datak] =zddscrash(cons,datak,gene,compo,phys,valmem,equiname,equicons)
%
% entree :
%
%      cons            =    parametre propre a la fonction (param.cons.mhd.limite)
%      datak           =    structure de donnees data au temps d'interet.
%      gene            =    parametres generaux (param.gne)
%      compo           =    composition du plasma: charge et masse des atomes ( param.compo)
%      phys            =    constantes physiques (param.phys)
%      valmem          =    valeurs au temps d'avant de datak.gene (attention a la decoupe 
%                           des intervalle de temps pour assurer la convergences)  
%      equiname        =    nom du module d'equilibre
%      equicons        =    cons du module d'equilibre  
% 
% sortie :
% 
%     cr               =  compte rendu d'execution (0 = ok)
%     datak            =  structure data a un temps modifiee
% 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 03/09/2004.
% 
% 
% liste des modifications :
%
%  * 24/04/2002 -> securite sur q0 apres le crash 
%  * 02/05/2002 -> correction forme des profils apres reconnexion 
%  * 17/05/2002 -> correction pour DDS centrale
%  * 17/05/2002 -> blindages pour profils creux
%  * 17/05/2002 -> correction calcul E//
%  * 15/10/2003 -> calcul de psistar en torique
%  * 23/10/2003 -> ajout du calcul de Kadomtsev
%  * 05/01/2004 -> modification Porcelli epsilon-q et q central
%  * 07/01/2004 -> modification caclcul q0 cas Porcelli
%  * 22/01/2004 -> ajout du module Kadomtsev incomplet
%  * 03/09/2004 -> nouveau zcaljmoy
%  * 24/11/2004 -> petite modif ligne 264 (psupra < 1e-9 au lieu de 0), mise a jour appel zprofdivers ligne 418
%  * 13/12/2004 -> blindage calcul indmix en cas de profil de q non monotone
%--------------------------------------------------------------
%
function [cr,datak,memoire] =zddscrash(cons,datak,gene,compo,phys,valmem,equiname,equicons,equimemoire,pfcur)

%keyboard 

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 
	valeur.w1         = 0.5;   % fraction de la zone reconnectee (en fraction de rmix)
	type.w1           = 'float';                % type reel
	borne.w1          = [0.1,1];               % valeurs possible 
	defaut.w1         = 0.1;                      % valeurs par defaut
	info.w1           = 'largeur de la zone de shear nul en fraction de psistar(q=1) [0.1,1] {0.5}';
	
	valeur.epsq         = 1e-2;   % pente de q dans la zone q = 1
	type.epsq           = 'float';                % type reel
	borne.epsq          = [0,0.1];               % valeurs possible 
	defaut.epsq         = 1e-2;                      % valeurs par defaut
	info.epsq           = 'pente de q dans la zone q = 1';

	valeur.mode         = 1;   % mode 0= Porcelli , 1 = Kadomtsev, 2 = Kadomtsev incomplet
	type.mode           = 'integer';                % type entier
	borne.mode          = {0,1,2};                  % valeurs possible 
	defaut.mode         = 1;                        % valeurs par defaut
	info.mode           = 'mode 0= Porcelli , 1 = Kadomtsev, 2 = Kadomtsev incomplet';
	
	valeur.invariant         = 2;  
	type.invariant           = 'integer';                % type entier
	borne.invariant          = {0,1,2};                  % valeurs possible 
	defaut.invariant         = 2;                        % valeurs par defaut
	info.invariant           = 'si invariant = 0, la valeur au centre du plasma du flux poloidal est conservee lors de la reconnexion;\n si invariant = 1, la valeur au bord du plasma du flux poloidal  est conservee lors de la reconnexion;\n si invariant = 2, le fonctionnement est choisi conformement au CL de l''equation de diffusion';
	
	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'Calcul du crash des DDS ';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	cr = sortie;
	return
end

% compatibilite ascendante
if ~isfield(cons,'invariant')
 	cons.invariant = 0;
        disp('undefined parameter "invariant"; forced to 0 (magnetic axis poloidal flux preserve during reconnexion)');
end

% si invariant =  2
if cons.invariant  == 2
	if datak.mode.cons.psi == 0
		cons.invariant = 1;
	else
		cons.invariant = 1;
	end
end
% initialisation des sortie
cr      = 0;
% si on ne passe pas dans le calcul de l'equilibre
memoire = equimemoire;
% pour les plots
profmem = datak.prof;
equimem = datak.equi;
ipmem = datak.gene.ip;

% calcul du flux helicoidal :
x         = gene.x;
psistar=cumtrapz(datak.prof.psi,datak.prof.q-1.);

% Calcul de la position de q=1
[psistarmax,ind1]=max(psistar);
ind1    = max(ind1);
x1      = gene.x(ind1);
aa = find(psistar<0);
aa = aa(aa>ind1);   % to avoid taking too small mixing lentgh indmix, in case of non monotonic q-profiles
indmix = aa(1);  % F.I. 13/12/2004
%indmix  = min(find(psistar<0)); % old version 
xmix    = x(indmix);

% Valeurs par defaut
indout = indmix - 1;
indin =1;

if ind1==1
	disp('no DDS')
	% pas de dds
	return
elseif x1 < 0.1
	disp('no DDS')
   % pas de DDS (trop petit)
	return
end

% Porcelli model
if cons.mode == 0
   d         = abs(psistar -psistar(1));
   d(1:ind1) = inf;

   % calcul de la position de la zone de shear nul
   psicons    = (1-cons.w1) .* psistarmax;
   dg         = abs(psistar(1:ind1) -psicons);
   indin      = max(find(min(dg) == dg));
   xin        = x(indin);
   indm       = ind1:indmix;
   dd         = abs(psistar(indm) -psicons);
   indout     = min(indm(min(find(min(dd) == dd))),indmix-1);
   xout       = x(indout);
   % profil de q
   if indin > 1
     nq0  = xin .^ 4 ./ 4 ./ trapz(x(1:indin),x(1:indin) .^ 3 ,2);
     q0   = xin .^ 4 ./ 4 ./ trapz(x(1:indin),x(1:indin) .^ 3 ./ datak.prof.q(1:indin),2)./nq0;
   else
     q0   = 1. ;
   end  
   dq        = cons.epsq .* linspace(-0.5,0.5,length(indin:indout)); 
   
   if indin == 1
          xq        = cat(2,x(indin:indout),x(indmix:end));
          qq        = cat(2,ones(size(indin:indout)) + dq,datak.prof.q(indmix:end));
          qnew      = spline(xq,qq,x);
   elseif indin == 2
          xq        = cat(2,x(1),x(indin:indout),x(indmix:end));
          qq        = cat(2,q0,ones(size(indin:indout)) + dq,datak.prof.q(indmix:end));
          qnew      = spline(xq,qq,x);
   else                    
          xq        = cat(2,x(1:(indin-1)),x(indin:indout),x(indmix:end));
          qq        = cat(2,q0 .*ones(1,indin-1),ones(size(indin:indout)) + dq,datak.prof.q(indmix:end));
          qnew      = spline(xq,qq,x);
   end
elseif cons.mode == 1
   % Kadomtsev model
   xxint=linspace(0,x1,50);
   psistarint=interp1(x,psistar,xxint);
   xxext=interp1(psistar(ind1:end),x(ind1:end),psistarint);
   xxf=sqrt(xxext.^2-xxint.^2);

   psistarf = psistar;
   % il faudrait imposer gradient(psistarf,x)=0 au centre
   psistarf(1:indmix-1) = interp1(xxf,psistarint,x(1:indmix-1));
   warning off 
   qf     = 1. + pdederive(datak.prof.psi,psistarf,0,2,2,1);
   warning on
   qf(1) = 1;
   if any(~isfinite(qf))
      % cas psi non monotone
      qf_alt = 1. + pdederive(datak.equi.psi,psistarf,0,2,2,1);  
      indbad_ = find(~isfinite(qf));
      qf(indbad_) = qf_alt(indbad_);
   end  
   qnew   = datak.prof.q;
   qnew(1:indmix-1) = qf(1:indmix -1);
   
else
    % Incomplete Kadomtsev model
   % calcul de la position de la zone Kadomtsev
   psicons    = (1-cons.w1) .* psistarmax;
   dg         = abs(psistar(1:ind1) -psicons);
   indin      = max(find(min(dg) == dg));
   xin        = x(indin);
   indm       = ind1:indmix;
   dd         = abs(psistar(indm) -psicons);
   indout     = min(indm(min(find(min(dd) == dd))),indmix-1);
   xout       = x(indout);
   % Zone Taylor : q = q0
   if indin > 1
     nq0       = xin .^ 4 ./ 4 ./ trapz(x(1:indin),x(1:indin) .^ 3 ,2);
     q0        = xin .^ 4 ./ 4 ./ trapz(x(1:indin),x(1:indin) .^ 3 ./ datak.prof.q(1:indin),2)./nq0;
   else
     q0        = 1. ;
   end
   
   % Zone Kadomtsev
   xxint=linspace(xin,x1,100);
   psistarint=interp1(x,psistar,xxint);
   xxext=interp1(psistar(ind1:end),x(ind1:end),psistarint);
   xxf=sqrt(xin^2+xxext.^2-xxint.^2);
   
   psistarf=psistar;
   psistarf(indin+1:indout-1) = interp1(xxf,psistarint,x(indin+1:indout-1));
   psistarf(indin)            = psistarmax ;
   
   if indin > 2
     psistartaylor = cumtrapz(datak.prof.psi(1:indin),q0*ones(size(x(1:indin)))-1.);
     psistarf(1:indin-1) = cumtrapz(datak.prof.psi(1:indin-1),q0*ones(size(x(1:indin-1)))-1.);
     % Continuity between Taylor relaxed zone and Kadomtsev zone
     psistarf(1:indin-1) = psistartaylor(1:indin-1)+psistarf(indin)-psistartaylor(indin);
   end
   warning off
   qf     = 1. + pdederive(datak.prof.psi,psistarf,0,2,2,1) ;
   warning on
   qf(1) = 1;
   if any(~isfinite(qf))
      % cas psi non monotone
      qf_alt = 1. + pdederive(datak.equi.psi,psistarf,0,2,2,1);  
      indbad_ = find(~isfinite(qf));
      qf(indbad_) = qf_alt(indbad_);
   end  
   qnew   = datak.prof.q ;
   if indin > 1
     qnew(1:indin-1) = q0*ones(size(x(1:indin-1))) ;
   end
   qnew(indin)            = 1. ;
   qnew(indin+1:indmix-1) = qf(indin+1:indmix -1);
   
end

% boucle de calcul  (convergence avec l'equilibre)
equi      = datak.equi;
prof      = datak.prof;
% on suppose un possible changement de geometrie pour le premier appel
change = 1;

% la boucle n'ameliore pas la precision: 2 appels suffisent
for k =1:2

	 % calcul des nouvaux profils
	 % nouvelle valeur de psi
	 inte      = - equi.F .* equi.vpr .* equi.r2i ./ (4*pi^2) ./ qnew;
	 psinew    = equi.rhomax .* cumtrapz(gene.x,inte,2);
	 % raccordement a rmix :
	 psinew(1:indmix)       = psinew(1:indmix) - psinew(indmix) + datak.prof.psi(indmix);
	 psinew((indmix+1):end) = datak.prof.psi((indmix+1):end);

	 % si le flux au centre est conserve
         if cons.invariant == 0
			psinew = psinew - psinew(1) + datak.prof.psi(1);
	 end

	 % nouveau profil de densite
	 dnedx                  = pdederive(gene.x,datak.prof.ne,0,2,2,1);
	 dnedx(1:indout)        = 0;
	 nenew                  = cumtrapz(gene.x,dnedx,2);
	 nenew                  = nenew - nenew(end) + datak.prof.ne(end);
	 % normalisation (conservation de la matiere)
	 ii                     = 1:indmix;
	 netot                  = equi.rhomax .* trapz(gene.x(ii),equi.vpr(ii) .* (datak.prof.ne(ii) - datak.prof.ne(indmix)),2);
	 netotnew               = equi.rhomax .* trapz(gene.x(ii),equi.vpr(ii) .* (nenew(ii) - nenew(indmix)),2);
	 nenew(ii)              = (nenew(ii) - nenew(indmix)) ./ max(1,netotnew) .* netot + nenew(indmix);
	 ninew                  = nenew .* datak.prof.ae;

	 % nouvelle valeur de la pression  supra
	 psupra                 = datak.prof.ptot - datak.prof.pe - datak.prof.pion;
	 if all(psupra <= 1e-9)      %Modif F.I. 24/11/04, le psupra n'est pas rigoureusement nul, mais negligeable
		psupranew              = zeros(size(psupra));
	 else
		porigine               = psupra;
		dpdx                   = rpdederive(gene.x,porigine,0,2,2,1);
		dpdx                   = dpdx .* (dpdx < 0);
		porigine               = cumtrapz(gene.x,dpdx,2);
		porigine               = porigine - porigine(end) + psupra(end);
		% valeur de dpdx minimal (non nulle)
		indm                   = find(dpdx < 0);
		dpdxok                 = dpdx(indm);
		%dpdxmin                = dpdxok(min(find(dpdx(indm) == max(dpdx(indm)))));
		dpdxmin                = 0;
		dpdx(1:indout)         = dpdxmin/3;
		pnew                   = cumtrapz(gene.x,dpdx,2);
		pnew                   = pnew - pnew(end) + porigine(end);
		% normalisation (conservation de l'energie)
		ii                     = 1:indmix;
		wtot                   = equi.rhomax .* trapz(gene.x(ii),equi.vpr(ii) .* (porigine(ii) - porigine(indmix)),2);
		wtotnew                = equi.rhomax .* trapz(gene.x(ii),equi.vpr(ii) .* (pnew(ii) - pnew(indmix)),2);
		rapp                   = max(1,wtot ./ max(1,wtotnew) );
		pnew(ii)               = (pnew(ii) - pnew(indmix)) .* rapp + pnew(indmix);
		psupranew              = pnew;
	 end
	 % nouvelle valeur de la pression eletronique
	 %porigine               = zmonotone(gene.x,datak.prof.pe);
	 porigine               = datak.prof.pe;
	 dpdx                   = rpdederive(gene.x,porigine,0,2,2,1);
	 dpdx                   = dpdx .* (dpdx < 0);
	 porigine               = cumtrapz(gene.x,dpdx,2);
	 porigine               = porigine - porigine(end) + datak.prof.pe(end);
	 % valeur de dpdx minimal (non nulle)
	 indm                   = find(dpdx < 0);
	 dpdxok                 = dpdx(indm);
	 %dpdxmin                = dpdxok(min(find(dpdx(indm) == max(dpdx(indm)))));
	 dpdxmin                = 0;
	 dpdx(1:indout)         = dpdxmin/3;
	 pnew                   = cumtrapz(gene.x,dpdx,2);
	 pnew                   = pnew - pnew(end) + porigine(end);
	 % normalisation (conservation de l'energie)
	 ii                     = 1:indmix;
	 wtot                   = equi.rhomax .* trapz(gene.x(ii),equi.vpr(ii) .* (porigine(ii) - porigine(indmix)),2);
	 wtotnew                = equi.rhomax .* trapz(gene.x(ii),equi.vpr(ii) .* (pnew(ii) - pnew(indmix)),2);
	 rapp                   = max(1,wtot ./ max(1,wtotnew) );
	 pnew(ii)               = (pnew(ii) - pnew(indmix)) .* rapp + pnew(indmix);
	 penew                  = pnew;
	 
	 % nouvelle valeur de la pression ionique
	 %porigine               = zmonotone(gene.x,datak.prof.pion);
	 porigine               = datak.prof.pion;
	 dpdx                   = rpdederive(gene.x,porigine,0,2,2,1);
	 dpdx                   = dpdx .* (dpdx < 0);
	 porigine               = cumtrapz(gene.x,dpdx,2);
	 porigine               = porigine - porigine(end) + datak.prof.pion(end);
	 % valeur de dpdx minimal (non nulle)
	 indm                   = find(dpdx < 0);
	 dpdxok                 = dpdx(indm);
	 %dpdxmin                = dpdxok(min(find(dpdx(indm) == max(dpdx(indm)))));
	 dpdxmin                = 0;
	 dpdx(1:indout)         = dpdxmin/3;
	 pnew                   = cumtrapz(gene.x,dpdx,2);
	 pnew                   = pnew - pnew(end) + porigine(end);
	 % normalisation (conservation de l'energie)
	 ii                     = 1:indmix;
	 wtot                   = equi.rhomax .* trapz(gene.x(ii),equi.vpr(ii) .* (porigine(ii) - porigine(indmix)),2);
	 wtotnew                = equi.rhomax .* trapz(gene.x(ii),equi.vpr(ii) .* (pnew(ii) - pnew(indmix)),2);
	 rapp                   = max(1,wtot ./ max(1,wtotnew));
	 pnew(ii)               = (pnew(ii) - pnew(indmix)) .* rapp + pnew(indmix);
	 pionnew                = pnew;
	 
	 % nouvelle structure de profil
	 prof.ne      = nenew;
	 prof.pe      = penew;
	 prof.pion    = pionnew;
	 prof.ptot    = psupranew + penew + pionnew;
	 prof.psi     = psinew;
	 prof.q       = qnew;
	 %prof.psid1   = pdederive(gene.x,prof.psi,0,2,2,1);
         % la valeur de psid1_bord est conservee lors de la reconnexion 
         sigma = [length(gene.x) ones(size(datak.prof.psi(2:end)))];
         %%%[yout,prof.psid1,prof.psid2] = interpos(gene.x,prof.psi,-1,[1 1],[0 datak.prof.psid1(end)],sigma);
         [prof.psid1,prof.psid2] = z0polyderive(gene.x,prof.psi,7,gene.x);
	 % psid1_1 est conserve lors de la reconnexion
	 ip        = - 1 ./ ( 2 * pi * phys.mu0 .* equi.rhomax) .* equi.c2c(end) .*  datak.prof.psid1(end); 
       	 prof.jmoy = zcaljmoy(gene.x,prof.psi,equi.vpr,equi.grho2r2,equi.rhomax,equi.ri,phys.mu0,gene.creux, ...
	                                             equi.rhoRZ,equi.c2c,prof.ptot,equi.F,equi.r2i,ip);

  	 % calcul du nouvel equilibre ( a ip donne)
         [equi,mhd_cd,memoire]  =  feval(equiname,equimemoire,equicons,pfcur,datak.geo,equi,prof,phys,ip,gene.x,change);
	 if equi.fail~= 0
	     fail = equi.fail
	     equi = datak.equi;
	     equi.fail  = fail;
	     zverbose('Probleme de convergence de l''equilibre avec les DDS');
             break
         end

         % flux libre au bord si equilibre a frontiere libre
	 if datak.geo.mode == 3
		prof.psi = prof.psi - prof.psi(end) + equi.psi(end);
	 end
     
end
equi.conv  = k;

% duree du crash
rc      = datak.equi.rhomax .* gene.x(indmix); 
taur    = phys.mu0 .* rc .^ 2 ./ mean(datak.coef.eta(1:indmix)); 
va      = datak.geo.b0 ./ sqrt(phys.mu0 .* phys.me .* datak.prof.ne(1));
taua    = rc ./ va;
dt      = sqrt(taur .* taua);

% calcul des derivees temporelles de l'equilibre
equi.drhomaxdt      = 0.5 .* (equi.rhomax - datak.equi.rhomax) ./ dt + 0.5 .* datak.equi.drhomaxdt;
equi.dvprdt         = 0.5 .* (equi.vpr - datak.equi.vpr) ./ dt + 0.5 .* datak.equi.dvprdt;
equi.dsprdt         = 0.5 .* (equi.spr - datak.equi.spr) ./ dt+ 0.5 .* datak.equi.dsprdt;
equi.dphidt         = 0.5 .* (equi.phi - datak.equi.phi) ./ dt+ 0.5 .* datak.equi.dphidt;
% calcul des derivees 
equi.phid1     = pdederive(gene.x,equi.phi,0,2,2,1);
equi.phid2     = pdederive(gene.x,equi.phi,1,2,2,2);

% calcul des derivees temporelles 
% cette formule est fausse en presence d'un ilot qui bouge
%prof.dpsidt  		 = (prof.psi - datak.prof.psi) ./ dt;
prof.dnedt			 = (prof.ne - datak.prof.ne) ./ dt;
prof.dpedt			 = (prof.pe - datak.prof.pe) ./ dt;
prof.dpiondt 		 = (prof.pion - datak.prof.pion) ./ dt;
% pas de changement dans la composition
%prof.daedt			 = (prof.ae - datak.prof.ae) ./ dt;
prof.dflucedt		 = (prof.fluce - datak.prof.fluce) ./ dt;
prof.dfluciondt 	 = (prof.flucion - datak.prof.flucion) ./ dt;
prof.drotdt  		 = (prof.rot - datak.prof.rot) ./ dt;

% calcul des derivees temporelles 3 points  
% dydt = zdydt3p(ynew,yold,dydtold,dt)
% cette formule est fausse en presence d'un ilot qui bouge
%prof.dpsidt3p			= zdydt3p(prof.psi,datak.prof.psi,datak.prof.dpsidt3p,dt);
prof.dnedt3p 			= zdydt3p(prof.ne,datak.prof.ne,datak.prof.dnedt3p,dt);
prof.dpedt3p 			= zdydt3p(prof.pe,datak.prof.pe,datak.prof.dpedt3p,dt);
prof.dpiondt3p  		= zdydt3p(prof.pion,datak.prof.pion,datak.prof.dpiondt3p,dt);
% pas de changement dans la composition
%prof.daedt3p 			= zdydt3p(prof.ae,datak.prof.ae,datak.prof.daedt3p,dt);
prof.dflucedt3p 		= zdydt3p(prof.fluce,datak.prof.fluce,datak.prof.dflucedt3p,dt);
prof.dfluciondt3p  	= zdydt3p(prof.flucion,datak.prof.flucion,datak.prof.dfluciondt3p,dt);
prof.drotdt3p			= zdydt3p(prof.rot,datak.prof.rot,datak.prof.drotdt3p,dt);

% memorisation des donnees
datak.prof = prof;
datak.equi = equi;

% nouvelle base temps avec securite pour la base temps
dt      = min(dt,gene.dt/2);
datak.gene.temps = datak.gene.temps + dt;
datak.gene.dt    = dt;

% calcul des divers profiles utiles (q,j,E//,Bpol ...)
datak = zprofequi(phys,gene,datak);

% calcul des divers profiles utiles (Te,Ti,  ...)
datak = zprofdivers(phys,gene,datak,compo);

% calcul des donnees generiques du plasma
datak = zcalcdivers(phys,gene,compo,datak,valmem,dt);

return


h = findobj(0,'type','figure','tag','ddscrash');
if isempty(h)
   h =figure('tag','ddscrash');
else
   figure(h);
end

clf
subplot(3,1,1)
plot(gene.x,datak.prof.pe,'r',gene.x,datak.prof.pion,'b',gene.x,datak.prof.ptot,'k', ...
gene.x,profmem.pe,'or',gene.x,profmem.pion,'bo',gene.x,profmem.ptot,'ko');
hold
plot([gene.x(ind1) gene.x(ind1)],[0,max(datak.prof.ptot)],'g');
plot([gene.x(indout) gene.x(indout)],[0,max(datak.prof.ptot)],'-.g');
plot([gene.x(indin) gene.x(indin)],[0,max(datak.prof.ptot)],'-.g');
plot([gene.x(indmix) gene.x(indmix)],[0,max(datak.prof.ptot)],':g');

title(num2str(datak.gene.temps));
subplot(3,1,2)
plot(gene.x,datak.prof.jmoy,'r',gene.x,profmem.jmoy,'or');
hold
plot([gene.x(ind1) gene.x(ind1)],[0,max(datak.prof.jmoy)],'g');
plot([gene.x(indout) gene.x(indout)],[0,max(datak.prof.jmoy)],'-.g');
plot([gene.x(indin) gene.x(indin)],[0,max(datak.prof.jmoy)],'-.g');
plot([gene.x(indmix) gene.x(indmix)],[0,max(datak.prof.jmoy)],':g');
title(sprintf('ip before = %g & ip after = %g',ipmem,datak.equi.ip));

subplot(3,1,3)
plot(gene.x,datak.prof.q,'r',gene.x,profmem.q,'or', ...
    gene.x,datak.prof.shear,'b',gene.x,profmem.shear,'bo');
hold
plot([gene.x(ind1) gene.x(ind1)],[min(datak.prof.shear),max(datak.prof.q)],'g');
plot([gene.x(indout) gene.x(indout)],[min(datak.prof.shear),max(datak.prof.q)],'-.g');
plot([gene.x(indin) gene.x(indin)],[min(datak.prof.shear),max(datak.prof.q)],'-.g');
plot([gene.x(indmix) gene.x(indmix)],[min(datak.prof.shear),max(datak.prof.q)],':g');
xlabel('x');
ylabel('q & s')
title(sprintf('li before = %g & li after = %g',equimem.li,equi.li));
pause(1)

%keyboard

