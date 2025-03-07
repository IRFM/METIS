% calcul de la puissance de fusion
function [pfus,salpha,ifus,xfus,jxfus,j0fus,taus,ecrit,pfus_nbi,pfus_loss,jfus,salf,palf,splustd,splusdt,splusff,splusicrh] = ...
	 zfus0tae(nDi,nTi,te,ne,zeff,tite,R,a,K,Bt,ane,ate,Vp,Sp, ...
         pnbi_th,taus_nbi,ecritnbi,e0nbi,e0nbi2,ftnbi,picrh_ion,taus_icrh,ecriticrh,e0icrh,temps,pnbi,d0,qa,qmin,te0, ...
         nebord,tebord,pped,ni,wth,taeonoff,nb_nbi,fspot,e_shielding,profli,fpolarized,forced_H_NBI)

% pour l'interaction faisceau plasma
tii = te .* tite;

% integartion su la section efficace avec effet de profils
if isfield(profli,'tep')
	x = profli.xli;
	ux  = 1 - x .^ 2;
	ve  = ones(size(x));
	vt  = ones(size(te));
	spr = profli.spr;
	vpr = profli.vpr;
	nep = profli.nep;
	tep = profli.tep;
	tip = profli.tip;
	nip = profli.nip;
	nD  = profli.n1p .* ((nDi./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
	nT  = profli.n1p .* ((nTi./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
	Raxe       = profli.Raxe;
	epsi       = profli.epsi;
	zeffp      = profli.zeff;
else
	x   = linspace(0,1,21);
	ux  = 1 - x .^ 2;
	ve  = ones(size(x));
	vt  = ones(size(te));
	% pour tenir compte des profils reels mode L/H
	spr   = (2 .* Sp) * x;
	vpr   = (2 .* Vp) * x;
	nep = ((ne .* (1 + ane)-nebord) * ve)  .* (vt * ux)  .^ (ane * ve) + nebord * ve;
	tep = ((te .* (1 + ate)-tebord) * ve)  .* (vt * ux)  .^ (ate *ve) + tebord * ve;
	nip = nep .* ((ni ./ ne) * ve);
	tip = tep .* (tite * ve);
	nD  = ((max(1,nDi) ./ ne) * ve) .* nep;
	nT  = ((max(1,nTi) ./ ne) * ve) .* nep;
	zeffp      = zeff * ve;
	Raxe       = R * ve + d0 * ux;
	epsi       = max(eps,(a * x)) ./ Raxe;
	% modeh /l
%  	ee           =   1.602176462e-19;
%  	teped            = pped ./ (nep(:,end-1) + nip(:,end-1).* tite ) ./ ee;
%  	ptotd            = ee .*(nep + nip .* (tite *ve)) .* tep;
%  	pp               = ptotd;
%  	pp(:,1:(end-1))  = pped * ve(1:(end-1));
%  	wpied            = (3/2) .* trapz(x,vpr .* pp,2);
%  	pc               = ptotd;
%  	pc               = pc - pc(:,end-1) * ve;
%  	pc(pc<0)         = 0;
%  	wcore            = (3/2) .* trapz(x,vpr .* pc,2);
%  	fact             = max(1,wcore) ./ max(1,wth - wpied);
%  	modehp            = (pped>1) .* ((wth - wpied) > 0);
%  	te21h            = max(tebord *ve,teped * ve  + (tep - tep(:,end -1) * ve) ./ (fact *ve));
%  	te21h(:,end)     = tebord;
%  	tep              = (modehp*ve) .* te21h +  ((~modehp)*ve) .* tep;
%  	tip              = tep .* (tite * ve);
end
if isfield(profli,'ftrap')
	%factrap = 1 - min(1,trapz(xp,profli.ftrap .* jdep .* abs(vt*xp),2) ./ max(1,trapz(xp,jdep .* abs(vt*xp),2)));
	ftrap = profli.ftrap;
else
	ftrap = 0.95 .* sqrt(vt*xp);
	%factrap = 0.5;
end

% source de particule alpha
%sigmavdt = zperes(max(tip,30));
%s        = nD .* nT .* sigmavdt;
[dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion(max(tip,30));
s        = nD .* nT .* dt.sv .* (1 + fpolarized / 2);


% calcul de l'elargissement du a la largeur d'orbite
% article L.G. Eriksson & Porcelli
% il y  a peu d'elargissement car les praticule passe leur temps pres des pointes de bananes.

% temps de thermalisation (plutot central, calibree avec spot)
taus  = 6.27e8 .* (te .* (1 + ate ./ 2)).^ (3/2) ./ (ne./1e6) ./ 17;
ecrit = 14.8 .* te .* (1 + ate ./ 2) .* (4 .* zeff) .^ (2/3);

% calcul de l'elargissement pour les TAE
% boucle de propagation
smem     = s;
%dsloss   = 0 .* vt;
if taeonoff == 1

  % boucle pour etalement du profil
  nbcount = 43;
  indinstable =vt;
  while (nbcount >0) & (any(indinstable))
      % calcul de l'effet des TAE sur la largeur de la source de alphas
      VA       = 2.18e16 .* (Bt * ve) ./ sqrt(2 .* nD + 3 .* nT);                 % vitesse de Alfven
      pres_alf = 3.56e6 .* 1.602176462e-19 .* s  ./ 2 .*  (taus * ve);      % pression des alpha en premiere approximation
      taun     = taus ./ 3 .* log(1 + (3.56e6 ./ ecrit) .^(3/2));           % temps d'arret des alpha
      n_alf    = s .* (taun * ve);                                          % nombre d'alpha dans la decharge nonthermalise
      emoy     = max(1e5,min(3.56e6,pres_alf ./ max(1,n_alf)./  1.602176462e-19));  % energy moyenne en eV pour F(x)
      v_alf    = sqrt(2 .* emoy .* 1.602176462e-19 ./ 4 ./ 1.66053873e-27); % vitesse des alpha a l'energie moyenne
      v_el     = sqrt(2 .* tep .* 1.602176462e-19 ./9.10938188e-31);       % vitesse des electrons
      % formulation Wesson p 404
      % parametre chosi selon ref : G.Y Fu and J.W Van Dam , phys. fluids B 1 (10) 1989.
      xpar     = VA ./ v_alf;
      fx       = xpar .* (1 + 2 .* xpar .^ 2 + 2 .* xpar .^ 4) .* exp(- xpar .^ 2);
      if isfield(profli,'qjli')
       	qp       = profli.qjli;
      else
      	qp       = z0qp(x,max(1,qmin),qa);
      end
      w0       = VA ./ 2 ./ qp ./ (R * ve);
      if isfield(profli,'rmx')
		rm  = profli.rmx(:,end);
      else
      		rm  = a .* (1+1./4.*(-1+K)-3./64.*(-1+K).^2+5./256.*(-1+K).^3-175./16384.*(-1+K).^4);
      end
      % pour le mode m = fix(q);
      wastar   = max(1,fix(qp)) .*  ((1 ./ 2 ./ Bt ./ rm .^ 2) * ve) .* emoy ./ max(vt * x,eps);
      wastar(:,1) = wastar(:,2);
      % la frequence du mode qui  est plus faible que la frequence propre du mode est prise pessimiste pour la stabilite
      % ref : M.N Rosenbluth et al , phys. fluids B 4 (7) 1992
      % ref : Vlad + MS
      cstar    = - max(wastar ./ w0,eps);
      cstar(:,1) = cstar(:,2);
      dstar    = 1.602176462e-19 .* nep .* tep ./ max(0.1,fx) .* VA ./ v_el;
      % les integrales sont compliquees a calculer ...
      %gpmax    = intgpmax(x,cstar,dstar);          % clacul du gradient qui declenche les TAEs (gamma == 0)
      % approximation petite correction a pression nulle
      gpmax    = (dstar + 0.5 .* pres_alf) ./ cstar;
      gpal     = pdederive(x,pres_alf,0,2,2,1);
      % reduction du gradient calculer a partir du point de forcage maximum
      forcage  = qp .^ 2 .* gpal;
      [mforce,indf] = min(forcage,[],2);
      xforce        = x(indf)';
      gpalforce     = gpal(:,indf);
      gpmaxforce    = gpmax(:,indf);
      indinstable   = (gpalforce(:,1) < gpmaxforce(:,1));
      % etalement de la source pour les temps instable
      s(indinstable,2:end) = (s(indinstable,2:end)  + s(indinstable,1:(end-1))) ./ 2;
      % conservation du nombre de particule
      s        = (s>1) .* s  .* ((trapz(x,smem .* vpr,2) ./ max(1,trapz(x,s .* vpr,2))) * ve);
      s        = s .* (s < (10 .* smem));

      if any(imag(s(:))) | any(s(:) <0) | any(~isfinite(s(:)))
          disp('pb in s')
          figure(45);clf
          subplot(2,2,1);
          plot(x,gpmax,'r',x,gpal,'b')
          subplot(2,2,2)
          plot(x,forcage,xforce,mforce,'o');
          subplot(2,2,3)
          plot(x,s,'r',x,smem,'bo');
          subplot(2,2,4)
          plot(x,s-smem,'g');
          drawnow;
      end
      % fin de la boucle
      nbcount = nbcount - 1;
  end
end

% valeur scalaire
salf     = s;
profli.alphashape = s;
palf     = 3.56e6 .* 1.602176462e-19 .* s;
salpha   = trapz(x,vpr .* s,2);
pfus     = 3.56e6 .* 1.602176462e-19 .* salpha;

%  % fspot  est un facteur de calibration calculer avec spot
%  % fspot   = 0.15;
%  % calcul du courant genere par les alpha (ref L-G Eriksson and F. Porcelli, Plasma Phys. Controlled Fus. 43 (2001) R145-)
%  % formule p R175, 9.4
%  valpha  = sqrt(2 .*  3.56e6 .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
%  wca     = 2 .* 1.602176462e-19 .* Bt ./ 1.66053873e-27 ./ 4;
%  dp0R    = (2 .* (max(1,qmin)*ve) .* valpha ./ (wca * ve) ./ (R*ve)) .^ (2/3);
%  % calcul de l'effet du courant de retour
%  %GZ      = (1.55 + 0.85./(zeff*ve)) .* sqrt((a./R) * x) - (0.2 + 1.55./(zeff*ve)) .* ((a./R) * x);
%  %jfus    = -fspot .* 2 .* 1.602176462e-19 .* valpha .* dp0R .* pdederive(x,s,2,0,2,1) .* (1 - 2 .* %(1-GZ) ./ (zeff*ve));
%  jfus    = -fspot .* 2 .* 1.602176462e-19 .* valpha .* dp0R .* pdederive(x,s,2,0,2,1);
%  xt = ftrap ./ (1 - ftrap);
%  D  = 1.414 .* zeffp + zeffp .^ 2 + xt .* (0.754 + 2.657 .* zeffp + 2 .* zeffp .^ 2) + ...
%       xt .^2 .* ( 0.348 + 1.243 .* zeffp + zeffp .^ 2);
%  GZ  = xt .* ((0.754 + 2.21 .* zeffp + zeffp .^2) + xt .* (0.348 + 1.243 .* zeffp + zeffp .^2)) ./ D;
%  jfus    =(1 - (1 - GZ) .* 2 ./ zeffp) .* jfus;
%  mask    = (jfus == (max(jfus,[],2)*ve));
%  xfus    = sum((vt*x) .* mask,2) ./max(1,sum(mask,2));
%  jxfus   = sum(jfus .* mask,2) ./max(1,sum(mask,2));
%  j0fus   = jfus(:,2);
%  spr     = (2.* Sp) * x;
%  ifus    = trapz(x,spr .* jfus,2);
%  profli.jfusshape = jfus;
%  % attention ifus doit etre multiplie par le temps de ralentissement pour obtenir un courant en A

% calcul de la puissance de perte 1ere orbite
% il faut ajouter le pertes ripples et de diffusion angulaire
rloc   = R * ve + a *x;
ral    = 4.55e-3.* sqrt(3.56e6 ./ 1e3) ./ (((Bt .* R)*ve) ./ rloc); % en m
%qp     = qmin * ve + (qa -qmin) * (x .^ 2);
if isfield(profli,'qjli')
       	qp       = profli.qjli;
else
	qp     = z0qp(x,max(1,qmin),qa);
end
% patato
dp1     = (R*ve) .* (2 .* qp .* ral ./ (R *ve)) .^ (2/3);
% banana
dp2     = sqrt((a * x) ./ (R * ve +  d0 * ux)) .* ral .* qp;
dp      = dp2 .* (dp2 < (a * x)) + dp1 .* (dp2 >= (a * x));
%figure(20);clf;plot(x,dp1,'r',x,dp2,'b',x,dp,'k',x,ral,'g');drawnow
% 0.05 = evite les artefacts petite machine
mask   = (ral + dp + d0 * ux + a *x - (a - 0.05) *ve) >0;
%sloss  = trapz(x,vpr .* s .* mask,2) + dsloss;
sloss  = trapz(x,vpr .* s .* mask,2);
frloss = max(0,min(1,sloss ./ (salpha + eps)));  % fraction d'alpha perdue

% fspot  est un facteur de calibration calculer avec spot
% fspot   = 0.15;
% calcul du courant genere par les alpha (ref L-G Eriksson and F. Porcelli, Plasma Phys. Controlled Fus. 43 (2001) R145-)
% formule p R175, 9.4
if fspot > 0
	dp0R    = (2 .* pi) .* dp ./ rloc;
	valpha  = sqrt(2 .*  3.56e6 .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
else
  	valpha  = sqrt(2 .*  3.56e6 .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
  	wca     = 2 .* 1.602176462e-19 .* Bt ./ 1.66053873e-27 ./ 4;
  	dp0R    = (2 .* (max(1,qmin)*ve) .* valpha ./ (wca * ve) ./ (R*ve)) .^ (2/3);
end
% calcul de l'effet du courant de retour
%GZ      = (1.55 + 0.85./(zeff*ve)) .* sqrt((a./R) * x) - (0.2 + 1.55./(zeff*ve)) .* ((a./R) * x);
%jfus    = -fspot .* 2 .* 1.602176462e-19 .* valpha .* dp0R .* pdederive(x,s,2,0,2,1) .* (1 - 2 .* %(1-GZ) ./ (zeff*ve));
if fspot > 0
	jfus    = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha .* dp0R .* pdederive(x,s,0,0,2,1);
        lnldei  = 15.2 - 0.5 .* log(nep./1e20) + log(tep ./1e3);
        taus_mat = 6.27e8 .* 4 .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei;
        taus_ave = trapz(x,vpr .* taus_mat,2) ./ trapz(x,vpr,2);
	jfus    = jfus .* (taus_mat ./ max(eps,taus_ave * ve));
else
	jfus    = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha .* dp0R .* pdederive(x,s,2,0,2,1);
end

switch e_shielding
    case 'Honda-NEO'
        % ref : M. Honda et al, Nucl. Fus. 52 (2012) p 023021
        % ref A. Redl et al, Phys. Plasmas 28, 022502 (2021); https://doi.org/10.1063/5.0012664
        %GZ = z0sauterL31(xp,tep,nep,qpr,zeffp,Raxe,ftrap,epsi);
        [~,~,GZ] = z0etaboot_neofit(x,tep,tep,nep,nep,qp,zeffp,1,Raxe,ftrap,epsi);   % only for DT plasma, so zion = 1.
    case 'Honda-Sauter'
        % ref : M. Honda et al, Nucl. Fus. 52 (2012) p 023021
        GZ = z0sauterL31(x,tep,nep,qp,zeffp,Raxe,ftrap,epsi);
    otherwise
        
        xt = ftrap ./ (1 - ftrap);
        D  = 1.414 .* zeffp + zeffp .^ 2 + xt .* (0.754 + 2.657 .* zeffp + 2 .* zeffp .^ 2) + ...
            xt .^2 .* ( 0.348 + 1.243 .* zeffp + zeffp .^ 2);
        GZ  = xt .* ((0.754 + 2.21 .* zeffp + zeffp .^2) + xt .* (0.348 + 1.243 .* zeffp + zeffp .^2)) ./ D;
end
	jfus    =(1 - (1 - GZ) .* 2 ./ zeffp) .* jfus;
if fspot > 0
	mask    = (jfus == (max(jfus,[],2)*ve)) .* mask;
else
	mask    = (jfus == (max(jfus,[],2)*ve));
end
xfus    = sum((vt*x) .* mask,2) ./max(1,sum(mask,2));
jxfus   = sum(jfus .* mask,2) ./max(1,sum(mask,2));
j0fus   = jfus(:,2);
spr     = (2.* Sp) * x;
ifus    = trapz(x,spr .* jfus,2);
profli.jfusshape = jfus;
% attention ifus doit etre multiplie par le temps de ralentissement pour obtenir un courant en A



% calcul de l'interraction faisceau-plasma (DT)
if (all(real(ftnbi == 0)) && all(imag(ftnbi == 0)))  || all( pnbi  < 1e3) || (forced_H_NBI ~= 0)
    esupranbi = zeros(size(temps));
    pTth      = zeros(size(temps));
    tauseff   = 1;
    splustd   = zeros(size(temps));
    emean     = e0nbi;
else
    %[esupranbi,pTth,tauseff] = zsupra0(temps,pnbi .* max(1e-6,ftnbi),taus_nbi,1e-6,ecritnbi,e0nbi,3);
    %[splustd,emean]   = zbpi0(nDi,tii,pTth,taus_nbi,e0nbi,ecritnbi,1);
    [splustd,emean]   = zbpi0(nDi,tii,real(pnbi_th) .* max(1e-6,real(ftnbi)),real(taus_nbi),e0nbi,real(ecritnbi),1);
    if nb_nbi > 1
        [splustd2,emean2]   = zbpi0(nDi,tii,imag(pnbi_th) .* max(1e-6,imag(ftnbi)),imag(taus_nbi),e0nbi2,imag(ecritnbi),1);
        splustd = splustd + splustd2;
        %emean   = (real(pnbi) .*emean   + imag(pnbi) .* emean2) ./ max(1,real(pnbi) + imag(pnbi));
        emean2(emean2 <=0) = 1;
        emean(emean <=0) = 1;
        emean   = emean   + sqrt(-1) .* emean2;
    end
end
% attention la matrice est T-> D, IL FAUDRA RAJOUTER 	D -> T
% maintenant c'est la bonne matrice
if	(all(real(ftnbi == 1)) && all(imag(ftnbi == 1)))  || all( pnbi  < 1e3) || (forced_H_NBI ~= 0)
    splusdt   = zeros(size(temps));
else
    splusdt   = zbpi0_dt(nTi,tii,real(pnbi_th),real(taus_nbi),e0nbi,real(ecritnbi),1-real(ftnbi));
    if nb_nbi > 1
        splusdt2   = zbpi0_dt(nTi,tii,imag(pnbi_th),imag(taus_nbi),e0nbi2,imag(ecritnbi),1-imag(ftnbi));
        splusdt    = splusdt + splusdt2;
    end
end
% approximation faisceau-faisceau avec la meme formule (c'est faux) cas ftnbi intermediaire (DT)
emean = max(max(30,te),real(emean)) + sqrt(-1) .* max(max(30,te),imag(emean));
if all(ftnbi == 0) || all( pnbi  < 1e3) || (forced_H_NBI ~= 0)
   splusff   = zeros(size(temps));
else
   [esupranbi,pTth,tauseff]   = zsupra0(temps,real(pnbi) .* max(1e-6,real(ftnbi)),real(taus_nbi),1e-6 .* ones(size(taus_nbi)),real(ecritnbi),e0nbi,3);
   esupranbi(~isfinite(esupranbi))= 0;
   splusff   = zbpi0((1-real(ftnbi)) .* real(pnbi_th) .* tauseff ./ 2 ./ Vp./ (1.602176462e-19 .* real(emean)), ...
               real(emean), real(pnbi) .* max(1e-6,real(ftnbi)),real(taus_nbi),e0nbi,real(ecritnbi), 1);
   if nb_nbi >1 
  	[esupranbi2,pTth2,tauseff2]   = zsupra0(temps,imag(pnbi) .* max(1e-6,imag(ftnbi)),imag(taus_nbi),1e-6 .* ones(size(taus_nbi)),imag(ecritnbi),e0nbi2,3);
   	esupranbi(~isfinite(esupranbi))= 0;
      	% 2 ->2
    	splusff22   = zbpi0((1-imag(ftnbi)) .* imag(pnbi_th) .* tauseff2 ./ 2 ./ Vp./ (1.602176462e-19 .* imag(emean)), ...
               		imag(emean), imag(pnbi) .* max(1e-6,imag(ftnbi)),imag(taus_nbi),e0nbi2,imag(ecritnbi), 1);

      	% 1 -> 2
    	splusff12   = zbpi0((1-real(ftnbi)) .* real(pnbi_th) .* tauseff ./ 2 ./ Vp./ (1.602176462e-19 .* real(emean)), ...
               		imag(emean), imag(pnbi) .* max(1e-6,imag(ftnbi)),imag(taus_nbi),e0nbi2,imag(ecritnbi), 1);

      	% 2 -> 1
        splusff21   = zbpi0((1-imag(ftnbi)) .* imag(pnbi_th) .* tauseff2 ./ 2 ./ Vp./ (1.602176462e-19 .* imag(emean)), ...
                      real(emean), real(pnbi) .* max(1e-6,real(ftnbi)),real(taus_nbi),e0nbi,real(ecritnbi), 1);

        splusff = splusff + splusff22 + splusff12 + splusff21;
   end
end
%  % acceleration de T par icrh (approximation type idn)
e0icrhm =zpmean(1:length(e0icrh),picrh_ion,e0icrh);
if all(picrh_ion < 1e3)
  	splusicrh = zeros(size(nDi));
elseif ~isempty(e0icrhm)
	if e0icrhm == 0
  		  	splusicrh = zeros(size(nDi));
 	else
  		rap 	  = zpmean(1:length(e0icrh),max(1,picrh_ion),nTi ./ max(1,nDi));
  		splusicrh = zbpi0(nDi,tii,picrh_ion,taus_icrh,e0icrhm,max(real(ecriticrh),imag(ecriticrh)),rap) ;
	end
else
  	splusicrh = zeros(size(nDi));
end

splus     = max(0,splustd) +  max(0,splusdt) +  max(0,splusff) + max(0,splusicrh) ;
%disp([salpha(end-2),splus(end-2) , nDi(end-2) .* Vp(end-2),nTi(end-2) .* Vp(end-2)])
salpha    = max(0,salpha) + max(0,splus);
pfus_nbi  = 3.56e6 .* 1.602176462e-19 .* splus; % approximation monocinetique
pfus      = max(0,pfus) + max(0,pfus_nbi);


% correction des pertes de premiere orbite
pfus_loss  = frloss .* pfus;
pfus       = pfus - pfus_loss;
salpha     = salpha .* (1-frloss);

%    figure(56);clf
%    tl = 1:length(splus);
%    semilogy(tl,salpha,'b',tl,splustd,'r',tl,splusdt,'m',tl,splusff,'c',tl,splusicrh,'k');
%    legend('alpha','T->D','D->T','BB','ICRH')
%  %plot(tl,pfus,'b',tl,pfus_nbi,'r',tl,pddfus,'k',tl,pddfus,'k');
%  %disp('in zfus0')
%  %keyboard
%  %plot(x,jfus,'b',x,z0j3(x,j0fus,xfus,jxfus),'r')
%      drawnow
%  


function  s  = zpmean(t,p,e)

indok = find(isfinite(p) & isfinite(e));
if length(indok) < 3
	indok = find(isfinite(e));
	if isempty(indok)
		s = NaN;
		disp('zpmean : invalid data ...');
	else
		s = mean(e(indok));
	end
else
	t = t(indok);
	p = p(indok);
	e = e(indok);
	indp = find(p > 0);
	if length(indp) > 2
		s = trapz(t(indp),p(indp) .* e(indp)) ./ trapz(t(indp),eps + p(indp));
	else
		s = mean(e);
	end
end


function gpmax = intgpmax(x,c,d)

ve = ones(size(x));
vt = ones(size(c,1),1);


int1   = cumtrapz(x, 1./ c,2) ./ 2;
imoins = exp(-int1);
iplus  = exp(int1);
iplus(~isfinite(imoins)) = 0;
imoins(~isfinite(imoins)) = 0;
imoins(~isfinite(iplus))   = 0;
iplus(~isfinite(iplus))   = 0;

gpmax  = (2 .* d + iplus .* cumtrapz(x,d ./ c .* imoins,2) -  ...
          2 .* (d(:,1) * ve) .* iplus) ./ 2 ./ c;

gpmax(~isfinite(gpmax)) = -Inf;
