% calcul de la puissance de fusion
function [pfus0,salpha0,pfus,salpha,ifus,xfus,jxfus,j0fus,taus,ecrit,pfus_nbi,pfus_loss,jfus,salf0,palf0,salf,palf,splusHe3D,splusD3He,splusDT,splusff,splusicrh] = ...
	 zfus0tae_dhe(nDi,nHe3i,nTi,te,ne,zeff,tite,R,a,K,Bt,ane,ate,Vp,Sp, ...
         pnbi_th,taus_nbi,ecritnbi,e0nbi,e0nbi2,fHe3nbi,picrh_ion,taus_icrh,ecriticrh,e0icrh,temps,pnbi,d0,qa,qmin,te0, ...
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
    nHe3  = profli.n1p .* ((nHe3i./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
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
    nHe3  = ((max(1,nHe3i) ./ ne) * ve) .* nep;
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
s_DHe3   = nD.*nHe3.*dhe3.sv.*(1 + fpolarized / 2);
s_DDp    = 0.5.*nD.^2.*dd_p.sv;
s_DDn    = 0.5.*nD.^2.*dd_n.sv;
s_DT     = nD .* nT .* dt.sv;
s = s_DHe3 + s_DDp + s_DDn + s_DT;

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


% % valeur scalaire
% salf     = s;
% profli.alphashape = s;
% palf     = 3.56e6 .* 1.602176462e-19 .* s;
% salpha   = trapz(x,vpr .* s,2);
% pfus     = 3.56e6 .* 1.602176462e-19 .* salpha;

salf0.he4_DHe3 = s_DHe3;
salf0.p_DHe3   = s_DHe3;

salf0.he3_DDn  = s_DDn;
salf0.n_DDn    = s_DDn;
salf0.t_DDp    = s_DDp;
salf0.p_DDp    = s_DDp;

salf0.he4_DT   = s_DT;
salf0.n_DT     = s_DT;

salf0.he4      = salf0.he4_DHe3 + salf0.he4_DT;
salf0.p        = salf0.p_DHe3 + salf0.p_DDp;
salf0.he3      = salf0.he3_DDn;
salf0.t        = salf0.t_DDp;
salf0.n        = salf0.n_DDn  + salf0.n_DT;

salf0.tot = salf0.he4 + salf0.p + salf0.he3 + salf0.t;
salf      = salf0.tot;

profli.aphashap = salf;

% power density of fast ions produced by thermal fusion reactions
palf0.he4_DHe3 = salf0.he4_DHe3 .* dhe3.he4 .* 1.602176462e-19;
palf0.p_DHe3   = salf0.p_DHe3   .* dhe3.p   .* 1.602176462e-19;

palf0.he3_DDn  = salf0.he3_DDn  .* dd_n.he3 .* 1.602176462e-19;
palf0.n_DDn    = salf0.n_DDn    .* dd_n.n   .* 1.602176462e-19;
palf0.t_DDp    = salf0.t_DDp    .* dd_p.t   .* 1.602176462e-19;
palf0.p_DDp    = salf0.p_DDp    .* dd_p.p   .* 1.602176462e-19;

palf0.he4_DT   = salf0.he4_DT   .* dt.he4   .* 1.602176462e-19;
palf0.n_DT     = salf0.n_DT     .* dt.n     .* 1.602176462e-19;

palf0.he4      = palf0.he4_DHe3 + palf0.he4_DT;
palf0.p        = palf0.p_DHe3 + palf0.p_DDp;
palf0.he3      = palf0.he3_DDn;
palf0.t        = palf0.t_DDp;
palf0.n        = palf0.n_DDn  + palf0.n_DT;

palf0.tot  = palf0.he4 + palf0.p + palf0.he3 + palf0.t;
palf       = palf0.tot;

% fast ion density produced by thermal fusion reactions
salpha0.he4_DHe3 = trapz(x,vpr .* salf0.he4_DHe3, 2);
salpha0.p_DHe3   = trapz(x,vpr .* salf0.p_DHe3,   2);

salpha0.he3_DDn  = trapz(x,vpr .* salf0.he3_DDn, 2);
salpha0.n_DDn    = trapz(x,vpr .* salf0.n_DDn,   2);
salpha0.t_DDp    = trapz(x,vpr .* salf0.t_DDp,    2);
salpha0.p_DDp    = trapz(x,vpr .* salf0.p_DDp,    2);

salpha0.he4_DT   = trapz(x,vpr .* salf0.he4_DT,   2);
salpha0.n_DT     = trapz(x,vpr .* salf0.n_DT,     2);

salpha0.he4      = salpha0.he4_DHe3 + salpha0.he4_DT;
salpha0.p        = salpha0.p_DHe3 + salpha0.p_DDp;
salpha0.he3      = salpha0.he3_DDn;
salpha0.t        = salpha0.t_DDp;
salpha0.n        = salpha0.n_DDn  + salpha0.n_DT;

salpha = salpha0.he4 + salpha0.p + salpha0.he3 + salpha0.t;

% volume intergrated power of fast ions produced by thermal fusion reactions
pfus0.he4_DHe3 = salpha0.he4_DHe3 .* dhe3.he4 .* 1.602176462e-19;
pfus0.p_DHe3   = salpha0.p_DHe3   .* dhe3.p   .* 1.602176462e-19;

pfus0.he3_DDn  = salpha0.he3_DDn  .* dd_n.he3 .* 1.602176462e-19;
pfus0.n_DDn    = salpha0.n_DDn    .* dd_n.n   .* 1.602176462e-19;
pfus0.t_DDp    = salpha0.t_DDp    .* dd_p.t   .* 1.602176462e-19;
pfus0.p_DDp    = salpha0.p_DDp    .* dd_p.p   .* 1.602176462e-19;

pfus0.he4_DT   = salpha0.he4_DT   .* dt.he4   .* 1.602176462e-19;
pfus0.n_DT     = salpha0.n_DT     .* dt.n     .* 1.602176462e-19;

pfus0.he4      = pfus0.he4_DHe3 + pfus0.he4_DT;
pfus0.p        = pfus0.p_DHe3 + pfus0.p_DDp;
pfus0.he3      = pfus0.he3_DDn;
pfus0.t        = pfus0.t_DDp;
pfus0.n        = pfus0.n_DDn  + pfus0.n_DT;

pfus = pfus0.he4 + pfus0.p + pfus0.he3 + pfus0.t;


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
sloss  = trapz(x,vpr .* salf .* mask,2);
frloss = max(0,min(1,sloss ./ (salpha + eps)));  % fraction d'alpha perdue

% fspot  est un facteur de calibration calculer avec spot
% fspot   = 0.15;
% calcul du courant genere par les alpha (ref L-G Eriksson and F. Porcelli, Plasma Phys. Controlled Fus. 43 (2001) R145-)
% formule p R175, 9.4

valpha.he4_DHe3 = sqrt(2 .*  dhe3.he4 .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
valpha.p_DHe3   = sqrt(2 .*  dhe3.p   .* 1.602176462e-19 ./  1.66053873e-27 ./ 1);
valpha.he3_DDn  = sqrt(2 .*  dd_n.he3 .* 1.602176462e-19 ./  1.66053873e-27 ./ 3);
valpha.t_DDp    = sqrt(2 .*  dd_p.t   .* 1.602176462e-19 ./  1.66053873e-27 ./ 3);
valpha.p_DDp    = sqrt(2 .*  dd_p.p   .* 1.602176462e-19 ./  1.66053873e-27 ./ 1);
valpha.he4_DT   = sqrt(2 .*  dt.he4   .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
if fspot > 0
	dp0R    = (2 .* pi) .* dp ./ rloc;
else
  	wca     = 2 .* 1.602176462e-19 .* Bt ./ 1.66053873e-27 ./ 4;
  	dp0R.he4_DHe3 = (2 .* (max(1,qmin)*ve) .* valpha.he4_DHe3 ./ (wca * ve) ./ (R*ve)) .^ (2/3);
    dp0R.p_DHe3   = (2 .* (max(1,qmin)*ve) .* valpha.p_DHe3   ./ (wca * ve) ./ (R*ve)) .^ (2/3);
    dp0R.he3_DDn  = (2 .* (max(1,qmin)*ve) .* valpha.he3_DDn  ./ (wca * ve) ./ (R*ve)) .^ (2/3);
    dp0R.t_DDp    = (2 .* (max(1,qmin)*ve) .* valpha.t_DDp    ./ (wca * ve) ./ (R*ve)) .^ (2/3);
    dp0R.p_DDp    = (2 .* (max(1,qmin)*ve) .* valpha.p_DDp    ./ (wca * ve) ./ (R*ve)) .^ (2/3);
    dp0R.he4_DT   = (2 .* (max(1,qmin)*ve) .* valpha.he4_DT   ./ (wca * ve) ./ (R*ve)) .^ (2/3);
end
% calcul de l'effet du courant de retour
%GZ      = (1.55 + 0.85./(zeff*ve)) .* sqrt((a./R) * x) - (0.2 + 1.55./(zeff*ve)) .* ((a./R) * x);
%jfus    = -fspot .* 2 .* 1.602176462e-19 .* valpha .* dp0R .* pdederive(x,s,2,0,2,1) .* (1 - 2 .* %(1-GZ) ./ (zeff*ve));
if fspot > 0
	jfus0.he4_DHe3 = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.he4_DHe3 .* dp0R .* pdederive(x,s,0,0,2,1);
    jfus0.p_DHe3   = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.p_DHe3   .* dp0R .* pdederive(x,s,0,0,2,1);
    jfus0.he3_DDn  = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.he3_DDn  .* dp0R .* pdederive(x,s,0,0,2,1);
    jfus0.t_DDp    = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.t_DDp    .* dp0R .* pdederive(x,s,0,0,2,1);
    jfus0.p_DDp    = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.p_DDp    .* dp0R .* pdederive(x,s,0,0,2,1);
    jfus0.he4_D    = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.he4_DT   .* dp0R .* pdederive(x,s,0,0,2,1);
    jfus0.tot = jfus0.he4_DHe3 + jfus0.p_DHe3 + jfus0.he3_DDn + jfus0.t_DDp + jfus0.p_DDp + jfus0.he4_D;
        lnldei  = 15.2 - 0.5 .* log(nep./1e20) + log(tep ./1e3);
        taus_mat = 6.27e8 .* 4 .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei;
        taus_ave = trapz(x,vpr .* taus_mat,2) ./ trapz(x,vpr,2);
	jfus0.tot    = jfus0.tot .* (taus_mat ./ max(eps,taus_ave * ve));
    jfus = jfus0.tot;
else
    jfus0.he4_DHe3 = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.he4_DHe3 .* dp0R.he4_DHe3 .* pdederive(x,s,0,0,2,1);
    jfus0.p_DHe3   = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.p_DHe3   .* dp0R.p_DHe3   .* pdederive(x,s,0,0,2,1);
    jfus0.he3_DDn  = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.he3_DDn  .* dp0R.he3_DDn  .* pdederive(x,s,0,0,2,1);
    jfus0.t_DDp    = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.t_DDp    .* dp0R.t_DDp    .* pdederive(x,s,0,0,2,1);
    jfus0.p_DDp    = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.p_DDp    .* dp0R.p_DDp    .* pdederive(x,s,0,0,2,1);
    jfus0.he4_D    = -abs(fspot) .* 2 .* 1.602176462e-19 .* valpha.he4_DT   .* dp0R.he4_DT   .* pdederive(x,s,0,0,2,1);
    jfus0.tot = jfus0.he4_DHe3 + jfus0.p_DHe3 + jfus0.he3_DDn + jfus0.t_DDp + jfus0.p_DDp + jfus0.he4_D;
    jfus = jfus0.tot;
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



% calcul de l'interraction faisceau-plasma (D3He)
if (all(real(fHe3nbi == 0)) && all(imag(fHe3nbi == 0)))  || all( pnbi  < 1e3) || (forced_H_NBI ~= 0)
    esupranbi = zeros(size(temps));
    pHe3th      = zeros(size(temps));
    tauseff   = 1;
    splusHe3D   = zeros(size(temps));
    emean     = e0nbi;
else
    %[esupranbi,pHe3th,tauseff] = zsupra0(temps,pnbi .* max(1e-6,fHe3nbi),taus_nbi,1e-6,ecritnbi,e0nbi,3);
    %[splusHe3D,emean]   = zbpi0(nDi,tii,pHe3th,taus_nbi,e0nbi,ecritnbi,1);
    [splusHe3D,emean]   = zbpi0_He3D(nDi,tii,real(pnbi_th), real(taus_nbi),e0nbi,real(ecritnbi), max(1e-6,real(fHe3nbi)) );
    if nb_nbi > 1
        [splusHe3D2,emean2]   = zbpi0_He3D(nDi,tii,imag(pnbi_th), imag(taus_nbi),e0nbi2,imag(ecritnbi), max(1e-6,imag(fHe3nbi)) );
        splusHe3D = splusHe3D + splusHe3D2;
        %emean   = (real(pnbi) .*emean   + imag(pnbi) .* emean2) ./ max(1,real(pnbi) + imag(pnbi));
        emean2(emean2 <=0) = 1;
        emean(emean <=0) = 1;
        emean   = emean   + sqrt(-1) .* emean2;
    end
end
% attention la matrice est 3He-> D, IL FAUDRA RAJOUTER 	D -> 3He
% add D beam injection to T ions in plasma,
% and do not consider the case of T beam injection for option.gaz=5.
% maintenant c'est la bonne matrice
if	(all(real(fHe3nbi == 1)) && all(imag(fHe3nbi == 1)))  || all( pnbi  < 1e3) || (forced_H_NBI ~= 0)
    splusD3He   = zeros(size(temps));
    splusDT   = zeros(size(temps));
else
    splusD3He   = zbpi0_D3He(nHe3i,tii,real(pnbi_th),real(taus_nbi),e0nbi,real(ecritnbi),1-real(fHe3nbi));
    splusDT     = zbpi0_dt(nTi,tii,real(pnbi_th),real(taus_nbi),e0nbi,real(ecritnbi),1-real(fHe3nbi));
    if nb_nbi > 1
        splusD3He2   = zbpi0_D3He(nHe3i,tii,imag(pnbi_th),imag(taus_nbi),e0nbi2,imag(ecritnbi),1-imag(fHe3nbi));
        splusD3He    = splusD3He + splusD3He2;
        splusDT2     = zbpi0_dt(nTi,tii,imag(pnbi_th),imag(taus_nbi),e0nbi2,imag(ecritnbi),1-imag(fHe3nbi));
        splusDT      = splusDT + splusDT2;
    end
end
% approximation faisceau-faisceau avec la meme formule (c'est faux) cas fHe3nbi intermediaire (D3He)
emean = max(max(30,te),real(emean)) + sqrt(-1) .* max(max(30,te),imag(emean));
if all(fHe3nbi == 0) || all( pnbi  < 1e3) || (forced_H_NBI ~= 0)
   splusff   = zeros(size(temps));
else
   [esupranbi,pHe3th,tauseff]   = zsupra0(temps,real(pnbi) .* max(1e-6,real(fHe3nbi)),real(taus_nbi),1e-6 .* ones(size(taus_nbi)),real(ecritnbi),e0nbi,3);
   esupranbi(~isfinite(esupranbi))= 0;
   splusff   = zbpi0_He3D((1-real(fHe3nbi)) .* real(pnbi_th) .* tauseff ./ 2 ./ Vp./ (1.602176462e-19 .* real(emean)), ...
               real(emean), real(pnbi), real(taus_nbi),e0nbi,real(ecritnbi), max(1e-6,real(fHe3nbi)) );
   if nb_nbi >1 
  	[esupranbi2,pHe3th2,tauseff2]   = zsupra0(temps,imag(pnbi) .* max(1e-6,imag(fHe3nbi)),imag(taus_nbi),1e-6 .* ones(size(taus_nbi)),imag(ecritnbi),e0nbi2,3);
   	esupranbi2(~isfinite(esupranbi2))= 0;
      	% 2 ->2
    	splusff22   = zbpi0_He3D((1-imag(fHe3nbi)) .* imag(pnbi_th) .* tauseff2 ./ 2 ./ Vp./ (1.602176462e-19 .* imag(emean)), ...
               		imag(emean), imag(pnbi), imag(taus_nbi),e0nbi2,imag(ecritnbi), max(1e-6,imag(fHe3nbi)));

      	% 1 -> 2
    	splusff12   = zbpi0_He3D((1-real(fHe3nbi)) .* real(pnbi_th) .* tauseff ./ 2 ./ Vp./ (1.602176462e-19 .* real(emean)), ...
               		imag(emean), imag(pnbi), imag(taus_nbi),e0nbi2,imag(ecritnbi), max(1e-6,imag(fHe3nbi)) );

      	% 2 -> 1
        splusff21   = zbpi0_He3D((1-imag(fHe3nbi)) .* imag(pnbi_th) .* tauseff2 ./ 2 ./ Vp./ (1.602176462e-19 .* imag(emean)), ...
                      real(emean), real(pnbi), real(taus_nbi),e0nbi,real(ecritnbi), max(1e-6,real(fHe3nbi)));

        splusff = splusff + splusff22 + splusff12 + splusff21;
   end
end
%  % acceleration de 3He par icrh (approximation type idn)
e0icrhm =zpmean(1:length(e0icrh),picrh_ion,e0icrh);
if all(picrh_ion < 1e3)
  	splusicrh = zeros(size(nDi));
elseif ~isempty(e0icrhm)
	if e0icrhm == 0
  		  	splusicrh = zeros(size(nDi));
 	else
  		rap 	  = zpmean(1:length(e0icrh),max(1,picrh_ion),nHe3i ./ max(1,nDi));
  		splusicrh = zbpi0_He3D(nDi,tii,picrh_ion,taus_icrh,e0icrhm,max(real(ecriticrh),imag(ecriticrh)),rap) ;
	end
else
  	splusicrh = zeros(size(nDi));
end

splus.he4_DHe3  = max(0,splusHe3D) + max(0,splusD3He) +  max(0,splusff) + max(0,splusicrh) ;
splus.p_DHe3    = splus.he4_DHe3;
splus.he4_DT    = max(0,splusDT);
splus.n_DT      = splus.he4_DT;
%disp([salpha(end-2),splus(end-2) , nDi(end-2) .* Vp(end-2),nTi(end-2) .* Vp(end-2)])
salpha    = max(0,salpha) + max(0,splus.he4_DHe3) + max(0,splus.he4_DT) + max(0,splus.p_DHe3);
pfus_nbi  =   dhe3.he4 .* 1.602176462e-19 .* splus.he4_DHe3 ...
            + dhe3.p   .* 1.602176462e-19 .* splus.p_DHe3  ...
            + dt.he4   .* 1.602176462e-19 .* splus.he4_DT ; % approximation monocinetique
pfus      = max(0,pfus) + max(0,pfus_nbi);
salpha0.n_DT = max(0,salpha0.n_DT) + max(0,splus.n_DT);

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
