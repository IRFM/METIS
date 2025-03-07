% compute fusion power for p-B11 reaction
% only for option.gaz == 11 
%splustd  become splusbp
%splusdt  become spluspb
function [pfus,salpha,ifus,xfus,jxfus,j0fus,taus,ecrit,pfus_nbi,pfus_loss,jfus,salf,palf,splusbp,spluspb,splusff,splusicrh] = ...
	     zfus0tae_pB11(nDi,nTi,te,ne,zeff,tite,R,a,K,Bt,ane,ate,Vp,Sp, ...
         pnbi_th,taus_nbi,ecritnbi,e0nbi,e0nbi2,ftnbi,picrh_ion,taus_icrh,ecriticrh,e0icrh,temps,pnbi,d0,qa,qmin,te0, ...
         nebord,tebord,pped,ni,wth,taeonoff,nb_nbi,fspot,e_shielding,profli,fpolarized,forced_H_NBI,mino)

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
    nH  = profli.n1p - nD;
    nHi = trapz(x,vpr .* abs(nH),2) ./ trapz(x,vpr,2);
	nB  = profli.n1p .* ((nTi./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
    nBi = trapz(x,vpr .* abs(nB),2) ./ trapz(x,vpr,2);
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
	nB  = ((max(1,nTi) ./ ne) * ve) .* nep;
    nBi = trapz(x,vpr .* abs(nB),2) ./ trapz(x,vpr,2);
    % coarse approximation for first evaluation
    nH  = max(0,nip - nD - 5 .* nB);
    nHi = trapz(x,vpr .* abs(nH),2) ./ trapz(x,vpr,2);
	zeffp      = zeff * ve;
	Raxe       = R * ve + d0 * ux;
	epsi       = max(eps,(a * x)) ./ Raxe;
end
if isfield(profli,'ftrap')
	ftrap = profli.ftrap;
else
	ftrap = 0.95 .* sqrt(epsi);
end

% source of alpha particles for thermal plasma
% p + B11 ->  3 .* alpha
rate_pB11 = pb11_fusion_rate_tentori(max(30,tip));
s         = 3 .*  nH .* nB .* rate_pB11 .* (1 + 0.6 .* fpolarized);
% mean energy of alpha particles
%e_alpha_mean = (1/3) * 4e6 + (2/3) * 2.3e6;
e_alpha_mean = 8.68e6 / 3;

% calcul de l'elargissement du a la largeur d'orbite
% article L.G. Eriksson & Porcelli
% il y  a peu d'elargissement car les particles passe leur temps pres des pointes de bananes.

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
    while (nbcount >0) && (any(indinstable))
        % calcul de l'effet des TAE sur la largeur de la source de alphas
        VA       = 2.18e16 .* (Bt * ve) ./ sqrt(1 .* nH + 11.0093054 .* nB);                 % vitesse de Alfven
        pres_alf = e_alpha_mean .* 1.602176462e-19 .* s  ./ 2 .*  (taus * ve);      % pression des alpha en premiere approximation
        taun     = (1/3) .* taus ./ 3 .* log(1 + (4e6 ./ ecrit) .^(3/2)) + ...
                   (2/3) .* taus ./ 3 .* log(1 + (2.3e6 ./ ecrit) .^(3/2));           % temps d'arret des alpha
        n_alf    = s .* (taun * ve);                                          % nombre d'alpha dans la decharge nonthermalise
        emoy     = max(1e5,min(e_alpha_mean,pres_alf ./ max(1,n_alf)./  1.602176462e-19));  % energy moyenne en eV pour F(x)
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
palf     = e_alpha_mean .* 1.602176462e-19 .* s;
salpha   = trapz(x,vpr .* s,2);
pfus     = e_alpha_mean .* 1.602176462e-19 .* salpha;


% calcul de la puissance de perte 1ere orbite
% il faut ajouter le pertes ripples et de diffusion angulaire
rloc    = R * ve + a *x;
ral_1   = 4.55e-3.* sqrt(4e6 ./ 1e3) ./ (((Bt .* R)*ve) ./ rloc); % en m
ral_23  = 4.55e-3.* sqrt(2.3e6 ./ 1e3) ./ (((Bt .* R)*ve) ./ rloc); % en m
%qp     = qmin * ve + (qa -qmin) * (x .^ 2);
if isfield(profli,'qjli')
       	qp       = profli.qjli;
else
	qp     = z0qp(x,max(1,qmin),qa);
end
% patato
dp1_1     = (R*ve) .* (2 .* qp .* ral_1 ./ (R *ve)) .^ (2/3);
dp1_23    = (R*ve) .* (2 .* qp .* ral_23 ./ (R *ve)) .^ (2/3);
% banana
dp2_1     = sqrt((a * x) ./ (R * ve +  d0 * ux)) .* ral_1 .* qp;
dp2_23    = sqrt((a * x) ./ (R * ve +  d0 * ux)) .* ral_23 .* qp;
dp_1      = dp2_1 .* (dp2_1 < (a * x)) + dp1_1 .* (dp2_1 >= (a * x));
dp_23     = dp2_23 .* (dp2_1 < (a * x)) + dp1_23 .* (dp2_23 >= (a * x));
%figure(20);clf;plot(x,dp1,'r',x,dp2,'b',x,dp,'k',x,ral,'g');drawnow
% 0.05 = evite les artefacts petite machine
mask_1   = (ral_1 + dp_1 + d0 * ux + a *x - (a - 0.05) *ve) >0;
mask_23   = (ral_23 + dp_23 + d0 * ux + a *x - (a - 0.05) *ve) >0;
%sloss  = trapz(x,vpr .* s .* mask,2) + dsloss;
sloss  = (1/3) .* trapz(x,vpr .* s .* mask_1,2) + (2/3) .* trapz(x,vpr .* s .* mask_23,2);
frloss = max(0,min(1,sloss ./ (salpha + eps)));  % fraction d'alpha perdue

% fspot  est un facteur de calibration calculer avec spot
% fspot   = 0.15;
% calcul du courant genere par les alpha (ref L-G Eriksson and F. Porcelli, Plasma Phys. Controlled Fus. 43 (2001) R145-)
% formule p R175, 9.4
if fspot > 0
    dp0R_1    = (2 .* pi) .* dp_1 ./ rloc;
    valpha_1  = sqrt(2 .*  4e6 .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
    dp0R_23    = (2 .* pi) .* dp_23 ./ rloc;
    valpha_23  = sqrt(2 .*  2.3e6 .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
else
    wca        = 2 .* 1.602176462e-19 .* Bt ./ 1.66053873e-27 ./ 4;
    valpha_1   = sqrt(2 .*  4e6 .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
    dp0R_1     = (2 .* (max(1,qmin)*ve) .* valpha_1 ./ (wca * ve) ./ (R*ve)) .^ (2/3);
    valpha_23  = sqrt(2 .*  2.3e6 .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
    dp0R_23    = (2 .* (max(1,qmin)*ve) .* valpha_23 ./ (wca * ve) ./ (R*ve)) .^ (2/3);
end
% calcul de l'effet du courant de retour
%GZ      = (1.55 + 0.85./(zeff*ve)) .* sqrt((a./R) * x) - (0.2 + 1.55./(zeff*ve)) .* ((a./R) * x);
%jfus    = -fspot .* 2 .* 1.602176462e-19 .* valpha .* dp0R .* pdederive(x,s,2,0,2,1) .* (1 - 2 .* %(1-GZ) ./ (zeff*ve));
fvd     = (1/3) .* valpha_1 .* dp0R_1 + (2/3) .* valpha_23 .* dp0R_23;
if fspot > 0
    jfus    = -abs(fspot) .* 2 .* 1.602176462e-19 .* fvd .* pdederive(x,s,0,0,2,1);
    lnldei  = 15.2 - 0.5 .* log(nep./1e20) + log(tep ./1e3);
    taus_mat = 6.27e8 .* 4 .* tep .^ (3/2) ./ (nep./ 1e6) ./ lnldei;
    taus_ave = trapz(x,vpr .* taus_mat,2) ./ trapz(x,vpr,2);
    jfus    = jfus .* (taus_mat ./ max(eps,taus_ave * ve));
else
    jfus    = -abs(fspot) .* 2 .* 1.602176462e-19 .* fvd .* pdederive(x,s,2,0,2,1);
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
    mask    = (1/3) .* (jfus == (max(jfus,[],2)*ve)) .* mask_1 + ...
              (2/3) .* (jfus == (max(jfus,[],2)*ve)) .* mask_23;
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



% calcul de l'interraction faisceau-plasma (B-> p)
if (all(real(ftnbi == 0)) && all(imag(ftnbi == 0)))  || all( pnbi  < 1e3) 
    splusbp   = zeros(size(temps));
    emeanb    = e0nbi;
else
    %[esupranbi,pTth,tauseff] = zsupra0(temps,pnbi .* max(1e-6,ftnbi),taus_nbi,1e-6,ecritnbi,e0nbi,3);
    %[splustd,emeanb]   = zbpi0(nDi,tii,pTth,taus_nbi,e0nbi,ecritnbi,1);
    [splusbp,emeanb]   = zbpi0_bp(nHi,tii,real(pnbi_th) .* max(1e-6,real(ftnbi)),real(taus_nbi),e0nbi,real(ecritnbi),1);
    if nb_nbi > 1
        [splusbp2,emeanb2]   = zbpi0_bp(nHi,tii,imag(pnbi_th) .* max(1e-6,imag(ftnbi)),imag(taus_nbi),e0nbi2,imag(ecritnbi),1);
        splusbp = splusbp + splusbp2;
        %emean   = (real(pnbi) .*emean   + imag(pnbi) .* emean2) ./ max(1,real(pnbi) + imag(pnbi));
        emeanb2(emeanb2 <=0) = 1;
        emeanb(emeanb <=0) = 1;
        emeanb   = emeanb   + sqrt(-1) .* emeanb2;
    end
end
% proton beam -> boron
if	(all(real(ftnbi == 1)) && all(imag(ftnbi == 1)))  || all( pnbi  < 1e3)
    spluspb   = zeros(size(temps));
    emeanh    = e0nbi;
else
    [spluspb,emeanh]   = zbpi0_pb(nBi,tii,real(pnbi_th),real(taus_nbi),e0nbi,real(ecritnbi),1-real(ftnbi));
    if nb_nbi > 1
        [spluspb2,emeanh2]   = zbpi0_pb(nBi,tii,imag(pnbi_th),imag(taus_nbi),e0nbi2,imag(ecritnbi),1-imag(ftnbi));
        spluspb    = spluspb + spluspb2;
        emeanh2(emeanh2 <=0) = 1;
        emeanh(emeanh <=0) = 1;
        emeanh   = emeanh   + sqrt(-1) .* emeanh2;
    end
end
%
emean = (real(emeanb) .* real(ftnbi) + real(emeanh) .* (1 - real(ftnbi))) + sqrt(-1) .* ...
        (imag(emeanb) .* imag(ftnbi) + imag(emeanh) .* (1 - imag(ftnbi)));
emean = max(max(30,te),real(emean)) + sqrt(-1) .* max(max(30,te),imag(emean));

% approximation faisceau-faisceau avec la meme formule (c'est faux) cas ftnbi intermediaire 
% Boron beam is slower than hydrogen beam at same energy and is considered
% to be the thermal species in this case with a effective temperature take
% on the averaged slowind down distribution function.
if all(ftnbi == 0) || all( pnbi  < 1e3)
    splusff   = zeros(size(temps));
else
    % compute boron equivalent background
    [esupranbiB,pBth,tauseffB]   = zsupra0(temps,real(pnbi) .* max(1e-6,real(ftnbi)),real(taus_nbi),1e-6 .* ones(size(taus_nbi)),real(ecritnbi),e0nbi,11.0093054);
    esupranbiB(~isfinite(esupranbiB))= 0;
    splusff   = zbpi0_pb(pBth .* tauseffB ./ 2 ./ Vp./ (1.602176462e-19 .* real(emeanb)), ...
                         real(emeanb), real(pnbi) .* max(1e-6,1-real(ftnbi)),real(taus_nbi),e0nbi,real(ecritnbi), 1);
    if nb_nbi >1
        [esupranbiB2,pBth2,tauseffB2]   = zsupra0(temps,imag(pnbi) .* max(1e-6,imag(ftnbi)),imag(taus_nbi),1e-6 .* ones(size(taus_nbi)),imag(ecritnbi),e0nbi2,11.0093054);
        supranbiB2(~isfinite(esupranbiB2))= 0;
        % 2 ->2
        splusff22   = zbpi0_pb(pBth2 .* tauseffB2 ./ 2 ./ Vp./ (1.602176462e-19 .* imag(emeanb)), ...
                      imag(emeanb), imag(pnbi) .* max(1e-6,1-imag(ftnbi)),imag(taus_nbi),e0nbi2,imag(ecritnbi), 1);
        
        % 1 -> 2
        splusff12   = zbpi0_pb(pBth .* tauseffB ./ 2 ./ Vp./ (1.602176462e-19 .* real(emeanb)), ...
                      imag(emeanb), imag(pnbi) .* max(1e-6,1-imag(ftnbi)),imag(taus_nbi),e0nbi2,imag(ecritnbi), 1);
        
        % 2 -> 1
        splusff21   = zbpi0_pb(pBth2 .* tauseffB2./ 2 ./ Vp./ (1.602176462e-19 .* imag(emeanb)), ...
                      real(emeanb), real(pnbi) .* max(1e-6,1 - real(ftnbi)),real(taus_nbi),e0nbi,real(ecritnbi), 1);
        
        splusff = splusff + splusff22 + splusff12 + splusff21;
    end
end

%  % acceleration de bore ou proton par icrh (approximation type idn)
e0icrhm =zpmean(1:length(e0icrh),picrh_ion,e0icrh);
if all(picrh_ion < 1e3)
    splusicrh = zeros(size(nBi));
elseif ~isempty(e0icrhm)
    if e0icrhm == 0
        splusicrh = zeros(size(nBi));
    else
        switch mino
            case 'H'
                rap 	  = zpmean(1:length(e0icrh),max(1,picrh_ion),nHi ./ max(1,nHi + nBi));
                splusicrh = zbpi0_pb(nBi,tii,picrh_ion,taus_icrh,e0icrhm,max(real(ecriticrh),imag(ecriticrh)),rap) ;
           case 'B'
                rap 	  = zpmean(1:length(e0icrh),max(1,picrh_ion),nBi ./ max(1,nHi + nBi));
                splusicrh = zbpi0_bp(nHi,tii,picrh_ion,taus_icrh,e0icrhm,max(real(ecriticrh),imag(ecriticrh)),rap) ;
            otherwise
                splusicrh = zeros(size(nBi));                
        end
    end
else
    splusicrh = zeros(size(nBi));
end

splus     = max(0,splusbp) +  max(0,spluspb) +  max(0,splusff) + max(0,splusicrh) ;
%disp([salpha(end-2),splus(end-2) , nDi(end-2) .* Vp(end-2),nTi(end-2) .* Vp(end-2)])
salpha    = max(0,salpha) + max(0,splus);
pfus_nbi  = e_alpha_mean .* 1.602176462e-19 .* splus; % approximation monocinetique
pfus      = max(0,pfus) + max(0,pfus_nbi);


% correction des pertes de premiere orbite
pfus_loss  = frloss .* pfus;
pfus       = pfus - pfus_loss;
salpha     = salpha .* (1-frloss);



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
