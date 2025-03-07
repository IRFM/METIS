function [ipout,iboot,ioh,poh,RR,vloop,qa,q95,qmin,q0,betap,piqj,wbp,dwbpdt,asser,tauj,lif,tauip,hitb,xdep,ate,aitb,hmhd, ...
                 te0,fwcorr,xlh,dlh,zeff,zmszl,qeff,efficiency,vmes,difcurconv,ipar,x,profil,phiplasma,indice_inv,poynting,kidds_evol]=  ...
                 zboot0diff(mode,modeh,vref,temps,ne,te,zeffin,R,a,K,d,Bt,ip,icd,li,wth, ...
                       ane,ate,Vp,Sp,Sext,hitb,transitoire,taue,tauip,peri,inbicd,xnbi,piqnbi,ieccd,xeccd,ifwcd, ...
                        ilh,xlh,dlh,ifus,xfus,jxfus,j0fus,runaway,tebord,nebord,swlh,pped,pfweh,picrh,plh,pecrh,pfus, ...
                        pnbi,tite,wtotal,ni,wlh,gaz,npar0,freqlh,pcyclo,pbrem,prad, ...
			vloopin,zu1,zu2,zimp,zmax,rimp,modezeff,nhem,frhe0, ...
			pion_icrh,pion_nbi,pion_fus,kishape,qdds,kidds,modeboot, ...
			Rsepa,Zsepa,amorti,directivity,friplh,drmdt,fracmino,xres_icrh,pioniz,irun,...
			difcurconv,laochange,evolution,edge_flux,breakdown,sitb,indice_inv_in,meff,mode_expo_inte,cronos_regul,...
                        t_switch,bootmul,ffit_ped,upshiftmode,fupshift,ddsmode,w1,epsq,npar_neg, ...
                        itb_sensitivity,itb_slope_max,faccu,xieorkie,berror,profil,flux_ip,L_ext,K_min,d0_ext, ...
                        tune_frac,width_ecrh,hollow,moments_mode,q0_dds_trig,betap1crit, ...
                        force_spitzer,neutral_friction,f_eta_turb,pioniz_i,Sn_fraction,cor_rel_spitzer,iso,nTm, ...
                        cmin,natural_nD_o_nH,collapse)

                    
                    
% compatibilite 
if all(size(bootmul) == 1)
	bootmul = bootmul * ones(size(te));
end
% decodage
freebie    = imag(evolution);
evolution  = real(evolution);
% decodage qdds
s1crit = imag(qdds);
qdds   = real(qdds);

% decode kidds
if any(imag(kidds))
    kidds_evol   = imag(kidds);
    kidds       = real(kidds(1));
else
    kidds_evol   = ones(size(temps));
end

% constants
mu0 = 4e-7 .* pi;
ee  = 1.602176462e-19;

% external equilibrium data  
transition_factor  = 1;                
if isappdata(0,'EQUILIBRIUM_EXP') || isappdata(0,'CURDIFF_EXP')
    equi_ext = getappdata(0,'METIS_EXTERNAL_CURDIF_EQUI');
    % get transition factor if defined
    if isappdata(0,'METIS_EXTERNAL_CURDIF_EQUI_TRANSITION_FACTOR')
        transition_factor = getappdata(0,'METIS_EXTERNAL_CURDIF_EQUI_TRANSITION_FACTOR');
    else
        transition_factor = 1;
    end
end

if isappdata(0,'EQUILIBRIUM_EXP') 
    li       = equi_ext.lif;
    piqj     = equi_ext.piqj;
    d95      = equi_ext.d95;
    K95      = equi_ext.K95;
    betap    = equi_ext.betap;
    d0       = equi_ext.d0;
    rm       = equi_ext.rm;
else
    %
    % le profil de courant
    %
    % piquage de j :
    lipj = li;
    lipj(~isfinite(lipj)) = 1;
    lipj = max(lipj,0.5);
    ip(ip == 0) =1;
    piqj = (exp(lipj)-1.65)./0.89;
    ind  = find(piqj <0.1);
    piqj(ind) = 0.1;
    ind  = find(piqj >10);
    piqj(ind) = 10;
    
    % d95  = d .* 0.95 .^ 2;         (de varie en x^ 2)
    % formule a partir de k0 = (K+1)/2 et K(x) = k0 + (K-k0) x^ 4
    K95 = 0.5 .* (K + 1)  + 0.5 .* (K - 1) .* 0.95 .^ 4;
    d95 = d .* (0.95 .^ 2);
    
    % betap (iter fdr documentaation N19 FDR 1 01-07-13 R0.1)
    ka = Vp ./ (2 .* pi .^ 2 .* R .* a.^ 2);

    if isfield(profil,'nep') && isfield(profil,'bpol')
        pep      = ee .* (profil.nep  .* profil.tep + profil.nip .* profil.tip);
        betap = trapz(profil.xli,pep .* profil.vpr,2) ./ trapz(profil.xli,profil.vpr,2) ./ ...
            (profil.bpol(:,end) .^ 2  ./ 2 ./ mu0);
    else
        betap = 8 .* wth ./ 3 ./ mu0 ./ ip .^ 2./ R .* (1 + K .^ 2) ./ 2 ./ ka;
    end
    % to prevent unconvergence during breakdown with 0 plasma current when there not yet a real confined plasma
    betap = min(max(5,2 .* R ./ a - 1),betap);
    
    if isfield(profil,'Raxe')
        d0   = profil.Raxe(:,1) - profil.Raxe(:,end);
        % 1ere approximation shafranov shift (formule plasma circulaire + correction Lao + Miller)
        d0lim   = a .^ 2  ./ 2 ./ R  .* (2.*(K .^ 2 + 1)./ (3 .* K .^ 2 + 1) .* ...
            ( (wtotal ./ max(1,wth) ) .* betap + li ./ 2) +  ...
            0.5 .* (K .^ 2 - 1)./ (3 .* K .^ 2 + 1)) .* ( 1 - a ./ R);
        d0lim   = max(0,min(a./4,d0lim));
    else
        % 1ere approximation shafranov shift (formule plasma circulaire + correction Lao + Miller)
        d0   = a .^ 2  ./ 2 ./ R  .* (2.*(K .^ 2 + 1)./ (3 .* K .^ 2 + 1) .* ...
            ( (wtotal ./ max(1,wth) ) .* betap + li ./ 2) +  ...
            0.5 .* (K .^ 2 - 1)./ (3 .* K .^ 2 + 1)) .* ( 1 - a ./ R);
        d0   = max(0,min(a./4,d0));
        d0lim = d0;
    end
    
    % coulage avec FREEBIE
    if ~isempty(d0_ext) && all(isfinite(d0_ext)) && all(d0_ext > 0)
        d0 = d0_ext;
    end
    
    % rayon moyen
    if isfield(profil,'rmx')
        rm  = profil.rmx(:,end);
    else
        rm  = peri ./ 2 ./ pi;
    end
end
% le bootstrap
iboot_sc = bootmul .* ip .*  0.45 .* sqrt( a ./ R) .* betap .* sqrt((1 + ane + ate) ./ (1 + piqj));
% securite breakdown
iboot_sc = min(2 .* ip, iboot_sc);

if (isfield(profil,'jboot') & (modeboot ~= 0))
	iboot = trapz(profil.xli,profil.jboot .* profil.spr,2);
else
	% de Th. Hoang ...
	% HOANG, G.T. et al, , Proc. of the 24nd EPS Conf., Berschtesgaden, Germany,
	%(1997), vol21A, part III, p965.
	iboot = iboot_sc;
end

if isappdata(0,'EQUILIBRIUM_EXP') 
    q95   = equi_ext.q95;
    qeff  = equi_ext.qjli(:,end);
    qa    = qeff;
    q0    = equi_ext.q0;
    j0    = equi_ext.jli(:,1);
else
    % profil de q (iter basis physics )
    q95 = 5 .* a .^ 2 .* Bt ./ (ip./1e6) ./ R .* (1 + K95 .^ 2 .*  ...
        (1 + 2 .* d95 .^ 2 - 1.2 .* d95 .^ 3) ./ 2) .* (1.17 - 0.65 .* a ./ R) ./  ...
        (1 - (a./R) .^ 2 ) .^ 2;
    qeff   = 5 .* a .^ 2 .* Bt ./ (ip./1e6) ./ R .* (1 + K .^ 2) ./ 2 .*  ( 1 + tanh((a ./ R) .^ 2 .* ...
        ( 1 + (betap + li ./ 2) .^ 2 ./ 2))) .* ...
        (1.24 - 0.54 .* K  + 0.3 .* (K .^ 2 + d .^ 2) + 0.13 .* d);
    
    % a ce  point
    qa  = qeff;
    q95 = min(q95 , qa ./ (1+ 0.25 .* modeh));
    
    j0   = ip .* (piqj + 1) ./ Sp;
    % formule Wesson standart sans correction
    % 1ere setimation
    q0   = max(0.8,q95 ./ (1 + piqj));
end
% flag commandant les donnees retournees selon le mode
flagout = 1;
% utilisation d'un forme de profil
if isfield(profil,'tep')
	x = profil.xli;
	ux  = (1  - x .^ 2);
	ve  = ones(size(x));
	vt  = ones(size(te));
	spr   = profil.spr;
	vpr   = profil.vpr;
	% le profil de ne est calculer ici
	%vppp = trapz(x,vpr,2);
	%fppp = trapz(x,vpr .* (vt * ux) .^ (ane * ve),2);
	%ne0 = (vppp .* ne  + nebord .* fppp - nebord .* vppp) ./ fppp;
	%nep = ((ne0 -nebord) * ve)  .* (vt * ux)  .^ (ane * ve) + nebord * ve;
	%nep(:,end) = nebord;

	% le profil de ne est calculer ici
	%ne0   = ne .* (1 + ane);
	%expo  = ane .* ne ./ max(1e13,ne - nebord);
	%nep = ((ne0 - nebord) * ve)  .* (vt * ux)  .^ (expo * ve) + nebord * ve;
	%nep(:,end) = nebord;
	
	nep = profil.nep;
	tep = profil.tep;
	nip = profil.nip;
	tip = profil.tip;
	
	modehp = (pped>1);
	wpied  = min(0.9 .* wth,(3/2) .* trapz(x,vpr .* (pped * ve),2));
	
	%wpied             = (3/2) .* trapz(x,vpr .* (pped * ve),2);
	%modehp            = (pped>1) .* ((wth - wpied) > 0);
	%flagout = 0;

else
	x   = linspace(0,1,21);
	ux  = (1  - x .^ 2);
	ve  = ones(size(x));
	vt  = ones(size(te));
	spr   = (2 .* Sp) * x;
	vpr   = (2 .* Vp) * x;

	% le profil de ne est calculer ici
	vppp = trapz(x,vpr,2);
	fppp = trapz(x,vpr .* (vt * ux) .^ (ane * ve),2);
    if any(~isfinite(fppp))
        indbad = find(~isfinite(fppp));
        fpppx = cumtrapz(x,vpr .* (vt * ux) .^ (ane * ve),2);
        fppp(indbad) = fpppx(indbad,end-1);
    end
	ne0 = max(nebord+1e13,(vppp .* ne  + nebord .* fppp - nebord .* vppp) ./ fppp);
	nep = ((ne0 -nebord) * ve)  .* (vt * ux)  .^ (ane * ve) + nebord * ve;
	nep(:,end) = nebord;

	tep = ((te .* (1 + ate)-tebord) * ve)  .* (vt * ux)  .^ (ate *ve) + tebord * ve;
	nip = nep .* ((ni ./ ne) * ve);
	tip = tep .* (tite * ve);
	
	wpied  = min(0.9 .* wth,(3/2) .* trapz(x,vpr .* (pped * ve),2));
	modehp = (pped>1);

	% modeh /l
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
  	profil.tep  = max(1,real(tep));
 	profil.tip  = max(1,real(tip));
end

%save('loc1');

% calcul de l'accumulation d'impurete du a l'itb
%  zacc = max(zimp,zmax);
%  dnzsnz = -(zacc / 2  + 1) .* pdederive(x,tip,0,2,2,1) ./ max(1,tip) + zacc .* pdederive(x,nip,0,2,2,1) ./ max(30,nip);
%  nzacc  = exp(cumtrapz(x,dnzsnz,2));
  
%  keyboard

% pour la densite de bord de He nhea  = nea  * (nhem/nem) * (taup / tauhe)
% le temps de confinement de la matiere est supposer 2 * taue dans metis
ftauhe   = min(100,max(1,((nhem .* Vp) ./ max(1,pfus ./ (3.56e6 .* 1.602176462e-19))) ./ max(1e-6, 2 .* taue)));
% calcul du profil de zeff
switch gaz
    case 5
        [zeff,zmszl,zeffp,n1p,nhep,nzp,nip] = z0zeff_DHe3(x,nep,tep,tip,profil.nwp,zu1,zu2,zimp,zmax,zeffin,modezeff,nhem,gaz,ne,te.*tite,frhe0,vpr,ftauhe,faccu,rimp,temps,Sn_fraction,iso,nTm,cmin);
    case 11
        [zeff,zmszl,zeffp,n1p,nhep,nzp,nip] = z0zeff_pB11(x,nep,tep,tip,profil.nwp,zu1,zu2,zimp,zmax,zeffin,modezeff,nhem,gaz,ne,te.*tite,frhe0,vpr,ftauhe,faccu,rimp,temps,Sn_fraction,iso,natural_nD_o_nH);
    otherwise
        [zeff,zmszl,zeffp,n1p,nhep,nzp] = z0zeff(x,nep,tep,tip,profil.nwp,zu1,zu2,zimp,zmax,zeffin,modezeff,nhem,gaz,ne,te.*tite,frhe0,vpr,ftauhe,faccu,rimp,temps,Sn_fraction);
end
% le profil nip est calcule ici
% for pB11 nTm is boron and has to be count separately
% for DHe3, He4 is encode separtely.
% for DHe3, the density of tritium is neglected at this level
switch gaz
    case 5
         % returned from z0zeff_DHe3
        nip  = max(1e13, n1p + nep*frhe0 + nhep + nzp .* ( 1 + rimp) + profil.nwp);
        warning('missing nTp !');
    case 11
         % returned from z0zeff_pB11
    otherwise
        nip  = max(1e13, n1p + nhep + nzp .* ( 1 + rimp) + profil.nwp);
end
pep          = ee .* (nep  .* tep + nip .* tip);
if isfield(profil,'ptot') && any(imag(profil.ptot(:)))
  ptot         = pep + imag(profil.ptot);
  profil.ptot  = real(profil.ptot);
else
  ptot         = pep .* max(1,((real(wtotal) ./ (3/2) ./ trapz(x,pep .* vpr,2)) * ve));
end
% cas des DDS
if any(indice_inv_in > 1) & (qdds < 0)
    % nouveau profil de pression supra
    ptot_mem                 = ptot;
    psupra_                  = ptot - pep;
    dptotdx                  = cat(2,0 * vt,diff(psupra_ - psupra_(:,end) * ve,1,2));
    mask                     = (vt * (1:length(x))) > (indice_inv_in * ve);
    dptotdx                  = dptotdx .* mask;
    ptotnew                  = cumsum(dptotdx,2);
    ptotnew                  = ptotnew - ptotnew(:,end) * ve + psupra_(:,end) *ve ;
    % normalisation (conservation de l'energie)
    dptot                  = trapz(x,(ptotnew - psupra_) .* vpr,2);
    maskc                   = (vt * (1:length(x))) <= (indice_inv_in * ve);
    maskc(:,1:2)           = 1;
    dvp                    = trapz(x,maskc .* vpr,2);
    dptot                  = min(0,dptot ./ max(eps,dvp));
    %ptot_mem                = ptot;
    ptot                    = ptotnew - (dptot* ve) .*  (~mask) + pep;
%     figure(23)
%     plot(x,ptot-ptot_mem)
%     drawnow
    
end

% amortissement ajouter sur ptot via la boucle externe
ptot_out = ptot;
if isfield(profil,'ptot')
	ptot = profil.ptot;
else
	%ptot = ptot_out;
end


% pour plus de precision
%fwcorr = trapz(x,ee .*(nep + nip .* (tite *ve)) .* (tep ./  (tep(:,1) * ve)) .* vpr,2);

%qp  = q0 * ve + (qa -q0) * (x .^ 2);
qpr   = z0qp(x,q0,qa);
if isfield(profil,'qjli')
	qp = profil.qjli;
else
	qp   = qpr;
end


if isappdata(0,'EQUILIBRIUM_EXP')
    kx = equi_ext.kx;
    dx = equi_ext.dx;

elseif moments_mode == 0
    % leading order
    % see for exemple dicussion in https://research.tue.nl/en/studentTheses/fast-2-d-equilibrium-reconstruction-for-tokamak-transport-calcula
    kx = K * ve;
    dx = d * x;

elseif isfield(profil,'rmx')  & isfield(profil,'psi')
	% H.J. de Blank, plasam equilibrium in tokamaks, lecture note @ www.rijnh.nl 
	% publier dans  Fusion Science and Technology
	% Volume 49 ?? Number 2T ?? February 2006 ?? Pages 111-117
        % Technical Paper ?? Plasma and Fusion Energy Physics - Equilibrium and Instabilities
	%
       	% dpsidx s'annule au centre
	psid1    = pdederive(x,profil.psi,0,2,2,1);
	% dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
	psid2    = pdederive(x,profil.psi,1,0,2,2);
	if cronos_regul == 4
	    psid2(:,1) = - Bt ./ profil.qjli(:,1) .* profil.rmx(:,end) .^ 2;
	end
	% r
	%rmx      = profil.rmx;
	%rmxd1    = pdederive(x,rmx,2,2,2,1);
	%rmxd2    = pdederive(x,rmx,2,2,2,2);
	r      = a * x;
	rd1    = a * ve;
	rd2    = zeros(size(r));
	% dpsidr et d2psidr2
	dpsidr   = psid1 ./ rd1;
	d2psidr2 = psid2 ./ rd1 .^ 2 - psid1 .* rd2 ./ rd1 .^ 3;
	% dr et dr2;
	dr       = rd1 .* mean(diff(x));
	dr2      = dr .^ 2;	
	% coef 
	ukp1     =  dpsidr .* r .^ 2 ./ dr2 + d2psidr2 .* r .^ 2 ./ dr + dpsidr .* r ./ 2 ./ dr;
	uk       =  - 2 .* dpsidr .* r .^ 2 ./ dr2 - 3 .* dpsidr;
	ukm1     =  dpsidr .* r .^ 2 ./ dr2 - d2psidr2 .* r .^ 2 ./ dr - dpsidr .* r ./ 2 ./ dr;
	% calcule de S2
	S2       = zeros(size(r));
	S2(:,1)  = 0;
	S2(:,2)  = r(:,2);
	for k = 2:(length(x)-1)
		S2(:,k+1) = - (uk(:,k) .* S2(:,k) + ukm1(:,k) .* S2(:,k-1)) ./ ukp1(:,k);
	end
	% normalisation par le bord
	S2a     = r(:,end) .* (K - 1) ./ (1 + K);
	%S2a      = (K-1) .* rmx(:,end) ./ 2;
	alpha   = S2a ./ S2(:,end);
	S2      = S2 .* (alpha * ones(size(x)));
	% calcul de kx
	kx      = min(max(3,K*ones(size(x))),max(min(1,K*ve),(r + S2) ./ max(eps,r - S2)));
	%kx       = 1 + 2 .* S2 ./ max(eps,rmx);
	%kx(:,1) = max(min(1,K),pchip(x(2:end),kx(:,2:end),0));
	kx(:,1) = kx(:,2);
	kx(:,end) = K;
	% 
	% forme de la triangularite
	%
	% coef 
	ukp1     =  dpsidr .* r .^ 2 ./ dr2 + d2psidr2 .* r .^ 2 ./ dr + dpsidr .* r ./ 2 ./ dr;
	uk       =  - 2 .* dpsidr .* r .^ 2 ./ dr2 - 8 .* dpsidr;
	ukm1     =  dpsidr .* r .^ 2 ./ dr2 - d2psidr2 .* r .^ 2 ./ dr - dpsidr .* r ./ 2 ./ dr;
	% calcule de S3
	S3       = zeros(size(r));
	S3(:,1)  = 0;
	S3(:,2)  = r(:,2).^2;
	for k = 2:(length(x)-1)
		S3(:,k+1) = - (uk(:,k) .* S3(:,k) + ukm1(:,k) .* S3(:,k-1)) ./ ukp1(:,k);
	end
	% normalisation par le bord
	S3a     = asin(d) .* r(:,end) ./ 4;
	alpha   = S3a ./ S3(:,end);
	S3      = S3 .* (alpha * ones(size(x)));
	% calcul de kx
	dx      = sin(min(1,max(-1,S3./max(eps,r).*4)));
	dx(:,1) = 0;
	dx(:,end) = d;
    dx      = max(0,min(d * ones(size(x)),dx));
else
    % just for initialisation
	dk    = max(0,(K-1) ./ 2);
	kx    =  max(1,1 + dk * (x .^ 4) + dk * ve);
	dx    =  d * (x .^ 3);
end

% pour le couplage avec FREEBIE
if isappdata(0,'EQUILIBRIUM_EXP')
      % nothing to do
elseif ~isempty(K_min) && all(isfinite(K_min)) && all(K_min > 0.5)
      %kx_mem =kx;
      kx = (kx - K * ve) ./ ((min(kx(:,3:end),[],2) - K) * ve) .* ((K_min - K) * ve) + K * ve;
      %figure(17);clf;plot(x,kx_mem,'b',x,kx,'r');drawnow
end

% rayon moyen des surface
if isappdata(0,'EQUILIBRIUM_EXP')
    rmx = equi_ext.rmx;
elseif isfield(profil,'rmx')
	rmx = profil.rmx;
else
	rmx       = (a * x) .* (1+1./4.*(-1+kx)-3./64.*(-1+kx).^2+5./256.*(-1+kx).^3-175./16384.*(-1+kx).^4);
	rmx       = rmx .* ((rm ./ rmx(:,end)) * ve);
end
if isappdata(0,'EQUILIBRIUM_EXP')
    Raxe     = equi_ext.Raxe;
    epsi     = equi_ext.epsi;
elseif isfield(profil,'Raxe')
	Raxe     = profil.Raxe;
	epsi     = profil.epsi;
else
	Raxe      = (R * ve  + d0 * (1 - x .^ 2));
	epsi      = max(eps,(a * x)) ./ Raxe;
end
flag_ftrap = false;
if isappdata(0,'EQUILIBRIUM_EXP')
    % tester si le champ ftrap exist, pas disponible dans tous les cas
    ftrap = equi_ext.ftrap;
elseif isfield(profil,'bpol') & isfield(profil,'fdia')
	b2m       = profil.bpol .^ 2  + profil.fdia .^ 2 .* profil.r2i;
	bm        = profil.fdia .* profil.ri .* (1 + 0.5 .* profil.bpol ./ sqrt(profil.grho2r2) ./ ...
	            profil.fdia .* sqrt(profil.grho2));
	% calcul de bmax
	btor         = (profil.fdia ./ (Raxe - a * x));
	grho         = abs((profil.rmx(:,end) * ve) ./ max(eps,abs(pdederive(x,profil.Raxe - a * x,0,2,2,1))));
	grho(:,1)    = 2 .* grho(:,2) - grho(:,3);
	bpol         = abs(pdederive(x,profil.psi,0,2,2,1))./ (profil.Raxe - a * x) .* ...
	                grho ./ (profil.rmx(:,end) * ve);
	bmax         = sqrt(btor .^ 2 + bpol .^ 2);
	
	% variable
	h  = min(1,bm ./ bmax);
	h2 = min(1,b2m ./ bmax .^ 2);
	 
	% expression de ftrap  Lin-Liu and Miller Phys. of Plasmas 2 (5) 1995
	ftu          = 1 - h2 ./ h .^ 2 .* ( 1 - sqrt(1 - h) .* (1 + 0.5 .* h));
	ftl          = 1 - h2 .* (1./ h2 - sqrt(1-sqrt(h2)) ./ h2  - sqrt(1 - sqrt(h2)) ./ 2 ./h);
	ftc          = 1 - (1-epsi) .^ 2 ./ sqrt(1 - epsi .^ 2) ./ (1 + 1.46 .* sqrt(epsi));  % Wesson , Tokamak
	fte          = 0.75 .* ftu + 0.25 .* ftl;
	ftrap        = pchip(cat(2,x(1:2),x(end-2:end)),cat(2,ftc(:,1:2),fte(:,end-2:end)),x);

else
	ftrap     = 1 - (1-epsi) .^ 2 ./ sqrt(1 - epsi .^ 2) ./ (1 + 1.46 .* sqrt(epsi));  % Wesson , Tokamak
end
% securite ftrap
ftrap = min(1- 2.* eps,max(0,ftrap));

% formule de Sauter calcul de jboot ete eta
Btr            = Bt * ve;
if isappdata(0,'CURDIFF_EXP') 
    if isfield(profil,'psi')
        psi_av     = transition_factor * equi_ext.psi   + (1 - transition_factor) * profil.psi;  
    else
        psi_av     = equi_ext.psi;          
    end
    dpsidx         = pdederive(x,psi_av,0,2,2,1);
    F              = equi_ext.fdia;
    if isfield(profil,'psi')	
        zion           = 1 + (profil.nhep >  profil.n1p);
    else
        zion = 1;
    end
elseif isfield(profil,'psi')
	dpsidx         = pdederive(x,profil.psi,0,2,2,1);
	F              = profil.fdia;
	zion           = 1 + (profil.nhep >  profil.n1p);
else
	zion           = 1;
	F              = ((Btr(:,1) .* Raxe(:,end)) * ve);
	bpol           = rmx .* Btr ./ Raxe ./ qp;
	dpsidx         = - (rmx(:,end) * ve) .*  bpol .* Raxe ;
	dpsidx(:,1)    = 0;
end

% modification for parametric studies
if force_spitzer == 1
    eta = z0eta_spitzer(nep,tep,zeffp); % we substitute zeff to zion in this case
else
    % correction pour le burn-through
    switch modeboot
        case {5,6}
            eta  = z0etaboot_neofit(x,tep,tip,nep,nip,qp,zeffp,zion,Raxe,ftrap,epsi);
            %             % ref : V. A. Belyakov et al , PhysCon 2003 Saint Petersburg Russia (IEEE)
            %             if (berror > 0) && isfield(profil,'n0m') && isfield(profil,'nep')
            %                 % la partie composante chaude (n0) est mal caclulee dans cette phase et introduit du bruit.
            %                 % de plus ne0 >> n0
            %                 eta =  eta  + 2.7e-6 .* profil.n0m ./ max(1,profil.nep);
            %             end
        otherwise
            %             % ref : V. A. Belyakov et al , PhysCon 2003 Saint Petersburg Russia (IEEE)
            %             if (berror > 0) && isfield(profil,'n0m') && isfield(profil,'nep')
            %                 % la partie composante chaude (n0) est mal caclulee dans cette phase et introduit du bruit.
            %                 % de plus ne0 >> n0
            %                 %eta = eta + 2.7e-6 .* profil.n0m ./ max(1e13,profil.nep);
            %                 eta = zeta0(tep,nep,qp,zeffp ,rmx ,Raxe,ftrap,epsi,2.7e-6 .* profil.n0m ./ max(1,profil.nep));
            %             else
            % resistivite
            eta = zeta0(tep,nep,qp,zeffp ,rmx ,Raxe,ftrap,epsi);
            %            end
    end
end

% appleid only to Spitzer resistivity and by extension to neoclassical one,
% even if there is no complete computation for neoclassical one.
switch cor_rel_spitzer
    case 'on'
        eta = resistivity_spitzer_rel_cor(nep,tep,zeffp) .* eta;
end

% effect of cold neutral on resistivity
% ref : V. A. Belyakov et al , PhysCon 2003 Saint Petersburg Russia (IEEE)
if (berror > 0) && isfield(profil,'n0m') && isfield(profil,'nep')
    if neutral_friction ~= 0
        eta =  eta  + neutral_friction .* 2.7e-6 .* profil.n0m ./ max(1,profil.nep);
    else
        eta =  eta  + 2.7e-6 .* profil.n0m ./ max(1,profil.nep);
    end
elseif (neutral_friction ~= 0) && isfield(profil,'n0m') && isfield(profil,'nep')
    %figure(21);clf;semilogy(x,eta,'r',x,neutral_friction .* 2.7e-6 .* profil.n0m ./ max(1,profil.nep),'b');drawnow
    eta =  eta  + neutral_friction .* 2.7e-6 .* profil.n0m ./ max(1,profil.nep);
end

% resistivity induced by turbulence
if any(f_eta_turb(:) ~= 0)
   eta_t = z0eta_turb(nep,tep,epsi,zeffp,abs(f_eta_turb));
   % eta can't become negative
   eta = eta + sign(f_eta_turb) .* min(0.9 .* eta,eta_t);
end

if isappdata(0,'CURDIFF_EXP')
    if isfield(profil,'jli')
        jmoy = transition_factor * equi_ext.jli   + (1 - transition_factor) * profil.jli;
    else
        jmoy = equi_ext.jli;        
    end
elseif isfield(profil,'jli')
	jmoy  = profil.jli;
else
	jmoy  = (j0 * ve)  .* (vt * ux) .^ (piqj *ve);
end
% normalisation de jmoy pour trouver ip dans ce cas, car ip consigne peu etre different de ip simulation
ipint = trapz(x,jmoy.* spr ,2);
rap   = ip ./ max(eps,ipint);
rap(ip==0) = 1;
rap(ipint == 0) = 1;
jmoy  = (rap * ve) .* jmoy;
etaref = eta;
ploc   = eta .* jmoy .^ 2;
pmax  = trapz(x,ploc .* vpr ,2);
RR    = pmax ./ ip  .^ 2 ; % Resistance du plasma en ohm
if evolution == 0
	RR(1) = mean(RR(1:2));
end
RR(~isfinite(RR)) = 1e-12;
% limite haute pour RR
RR    = max(1e-12,min(RR,1e6 .* Vp ./max(1,ip)  .^ 2));

if isappdata(0,'CURDIFF_EXP')
    if transition_factor == 1
        wbp     = equi_ext.wbp;
        dwbpdt  = equi_ext.dwbpdt;
        drmdt   = equi_ext.drmdt;
        dphidt  = equi_ext.dphidt;
    else
        % energie magnetique
        wbp     = mu0 .* ip .^ 2 .* R ./ 4 .* li;
        wbp      = transition_factor * equi_ext.wbp   + (1 - transition_factor) * wbp;       
        if (evolution == 1) && (freebie == 1)
            dwbpdt  = z0dxdt_freebie(wbp,temps);
        else
            dwbpdt  = z0dxdt(wbp,temps);
        end
        dwbpdt     = transition_factor * equi_ext.dwbpdt   + (1 - transition_factor) * dwbpdt;       
        if (evolution == 1) && (freebie == 1)
            drmdt  = z0dxdt_freebie(rm,temps);
        else
            drmdt  = z0dxdt(rm,temps);
        end
        drmdt     = transition_factor * equi_ext.drmdt   + (1 - transition_factor) * drmdt;       
        % derivee temporelle de dphi/dt
        if (evolution == 1) && (freebie == 1)
            dphidt    = 2 .* pi .* rm .* Bt .* drmdt +  pi .* rm .^ 2 .* z0dxdt_freebie(Bt,temps);
        else
            dphidt    = 2 .* pi .* rm .* Bt .* drmdt +  pi .* rm .^ 2 .* z0dxdt(Bt,temps);
        end
        dphidt     = transition_factor * equi_ext.dphidt   + (1 - transition_factor) * dphidt;       
    end
else
    % energie magnetique
    wbp     = mu0 .* ip .^ 2 .* R ./ 4 .* li;
    %wbpf    = wbp;
    %wbpf(~isfinite(wbpf))= 0;
    %dwbpdt  = zdxdt(wbpf,temps);
    % wbp est impose par li et ip (donnees externes)
    % l'operateur est ici utilise comme filtre passe bas
    % c'est la constantes de temps de Mikkelsen pour les courant constant
    %taufilt  =  mu0 .* R ./ RR ./ 6;
    %[wbp_void,dwbpdt]  = zdwdt0(wbp, wbp ./ taufilt,taufilt,temps,Vp);

    % on repasse au calcul direct
    if (evolution == 1) && (freebie == 1)
	    %dwbpdt  = z0dxdt(wbp,temps);
	    %dwbpdt(3:end) = dwbpdt(2);
	    %dwbpdt = ones(size(temps)) * polyder(polyfit(temps,wbp,1));
	    %dwbpdt(:)   =  (wbp(end) - wbp(1)) ./ (temps(end) - temps(1));
	    dwbpdt  = z0dxdt_freebie(wbp,temps);

    else
	    dwbpdt  = z0dxdt(wbp,temps);
    end
    % derivee temporellle de rm
    %[rm_void,drmdt]  = zdwdt0(rm, rm ./ taufilt,taufilt,temps,NaN);
    if (evolution == 1) && (freebie == 1)
	    %drmdt  = z0dxdt(rm,temps);
	    %drmdt(3:end) = drmdt(2);		
	    %drmdt = ones(size(temps)) * polyder(polyfit(temps,rm,1));
	    drmdt  = z0dxdt_freebie(rm,temps);
	    %drmdt(:)   =  (rm(end) - rm(1)) ./ (temps(end) - temps(1));
    else
	    drmdt  = z0dxdt(rm,temps);
    end
    % derivee temporelle de dphi/dt
    if (evolution == 1) && (freebie == 1)
	    %dphidt    = 2 .* pi .* rm .* Bt .* drmdt +  pi .* rm .^ 2 .* z0dxdt(Bt,temps);
	    %dphidt = 2 .* pi .* rm .* Bt .* drmdt +  pi .* rm .^ 2 .* polyder(polyfit(temps,Bt,1));
	    dphidt    = 2 .* pi .* rm .* Bt .* drmdt +  pi .* rm .^ 2 .* z0dxdt_freebie(Bt,temps);
    else
	    dphidt    = 2 .* pi .* rm .* Bt .* drmdt +  pi .* rm .^ 2 .* z0dxdt(Bt,temps);
    end
end
%figure(22);clf;plot(temps,drmdt,'b',temps,dphidt,'r');drawnow
% securite anti derive
indbad = find(iboot >= (0.999 .* ip));
if ~isempty(indbad)
	iboot(indbad) = min(iboot_sc(indbad),ip(indbad));
end
indbad = find(iboot <= 0);
if ~isempty(indbad)
	iboot(indbad) = max(1,iboot_sc(indbad));
end

%save('loc2');

% debut de la boucle de limitation de vloop
% poh et joh
% fonctionnement a ip donne
ioh   = ip - iboot - icd;
poh   = RR .* ioh .* ioh;
if transitoire == 0
    vloop = sign(ioh) .* sqrt(RR .* poh);
else
    vloop = vloopin;
end
% expression complete de vloop pour la valuer initiale
%  vmes_an = sign(ioh) .* sqrt(RR .* poh) + dwbpdt ./ ip +  ...
%             drmdt .* (R .* mu0 .* ip) ./ (4 .* pi .^ 2 .* rm .^ 2) + ...
%             dphidt ./ qa;
%  vloop_an = sign(ioh) .* sqrt(RR .* poh) + dwbpdt ./ ip;
%  vloop(1) = 0.5 .* (vloop_an(1) + vloop(2));
	  
	  
asser       = zeros(size(vloop));
if (mode == 5) && (evolution == 1)
    asser(:)  =  -2;
elseif (mode == 5) && (evolution == 0)
    error('METIS internal error : the mode option.vloop = 5 can be used only in evolution run of METIS');
elseif (mode == 4) && (transitoire == 1)
    %  asser a 1 tous les temps
    if isfinite(t_switch)
        asser(temps >= t_switch) = -1;
    else
        asser(:)  =  -1;
    end
elseif (mode > 0) && (transitoire == 1)
    % temps ou c'est possible
    mode0       = (vloop <= vref)  & ((ioh ./ ip) < 0.5);
    mode1       = (vloop <= (1.2 .* vref + 0.05)) & ((ioh ./ ip) < 0.5);
    if evolution == 0
        if isfinite(t_switch)
            mode0 = mode0 & (temps >= t_switch);
            mode1 = mode1 & (temps >= t_switch);
        end
        mode0(1)    = 0;
        mode1(1)    = 0;
    end
    mode1       = cat(1,0,mode1(2:end));
    ind         = find( mode0 | mode1);
    ind(ind == 1) = [];
    asser(ind)  =  1;
    % fonctionnement a psi fixe au bord
    if (mode == 2) | (mode == 3) | (mode == 6)
        ioh(ind)   =  vref ./ RR(ind);
        vloop(ind) =  vref;
        poh(ind)   =  vref .^ 2 ./ RR(ind);
    else
        ioh(ind)    =  0;
        vloop(ind)  =  0;
        poh(ind)    =  0;
    end
end

%  % courant reel quel que soit le mode
%  ipv          = iboot + icd + ioh;
%  ipv(1)       = ip(1);
%  ipv(ipv <=0) = ip(ipv<=0);
% pour eviter un trop fort over-shoot au debut
%vloop(1) = ipv(1) .* RR(1);

% temps de diffusion du courant  (definition ref : D;R. Mikkelsen, Phys. of Fluids B1 (2), 1989, p 333-
tauj     = mu0 .* rm .^2 .* trapz(x,spr./eta,2) ./ trapz(x,spr,2);
tauj(1)  = tauip(1);
tauj = min(1e6,max(1e-6,tauj));
% temps de diffusion interne du courant (securite)
%tauli = min(1e6,max(1e-6,tauip));


if isappdata(0,'EQUILIBRIUM_EXP')
    r2i = equi_ext.r2i;
elseif isfield(profil,'r2i')
    r2i = profil.r2i;
else
    r2i = 1 ./ Raxe .^ 2;
end
if isappdata(0,'EQUILIBRIUM_EXP')
    grho2r2 = equi_ext.grho2r2;
elseif isfield(profil,'grho2r2')
    grho2r2 = profil.grho2r2;
else
    grho2r2 = 0.5 .* (1 + kx .^ 2) ./ Raxe .^ 2;
end

switch modeboot
case {3,4}
    jboot      = zsauter0d_HC(x,tep,tip,nep,nip,qp,zeffp,zion,Raxe,ftrap,epsi,dpsidx,F,modehp,ffit_ped,1,rmx,r2i,grho2r2,0);
    %jboots      = zsauter0d(x,tep,tip,nep,nip,qp,zeffp,zion,Raxe,ftrap,epsi,dpsidx,F,modehp,ffit_ped);
    %figure(21);clf;plot(x,jboots,'b',x,jboot,'r');drawnow
    
case {5,6}
   [~,jboot]      = z0etaboot_neofit(x,tep,tip,nep,nip,qp,zeffp,zion,Raxe,ftrap,epsi,dpsidx,F,modehp,ffit_ped);   
otherwise 
    jboot      = zsauter0d(x,tep,tip,nep,nip,qp,zeffp,zion,Raxe,ftrap,epsi,dpsidx,F,modehp,ffit_ped);
end
jboot          =    (bootmul * ones(size(x))) .* real(jboot) ./ Btr;

% ajout du courant asymetric si demander
if ((modeboot == 2) || (modeboot == 4) || (modeboot == 6))  && isfield(profil,'psi') && isfield(profil,'nep') && isfield(profil,'tep') && isfield(profil,'bpol')
    btot      = sqrt(profil.bpol .^ 2 + (profil.fdia ./ profil.Raxe) .^ 2);
    [jae,jbe] = zbootasymetric(x,profil.qjli,profil.epsi,profil.Raxe,btot,profil.nep,profil.tep,0,-1,Bt);
    [jai,jbi] = zbootasymetric(x,profil.qjli,profil.epsi,profil.Raxe,btot,profil.nip,profil.tip,meff,zeff,Bt);
    ja = jae + jai;
    jb = jbe + jbi;

    % il manque les effet torique dans la formule 
    % onsuppose qu'il sont les memes pour les 2 courants 
    %fna = trapz(x,spr .* jboot,2) ./ max(1,trapz(x,spr .* jb,2));
    %ja  = ja .* (fna * ve);
    warning off
    fna = jboot ./ max(1,jb);
    fna(jb < 1) = 0;
    fna(:,1) = fna(:,3);
    fna(:,2) = fna(:,3);
    ja       = ja .* fna;
    %figure(21);clf;plot(x,jboot,'b',x,ja,'r',x,jboot+ja,'k');drawnow
    jboot = jboot + ja;

end
% probleme de precision dans les gradients 
jboot(:,2)     = pchip(x([1,3:7]),jboot(:,[1,3:7]),x(2));
% normalisation  de jboot
if (~isfield(profil,'jboot') | (modeboot == 0))
    inteboot = trapz(x,spr .* jboot,2);
    inteboot(inteboot == 0) = 1;
    jboot    = jboot .* ((iboot ./  inteboot) * ve);
else
    rapis = trapz(x,spr .* abs(jboot),2) ./ max(1,ip);
    indbad = find((rapis > 1) | (rapis <= 0));
    if ~isempty(indbad)
        inteboot = trapz(x,spr .* jboot,2);
        inteboot(inteboot == 0) = 1;
        jboot_sc    = jboot .* ((iboot_sc ./  inteboot) * ve);
        jboot(indbad,:) = jboot_sc(indbad,:);
    end
    
end
jboot = real(jboot);
% pour la precision
iboot = trapz(x,spr .* jboot,2);

switch gaz
    case 1
        agaz    = 1 .* vt;
        zgaz    = 1 .*vt;
    case 2
        agaz    = 2 .* vt;
        zgaz    = 1 .* vt;
    case 3
        agaz   = 2.5 .* vt;
        zgaz    = 1 .* vt;
    case 5
        agaz   = 2.5 .* vt;
        zgaz    = 1.5 .* vt;
    case 11
        agaz   = 1 .* vt;
        zgaz    = 1 .* vt;
    otherwise
        agaz    = 4 .* vt;
        zgaz    = 2 .* vt;
end
   


% selon le mode LH
efficiency = NaN .* ip;
lh_numeric = 1;
if isappdata(0,'LH_SHAPE_EXP') & isfield(profil,'qjli');
	lhexp = getappdata(0,'LH_SHAPE_EXP');
	fplh  = interp1_ex(lhexp.temps,lhexp.jlh,temps,'nearest','extrap');
	fplh(abs(fplh) < eps) = sign(fplh(abs(fplh) < eps)) * eps; 
	fplh  = pchip(lhexp.x,fplh,profil.xli);
   	indnok = find(any(~isfinite(fplh),2));
   	indok  = find(all(isfinite(fplh),2));
   	fplh(indnok,:) = ones(length(indnok),1) * mean(fplh(indok,:),1);
   	jlh   = fplh .* ((ilh./  max(1,trapz(x,spr .* fplh,2))) * ve);

elseif isfield(profil,'fxdurplh') && all(dlh == 0) && (wlh == 0) && all(isfinite(profil.fxdurplh(:)))
   % cas TS experimental
   fplh  = profil.fxdurplh;
   indnok = find(any(~isfinite(fplh),2));
   indok  = find(all(isfinite(fplh),2));
   fplh(indnok,:) = ones(length(indnok),1) * mean(fplh(indok,:),1);
   jlh   = fplh .* ((ilh./  max(1,trapz(x,spr .* fplh,2))) * ve);
   %figure(31);clf;plot(x,jlh);drawnow;
elseif wlh == 0
   fplh    = exp(-(vt*x -(xlh * ve)).^ 2 ./ ((max(dlh,0.05)*ve) .^ 2));
   % c'est un probleme de definition entre ici et zerod_init
   %fplh    = exp(-(vt*x -(xlh * ve)).^ 2 ./ ((max(dlh,0.05)*ve) .^ 2) ./ 2);
   jlh    = fplh .* ((ilh./  trapz(x,spr .* fplh,2)) * ve);
   lh_numeric = 0;
else
	% facteur de propagation LHCD
	%qcyl          =  5 .* a .^ 2 .* Bt ./ (ip./1e6) ./ R .* ( 1 + (a ./ R) .^ 2);
	qcyl          =  5 .* a .^ 2 .* Bt ./ (ip./1e6) ./ R;
	if isfield(profil,'qjli')
		qbord = profil.qjli(:,end);
	else
		qbord = qeff;
	end
	switch upshiftmode
	case {'newmodel','newmodel + tail'}
	  upshift       = fupshift .* vt;
	otherwise
	  upshift       = fupshift .* max(eps,qbord ./ qcyl - 1);
	end
	%   factproplhcd  = max(eps,erf(upshift));
	%    figure(21);clf;
	%     subplot(2,1,1);
	%plot(temps,factproplhcd);
	%     subplot(2,1,2);
	%      plot(temps,upshift);  
	%     set(gca,'ylim',[0 3]); 
	%     drawnow
	[xvoid,fplh,xlh,dlh,efficiency,rapnegpos] = z0lhacc2lobes(freqlh.*1e9,npar0,wlh,agaz,zgaz,temps,x,nep,tep,...
					qp,Raxe,rmx,spr,vpr,Bt,max(1,plh),xlh,dlh,transitoire,directivity, ...
					friplh,upshift,0,upshiftmode,npar_neg);
        if isfield(profil,'zeff')
		fjlh = imag(fplh) ./ (5 + profil.zeff) .* (1 - profil.epsi .^ ((5 + profil.zeff) ./ 2 ./ (1 + profil.zeff)));
	else
		fjlh = imag(fplh);
	end
	fplh = real(fplh);

	if (transitoire == 0) | (evolution == 1)
		fplh = ones(size(temps)) * mean(fplh,1);
		efficiency = ones(size(temps)) *  mean(efficiency);
		rapnegpos = ones(size(temps))  *  mean(rapnegpos);
		xlh       = ones(size(temps))  *  mean(xlh);
		dlh       = ones(size(temps))  *  mean(dlh);
	end
	%
	%  profil de densite de courant
	%  cette facon de le decrire le rend compatible avec le module sans lobe negatif
	%  
	jlh    = fplh .* ((fplh > 0) + (rapnegpos * ve) .* (fplh <=0));
        indjlh = find(any(fjlh~=0,2));
	if ~isempty(indjlh)
		%figure(21);clf;subplot(2,1,1);plot(x,fjlh,'b',x,jlh,'r');
		jlh(indjlh,:) = fjlh(indjlh,:);
	end
	% normalisation sur le courant LH
	indbadlh = find(trapz(x,spr .* jlh,2) < eps);
	if ~isempty(indbadlh)		
		jlh(indbadlh,:) = ones(length(indbadlh),1) *mean(fplh .* (fplh > 0),1);
	end
	jlh    = jlh  .* ((ilh./  max(eps,trapz(x,spr .* jlh,2))) * ve);
	%
	% ajout de l'effet du Zeff
	efficiency = efficiency ./ (5 + zeff) .* (1 - (xlh .* a ./ R) .^ ((5 + zeff) ./ 2 ./ (1 + zeff))); 
    if isfield(profil,'zeff') && ~isempty(indjlh)
		%efficiency_full    =  trapz(x,spr .* fjlh ,2) ./ max(1,trapz(x,vpr .* profil.plh,2)) .* trapz(x, 2 .* pi .* profil.nep .* profil.Raxe,2);
		efficiency_full    =  trapz(x,spr .* fjlh ,2) ./ max(1,trapz(x,spr .* profil.plh,2)) .* trapz(x, profil.nep,2);
		%figure(21);subplot(2,1,2);plot(temps,efficiency_full,temps,efficiency,'r');drawnow
		efficiency(indjlh) =  efficiency_full(indjlh);
	end

	%
	% profil de source de chaleur eletronique LH
	%
	fplh   = abs(fplh);
	%     figure(161);clf
	%     plot(efficiency)
	%     drawnow
	%    keyboard
end
if any(~isfinite(efficiency))
    fplh_loc = fplh .* ((max(1,plh) ./ max(eps,trapz(x,vpr .* fplh,2))) * ve);
    efficiency_nonan = max(1,min(1e21,trapz(x,spr .* jlh,2) ./ trapz(x,spr .* fplh_loc,2) .* ne .* R));
    efficiency(~isfinite(efficiency)) = efficiency_nonan(~isfinite(efficiency));
end
% ECCD
if isappdata(0,'ECCD_SHAPE_EXP') & isfield(profil,'qjli');
	ecexp = getappdata(0,'ECCD_SHAPE_EXP');
	jeccd  = interp1_ex(ecexp.temps,ecexp.jeccd,temps,'nearest','extrap');
	jeccd  = pchip(ecexp.x,jeccd,profil.xli);
   	indnok = find(any(~isfinite(jeccd),2));
   	indok  = find(all(isfinite(jeccd),2));
   	jeccd(indnok,:) = ones(length(indnok),1) * mean(jeccd(indok,:),1);
	jeccd(jeccd == 0) = eps;
	%figure(21);plot(profil.xli,jeccd,'b',ecexp.x,ecexp.jeccd,'r');drawnow
	fpecrh  = max(1,interp1_ex(ecexp.temps,ecexp.peccd,temps,'nearest','extrap'));
	fpecrh  = pchip(ecexp.x,fpecrh,profil.xli);
   	indnok = find(any(~isfinite(fpecrh),2));
   	indok  = find(all(isfinite(fpecrh),2));
   	fpecrh(indnok,:) = ones(length(indnok),1) * mean(fpecrh(indok,:),1);	
	
elseif isfield(profil,'tep')
	% d'apres le calcul simple de largeur dans le Wesson
	% verifier sur de cas avec rema.
	dxm      = abs(vt*x - (abs(xeccd) * ve));
	mask     = (dxm == (min(dxm,[],2) * ve));
	tex      = sum(profil.tep .* mask,2) ./ max(1,sum(mask,2));
	vth2     = 2 .* 1.602176462e-19 .* tex ./9.10938188e-31;
	if width_ecrh > 0
	    deccd    = abs(width_ecrh / sqrt(2)) * vt;
	elseif width_ecrh < 0
	    deccd    = abs(width_ecrh) .* sqrt(vth2 .^ 2 ./ 4 ./ 2.99792458e8 .^ 4  + ...
		       min(1,abs(ieccd./max(1,pecrh))) .* vth2 ./ 2.99792458e8 .^ 2);	
	else
	    deccd    = sqrt(vth2 .^ 2 ./ 4 ./ 2.99792458e8 .^ 4  + ...
			min(1,abs(ieccd./max(1,pecrh))) .* vth2 ./ 2.99792458e8 .^ 2);	
	end
	deccd    = max(0.5 /length(x),deccd);
	jeccd    = exp(-(vt*x -(abs(xeccd) * ve)).^ 2 ./2 ./ (deccd * ve) .^ 2);  
        fpecrh   = jeccd;

else
	jeccd    = exp(-(vt*x -(abs(xeccd) * ve)).^ 2 .* 100);
        fpecrh   = jeccd;
end
jeccd    = jeccd .* ((ieccd./  trapz(x,spr .* jeccd,2)) * ve);

% FWCD
if isappdata(0,'ICRH_SHAPE_EXP') & isfield(profil,'qjli');
	icexp = getappdata(0,'ICRH_SHAPE_EXP');
	jfwcd  = max(1,abs(interp1_ex(icexp.temps,icexp.jfwcd,temps,'nearest','extrap')));
	jfwcd  = pchip(icexp.x,jfwcd,profil.xli);
   	indnok = find(any(~isfinite(jfwcd),2));
   	indok  = find(all(isfinite(jfwcd),2));
   	jfwcd(indnok,:) = ones(length(indnok),1) * mean(jfwcd(indok,:),1);
else
	% ref : F. Meo et al  (F. NGuyen)
	%jfwcd     = exp(-(vt*x - 0.1 ).^ 2 .* 25);
	bp1d       = rmx .* (Bt * ve) ./ Raxe ./ qp;
	b1d        = sqrt(bp1d.^ 2 + (((Bt .* Raxe(:,end)) * ve) ./ (Raxe + a * x)) .^ 2);
	% E^2 = P, J = E /eta
	jfwcd      = sqrt(tep .* nep  ./ b1d .^ 3  ./ (vpr+ eps)) ./ eta;
	jfwcd(:,1) = 2 .* jfwcd(:,2) - jfwcd(:,3);
end
jfwcd      =  jfwcd .* ((ifwcd./  trapz(x,spr .* jfwcd,2)) * ve);

%NBI
% external data take into account in zicd0
if isfield(profil,'jnbishape')
	jnbicd   = real(profil.jnbishape);
	jnbicd   = jnbicd .* ((real(inbicd)./  max(1,trapz(x,spr .* jnbicd,2))) * ve);
        if any(imag(inbicd))
	  jnbicd2   = imag(profil.jnbishape);
	  jnbicd2   = jnbicd2 .* ((imag(inbicd)./  max(1,trapz(x,spr .* jnbicd2,2))) * ve);
	  jnbicd   = jnbicd  + jnbicd2;
	end       
else
	jnbicd   = exp(-(vt*x -(xnbi * ve)).^ 2 ./ ((piqnbi * ve) .^ 2) ./ 2);
	jnbicd   = jnbicd .* (((real(inbicd) + imag(inbicd))./  max(1,trapz(x,spr .* jnbicd,2))) * ve);
end

% fusion
if (tune_frac == 0) && isfield(profil,'jfusshape');
	jfus     = profil.jfusshape;
	ifus_int = trapz(x,spr .* jfus,2);
	indok    = find(abs(ifus_int) > sqrt(eps));
	indnok   = find(abs(ifus_int) <= sqrt(eps));
	if ~isempty(indok)
	    jfus(indok,:)    = jfus(indok,:) .* ((ifus(indok)./ ifus_int(indok)) * ve);
	end
	if ~isempty(indnok)	
	  jfus(indnok,:)    =  0;
	end
else
  if isfield(profil,'jfusshape')
	jfus    = profil.jfusshape;
  else
	jfus    = z0j3(x,j0fus,xfus,jxfus);
  end
  jfus    = jfus .* (jfus > 0);
  jfus    = jfus .* ((ifus./  max(1,trapz(x,spr .* jfus,2))) * ve);
end

% courant runaway
switch runaway
case  0
	jrun = zeros(size(jfus));
	irun = zeros(size(ip));
otherwise
	% selon le signe du parametre breakdown
	if breakdown < 0
		epr_in = breakdown ./ 2 ./ pi ./ Raxe(1,end);		 
	else
		epr_in = breakdown;
	end
 
	if isappdata(0,'RUNAWAY_EXP')
	    runexp = getappdata(0,'RUNAWAY_EXP');
	    jrun  = interp1_ex(runexp.temps,runexp.jrun,temps,'nearest','extrap');
	    jrun  = pchip(runexp.x,jrun,x);
	    indnok = find(any(~isfinite(jrun),2));
	    jrun(indnok,:) = 0;
	elseif isfield(profil,'epar')
		% si epar est defini
		jrun = zjrun(profil.epar,nep,tep,zeffp,ftrap,tauj,epr_in);
	else
		jrun = zjrun(((vloop ./ 2 ./ pi) * ve) ./ Raxe,nep,tep,zeffp,ftrap,tauj,epr_in);
	end
   	% effet des electrons pieges (Y.Peysson)
   	%frp   = max(0,1 - 2.5 .* (a./R) * x);
	%jrun  = max(1e-9,jrun .* frp);
	irunx  = trapz(x,spr .* jrun,2);
	warning off
	jrun  = jrun .* ((irun ./ irunx) * ve);
	warning on
	jrun(irunx == 0,:) = 0;
	
	%figure(21);clf;plot(x,jrun);drawnow
end

%  % courant resistif
%  if mode == 0
%        ires     = ip  - iboot - icd;
%  else
%        ires     = ipv - iboot - icd;
%  end
%  
%  % forme jres
%  jres     = 1 ./ eta ./ Raxe;
%  jres     = jres .* (((ires-irun) ./  max(1,trapz(x,spr .* jres,2))) * ve);
%  
%  
% total
jni      = jboot + jlh + jfwcd + jeccd + jnbicd + jfus+jrun;
%  jli      = jres  + jni;
%  
%  
%  % trou de courant
%  jli      = jli .* (jli >0) + eps;
%  if mode == 0
%     jli      = jli .* ((ip./  trapz(x,spr .* jli,2)) * ve);
%  else
%     jli      = jli .* ((ipv./  trapz(x,spr .* jli,2)) * ve);
%  end

%save('loc3');

% calcul d'une aproximation de bpol like pour li
% correction de l'elongation ITER FDR p 116
% la pression
%  pep      = ee .* (nep  .* tep + nip .* tip);
%  ptot     = pep .* max(1,((real(wtotal) ./ (3/2) ./ trapz(x,pep .* vpr,2)) * ve));
%dpdx     = min(-eps,pdederive(vt*x, ptot,0,2,2,1));
%dpdr     = dpdx ./ (a*ve);  % dans le plan median
% nouvelle formulation plus proche de la realite
if isappdata(0,'EQUILIBRIUM_EXP') 
    Raxe        = equi_ext.Raxe;
	epsi        = equi_ext.epsi;
    
elseif isfield(profil,'bpol')
	% calcul du decentrement de Shafranov (formule dans Lao Phys of fluids vol 24 (8) 1981 ,pp 1436)
        xs          = vt *x;
	% kx et ka donne des resultats similaira, on garde kx.
	kax         = kx;
	% il faut utiliser cette definition de bpol qui est dans l'article ? non
	bpolm       = profil.bpol;
	bpolm(:,1)  = 0;
	ep          = (a./R) * ve;
        % self interne glissante
	lim         = 4 .* (kax .^ 2 + 1) ./ (3.* kax .^ 2 + 1) ./ max(eps,xs .^ 2) ./ max(eps,bpolm .^ 2) .* ...
	              cumtrapz(x, xs .* bpolm .^ 2,2);
        lim(:,1)    = 2 .* lim(:,2) - lim(:,3);
	% beta moyen glissant
	ptm          = 2 ./ max(eps,xs .^ 2) .* cumtrapz(x, xs .* ptot,2);
        ptm(:,1)     = 2 .* ptm(:,2) - ptm(:,3);
        % lambda
	inter        = 2 .* (kax .^ 2 + 1) ./ (3.* kax .^ 2 + 1) .* 2 .* mu0 .* (ptm - ptot) ./ max(eps,bpolm .^ 2);
	inter(:,1)   = 2 .* inter(:,2) - inter(:,3);
	lam          = inter + 0.5 .* lim + 0.5 .* (kax .^ 2 - 1) ./ (3.* kax .^ 2 + 1) - 1;
	%deltax       = ep .* cumtrapz(x, xs .*  (1 + lam),2);
 	deltax       = cumtrapz(x, max(eps,min(pi,ep .* xs .*  (1 + lam))),2);
        deltax       = deltax(:,end) * ve -deltax;
	% correction final (la formule de ep dans l'article fait intervenir Raxe)
	%deltax       = min((a*ve)/2,max(0,(a*ve) .* ((R./(R + deltax(:,1))) * ve) .*  deltax .* (1 - ep)));
	deltax       = max(0,(a*ve) .* ((R./(R + deltax(:,1))) * ve) .*  deltax ./ (1 + dx(:,end) * ve));
	% recherche des temps ou le calcul n'est pas valide
	indnok = find(any(deltax > (a*ve)/4,2));
	if ~isempty(indnok)
		deltax(indnok,:) = d0lim(indnok) * (1 - x .^2);
	end
        % couplage avec FREEBIE
	if ~isempty(d0_ext) && all(isfinite(d0_ext)) && all(d0_ext > 0)
	      %deltax_mem = deltax;
	      deltax = deltax .* ((d0_ext ./ max(deltax,[],2)) * ve);
	      %figure(19);clf;plot(x,deltax_mem,'b',x,deltax,'r');drawnow
	end
  
	% nouvelle donnees
	Raxe        = (R * ve  + deltax);
	epsi        = max(eps,(a * x)) ./ Raxe;


end

% appel du module de diffusion du courant rapide
if isfield(profil,'fdia')
	Flast = profil.fdia;
else
	Flast = ((Bt .* R) * ones(1,size(Raxe,2)));
end
if isfield(profil,'kx')
	kx_up = profil.kx;
else
	kx_up = kx;
end
if isfield(profil,'dx')
	dx_up = profil.dx;
else
	dx_up = dx;
end
if isfield(profil,'Raxe')
	Raxe_up = profil.Raxe;
else
	Raxe_up = Raxe;
end

%save('loc4');
% debut probleme de precision numerique
if isappdata(0,'CURDIFF_EXP') 
                [psi,dpsidt,qjli,jli,jres,epar,fdia,bpol,lif,rmx,vpr_tor, ...
		grho2r2,r2i,ri,C2,C3,ej,ipout,grho,grho2,jeff,spr_tor,phi_tor, ...
                dphidx_tor,difcurconv,phiplasma,indice_inv,poynting,jgs,df2dpsi,dptotdpsi] =  ...
		z0curdiff_external(temps,x,eta,jni,ptot,evolution,qdds,s1crit,ddsmode,w1,epsq,q0_dds_trig,betap1crit);
    
    if transition_factor < 1
        [psi_,dpsidt_,qjli_,jli_,jres_,epar_,fdia_,bpol_,lif_,rmx_,vpr_tor_, ...
            grho2r2_,r2i_,ri_,C2_,C3_,ej_,ipout_,grho_,grho2_,jeff_,spr_tor_,phi_tor_, ...
            dphidx_tor_,difcurconv_,phiplasma_,indice_inv_,poynting_,jgs_,df2dpsi_,dptotdpsi_] =  ...
            z0curdiff(temps,x,eta,jni,ptot + sqrt(-1) .* hollow,Sp,kx_up,Raxe_up,Bt,a,dx_up, ...
            ip,vloop,li,asser,Vp,Flast,qdds + sqrt(-1) .*  s1crit,peri,Rsepa,Zsepa,amorti,drmdt,qeff,profil.psi,profil.phi,difcurconv,laochange,evolution + sqrt(-1) * freebie, ...
            edge_flux,mode_expo_inte,cronos_regul,ddsmode,w1,epsq,[],[],q0_dds_trig,betap1crit);
        
        psi     = transition_factor * psi   + (1 - transition_factor) * psi_;  
        dpsidt  = transition_factor * dpsidt   + (1 - transition_factor) * dpsidt_;  
        qjli    = transition_factor * qjli   + (1 - transition_factor) * qjli_;  
        jli     = transition_factor * jli   + (1 - transition_factor) * jli_;  
        jres    = transition_factor * jres   + (1 - transition_factor) * jres_;  
        epar    = transition_factor * epar   + (1 - transition_factor) * epar_;  
        fdia    = transition_factor * fdia   + (1 - transition_factor) * fdia_;  
        bpol    = transition_factor * bpol   + (1 - transition_factor) * bpol_;  
        lif     = transition_factor * lif   + (1 - transition_factor) * lif_;  
        rmx     = transition_factor * rmx   + (1 - transition_factor) * rmx_;  
        vpr_tor = transition_factor * vpr_tor   + (1 - transition_factor) * vpr_tor_;  
        grho2r2 = transition_factor * grho2r2   + (1 - transition_factor) * grho2r2_;  
        r2i     = transition_factor * r2i   + (1 - transition_factor) * r2i_;  
        ri      = transition_factor * ri   + (1 - transition_factor) * ri_;  
        C2      = transition_factor * C2   + (1 - transition_factor) * C2_;  
        C3      = transition_factor * C3   + (1 - transition_factor) * C3_;  
        ej      = transition_factor * ej   + (1 - transition_factor) * ej_;  
        ipout   = transition_factor * ipout   + (1 - transition_factor) * ipout_;  
        grho    = transition_factor * grho   + (1 - transition_factor) * grho_;  
        grho2   = transition_factor * grho2   + (1 - transition_factor) * grho2_;  
        jeff    = transition_factor * jeff   + (1 - transition_factor) * jeff_;
        spr_tor = transition_factor * spr_tor   + (1 - transition_factor) * spr_tor_;  
        phi_tor = transition_factor * phi_tor   + (1 - transition_factor) * phi_tor_;  
        dphidx_tor = transition_factor * dphidx_tor   + (1 - transition_factor) * dphidx_tor_;
        difcurconv = ceil(transition_factor * difcurconv   + (1 - transition_factor) * difcurconv_);
        phiplasma  = transition_factor * phiplasma  + (1 - transition_factor) * phiplasma_;  
        indice_inv = ceil(transition_factor * indice_inv   + (1 - transition_factor) * indice_inv_);
        poynting   = transition_factor * poynting   + (1 - transition_factor) * poynting_;
        jgs        = transition_factor * jgs   + (1 - transition_factor) * jgs_;
        df2dpsi    = transition_factor * df2dpsi   + (1 - transition_factor) * df2dpsi_;
        dptotdpsi  = transition_factor * dptotdpsi   + (1 - transition_factor) * dptotdpsi_;
        
     end
elseif isfield(profil,'psi') 
        if (evolution == 1) && any(asser == -2)
		[psi,dpsidt,qjli,jli,jres,epar,fdia,bpol,lif,rmx,vpr_tor, ...
		grho2r2,r2i,ri,C2,C3,ej,ipout,grho,grho2,jeff,spr_tor,phi_tor, ...
                dphidx_tor,difcurconv,phiplasma,indice_inv,poynting,jgs,df2dpsi,dptotdpsi] =  ...
		z0curdiff(temps,x,eta,jni,ptot + sqrt(-1) .* hollow,Sp,kx_up,Raxe_up,Bt,a,dx_up, ...
		ip,vloop,li,asser,Vp,Flast,qdds + sqrt(-1) .*  s1crit,peri,Rsepa,Zsepa,amorti,drmdt,qeff,profil.psi,profil.phi,difcurconv,laochange,evolution + sqrt(-1) * freebie, ...
                edge_flux,mode_expo_inte,cronos_regul,ddsmode,w1,epsq,flux_ip,L_ext,q0_dds_trig,betap1crit);
	elseif evolution == 1
		[psi,dpsidt,qjli,jli,jres,epar,fdia,bpol,lif,rmx,vpr_tor, ...
		grho2r2,r2i,ri,C2,C3,ej,ipout,grho,grho2,jeff,spr_tor,phi_tor, ...
                dphidx_tor,difcurconv,phiplasma,indice_inv,poynting,jgs,df2dpsi,dptotdpsi] =  ...
		z0curdiff(temps,x,eta,jni,ptot + sqrt(-1) .* hollow,Sp,kx_up,Raxe_up,Bt,a,dx_up, ...
		ip,vloop,li,asser,Vp,Flast,qdds + sqrt(-1) .*  s1crit,peri,Rsepa,Zsepa,amorti,drmdt,qeff,profil.psi,profil.phi,difcurconv,laochange,evolution + sqrt(-1) * freebie, ...
                edge_flux,mode_expo_inte,cronos_regul,ddsmode,w1,epsq,[],[],q0_dds_trig,betap1crit);
	else
		[psi,dpsidt,qjli,jli,jres,epar,fdia,bpol,lif,rmx,vpr_tor, ...
		grho2r2,r2i,ri,C2,C3,ej,ipout,grho,grho2,jeff,spr_tor,phi_tor, ...
                dphidx_tor,difcurconv,phiplasma,indice_inv,poynting,jgs,df2dpsi,dptotdpsi] =  ...
		z0curdiff(temps,x,eta,jni,ptot + sqrt(-1) .* hollow,Sp,kx_up,Raxe_up,Bt,a,dx_up, ...
		ip,vloop,li,asser,Vp,Flast,qdds + sqrt(-1) .*  s1crit,peri,Rsepa,Zsepa,amorti,drmdt,qeff,[],profil.phi,difcurconv,laochange,evolution + sqrt(-1) * freebie, ...
                edge_flux,mode_expo_inte,cronos_regul,ddsmode,w1,epsq,[],[],q0_dds_trig,betap1crit);	
	end
else
	[psi,dpsidt,qjli,jli,jres,epar,fdia,bpol,lif,rmx,vpr_tor, ...
	grho2r2,r2i,ri,C2,C3,ej,ipout,grho,grho2,jeff,spr_tor,phi_tor,...
        dphidx_tor,difcurconv,phiplasma,indice_inv,poynting,jgs,df2dpsi,dptotdpsi] =  ...
	z0curdiff(temps,x,eta,jni,ptot + sqrt(-1) .* hollow,Sp,kx_up,Raxe_up,Bt,a,dx_up, ...
		ip,vloop,li,asser,Vp,Flast,qdds + sqrt(-1) .*  s1crit,peri,Rsepa,Zsepa,amorti,drmdt,qeff,[], ...
		rmx .^ 2 .* pi .* ((Bt .* R) * ones(1,size(rmx,2))) ,difcurconv,laochange,evolution + sqrt(-1) * freebie, ...
                edge_flux,mode_expo_inte,cronos_regul,ddsmode,w1,epsq,[],[],q0_dds_trig,betap1crit);
end

% conservation de la valeur initiale
if evolution == 0
	lif(1) = li(1);
end
lif    = max(0.1,min(10,lif));

% changement de coordonnee + normalisation
vpr  =  vpr_tor .* pdederive(x,rmx,0,2,2,1);
spr  =  spr_tor .* pdederive(x,rmx,0,2,2,1);
fact = (Vp ./ max(eps,trapz(x,vpr,2))) * ve;
vpr  = vpr .* fact;
fact = (Sp ./ max(eps,trapz(x,spr,2))) * ve;
spr  = spr .* fact;

if isappdata(0,'CURDIFF_EXP')
    vloop = -(2*pi) .* equi_ext.dpsidt(:,end);
    poh   = equi_ext.pohm;
    ioh   = min(2* equi_ext.ipout,max(-equi_ext.ipout,trapz(x,spr .* jres,2)));
    ipar  = trapz(x,spr .* jeff,2);
    if transition_factor < 1
        vlmax    = 1e2;
        vloop    = transition_factor * vloop + ...
                   (1 - transition_factor) * max(-vlmax,min(vlmax,- 2 .* pi .* dpsidt(:,end)));
        % ioh final
        iplimite = max(max(abs(ip),abs(ipout)),abs(vloop ./ RR));
        iplimite_alt = max(abs(ip),abs(vloop ./ RR));
        iplimite(~isfinite(iplimite)) = iplimite_alt(~isfinite(iplimite));
        iplimite(~isfinite(iplimite)) = abs(ip(~isfinite(iplimite)));
        pohlimite = iplimite .^ 2 .* RR;
        ioh     = transition_factor * ioh + ...
                  (1 - transition_factor) * max(-iplimite,min(2 .* iplimite,trapz(x,spr .* jres,2)));
        % courant // a B
        ipar    = transition_factor * ipar +  (1 - transition_factor) * trapz(x,spr .* jeff,2);
        
        % puissance ohmique final
        poh        = transition_factor * poh + (1 - transition_factor) * min(pohlimite,max(-pohlimite,trapz(x,vpr .* ej,2)));
    end
    
elseif evolution == 1
    vlmax    = 1e2;
    vloop    = max(-vlmax,min(vlmax,- 2 .* pi .* dpsidt(:,end))); 
    % ioh final
    iplimite = max(max(abs(ip),abs(ipout)),abs(vloop ./ RR));
    iplimite_alt = max(abs(ip),abs(vloop ./ RR));
    iplimite(~isfinite(iplimite)) = iplimite_alt(~isfinite(iplimite));
    iplimite(~isfinite(iplimite)) = abs(ip(~isfinite(iplimite)));
    pohlimite = iplimite .^ 2 .* RR; 
    ioh     = max(-iplimite,min(2 .* iplimite,trapz(x,spr .* jres,2)));
    % courant // a B
    ipar    = trapz(x,spr .* jeff,2);

    % puissance ohmique final
    %poh        = min(pohlimite,max(1,trapz(x,vpr .* ej,2)));
    poh        = min(pohlimite,max(-pohlimite,trapz(x,vpr .* ej,2)));
else
    % ioh final
    ioh     = max(-abs(ip),min(2 .* abs(ip),trapz(x,spr .* jres,2)));
    % courant // a B
    ipar    = trapz(x,spr .* jeff,2);

    % puissance ohmique final
    vlmax    = 1e2;
    pohlimite = min(abs(vlmax) .* ip,2 .* RR .* ip .* ip);
    poh       = min(pohlimite,max(-pohlimite,trapz(x,vpr .* ej,2)));
    % vloop final :
    vloop    = max(-vlmax,min(vlmax,- 2 .* pi .* dpsidt(:,end))); 
end


if evolution == 0
	switch runaway
	case 0
		vloop(1) = pchip(cat(1,2 .* temps(1) - temps(2),temps(2:end)), ...
                 	   cat(1,sign(ioh(1)) .* sqrt(RR(1) .* abs(poh(1))),vloop(2:end)),temps(1));
	otherwise
		if transitoire == 0
                	vloop(1) = pchip(cat(1,2 .* temps(1) - temps(2),temps(2:end)), ...
                 	   cat(1,sign(ioh(1)) .* sqrt(RR(1) .* abs(poh(1))),vloop(2:end)),temps(1));
		else
			vloop(1) = vloopin(1);
		end
	end
end
if (evolution == 1)
    vmes     = vloop + drmdt .* (R .* mu0 .* ipout) ./ (4 .* pi .^ 2 .* rm .^ 2) + dphidt ./ qa;
else
    vmes     = vloop + drmdt .* (R .* mu0 .* ip) ./ (4 .* pi .^ 2 .* rm .^ 2) + dphidt ./ qa;    
end
%save('loc5');
% le probleme de precision numerique est avant

%verification
%  psiedge = - cumtrapz(temps,vloop ./ 2 ./ pi);
%  psiedge = - cat(1,0,cumsum(vloop(2:end) .* diff(temps))) ./ 2 ./ pi;
%  psiedge2 = cumtrapz(temps,dpsidt(:,end));
%  psiedge2 = cat(1,0,cumsum(dpsidt(2:end,end) .* diff(temps)));
%  figure(21);clf
%  plot(temps,psiedge,'b',temps,psi(:,end) - psi(1,end),'or',temps,psiedge2,'+m')
%  drawnow

% energie magnetique
if isappdata(0,'CURDIFF_EXP')
  wbp    = equi_ext.wbp;
  dwbpdt = equi_ext.dwbpdt;
  if transition_factor < 1
      wbp        = transition_factor * wbp + (1 - transition_factor) * mu0 .* ipout .^ 2 .* R ./ 4 .* li;
      if freebie == 1
          dwbpdt  = transition_factor * dwbpdt + (1 - transition_factor) * z0dxdt_freebie(wbp,temps);
      else
          dwbpdt  = z0dxdt(wbp,temps);
      end      
  end  
elseif (evolution == 1)
   wbp     = mu0 .* ipout .^ 2 .* R ./ 4 .* li;
   if freebie == 1
	dwbpdt  = z0dxdt_freebie(wbp,temps);
   else
	dwbpdt  = z0dxdt(wbp,temps);
   end
end


% constante de temps de Mikkelsen pour dpsidt donn??
lext     = log(8 .* R ./ a) - 2;
taufree  = mu0 .* R .* (lext + lif ./ 2)  ./ RR ./ 2;
tauvloop = mu0 .* R .* lif ./ RR ./ 1.8;
tauipc   = mu0 .* R  ./ RR / 6;
modetau  = (vloop - RR .* ioh) ./ ip ./ RR;
neutral  = 0.5;
modefree  = modetau <= - neutral;
tauip    =  (~modefree) .* tauvloop + modefree .* taufree;
tauip    = min(1e6,max(1e-6,tauip));

% securite
indnok  = find(qjli(:,end) == min(qjli,[],2));
if ~ isempty(indnok)
    qjli(indnok,:) = qpr(indnok,:);
end
indnok  = find(~isfinite(qjli));
if ~ isempty(indnok)
    %warning('NaN in q profile !');
    fprintf('Q');
    qjli(indnok) = qpr(indnok);
end
qa          = qjli(:,end);
psin = (psi - psi(:,1) * ve) ./ ((psi(:,end) - psi(:,1)) * ve);
for k=1:size(qjli,1)
    if any(diff(psin(k,:))<=0)
        psin_ = cat(2,0,cumsum(max(eps*abs(max(psin(:)) - min(psin(:))),abs(diff(psin(k,:))))));
        psin_ = psin_ - max(psin_) + 1;
        q95(k)         = interp1(psin_,qjli(k,:),0.95,'linear','extrap');
        %warning('non monotonic Psi_n !');
    else
        q95(k)         = interp1(psin(k,:),qjli(k,:),0.95,'linear','extrap');
    end
end

if cronos_regul > 2
    q0          = min(qa - 0.5,max(sqrt(eps),qjli(:,1)));
    qmin        = min(q0,max(sqrt(eps),min(qjli,[],2)));
    % securite pour les point trop bas
    qjli        =  max(sqrt(eps),qjli);
else
    q0          = min(qa - 0.5,max(0.5,qjli(:,1)));
    qmin        = min(q0,max(0.5,min(qjli,[],2)));
    % securite pour les point trop bas
    qjli        =  max(0.5,qjli);
end

if isappdata(0,'CURDIFF_EXP')
    if transition_factor < 1
        qjli = transition_factor * equi_ext.qjli + (1 - transition_factor) * qjli;
        qa   = qjli(:,end);
        q0   = transition_factor * equi_ext.q0 + (1 - transition_factor) * q0;
        qmin = transition_factor * equi_ext.qmin + (1 - transition_factor) * qmin;
        q95  = transition_factor * equi_ext.q95 + (1 - transition_factor) * q95;
    else
        qjli = equi_ext.qjli;
        qa   = qjli(:,end);
        q0   = equi_ext.q0;
        qmin = equi_ext.qmin;
        q95  = equi_ext.q95;
    end
end

% calcul de hitb (nouvelle version inspiree de GLF23, ce calcul relatif n'est pas sensible a l'effet de  vei)
% ref : R.E. Waltz et al,Phys. Plasma. 4 (7) 1997
% ref : T. Tala et al, NF 46 (2006) p 548-561

% effet de la rotation
if isfield(profil,'web')
    %gitg    = sqrt(2 .* 1.602176462e-19 .* profil.tip ./ (agaz * ve) ./ 1.6726485e-27) ./ (profil.Raxe + a * x);
    % expression alternative (X. Garbet) pour la normalisation
    % ref : R.E. Waltz PoP 1 (1994), p 2229
    gitg    = sqrt(1.602176462e-19 .* (profil.tip + profil.zeff .* profil.tep)./1.6726485e-27 ./ (agaz * ve)) ./ ...
        (profil.rmx(:,end) * ve);
    srotg   = profil.web ./ gitg;
else
    srotg   = zeros(size(qjli));
end
% calcul du shear magnetique
s        = pdederive(vt*x, qjli,0,2,2,1) ./ qjli .* (vt*x);
% alpha
if hollow == 1
    dpdr     = pdederive(vt*x, ptot,0,2,2,1) ./ (rmx(:,end) * ve);
else
    dpdr     = min(-eps,pdederive(vt*x, ptot,0,2,2,1)) ./ (rmx(:,end) * ve);
end
al       = - 2 .* mu0 .* (R *ve) .* qjli .^ 2 .* dpdr ./ (Bt*ve) .^ 2;

% profil de q sans itb
qmx        = min(qjli,[],2);
alx        = (qjli - qmx * ve) ./ max(0.5,(qjli(:,end) - qmx) * ve);
indok      = find((x > 0) & (x < 1));
warning off
alq        = log(alx(:,indok)) ./ log(vt * x(indok));
warning on
mask       = isfinite(alq);
alq(~isfinite(alq)) = 0;
alq        = max(1,sum(alq,2) ./ max(1,sum(mask,2)));
qdb        = qmx * ve  + (vt * x) .^ (alq  * ve) .* ((qjli(:,end) - qmx) * ve);
sdb        = pdederive(vt*x, qdb,0,2,2,1) ./ qdb .* (vt*x);


% facteur de reduction max (rapport de gradient critique ITM/TEM avec s et s null ref :Fourment)
% le shear moyen est suppose = 1
% prendre le vrai profil de q met une contre reaction sur les fort q au centre en reverse shear
rapk       = 1.1 ./(1.1 + 1.4  + 1.9 ./ qjli);
% limite par mode ballooning (sans effet de shear) 
if isfield(profil,'fdia')
      % ref :Onjun 2002 (valide pour s<3 )
      % (Wesson montre que en dessous de s = 0.5 alpha_c ne change plus)
      % normalement cette formulation n'est pas valide pour q < 1
      % il s'agit d'une simple limitation du gradient dans la barriere
      % cette modification de la formule permet de limiter fortement l'effet des barriere en dessous de la valeur optimale de q=2
      balcrit   =  0.8 .* max(0.5,min(3,s)) .* (profil.fdia .* profil.ri) .^ 2 ./ mu0 .* profil.ri ./ (qjli .^ 2 + 1) .* (5/4);
      %balcrit   =  0.8 .* max(0.5,min(3,s)) .*(1 + kx .^ 2 .* (1 + 5.* dx .^ 2)) ./ 2 .* (profil.fdia .* profil.ri) .^ 2 ./ mu0 .* profil.ri ./ (qjli .^ 2 + 1) .* (5/4);
      %balcrit   =  (profil.fdia .* profil.ri) .^ 2 ./ mu0 .* profil.ri ./ qjli .^ 2; % pour alpha_c =1
      %dpthdr     = min(-eps,pdederive(vt*x, pep,0,2,2,1)) ./ (rmx(:,end) * ve);
      if hollow == 1
            dpthdr     = pdederive(vt*x, profil.ptot,0,2,2,1) ./ (rmx(:,end) * ve);
     else
            dpthdr     = min(-eps,pdederive(vt*x, profil.ptot,0,2,2,1)) ./ (rmx(:,end) * ve);
      end
      %rapk      =  max(0.1,tanh(abs(dpthdr ./ balcrit)));
      %rapk      =  min(1,max(0.1,abs(dpthdr ./ balcrit) >= 1));
      %rapk      =  max(0.1,tanh(30 .* (abs(dpthdr ./ balcrit)-1)));
      rapk      =  max(0.01,abs(dpthdr ./ max(eps,balcrit)) - (itb_slope_max - 1));
end


% fonction de shear generalisee
% coef defavorable de Tala
% fsg      = 0.2 + s - min(0,(3/5)  .* al - 3/8) - (sitb > 1) .* 0.5 .* srotg;
% coef de Tala (les 2 references + experience TS)
fsg      = -0.14 + s - min(0,(3/5)  .* al - 3/8) - (sitb > 1) .* srotg;
% fonction de modification du transport
fs      = max(rapk,1 + itb_sensitivity .* fsg);
% fonction limite seuil
fsdb    = max(rapk,1 + itb_sensitivity .* max(0,fsg));
fs(:,end-1:end) = fsdb(:,end-1:end);  % securite au bord

%figure(21);clf;semilogy(x,fsdb,'b',x,fs,'r',x,rapk,'g');drawnow
%figure(21);clf;plot(x,s,'r',x,sdb,'b');drawnow
%  	figure(21);clf
%  	subplot(2,2,1)
%  	plot(x,rapk);
%  	subplot(2,2,2)
%  	plot(x,fs ./ fsdb);
%  	subplot(2,2,3)
%  	plot(x,min(0,(3/5)  .* al - 3/8),'b',x,0.5.* srotg,'r');
%  	subplot(2,2,4)
%  	plot(x,fsdb,'b',x,fs,'r');
%   	drawnow

% calcul du profil de depot dont depent le piquage
if isappdata(0,'LH_SHAPE_EXP') & isfield(profil,'qjli');
	lhexp = getappdata(0,'LH_SHAPE_EXP');
	fplh  = max(1,interp1_ex(lhexp.temps,lhexp.plh,temps,'nearest','extrap'));
	fplh  = pchip(lhexp.x,fplh,profil.xli);
   	indnok = find(any(~isfinite(fplh),2));
   	indok  = find(all(isfinite(fplh),2));
   	fplh(indnok,:) = ones(length(indnok),1) * mean(fplh(indok,:),1);
end



% puissance lh
if lh_numeric == 1
  fplh    = abs(fplh) .* ((plh./  max(1,trapz(x,vpr .* abs(fplh),2))) * ve);
else
  fplh    = abs(fplh) .* ((plh./  trapz(x,vpr .* abs(fplh),2)) * ve);
end
%  if isfield(profil,'plh')
%  	fplh    = (profil.plh + fplh) ./ 2;
%  end



if isappdata(0,'ECCD_SHAPE_EXP') & isfield(profil,'qjli');
	ecexp = getappdata(0,'ECCD_SHAPE_EXP');
	fpecrh  = max(1,interp1_ex(ecexp.temps,ecexp.peccd,temps,'nearest','extrap'));
	fpecrh  = pchip(ecexp.x,fpecrh,profil.xli);
   	indnok = find(any(~isfinite(fpecrh),2));
   	indok  = find(all(isfinite(fpecrh),2));
   	fpecrh(indnok,:) = ones(length(indnok),1) * mean(fpecrh(indok,:),1);
end
% puissance ecrh
fpecrh   = fpecrh .* ((pecrh./  trapz(x,vpr .* fpecrh,2)) * ve);
%figure(21); plot(x,fpecrh,'b',x,fplh,'.r');drawnow

% ref : F. Meo et al  (F. NGuyen)
%fpfweh   = exp(-(vt*x - 0.05 ).^ 2 .* 33);
% absorbtion proportionnelle a beta et propagation 1D et volume chauffer
b1d      = sqrt(bpol .^ 2 + (fdia ./ (Raxe + a * x)) .^ 2);
fpfweh   = pep ./ b1d .^ 3  ./ (vpr+ eps);
fpfweh(:,1) = 2 .* fpfweh(:,2) - fpfweh(:,3);
fpfweh   = fpfweh .* ((pfweh./  trapz(x,vpr .* fpfweh,2)) * ve);

% NBI
% external data take into account in zicd0
if isfield(profil,'nbishape_el')
    fpnbi_el    = real(profil.nbishape_el);
    fpnbi_ion   = real(profil.nbishape_ion);
    fpnbi_el    = fpnbi_el .* ((max(0,real(pnbi) - real(pion_nbi))./  max(1,trapz(x,vpr .* fpnbi_el,2))) * ve);
    fpnbi_ion   = fpnbi_ion .* ((max(0,real(pion_nbi))./  max(1,trapz(x,vpr .* fpnbi_ion,2))) * ve);
    fpnbi       = fpnbi_ion + fpnbi_el;
    fpnbi       = fpnbi .* ((real(pnbi)./  max(1,trapz(x,vpr .* fpnbi,2))) * ve);
    
    if any(imag(pnbi))
        fpnbi_el2  = imag(profil.nbishape_el);
        fpnbi_ion2 = imag(profil.nbishape_ion);
        fpnbi_el2    = fpnbi_el2 .* ((max(0,imag(pnbi) - imag(pion_nbi))./  max(1,trapz(x,vpr .* fpnbi_el2,2))) * ve);
        fpnbi_ion2   = fpnbi_ion2 .* ((max(0,imag(pion_nbi))./  max(1,trapz(x,vpr .* fpnbi_ion2,2))) * ve);
        fpnbi2       = fpnbi_ion2 + fpnbi_el2;
        fpnbi2       = fpnbi2 .* ((imag(pnbi)./  max(1,trapz(x,vpr .* fpnbi2,2))) * ve);
        
        fpnbi_el     = fpnbi_el  + fpnbi_el2;
        fpnbi_ion    = fpnbi_ion + fpnbi_ion2;
        fpnbi        = fpnbi     + fpnbi2;
        
    end
else
    fpnbi_el    = exp(-(vt*x -(xnbi * ve)).^ 2 ./ ((piqnbi * ve) .^ 2) ./ 2);
    fpnbi_ion   = fpnbi_el;
    fpnbi_el    = fpnbi_el .* ((max(0,real(pnbi) + imag(pnbi) - real(pion_nbi) - imag(pion_nbi))./   ...
        max(1,trapz(x,vpr .* fpnbi_el,2))) * ve);
    fpnbi_ion   = fpnbi_ion .* ((max(0,real(pion_nbi) + imag(pion_nbi))./  max(1,trapz(x,vpr .* fpnbi_ion,2))) * ve);
    fpnbi       = fpnbi_ion + fpnbi_el;
    fpnbi       = fpnbi .* (((real(pnbi) + imag(pnbi))./  max(1,trapz(x,vpr .* fpnbi,2))) * ve);
end

% fusion
if isfield(profil,'salf')
	fpfus    = profil.salf;
else
	sigmavdt = zperes(tep .* (tite * ve) + 13.6);
	fpfus    = nep .^ 2   .* sigmavdt;
end
fpfus    = fpfus .* ((pfus./  max(1,trapz(x,vpr .* fpfus,2))) * ve);



% icrh
if isappdata(0,'ICRH_SHAPE_EXP') & isfield(profil,'qjli');
	% il faut aussi modifier dans zero1t
	icrhexp = getappdata(0,'ICRH_SHAPE_EXP');

	fpicrh_el  = max(1,interp1_ex(icrhexp.temps,icrhexp.pel,temps,'nearest','extrap'));
	fpicrh_el  = pchip(icrhexp.x,fpicrh_el,profil.xli);
   	indnok = find(any(~isfinite(fpicrh_el),2));
   	indok  = find(all(isfinite(fpicrh_el),2));
   	fpicrh_el(indnok,:) = ones(length(indnok),1) * mean(fpicrh_el(indok,:),1);

	fpicrh_ion  = max(1,interp1_ex(icrhexp.temps,icrhexp.pion,temps,'nearest','extrap'));
	fpicrh_ion  = pchip(icrhexp.x,fpicrh_ion,profil.xli);
   	indnok = find(any(~isfinite(fpicrh_ion),2));
   	indok  = find(all(isfinite(fpicrh_ion),2));
   	fpicrh_ion(indnok,:) = ones(length(indnok),1) * mean(fpicrh_ion(indok,:),1);

	fpfweh  = max(1,interp1_ex(icrhexp.temps,icrhexp.pfw,temps,'nearest','extrap'));
	fpfweh  = pchip(icrhexp.x,fpfweh,profil.xli);
   	indnok = find(any(~isfinite(fpfweh),2));
   	indok  = find(all(isfinite(fpfweh),2));
   	fpfweh(indnok,:) = ones(length(indnok),1) * mean(fpfweh(indok,:),1);

	fpicrh = fpicrh_ion + fpicrh_el + fpfweh;

%figure(51);clf;plot(temps,trapz(x,vpr .* fpicrh,2),'or',temps,picrh+pfweh,'+b');drawnow

	fact_icrh  = ((picrh+pfweh) ./  max(1,trapz(x,vpr .* fpicrh,2))) * ve;	
	fpicrh     = (fpicrh_ion + fpicrh_el) .* fact_icrh;
	fpicrh_ion = fpicrh_ion .* fact_icrh;	
	fpfweh     = fpfweh .* fact_icrh;


elseif tune_frac ~= 0

	% ce piquage est calibre sur pion
	%tune_frac = 0.5;
	piqicrh    = tune_frac .* fracmino;
	fpicrh     = exp(-(vt*x-xres_icrh*ve).^ 2 ./ ((piqicrh * ve) .^ 2) );
	fpicrh     = fpicrh .* ((picrh./  trapz(x,vpr .* fpicrh,2)) * ve);
	fpicrh_ion = ((pion_icrh ./ max(1e3,picrh)) * ve) .* fpicrh;
elseif isfield(profil,'picrh')
	% profiles are already computed in new ICRH model
	fpicrh     = profil.picrh;
	fpicrh_ion = profil.picrh_ion;
	fpicrh     = fpicrh .* ((picrh./  max(eps,trapz(x,vpr .* fpicrh,2))) * ve);
	fpicrh_ion = fpicrh_ion .* ((pion_icrh./  max(eps,trapz(x,vpr .* fpicrh_ion,2))) * ve);
else
	piqicrh    = 0.5 .* fracmino;
	fpicrh     = exp(-(vt*x-xres_icrh*ve).^ 2 ./ ((piqicrh * ve) .^ 2) );
	fpicrh     = fpicrh .* ((picrh./  trapz(x,vpr .* fpicrh,2)) * ve);
	fpicrh_ion = ((pion_icrh ./ max(1e3,picrh)) * ve) .* fpicrh;
end

% Poh  (attention il faut ajouter une contre reaction pour simuler le piquage)
%fpoh    = abs(((vloop ./ 2 ./ pi) * ve) ./ Raxe .* jres);
fpoh     = ej .* ((abs(poh) ./ max(1,abs(trapz(x,vpr .* ej,2)))) * ve);
%poh    = min(min(abs(vlmax) .* ip,2 .* RR .* ip .* ip),max(1,rm .* trapz(x,vpr_tor .* ej,2)));

% Pbrem
if isfield(profil,'pbrem')
	fpbrem   = profil.pbrem;
else
	fpbrem    = nep .^  2 .* sqrt(tep);
	fpbrem   = fpbrem .* ((pbrem./  trapz(x,vpr .* fpbrem,2)) * ve);
end

% Pcyclo
%fpcyclo    = nep .* tep;
%fpcyclo     = real(Raxe .* kx .^ 0.79 .* (fdia./Raxe) .^ 2.62 .* nep .^ 0.38 .* tep .* (16 + tep) .^ 2.61);
fpcyclo     = real(Raxe .* kx .^ 0.79 .* (fdia./Raxe) .^ 2.62 .* nep .^ 0.38 .* tep .^ 3.61);
fpcyclo     = fpcyclo .* ((pcyclo./  trapz(x,vpr .* fpcyclo,2)) * ve);


% Pline ou prad
if isfield(profil,'fprad')
	fprad      = profil.fprad;
else
	fprad      = vt*exp(-(x - 1) .^ 2 ./ 0.01);
end
% dans le plasma de coeur
fprad    = fprad .* ((prad./  max(1,trapz(x,vpr .* fprad,2))) * ve);


% ionization
if isfield(profil,'s0m')
	fpioniz   = max(1,profil.s0m) + max(1,profil.s0) + max(1,profil.spellet);  
	fpioniz   = fpioniz .* ((pioniz./  max(1,trapz(x,vpr .* fpioniz,2))) * ve);
	fpioniz_i = max(1,profil.s0m) + max(1,profil.s0);
	fpioniz_i = fpioniz_i .* ((pioniz_i./  max(1,trapz(x,vpr .* fpioniz_i,2))) * ve);
else
	fpioniz   = zeros(size(fprad));
	fpioniz_i = zeros(size(fprad));
end


% fptot
if collapse == 0
    fptot = max(1./(Vp*ve),fplh + fpecrh + fpfweh + fpnbi + fpfus + fpicrh + fpoh - fprad - fpcyclo - fpbrem - fpioniz);
else
    fptot = fplh + fpecrh + fpfweh + fpnbi + fpfus + fpicrh + fpoh - fprad - fpcyclo - fpbrem - fpioniz;
end

% calcul du piquage de te et du facteur d'amelioration du confinement

% dependance du transport en q^2
%save('loc6'); 
[hitb,ateloc,aitb,teshape,xieshape,xieshape_itb]       = z0piqpe(x,fs,fsdb,nep,qjli,fptot,tep,kishape,qdds,kidds,xieorkie); % il faut le meme profil de q dans les 2 appels

% ajout de l'effet des dents de scie resolu dans le temps
indice_inv_min = 1.5;
if (qdds < 0) && any(fix(indice_inv) > indice_inv_min) && (evolution == 0)
    indice_crash = find(indice_inv > indice_inv_min);
    indice_bad   = find(diff(indice_crash) < 3);
    if ~isempty(indice_bad)
        indice_crash(indice_bad+1) = [];
    end
    indice_crash(indice_crash == 1) = [];
    indice_crash(indice_crash == length(temps)) = [];
    kidds_evol = ones(size(qjli,1),1);
    taue_dds    = max(taue ./ 2,1e-6) ;
    for ln=1:length(indice_crash)
        kidds_evol(indice_crash(ln)) = kidds;
        taue_st = max(taue_dds(indice_crash(ln)) .* (indice_inv(indice_crash(ln)) ./ 21) .^ 2,1e-6);
        if ln <length(indice_crash)
            taue_st = min(taue_st,(temps(indice_crash(ln+1)) - temps(indice_crash(ln))) / 3);
        end
        %%%%skidds_evol = kidds_evol ./ taue_st;
        skidds_evol = ones(size(kidds_evol)) ./ taue_st;
        segment_time = indice_crash(ln):length(temps);
        [tntot,kidds_inter] = z0ode(temps(segment_time),skidds_evol(segment_time),taue_st * ones(length(segment_time),1),kidds);
        kidds_evol(segment_time) = max(1,min(kidds,kidds_inter));
    end
    for ln=1:length(kidds_evol)
        indice_modif = max(find(qjli(ln,:) < abs(qdds)));
        if indice_inv(ln) > indice_inv_min
            xieshape(ln,1:indice_inv(ln))      = max(xieshape(ln,1:indice_inv(ln))      * kidds);
            xieshape_itb(ln,1:indice_inv(ln))  = max(xieshape_itb(ln,1:indice_inv(ln))  * kidds);            
        elseif ~isempty(indice_modif)
            xieshape(ln,1:indice_modif)      = xieshape(ln,1:indice_modif)      * kidds_evol(ln);
            xieshape_itb(ln,1:indice_modif)  = xieshape_itb(ln,1:indice_modif)  * kidds_evol(ln);
        end
    end
    %figure(21);plot(temps,kidds_evol);drawnow
    %figure(24);plot(x,xieshape);drawnow
    
elseif (qdds < 0) && (evolution ~= 0)
    indice_inv = max(1,fix(indice_inv));
    if (indice_inv(3) > indice_inv_min) &&  ...
            (indice_inv(1) <= indice_inv_min)&&  ...
            (indice_inv(2) <= indice_inv_min)
        kidds_evol(2) = kidds;
        taue_dds    = max(taue .* (indice_inv(3) ./ 21) .^ 2,1e-6);
        crash_indice = 3;
    elseif (indice_inv(1) <= indice_inv_min)&&  ...
            (indice_inv(2) <= indice_inv_min)
        crash_indice = [];
        taue_dds    = max(min(min(taue),min(diff(temps))),1e-6) .* ones(size(temps));
        for ln=1:length(temps)
            indice_modif = max(find(qjli(ln,:) < abs(qdds)));
            if ~isempty(indice_modif)
                taue_dds(ln)    = max(kidds_evol(ln) ./ kidds .* taue(ln) .* (indice_modif ./ 21) .^ 2 ,1e-6);
            end
        end
    else
        % antialiasing 
        crash_indice = [];
        taue_dds    = max(min(min(taue),min(diff(temps))),1e-6) .* ones(size(temps));        
    end
    skidds_piq  = ones(size(kidds_evol)) ./ taue_dds;
    [tntot,kidds_inter] = z0ode(temps(2:end),skidds_piq(2:end),taue_dds(2:end),kidds_evol(2));
    kidds_evol(3:end)   = max(1,min(kidds,kidds_inter(2:end)));
%     if isappdata(0,'debug_st_evol')
%         debug_st_evol = getappdata(0,'debug_st_evol');
%         debug_st_evol.temps      = cat(1,debug_st_evol.temps,temps(:));
%         debug_st_evol.kidds_evol = cat(1,debug_st_evol.kidds_evol,kidds_evol(:));
%         debug_st_evol.skidds_piq = cat(1,debug_st_evol.skidds_piq,skidds_piq(:));
%         debug_st_evol.taue_dds   = cat(1,debug_st_evol.taue_dds,taue_dds(:));
%         debug_st_evol.indice_inv = cat(1,debug_st_evol.indice_inv,indice_inv(:));
%         crash_indice_debug       = zeros(size(temps(:)));
%         if ~isempty(crash_indice)
%             crash_indice_debug(crash_indice) = 1;
%         end
%         debug_st_evol.crash_indice_debug = cat(1,debug_st_evol.crash_indice_debug,crash_indice_debug(:));
%         debug_st_evol.tntot = cat(1,debug_st_evol.tntot,tntot(:));
%         debug_st_evol.kidds_inter = cat(1,debug_st_evol.kidds_inter,kidds_inter(:));
%         setappdata(0,'debug_st_evol',debug_st_evol);
%     else
%         debug_st_evol.temps      = temps(:);
%         debug_st_evol.kidds_evol = kidds_evol(:);
%         debug_st_evol.skidds_piq = skidds_piq(:);
%         debug_st_evol.taue_dds   = taue_dds(:);
%         debug_st_evol.indice_inv = indice_inv(:);
%         crash_indice_debug       = zeros(size(temps(:)));
%         if ~isempty(crash_indice)
%             crash_indice_debug(crash_indice) = 1;
%         end
%         debug_st_evol.crash_indice_debug = crash_indice_debug(:);
%         debug_st_evol.tntot = tntot(:);
%         debug_st_evol.kidds_inter = kidds_inter(:);       
%         setappdata(0,'debug_st_evol',debug_st_evol);
%     end
    
    for ln=1:length(kidds_evol)
        if ln == crash_indice 
            xieshape(ln,1:indice_inv(crash_indice))      = max(xieshape(ln,1:indice_inv(crash_indice)) .* kidds);
            xieshape_itb(ln,1:indice_inv(crash_indice))  = max(xieshape_itb(ln,1:indice_inv(crash_indice)) .* kidds);
        else
            indice_modif = max(find(qjli(ln,:) < abs(qdds)));
            if ~isempty(indice_modif)
                if kidds_evol(ln) == kidds
                    xieshape(ln,1:indice_modif)      = max(xieshape(ln,1:indice_modif)      .* kidds_evol(ln));
                    xieshape_itb(ln,1:indice_modif)  = max(xieshape_itb(ln,1:indice_modif)  .* kidds_evol(ln));
                else
                    xieshape(ln,1:indice_modif)      = xieshape(ln,1:indice_modif)      .* kidds_evol(ln);
                    xieshape_itb(ln,1:indice_modif)  = xieshape_itb(ln,1:indice_modif)  .* kidds_evol(ln);
                end
            end
        end
    end
    %figure(21);plot(temps,kidds_evol);drawnow
    %figure(24);plot(x,xieshape);drawnow

end

% test de barriere imposee
%ix = find(temps >= 10);
%hitb(ix) = 1.5;
%ipositb  = find(x >= 0.6  & x <= 0.7);
%xieshape_itb(ix,ipositb) = xieshape(ix,ipositb)/10;
%save('loc7');

% pour plus de precision
fwcorr = trapz(x,ee .*(nep .* tep + nip .* tip) .* vpr,2) ./ tep(:,1);

% position de la bariere
xdep1       = max((vt*x) .* (qjli == (min(qjli,[],2)*ve)),[],2); % profil creux
xdep2       = max((vt*x) .* ((jli <= 1) & ((vt *x ) < 0.7)),[],2); % trou de courant
xdep        = max(xdep2,xdep1);

% autres effets
% xq1,xq2,xq3
dqx   = abs(qjli - 1) + ((vt*x) == 0) ./ eps;
xq1  = max((vt*x) .* (dqx == (min(dqx,[],2)*ve)),[],2);
xq1(qmin>1.2) = 0;
dqx   = abs(qjli - 1.5) + ((vt*x) == 0) ./ eps;
xq15  = max((vt*x) .* (dqx == (min(dqx,[],2)*ve)),[],2);
xq15(qmin > 1.7) = 0;
dqx   = abs(qjli - 2) + ((vt*x) == 0) ./ eps;
xq2  = max((vt*x) .* (dqx == (min(dqx,[],2)*ve)),[],2);
xq2(qmin >2.2) = 0;
dqmin   = abs(qjli - min(qjli,[],2) * ve) + ((vt*x) == 0) ./ eps;
xqmin  = max((vt*x) .* (dqmin == (min(dqmin,[],2)*ve)),[],2);

% larmor radius
rhoi    = 6.46e-3 .* sqrt(tep ./1e3) ./ (Bt * ve);
qd2     = abs(pdederive(vt*x, qjli,2,2,2,2));
qming   = min(qjli,[],2);
maskg   = (qjli == (qming * ve));
qd2ming = sum(qd2 .* maskg,2) ./ sum(maskg,2);
rhoig   = sum(rhoi .* maskg,2) ./ sum(maskg,2);
% largeur du gap normalise
warning off
dgap1  = sqrt(2 .* qming .* rhoig ./ xq1 ./ qd2ming);
dgap15 = sqrt(2 .* qming .* rhoig ./ xq15 ./ qd2ming ./ 2);
dgap2 = sqrt(2 .* qming .* rhoig ./ xq2 ./ qd2ming ./ 2);
warning on



% aide au declenchement sur surface rationnelle
precq   = length(x) ./ max(1,qa);
%ratq =  exp(-abs(qmin - 2) .* precq)  .* dgap2 .* xq2 +  ...
%        exp(-abs(qmin - 1.5) .* precq) .* dgap15 .* xq15+  ...
%        exp(-abs(qmin - 1) .* precq) .* dgap1 .* xq1;
ratq =  exp(-abs(qmin - 2) .* precq)  .* xq2 .^2 +  ...
        exp(-abs(qmin - 1.5) .* precq) .*  xq15 .^2+  ...
        exp(-abs(qmin - 1) .* precq) .* xq1 .^2 .* (qmin >= 1);


%figure(61);clf;plot(temps,hitb,temps,ratq);drawnow


% ajout a hitb + securite (on peu mettre plus)
hitb = min(2,hitb + min(0.1,max(0,ratq)));

% mhd simulee : test
hmhd  = 1 - (q0 >  (1.1 .* qmin)) .* (exp(-abs(q0 - 2) .* precq) .* xq2 .^ 2 +  ...
            exp(-abs(q0 - 1.5) .*precq ) .* xq15 .^ 2) - ...
            max(0,min(1,(1 - q0))) .* xq1 .^ 2  - ...
	    max(0,tanh(q0 -  (3 + qmin))) .* xqmin .^ 2;

%  figure(23);clf;plot(hmhd);drawnow	     

	     
% dds (effet moyen = 0.5 )
% securite
hmhd  = max(0.1,hmhd);
%save('loc8');


% les profils pour z0dacces
profil.jli  = jli;
profil.jeff  = jeff;
profil.qjli = qjli;
profil.jboot = jboot;
profil.eta   = eta;
profil.jnbicd = jnbicd;
profil.jlh    = jlh;
profil.jeccd  = jeccd;
profil.jfwcd = jfwcd;
profil.jfus  = jfus;
profil.jrun  = jrun;
profil.plh   = fplh;
profil.pnbi  = fpnbi;
profil.pnbi_ion = fpnbi_ion;
profil.pecrh = fpecrh;
profil.pfweh = fpfweh;
profil.picrh = fpicrh;
profil.picrh_ion = fpicrh_ion;
if isfield(profil,'tep') && isfield(profil,'salf') && isfield(profil,'pfus')
  % nothing, it is computed before
else
   profil.pfus_ion = ((pion_fus ./ max(1e3,pfus)) * ve) .* fpfus;
   profil.pfus  = fpfus;
end
profil.pbrem = fpbrem;
profil.prad = fprad;
profil.pioniz = fpioniz;
profil.pioniz_i = fpioniz_i;
profil.pcyclo = fpcyclo;
profil.pohm = fpoh;
profil.nep  = nep;
profil.nip  = nip;
profil.vpr  = vpr;
profil.vpr_tor  = vpr_tor;
profil.spr  = spr;
profil.grho2r2  = grho2r2;
profil.r2i  = r2i;
profil.ri  = ri;
profil.C2  = C2;
profil.C3  = C3;
profil.grho  = grho;
profil.grho2  = grho2;
profil.kx  = kx;
profil.dx  = dx;
profil.Raxe  = Raxe;
profil.epsi = epsi;
profil.rmx = rmx;
profil.bpol = bpol;
profil.fdia = fdia;
profil.psi = psi;
profil.phi = phi_tor;
profil.dphidx = dphidx_tor;
profil.dpsidt = dpsidt;
profil.epar = epar;
profil.zeff = zeffp;
profil.n1p = n1p;
profil.nhep = nhep;
profil.nzp = nzp;
profil.xieshape     = xieshape;
profil.xieshape_itb = xieshape_itb;
profil.source_ion = profil.pnbi_ion + profil.picrh_ion + profil.pfus_ion - profil.pioniz_i;
profil.source_el  = fptot - profil.source_ion;
profil.jni        =jni;
profil.ftrap     = ftrap;
profil.ptot      = ptot_out;
profil.ej        = ej;
profil.xli = x;

% profile for coupling to FREEBIE
profil.jgs       = jgs;
profil.df2dpsi   = df2dpsi;
profil.dptotdpsi = dptotdpsi;

% te0 reel
te0 = tep(:,1);
% keyboard
%save('loc9');
%keyboard


% fonction de limitation soft du shear (saturation)
function o = slim(e)

o = pi .* tanh(e ./pi);

