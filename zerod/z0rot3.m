% dans metis la dirction toroidal est dans le sens du courant plasma, 
% de telle sorte que Btheta est toujours positif
% le moment injecte par l'IDN est posistif si l'injection est co courant
% le champs toroidal est positif s'il est dans le sens du courant
% calcul de la vitesse de rotation moyenne
function [wrad,snbi,sicrh,sfus,sripth,sriplh,sripicrh,sturb,fact,wrot,slh,profli] = z0rot3(zs,option,cons,geo,profli)

% constante physique (phys)
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

% case of Hydrogen NBI in DT plasma
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
   gas_nbi = -1; 
else
   gas_nbi = option.gaz; 
end

% isotopic composition for option.gaz == 5
if option.gaz == 5
    nHe3onD = real(cons.iso);
    nTonD   = imag(cons.iso);
    warning('nHe3onD & nTonD not yet used !');
else
    nHe3onD = zeros(size(cons.iso));
    nTonD   = real(cons.iso);
end
cons.iso = real(cons.iso);



% improve precision
[A_el,Z_el,name_el] = chargemasse;
dd   = abs(Z_el - option.zimp);
mask = (dd == min(dd));
aimp = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(aimp)
    aimp = 7/3 .* option.zimp;
end
dd   = abs(Z_el - option.zmax);
mask = (dd == min(dd));
amax = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(amax)
    amax = 7/3 .* option.zmax;
end



% compatibilite mode fast pour plot
if isfield(profli,'temps') && (length(cons.temps) ~= length(profli.temps))
        times = cons.temps;
        temps = profli.temps;
	% modification des donnees
	noms = fieldnames(profli);
	for l=1:length(noms)
		nomc = noms{l};
		val  = getfield(profli,nomc);
		if all(size(val) > 1)
			warning off
			valn  = interp1(temps,val,times,'linear');
			warning on
			indbad      = find(any(~isfinite(valn),2));
			if ~isempty(indbad)
				valn(indbad,:) = ones(length(indbad),1) * val(end,:);
			end
			profli = setfield(profli,nomc,valn);
		end
	end
	profli.temps = times;
	
end

% parametre pour jouer sur la rotation intrinseque
if ~isfield(option,'fintrinsic')
	option.fintrinsic = 1;
end
% pour une discussion general sur la rotation d'ensemble du plasma :
% ref : NF 2011 013006 (7p) V.D. Pustovitov : Integral torque balance in tokamaks

% sacling pour taurotmul
if option.taurotmul  == 0
	% il semble que le temps de confinement de la rotation ne depasse jamais le temps de confinement de l'energie
	%option.taurotmul = min(1,max(0.1,zs.nbar ./ zs.nsat));
	%figure(21);clf;plot(zs.temps,option.taurotmul);drawnow
	% reference : P.C. de Vries PPCF 48 (2006) p 1693
	% reference : P.C. de Vries NF 48 (2008) 
	tau_rot = min(zs.taue,zs.tauii);
else
	tau_rot = option.taurotmul .* zs.taue;
end


% profils utiles
x   = profli.xli;
ux  = (1  - x .^ 2);
ve  = ones(size(x));
vt  = ones(size(geo.R));


% facteur pour passer de V au moment total : V  * fact = M
% equation resolue : dV/dt = -V/tau + S;
% calcul plus precis (le profil de rotation est suppose homothetique a Ti) <m R> Vp
mass = ((((zs.n1m - zs.nDm - zs.nTm) + 2 .* zs.nDm  + 3 .* zs.nTm) ./ zs.n1m) * ve)  .* profli.n1p + ...
       4 .* profli.nhep  + ((option.zimp+1) * 2 + option.rimp .* (option.zmax+1) .* 2) .* profli.nzp;
if isfield(profli,'nwp')
    if option.Sn_fraction > 0
        mass = mass + (1 - option.Sn_fraction) .* 183.84 .* profli.nwp + ...
                           option.Sn_fraction  .* 118.71 .* profli.nwp;
    else
        mass = mass + 183.84 .* profli.nwp;
    end
end
if isfield(profli,'omega')
	poids = max(eps,abs(profli.omega));
	fact = phys.ua .* trapz(x,profli.vpr .* mass .* poids .* profli.Raxe,2) ./  ...
       		trapz(x,profli.vpr .* poids,2) .* trapz(x,profli.vpr,2) ;
else
	fact = phys.ua .* trapz(x,profli.vpr .* mass .* profli.tip .* profli.Raxe,2) ./  ...
       		trapz(x,profli.vpr .* profli.tip,2) .* trapz(x,profli.vpr,2) ;
end

% coefficient de friction sur le neutre de bord
warning off
sn0fr_core      = zs.sn0fr ./ zs.wrad;
sn0fr_core(zs.wrad == 0) = 0;
warning on
%
if isfield(profli,'omega')
	warning off
	sn0fr_edge      = zs.sn0fr ./ profli.omega(:,end);
	sn0fr_edge(profli.omega(:,end)== 0) = 0;	
	warning on
else
        sn0fr_edge = sn0fr_core;
end



% temps de confinement (prise en compte du multiplicateur et de la friction)
invtau   =  1 ./ max(1e-6,tau_rot);
invtaufr =  invtau - sn0fr_core ./ fact ./ geo.R;
if any(invtaufr <= 0)
	invtaufr(invtaufr <= 0) = invtau(invtaufr <= 0);
end
taufr    =  1./ invtaufr; 


% ces donnnees ne sont pas calculer dans cette version (il faut utiliser la version complete zrot2)
sicrh    = NaN .* vt;
sfus     = NaN .* vt;
sripth   = NaN .* vt;
sriplh   = NaN .* vt;
sripicrh = NaN .* vt;
s_epar   = NaN .* vt;


% source du aux pertes ripples
% For bulk plasma (JxB torque);
% ion loss =  dV_phi counter-current on T-S
% electron oss = dV_phi co-current on T-S
% torque ~ I_return x B 
if isfield(profli,'bpol')
    bpol     = profli.bpol(:,end);
else
    bpol     = (geo.a ./ geo.R) ./ zs.q95 .* geo.b0;
end
% for edge effect on E_r :
% Ref T.M. Wilks, Pop 24 (2017) 012505 and Pop 23 (2016) 122505
% Ref J. S. deGrassie et al, PoP 22 (2015) 080701
% sign opposite to bulk
lne            = 31.3 - log(sqrt(zs.nebord)./(zs.tebord));
eta_per        = 1.03e-4 .* zs.zeff .* lne .* zs.tebord .^ -(3/2);
edge_t2        =  fact  ./ max(1e-6,tau_rot) .* eta_per ./ bpol .^ 2  ./ ( zs.dsol .* 2 .* pi .* geo.R) ./ geo.R;
%figure(21);clf;plot(cons.temps,edge_t2);drawnow
%
switch option.rip
case 0
	sripth   = 0 .* vt;
	sriplh   = 0 .* vt;
	sripicrh = 0 .* vt;
	
case 1
	ilhrip       = max(0,cons.plh - zs.plh) ./ max(1,zs.einj_lh);
	sriplh       = ilhrip .* bpol .* geo.R;
    switch option.mino
        case 'He3'
            ag = 3;
            zg = 2;
        case 'T'
            ag = 3;
            zg = 1;
        case 'He4'
            ag = 4;
            zg = 2;
        case 'D'
            ag = 2;
            zg = 1;
        case 'B'
            ag = 11;
            zg = 5;
        otherwise
            ag = 1;
            zg = 1;
    end
        iicrhrip   = -double(option.fwcd == 0) .* zg .* max(0,cons.picrh - zs.picrh.* (1 - zs.frloss_icrh)) ./ max(1,zs.einj_icrh);
	sripicrh   = iicrhrip .* bpol .* geo.R;

	% effet du ripple :
	warning off
	z1   =  (1+(zs.meff >3));
	qeff = 1./(1./zs.qa + 1./zs.qmin);
	% Tokamaks, Wesson p 661
	ras    = geo.a ./ geo.R;
	lnii   = 17.3 - 0.5.*log(zs.nem ./ 1e20) + 3/2 .* log(zs.tem .* zs.tite ./ 1e3); % pour H, D et T uniquement
	% pour l'espece principale
	% Tokamaks, Wesson p 663
	taui   =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* zs.meff)  .* ...
				(phys.e .* zs.tem .* zs.tite) .^ (3/2) ./ zs.nim ./ lnii ./ z1 .^ 4 ;
	% Plasma rotation in Tokamaks, V. Rozhansky and M. Tendler,Review of plasma physics, tome 19, p 163
	nui    = 1./ min(taui,1);
	vthi   = sqrt(2 .* phys.e .*zs.tem .* zs.tite ./ zs.meff ./ phys.mp);
	nuis   = nui .* geo.R .* qeff .* ras .^ (-3/2) ./ vthi;
	fk     = (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* ras .^ 2) ./ ...
			(1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* ras .^ 3);
	fk     = max(-2.1,min(1.7,fk));
	frot   = (fk + (7/2) - 1) ;
	warning on
	% terme ohmique du au champ electrique // (C. Bourdelle these  5-15)
	% EXPRESSION THESE C. BOURDELLE 5-19 (il y a toujours un peu de ripple dans un tokamak)
	sripth    =  - frot   .* zs.tem .* zs.tite ./ geo.a ./ z1 ./ ( phys.mu0 .* max(1e3,zs.ip ./ zs.peri)) .* ...
                       fact ./ max(1e-6,tau_rot);
                      
end

switch option.rot_jr_loss
case 'on'
    switch option.mino
        case 'He3'
            ag = 3;
            zg = 2;
        case 'T'
            ag = 3;
            zg = 1;
        case 'He4'
            ag = 4;
            zg = 2;
        case 'D'
            ag = 2;
            zg = 1;
        case 'B'
            ag = 11;
            zg = 5;
        otherwise
            ag = 1;
            zg = 1;
    end
    iicrhloss   = - double(option.fwcd == 0) .* zg .*  zs.picrh .* zs.frloss_icrh ./ max(1,zs.einj_icrh);
	slossicrh   = iicrhloss .* bpol .* geo.R;
        inbiloss    = - real(zs.pnbi) .* real(zs.firstorb_nbi) ./ max(1,option.einj) - ...
                        imag(zs.pnbi) .* imag(zs.firstorb_nbi) ./ max(1,option.einj2);
	slossnbi    = inbiloss .* bpol .* geo.R;
        ifusloss   = - 2 .* zs.pfus_loss ./ max(1,3.56e6);
	slossfus   = iicrhloss .* bpol .* geo.R;
	% thermal losses appears at the very edge. METIS grid step is to large to take it into account

otherwise
	slossfus = 0 .* vt;
	slossnbi   = 0 .* vt;
	slossicrh  = 0 .* vt;
end
% summation of term due to particles loss
slossandrip = slossicrh + slossnbi + slossfus + sripth + sriplh + sripicrh;
% edge value (JxB - Ohm radial current effect)
% very limited effect up to edge become cold
slossandrip_edge = (1 - edge_t2) .* slossandrip;
%figure(21);plot(cons.temps,slossandrip,'r',cons.temps,slossandrip,'b');drawnow

% calcul des sources
% idn -> injection de moment
rnbi      = zs.mu0_nbi;
%snbi      = geo.R .* (2 .* (1-cons.ftnbi) + 3 .* cons.ftnbi) .* phys.ua .* (zs.pnbi_th ./ option.einj./phys.e)  .* ...
%            sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ (2 .* (1-cons.ftnbi) + 3 .* cons.ftnbi)) .* rnbi;
% le sinus qui sert a moduler la fraction // et perpendiculaire de l'injection de neutre est deja inclus dans rnbi
if isfield(profli,'pnbi')
	% on suppose le depot cote champ faible
    if gas_nbi ==  -1
        cmb1 = (ones(size(cons.ftnbi)) .* phys.ua) * ve;
        vmb1 = sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ ones(size(cons.ftnbi))) * ve;
        cmb2 = (ones(size(cons.ftnbi)) .* phys.ua) * ve;
        vmb2 = sqrt(2 .*option.einj2 .* phys.e ./ phys.ua ./ ones(size(cons.ftnbi))) * ve;
    elseif gas_nbi == 3
        cmb1 = (((2 .* (1-real(cons.ftnbi)) + 3 .* real(cons.ftnbi)) .* phys.ua) * ve);
        vmb1 = sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ (2 .* (1-real(cons.ftnbi)) + 3 .* real(cons.ftnbi))) * ve;
        cmb2 = (((2 .* (1-imag(cons.ftnbi)) + 3 .* imag(cons.ftnbi)) .* phys.ua) * ve);
        vmb2 = sqrt(2 .*option.einj2 .* phys.e ./ phys.ua ./ (2 .* (1-imag(cons.ftnbi)) + 3 .* imag(cons.ftnbi))) * ve;
    else
        % still OK for 5 and 11; no boron or He3 NBI injector planned
        cmb1 = (((2 .* (1-real(cons.ftnbi)) + 1 .* real(cons.ftnbi)) .* phys.ua) * ve);
        vmb1 = sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ (2 .* (1-real(cons.ftnbi)) + 1 .* real(cons.ftnbi))) * ve;
        cmb2 = (((2 .* (1-imag(cons.ftnbi)) + 1 .* imag(cons.ftnbi)) .* phys.ua) * ve);
        vmb2 = sqrt(2 .*option.einj2 .* phys.e ./ phys.ua ./ (2 .* (1-imag(cons.ftnbi)) + 1 .* imag(cons.ftnbi))) * ve;
    end
    if isfield(profli,'pitch')
		pitch = profli.pitch;
	else	
		pitch = vt * ve;
	end
	psnbi1 = profli.Raxe .* cmb1 .* (real(profli.pnbi) ./ option.einj./phys.e) .* vmb1 .* real(pitch);
	snbi1  = sin(option.angle_nbi ./ 180 .* pi) .* trapz(x,profli.vpr .* psnbi1,2);
	psnbi2 = profli.Raxe .* cmb2 .* (imag(profli.pnbi) ./ option.einj2./phys.e) .* vmb2 .* imag(pitch);
	snbi2  = sin(option.angle_nbi2 ./ 180 .* pi) .* trapz(x,profli.vpr .* psnbi2,2);
        psnbi  = psnbi1 + sqrt(-1) .* psnbi2;
        snbi   = snbi1  + sqrt(-1) .* snbi2;

else
    if gas_nbi ==  -1
        snbi = geo.R .* 1 .* phys.ua .* (real(zs.pnbi_th) ./ option.einj./phys.e)  .* ...
            sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ 1) .* real(rnbi);
        snbi = snbi + sqrt(-1) .* ...
            geo.R .* 1 .* phys.ua .* (imag(zs.pnbi_th) ./ option.einj2./phys.e)  .* ...
            sqrt(2 .*option.einj2 .* phys.e ./ phys.ua ./ 1) .* imag(rnbi);
    elseif gas_nbi == 3
        snbi = geo.R .* (2 .* (1-real(cons.ftnbi)) + 3 .* real(cons.ftnbi)) .* phys.ua .* (real(zs.pnbi_th) ./ option.einj./phys.e)  .* ...
            sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ (2 .* (1-real(cons.ftnbi)) + 3 .* real(cons.ftnbi))) .* real(rnbi);
        snbi = snbi + sqrt(-1) .* ...
            geo.R .* (2 .* (1-imag(cons.ftnbi)) + 3 .* imag(cons.ftnbi)) .* phys.ua .* (imag(zs.pnbi_th) ./ option.einj2./phys.e)  .* ...
            sqrt(2 .*option.einj2 .* phys.e ./ phys.ua ./ (2 .* (1-imag(cons.ftnbi)) + 3 .* imag(cons.ftnbi))) .* imag(rnbi);
        
    else
        % still OK for 5 and 11; no boron or He3 NBI injector planned
        snbi = geo.R .* (2 .* (1-real(cons.ftnbi)) + 1 .* real(cons.ftnbi)) .* phys.ua .* (real(zs.pnbi_th) ./ option.einj./phys.e)  .* ...
            sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ (2 .* (1-real(cons.ftnbi)) + 1 .* real(cons.ftnbi))) .* real(rnbi);
        snbi = snbi + sqrt(-1) .* ...
            geo.R .* (2 .* (1-imag(cons.ftnbi)) + 1 .* imag(cons.ftnbi)) .* phys.ua .* (imag(zs.pnbi_th) ./ option.einj2./phys.e)  .* ...
            sqrt(2 .*option.einj2 .* phys.e ./ phys.ua ./ (2 .* (1-imag(cons.ftnbi)) + 1 .* imag(cons.ftnbi))) .* imag(rnbi);
    end
	if isfield(profli,profli.nbishape_ion)
		psnbi = real(profli.nbishape_ion) + sqrt(-1) .* imag(profli.nbishape_ion);
	else
		psnbi = 0 .* (vt * ve);
	end
end
% unite de nbi : m kg  s^-1 m s^-1 = N.m


% source de moment du a LH
% reference : L.G. Eriksson PPCF 51 (2009) p 044008.
directivite = abs(option.etalh);
% pour les modes ou elle n'est pas donnees
if directivite > 1
	directivite = 0.75;
end
% le  lobe > 0 injecte un moment a contre courant (si lh est co courant) et le lobe < 0 dans le sens du courant
slh = - sign(option.etalh) .* geo.R .* option.npar0 ./ phys.c .* zs.plh_th .* (2 .* directivite - 1);
% on ajoute la contribution de ECCD (le N// de ECCD est pris a 0.5 par defaut ; proche de l'optimal )
slh = slh - sign(option.sens) .*  geo.R .* 0.5 ./ phys.c .* zs.pecrh;


% source du au champ electrique // induit. cette source est negligeable mais brise la symetrie co et contre courant (ref: Kim 1993)
vion_epar = zs.meff .* (phys.me ./ phys.mp) .*  (zs.iohm + zs.irun) ./ (zs.nem .* zs.sp) ./ (phys.e .* zs.zeff);
s_epar = fact  ./ max(1e-6,tau_rot) .* vion_epar;

% calcul de la vitesse au bord en divertor (calcul Vassili Parail)
%wrad_edge  = (snbi + zs.sn0fr) ./ max(1,zs.meff) ./ phys.ua ./ max(1,zs.n0a ./ max(0.1,1 - zs.frac_pellet) + (zs.pnbi ./ option.einj ./ phys.e)) ./ geo.R .^ 2;
% rapport masse sur charge du plasma sortant	
nDa  = profli.n1p(:,end) .* (zs.nDm ./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2));
nTa  = profli.n1p(:,end) .* (zs.nTm ./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2));
nHa  = max(0,profli.n1p(:,end)) - nDa -nTa;
if isfield(profli,'nwp')
    
    if option.Sn_fraction > 0
        factor_w_sn = (1 - option.Sn_fraction) .* z0wavez(profli.tep(:,end)) .* profli.nwp(:,end) + ...
                      option.Sn_fraction .* z0snavez(profli.tep(:,end)) .* profli.nwp(:,end);
        msz  = (nHa + 2 .* nDa + 3 .* nTa + 4 .* profli.nhep(:,end) + ceil(aimp) .* profli.nzp(:,end) +  ...
            option.rimp .*  ceil(amax) .* profli.nzp(:,end) + 183.84 .* profli.nwp(:,end)) ./ ...
            (nHa + nDa + nTa + 2 .* profli.nhep(:,end) + option.zimp .* profli.nzp(:,end) +  ...
            option.rimp .*  option.zmax .* profli.nzp(:,end) + factor_w_sn);
    else
        msz  = (nHa + 2 .* nDa + 3 .* nTa + 4 .* profli.nhep(:,end) + ceil(aimp) .* profli.nzp(:,end) +  ...
            option.rimp .*  ceil(amax) .* profli.nzp(:,end) + 183.84 .* profli.nwp(:,end)) ./ ...
            (nHa + nDa + nTa + 2 .* profli.nhep(:,end) + option.zimp .* profli.nzp(:,end) +  ...
            option.rimp .*  option.zmax .* profli.nzp(:,end) + z0wavez(profli.tep(:,end)) .* profli.nwp(:,end));
    end
else
	msz  = (nHa + 2 .* nDa + 3 .* nTa + 4 .* profli.nhep(:,end) + ceil(aimp) .* profli.nzp(:,end) +  ...
        	option.rimp .*  ceil(amax) .* profli.nzp(:,end)) ./ ...
       		(nHa + nDa + nTa + 2 .* profli.nhep(:,end) + option.zimp .* profli.nzp(:,end) +  ...
        	option.rimp .*  option.zmax .* profli.nzp(:,end));
end
% calcul du flux de particule traversant la LCMS, y compris celle qui sont echangees contre des froides (interchange)
% pas seulement le flux net de matiere 
% dans JET taup (sans recyclage) est de l'ordre de 0.1 tauE (Stangeby)
% a partir du transport d'energie
% interchange chaleur (ion) + flux de matiere net  + NBI
nout = (max(1,zs.pion + zs.pei) ./ max(13.6,profli.tip(:,end)) ./ 1.602176462e-19 ./ 5 .* 2) +  ...
       profli.nip(:,end) ./  max(1e13,profli.nep(:,end)) .* zs.n0a ./ max(0.1,1 - zs.frac_pellet) +  ...
       (real(zs.pnbi) ./ option.einj ./ phys.e) + (imag(zs.pnbi) ./ option.einj2 ./ phys.e);


%  figure(21)
%  plot(cons.temps,max(1,zs.pion + zs.pei) ./ max(13.6,profli.tip(:,end)) ./ 1.602176462e-19 ./ 5 .* 2,'r', ...
%       cons.temps,profli.nip(:,end) ./  max(1e13,profli.nep(:,end)) .* zs.n0a,'b', ...
%       cons.temps,(real(zs.pnbi) ./ option.einj ./ phys.e) + (imag(zs.pnbi) ./ option.einj2 ./ phys.e),'c', ...
%       cons.temps,nout,'k.');
%  drawnow

% on ajoute l'effet du freinage du aux neutres froid au bord       
rap  = max(1,msz) .* phys.ua .* max(1,nout) .* geo.R .^ 2 - sn0fr_edge;
warning off
%wrad_edge =  (real(snbi) + imag(snbi) + slh + sripicrh + sriplh + sripth + s_epar)  ./ rap;
wrad_edge =  (real(snbi) + imag(snbi) + slh + slossandrip_edge + s_epar)  ./ rap;
wrad_edge(rap == 0) = 0;
warning on
wrad_edge(nout <= 1) = 0;
wrad_edge = min(zs.wrad,wrad_edge) .* (zs.wrad >= 0) + max(zs.wrad,wrad_edge) .* (zs.wrad < 0);
%wrad_edge = min(abs(zs.wrad),max(- abs(zs.wrad),wrad_edge));

% calcul de fintrinsic dependant de la collisionalite
if  option.fintrinsic == 0
    % from the equivalent intrinsic source shape model from reference:
    % J. C. Hillesheim et al, arxiv:1407.2121v1 physics.plasma-ph 8 Jul 2014
    if isfield(profli,'qjli') && isfield(profli,'tip') && isfield(profli,'omega')
	z1   =  (1+(zs.meff >3)) * ve;
	% Tokamaks, Wesson p 661
	lnii   = 17.3 - 0.5.*log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3); % pour H, D et T uniquement
	%
	% Tokamaks, Wesson p 663
	taui   =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* (zs.meff * ve))  .* ...
		  (phys.e .* profli.tip) .^ (3/2) ./ profli.nip ./ lnii ./ z1 .^ 2 ./ profli.zeff ;
	% attention sqrt(2) par rapport a la definition habituelle
	vthi   = sqrt(2 .* phys.e .* profli.tip ./ (zs.meff * ve) ./ phys.mp);
	nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ (vthi .* taui);
	nuc    = 1.7;
        fintrinsic_prof =  0.3 .* (nuis ./ nuc -1) ./ (1 + nuis ./ nuc .* 0.3);
        % volume averaged 
        % Banana /  low nuis = co current + normalisation to retrieve  scaling
        fintrinsic =  - min(1,max(-0.3,trapz(x(2:end),profli.vpr(:,2:end)  .* fintrinsic_prof(:,2:end),2) ./ max(eps, trapz(x(2:end),profli.vpr(:,2:end),2)))) ./ 0.3;
        %figure(21);clf;plot(cons.temps,fintrinsic);drawnow;
        
    else
	z1   =  (1+(zs.meff >3));
	qeff = 1./(1./zs.qa + 1./zs.qmin);
	% Tokamaks, Wesson p 661
	ras    = geo.a ./ geo.R;
	lnii   = 17.3 - 0.5.*log(zs.nem ./ 1e20) + 3/2 .* log(zs.tem .* zs.tite ./ 1e3); % pour H, D et T uniquement
	% pour l'espece principale
	% Tokamaks, Wesson p 663
	taui   =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* zs.meff)  .* ...
				(phys.e .* zs.tem .* zs.tite) .^ (3/2) ./ zs.nim ./ lnii ./ z1 .^ 2 ./ zs.zeff;
	% Plasma rotation in Tokamaks, V. Rozhansky and M. Tendler,Review of plasma physics, tome 19, p 163
	vthi   = sqrt(2 .* phys.e .*zs.tem .* zs.tite ./ zs.meff ./ phys.mp);
	nuis   = geo.R .* qeff .* ras .^ (-3/2) ./ (vthi .* taui);
	nuc    = 1.7;
        fintrinsic = - (nuis ./ nuc -1) ./ (1 + nuis ./ nuc .* 0.3);
        %figure(21);clf;subplot(3,1,1);plot(cons.temps,fintrinsic);
        %subplot(3,1,2);plot(cons.temps,nuis);
        %subplot(3,1,3);plot(nuis,fintrinsic,'.');drawnow;
    end
else
    fintrinsic = option.fintrinsic;
end

% backward compatibility
if ~isfield(option,'rotation_scl')
    if option.fintrinsic < 0
        option.rotation_scl = 'Barnes & Parra';
    else
        option.rotation_scl = 'Rice';
    end
end
% equation sur le moment total M = fact* V ~  m R V
% vitesse induite
if option.fintrinsic >=0
 switch option.rotation_scl
 case 'Rice'
      % source de rotation spontannee (expose ttf 2006 Marseille, Rice , ...)
      % utilisation de wion uniquement pour TS et TCV
      % nouveau scaling Rice NF  47, (2007) p 1618-1624
      % utilisation de Pion uniquement pour TS et TCV
      % on ajoute la pression suprathermique pour prendre en compte les effets des ion rapides.
      press  = min(trapz(x,(profli.nip .* profli.tip + profli.nep .* profli.tep) .* phys.e  .* profli.vpr,2) ./ zs.vp, ...
		  2 .* trapz(x,profli.nip .* profli.tip .* phys.e  .* profli.vpr,2) ./ zs.vp) + ...
		  (2/3) .* max(0,zs.w - zs.wth - zs.esup_lh) ./ zs.vp;
      % nouveau scaling Rice NF  47, (2007) p 1618-1624
      vturb  = fintrinsic .* 0.46e11 .* geo.b0 .^ 1.1 .* press .* zs.ip .^ (-1.9) .* geo.R .^ 2.2;
      
 case 'Barnes & Parra'
      % Intrinsic rotation driven by non-Maxwellian equilibria in tokamak plasmas
      % M. Barnes, F. I. Parra, J. P. Lee, E. A. Belli, M. F. F. Nave, and A. E. White
      % PRL
      % Accepted Thursday Jun 27, 2013 
      %
      % F.I. Parra et al, arXiv:1108.6106v2  22 Mar 2012
      %
      delta_ti = trapz(x,(profli.nip .* profli.tip)  .* profli.vpr,2) ./ max(1,trapz(x,profli.nip .* profli.vpr,2));
      vturb    = abs(option.fintrinsic) .* 18 .* delta_ti  ./ (max(1,zs.ip) ./ 1e6);
      
 case 'deGrassie'
      % reference:  J.S. deGrassie et all, PoP 23 (2016) 082501
      delta_ti = trapz(x,(profli.nip .* profli.tip)  .* profli.vpr,2) ./ max(1,trapz(x,profli.nip .* profli.vpr,2));
      vbar_ti  = sqrt(phys.e .*  delta_ti ./ phys.ua ./ zs.meff);
      V_A      = profli.fdia(:,end) ./ profli.Raxe(:,end) ./ sqrt(phys.mu0 .*  zs.meff .* phys.ua .* zs.nim);
      % no information for other mixture; keep 1 for option.gaz 5 and 11
      omega_c  = (1 + double(option.gaz == 4)) .* phys.e .* profli.fdia(:,end) ./ profli.Raxe(:,end) ./ phys.ua ./ zs.meff;
      rho_star = sqrt(2) .* vbar_ti ./ (profli.Raxe(:,end) .* profli.epsi(:,end)) ./ omega_c;
      M_A      = 1.25 .* (100 .* zs.betan) .* rho_star;
      vturb    = M_A .* V_A;
           
 otherwise
      error(sprintf('Z0ROT3: %s scaling is not yet implemanted',option.rotation_scl));
 end
else
  % Intrinsic rotation driven by non-Maxwellian equilibria in tokamak plasmas
  % M. Barnes, F. I. Parra, J. P. Lee, E. A. Belli, M. F. F. Nave, and A. E. White
  % PRL
  % Accepted Thursday Jun 27, 2013 
  %
  % F.I. Parra et al, arXiv:1108.6106v2  22 Mar 2012
  %
  delta_ti = trapz(x,(profli.nip .* profli.tip)  .* profli.vpr,2) ./ max(1,trapz(x,profli.nip .* profli.vpr,2));
  vturb    = abs(option.fintrinsic) .* 18 .* delta_ti  ./ (max(1,zs.ip) ./ 1e6);
end  
% source associe
sturb  = fact  ./ max(1e-6,tau_rot) .* vturb;
% c'est un pinch, le signe depend  de la rotation au bord => les mesures TS ne verifient pas cet effet
% et le ripple fort ne block pas cet effet (mesure TS)
%sturb  = sign(option.rip + sign(wrad_edge)) .* fact  ./ max(1e-6,tau_rot) .* vturb;
%sturb  = sign(wrad_edge) .* fact  ./ max(1e-6,tau_rot) .* vturb;

% equa diff (la friction est prise en compte dans la constante de temps)
%srot = (sturb + real(snbi) + imag(snbi) + slh  + sripicrh + sriplh + sripth + s_epar);
srot = (sturb + real(snbi) + imag(snbi) + slh  + slossandrip + s_epar);
srot(~isfinite(srot)) = 0;


if option.evolution == 1
	sini = zs.wrad(1) .* fact .* geo.R;
else
	sini           = srot(1) .* taufr(1);
	%sini           = srot(1) .* max(1e-6,tau_rot);
	if ~isfinite(sini)
   		sini = 0;
	end
end
[tvtorm,storm] = z0ode(cons.temps,srot,taufr,sini);
wrad           = storm ./ fact ./ geo.R;

% limite type viriel pour eviter les divergence numerique
wrot      = 0.5 .* fact .*  geo.R .^ 2 .* wrad .^ 2;
wmax      = sqrt(zs.wth ./ fact ./ geo.R .^ 2);
indbad    = find(abs(wrad) > wmax);
if ~isempty(indbad)
	wrad(indbad) = sign(wrad(indbad)) .* wmax(indbad);
end

% calcul du profil de rotation toroidal
% palsma de fond
switch option.gaz
    case 1
        zj = 1;
        aj = 1;
    case 2
        zj = 1;
        aj = 2;
    case 3
        zj = 1;
        aj = mean((2  + 3 .* cons.iso)  ./  (1+ cons.iso));
    case 4
        zj = 2;
        aj = 4;
    case 5
        zj = mean((1  + 4 .* cons.iso)   ./  (1 + 2 .* cons.iso));
        aj = mean((2  + 3 .* cons.iso)   ./  (1 + cons.iso));
    case 11
        zj = mean((1  + 25 .* cons.iso)   ./ (1 + 5.* cons.iso));
        aj = mean((1  + 11 .* cons.iso)  ./  (1 + cons.iso));
end

% choix de l'imprurete pour la rotation
switch option.impur_rot
case 'max'

  % impurete principale
  zimp = option.zmax;
  %aimp = ceil(option.zmax .* (7/3));

  % 2ieme impurete
  zmax = option.zimp;
  %amax = ceil(option.zimp .* (7/3));

otherwise

  % impurete principale
  zimp = option.zimp;
  %aimp = ceil(option.zimp .* (7/3));

  % 2ieme impurete
  zmax = option.zmax;
  %amax = ceil(option.zmax .* (7/3));

end

% improve precision
[A_el,Z_el,name_el] = chargemasse;
dd   = abs(Z_el - zimp);
mask = (dd == min(dd));
aimp = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(aimp)
    aimp = 7/3 .* zimp;
end
dd   = abs(Z_el - zmax);
mask = (dd == min(dd));
amax = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(amax)
    amax = 7/3 .* zmax;
end

% pour chaque espece d'ions
nDp   = max(1e13,profli.n1p .* ((zs.nDm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
switch option.gaz
    case 11
        nbp   = max(1e13,profli.n1p .* ((zs.nTm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
        nTp   = zeros(size(nbp));
    otherwise
        nTp   = max(1e13,profli.n1p .* ((zs.nTm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
        nbp   = zeros(size(nTp));
end
nHp   = max(1e13,profli.n1p - nTp - nDp);
%nhep  = max(1e13,profli.nhep);
switch option.gaz
    case 5
        nhe3p  = profli.nhep;
        nhep   = option.frhe0 .* profli.nep;
    otherwise
        nhep  = profli.nhep;
        nhe3p = zeros(size(profli.nhep));
end
if isfield(profli,'nwp')
	nwp  = profli.nwp;
else
	nwp = 0 .* nhep;
end
% choix de l'imprurete pour la rotation
switch option.impur_rot
case 'max'

  nz2p  = max(1e13,profli.nzp);
  nz1p  = max(1e11,profli.nzp .* option.rimp);   

otherwise

  nz1p  = max(1e13,profli.nzp);
  nz2p  = max(1e11,profli.nzp .* option.rimp);   

end

% masse
if option.Sn_fraction > 0
    Mtor    = phys.mp .*  max(1e13,nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + aimp .* nz1p + amax .* nz2p +  ...
              (1 - option.Sn_fraction) .* 183.84 .* nwp + option.Sn_fraction .* 118.71 .* nwp + 3.02 .* nhe3p + 11 .* nbp);
else
    Mtor    = phys.mp .*  max(1e13,nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + aimp .* nz1p + amax .* nz2p + 183.84 .* nwp + 3.02 .* nhe3p + 11 .* nbp);
end
% omega homothetic a Ti (cas Jet avec NBI)
if (option.omega_shape == 1) && isfield(profli,'omega')
    % utilisation scaling 
    % ref : NF 52 (2012) 042001, H. Weisen et al
    % nous supposons pour des problemes de stabilite numerique que xii est 
    xphi_loc = max(0.1,profli.xii ./ ((max(1e-6,tau_rot) ./ max(1e-6,zs.taue)) * ve));
    rot_epar = abs(profli.epar) .* ((s_epar ./ max(1,trapz(x,profli.vpr .* abs(profli.epar),2))) * ve);
    tstar_ltor_1 = (profli.rot_nbi + profli.rot_n0 + profli.rot_lh + rot_epar) ./ Mtor ./ xphi_loc;
    % le signe de ge n'est pas defini dans le papier !
    % utilisation de ge stationnaire pour de raison de bruit numerique
    % le signe semble ok
    einj_w = ((option.einj .* max(1,real(zs.pnbi)) + option.einj2 .* imag(zs.pnbi)) ./ max(1, real(zs.pnbi)+ imag(zs.pnbi))) * ve;
    sntot = profli.s0 + profli. s0m + profli.spellet + profli.pnbi ./ (einj_w .* phys.e);
    ge_ =   cumtrapz(profli.xli, profli.vpr .* sntot,2) ./ max(eps,profli.vpr_tor .* profli.grho2);
    ge_(:,1) = 0;
    tstar_2      = - ge_ ./ profli.nep ./ xphi_loc ./ profli.ri;

    lmin_   = (profli.rmx(:,end) ./ 101) * ve;
    lmax_   = (2 .* pi .* profli.qjli(:,end) .* profli.Raxe(:,end)) * ve;
    nepd1_  = pdederive(profli.xli,profli.nep,0,2,2,1);
    rhomax_ = profli.rmx(:,end) * ve;
    rolne_  = min(1./lmin_,max(1./lmax_,abs(nepd1_) ./ max(profli.nep,1e13) ./ rhomax_)) .* profli.Raxe;
    
    rolomega_2 = 1.2 .* tstar_2 +  0.41 .* rolne_ + 12 .* profli.ftrap ./ 1.46  + 0.41 .* profli.qjli  -  ...
                 1.9 .* (1 + profli.tip ./ max(13.6,profli.tep));
    terme_1    = tstar_ltor_1 .* profli.ri .*rhomax_;
    terme_2    = rolomega_2 .* profli.ri .*rhomax_ ;
    % resolution de proche en proche
    omega      = wrad_edge * ve;
    dx_ = mean(diff(profli.xli));
    for kz=(length(profli.xli)-1):-1:1
	omega(:,kz) = omega(:,kz + 1) + 1.2 .* (terme_1(:,kz) + terme_1(:,kz + 1)) ./ 2 .* dx_ + ...
                      (terme_2(:,kz) .* profli.omega(:,kz) + terme_2(:,kz + 1) .* profli.omega(:,kz + 1)) ./ 2 .* dx_;
    end

%      figure(21);
%      clf
%      plot(profli.xli,omega);
%      drawnow

    % pour la normalisation
    omega  = omega - omega(:,end) * ve;
elseif (option.omega_shape == 2) && isfield(profli,'nip')
    omega  = profli.tip .* profli.nip;
    omega  = abs(omega - omega(:,end) * ve);   
elseif (option.omega_shape == 3) && isfield(profli,'ptot')
    omega  = abs(profli.ptot - profli.ptot(:,end) * ve);
elseif (option.omega_shape == 4) && isfield(profli,'nip')
    omega  = abs(profli.nip - profli.nip(:,end) * ve);
elseif (option.omega_shape == 5) && isfield(profli,'nip') 
    % used of deGrassie scaling (as it is a local scaling)
    % reference:  J.S. deGrassie et all, PoP 23 (2016) 082501
    vbar_ti  = sqrt(profli.tip);
    V_A      = 1 ./ sqrt(profli.nip);
    omega_c  = sqrt((profli.fdia ./ profli.Raxe) .^ 2 + profli.bpol .^ 2);
    rho_star = vbar_ti ./ omega_c;
    M_A      = rho_star;
    omega    = M_A .* V_A ./ (profli.Raxe .* (1 + profli.epsi));
    omega(:,1) = omega(:,2);
    omega    = abs(omega - omega(:,end) * ve);
elseif (option.omega_shape == 6) && isfield(profli,'nip') 
    fshape = max(0,min(1,2 .* (abs(sturb) ./ abs(max(eps,srot)) - 0.5))) * ve;
    omega_ti  = abs(profli.tip - profli.tip(:,end) * ve);
    omega_ti  = omega_ti ./ max(eps,max(omega_ti,[],2) * ve);
    % used of deGrassie scaling (as it is a local scaling)
    % reference:  J.S. deGrassie et all, PoP 23 (2016) 082501
    vbar_ti  = sqrt(profli.tip);
    V_A      = 1 ./ sqrt(profli.nip);
    omega_c  = sqrt((profli.fdia ./ profli.Raxe) .^ 2 + profli.bpol .^ 2);
    rho_star = vbar_ti ./ omega_c;
    M_A      = rho_star;
    omega_dg = M_A .* V_A ./ (profli.Raxe .* (1 + profli.epsi));
    omega_dg(:,1) = omega_dg(:,2);
    omega_dg    = abs(omega_dg - omega_dg(:,end) * ve);
    omega_dg  = omega_dg ./ max(eps,max(omega_dg,[],2) * ve);
    omega     = fshape .* omega_dg + (1 - fshape) .* omega_ti;
    %figure(21);subplot(2,2,1);plot(profli.xli,omega,'r',profli.xli,omega_dg,'b',profli.xli,omega_ti,'g');
    %subplot(2,2,2);plot(cons.temps,fshape(:,1));
    %subplot(2,2,3);plot(cons.temps,sturb,'r',cons.temps,srot,'b',cons.temps,snbi,'g');drawnow
    
else
    omega  = abs(profli.tip - profli.tip(:,end) * ve);
end

% extranal shape for rotation
if isappdata(0,'VTOR_SHAPE_EXP') & isfield(profli,'qjli');
      % il faut aussi modifier dans zicd0
      vtor_shape = getappdata(0,'VTOR_SHAPE_EXP');
      omega_  = max(eps,interp1_ex(vtor_shape.temps,vtor_shape.omega,cons.temps,'nearest','extrap'));
      omega_  = pchip(vtor_shape.x,omega_,profli.xli);
      indok  = find(all(isfinite(omega_),2));
      omega(indok,:) = abs(omega_(indok,:) - omega_(indok,end) * ve);
      %figure(31);plot(profli.xli,omega)
end
% effet des dent de scie
if abs(option.qdds) > 0
	mask1  = (profli.qjli <= max(1,abs(option.qdds)));
	omega1 = max((~mask1) .* omega,[],2) * ones(size(x));
	omega  = (~mask1) .* omega + mask1 .* omega1;

end

% renormalisation de la quantite de mouvement
normv = trapz(x,profli.vpr .* Mtor .* omega .* profli.Raxe,2);
normf = trapz(x,profli.vpr .* Mtor .* profli.Raxe,2);
%profli.omega  = omega .* ((((wrad - wrad_edge) .* normf) ./ max(eps,normv)) * ve) + wrad_edge * ve;
omega  = (sign(normv) * ve).* omega .* ((((wrad - wrad_edge) .* normf) ./ max(eps,abs(normv))) * ve) + wrad_edge * ve;
if ~isfield(profli,'omega')
	profli.omega  = omega;
end

% les champ en Rmax
rmax         = profli.Raxe + geo.a * x;
btor         = option.signe .* (profli.fdia ./rmax);
grho         = abs((profli.rmx(:,end) * ve) ./ max(eps,pdederive(x,rmax,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,profli.psi,0,2,2,1)./ rmax .* grho ./ (profli.rmx(:,end) * ve);
btot         = sqrt(btor .^ 2 + bpol .^ 2);

% calcul de la rotation poloidal pour l'impurete principale
% Y. B. Kim et all Phys. Fluids. B 3  (8) 1991  p 2050-
switch option.gaz
    case 4
        alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.nhep) .* 4);
        nii   = max(1e13,profli.nhep);
    case 5
        alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.n1p) + max(1e13,profli.nhep) .* 4);
        nii   = max(1e13,profli.n1p + profli.nhep);        
    case 11
        alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.n1p) .* (1 + 25 .* cons.iso * ones(1,size(profli.n1p,2))));
        nii   = max(1e13,profli.n1p .* (1 + cons.iso * ones(1,size(profli.n1p,2))));
   otherwise
        alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.n1p) .* 1);
        nii   = max(1e13,profli.n1p);
end
% on traite les impuretes a l'etat de trace ...
alpha  = min(zimp,alpha);
beta   = (27/4) .^ 2 .* (aj ./ aimp) .^ 2 ./ (15/2 + sqrt(2 .* alpha) .* sqrt(aimp ./ aj));
g      = profli.ftrap ./ max(0.01,1 - profli.ftrap);
mui_00 = g .* (          alpha + sqrt(2)             -           log(1 + sqrt(2)));
mui_10 = g .* ((3/2)  .* alpha + 4  ./ sqrt(2)       - (5/2)  .* log(1 + sqrt(2)));
mui_01 = mui_10;
mui_11 = g .* ((13/4) .* alpha + 39 ./ (4 * sqrt(2)) - (25/4) .* log(1 + sqrt(2)));
D      = mui_00 .* (mui_11 + sqrt(2) + alpha - alpha .* beta) - mui_01 .^ 2;
D(D==0) = 1e38;
K1     = mui_01 ./ D .* (sqrt(2) + alpha - alpha .* beta);
K2     = (mui_00 .* mui_11 - mui_01 .* mui_10) ./ D;
vth    = sqrt(2 .* profli.tip .* phys.e ./ (phys.mp .* aj));
b2     = sqrt(profli.bpol .^ 2 + (profli.fdia .* profli.ri) .^ 2);
rhoi   = 4.57e-3 .* sqrt(aj .* profli.tip ./ 1e3 ./ b2);
ltim1  = pdederive(x,profli.tip,0,2,2,1) ./ profli.tip ./ (profli.rmx(:,end) * ve); 
lpiim1 = pdederive(x,profli.tip .* nii,0,2,2,1)  ./ (profli.tip .* nii)  ./ (profli.rmx(:,end) * ve); 
lpiIm1 = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.tip .* nz1p) ./ (profli.rmx(:,end) * ve); 
%figure(21);plot(x,1./ltim1,'b',x,1./lpiim1,'r',x,1./lpiIm1,'g');drawnow
%
% formule 34
vtheta_z1 = option.signe .*  0.5 .* vth .* rhoi .* ((K1 + (3/2) .* K2) .* ltim1 - lpiim1 + (zj ./ zimp) .* 1 .* lpiIm1) .* (profli.fdia .* profli.ri) ./ sqrt(b2);
utheta_z1 = vtheta_z1 ./ max(eps,profli.bpol);
utheta_z1(:,1) = 0;
% formule 33
vtheta_main = option.signe .*  0.5 .* vth .* rhoi .* (K1 .* ltim1) .* (profli.fdia .* profli.ri) ./ sqrt(b2);
utheta_main = vtheta_z1 ./ max(eps,profli.bpol);
utheta_main(:,1) = 0;
% calcul du moment poloidal injecte par IdN
% pour utiliser cette formulation il faut connaitre le courant radial
% ref F.L. Hinton et al , Physics Letters A 259 (1999) 267-275
% ce n'est pas possible ici
% ce terme est generalement negligeable ...
% sauf lors de la msie en marchedu chauffage
profli.utheta  = utheta_z1 + 0;
profli.vtheta  = profli.utheta .* bpol;

% changement de repere
dpsidrho =  pdederive(x,profli.psi,0,2,2,1) ./ (profli.rmx(:,end) * ve);
dpsidrho(:,1) = 0;

% computation of U_theta (term full E_r equation)
if option.solid_rotation <= 0
    if option.solid_rotation <  0
        % formulaire ORNL
        lnii   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(1);
        lnhe   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(2);
        lnhe3  = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(2);
        lnz1   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(zimp);
        lnz2   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(zmax);
        lnw    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(z0wavez(profli.tep));
        lnb    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(5);
        
        % pour l'espece principale
        % Tokamaks, Wesson p 663
        taui_h     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 1)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nHp ./ lnii ./ 1 .^ 4 ;
        taui_d     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 2)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nDp ./ lnii ./ 1 .^ 4 ;
        taui_t     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 3)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nTp ./ lnii ./ 1 .^ 4 ;
        taui_he     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 4)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nhep ./ lnhe ./ 2 .^ 4 ;
        taui_he3     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 3.02)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nhe3p ./ lnhe3 ./ 2 .^ 4 ;
        taui_z1     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* aimp)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nz1p ./ lnz1 ./ zimp .^ 4 ;
        taui_z2     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* amax)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nz2p ./ lnz2 ./ zmax .^ 4 ;
        taui_w      =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 183.84)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nwp ./ lnw ./ max(1,z0wavez(profli.tep)) .^ 4 ;
        taui_b      =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 11)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nbp ./ lnb ./ 5 .^ 4 ;
        if option.Sn_fraction > 0
            lnsn    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(z0snavez(profli.tep));
            taui_sn     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 118.71)  .* ...
                (phys.e .* profli.tip) .^ (3/2) ./ (option.Sn_fraction .* nwp) ./ lnsn ./ max(1,z0snavez(profli.tep)) .^ 4 ;
            if option.Sn_fraction < 1
                taui_w      =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 183.84)  .* ...
                    (phys.e .* profli.tip) .^ (3/2) ./ ((1 - option.Sn_fraction) .* nwp) ./ lnw ./ max(1,z0wavez(profli.tep)) .^ 4 ;
            end
        end
        
        % Plasma rotation in Tokamaks, V. Rozhansky and M. Tendler,Review of plasma physics, tome 19, p 163
        vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 1 ./ phys.mp);
        nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_h;
        nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
        fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
        fkh     = max(-2.1,min(1.7,fk));
        
        vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 2 ./ phys.mp);
        nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi./ taui_d;
        nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
        fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
        fkd     = max(-2.1,min(1.7,fk));
        
        vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 3 ./ phys.mp);
        nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_t;
        nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
        fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
        fkt     = max(-2.1,min(1.7,fk));
        
        vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 4 ./ phys.mp);
        nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_he;
        nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
        fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
        fkhe     = max(-2.1,min(1.7,fk));
        
        vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 3.02 ./ phys.mp);
        nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_he3;
        nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
        fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
        fkhe3 = max(-2.1,min(1.7,fk));
        fkhe3(~isfinite(fkhe3)) = 0;
        
        vthi   = sqrt(2 .* phys.e .* profli.tip  ./ aimp ./ phys.mp);
        nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_z1;
        nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
        fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
        fkz1     = max(-2.1,min(1.7,fk));
        
        vthi   = sqrt(2 .* phys.e .* profli.tip  ./ amax ./ phys.mp);
        nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_z2;
        nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
        fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
        fkz2     = max(-2.1,min(1.7,fk));
        
        vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 11 ./ phys.mp);
        nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_b;
        nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
        fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
        fkb    = max(-2.1,min(1.7,fk));
        fkb(~isfinite(fkb)) = 0;
        
        if option.Sn_fraction  == 1
            vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 183.84 ./ phys.mp);
            nuis   = zeros(size(profli.Raxe));
        else
            vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 183.84 ./ phys.mp);
            nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_w;
            nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
        end
        fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
        fkw     = max(-2.1,min(1.7,fk));
        
        if option.Sn_fraction > 0
            vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 118.71 ./ phys.mp);
            nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_sn;
            nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
            fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
                (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
            fksn     = max(-2.1,min(1.7,fk));
        end
        
        warning on
        % vitesse poloidal
        b2 = (geo.b0 .^ 2) * ve;
        % standard  + orbit squeezing
        gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
        gtheta2 = pdederive(x,phys.e .* profli.tip .*nHp,0,2,2,1) ./ nHp;
        stheta  = min(sign(gtheta1),sign(gtheta2));
        gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
        utheta_h  = - fkh .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho;
        utheta_h(:,1) = 2 .* utheta_h(:,2) - utheta_h(:,3);
        
        gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
        gtheta2 = pdederive(x,phys.e .* profli.tip .* nDp,0,2,2,1) ./ nDp;
        stheta  = min(sign(gtheta1),sign(gtheta2));
        gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
        utheta_d  = - fkd .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho;
        utheta_d(:,1) = 2 .* utheta_d(:,2) - utheta_d(:,3);
        
        gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
        gtheta2 = pdederive(x,phys.e .* profli.tip .* nTp,0,2,2,1) ./ nTp;
        stheta  = min(sign(gtheta1),sign(gtheta2));
        gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
        utheta_t  = - fkt .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho;
        utheta_t(:,1) = 2 .* utheta_t(:,2) - utheta_t(:,3);
        
        gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
        gtheta2 = pdederive(x,phys.e .* profli.tip .* nhep,0,2,2,1) ./ nhep;
        stheta  = min(sign(gtheta1),sign(gtheta2));
        gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
        utheta_he  = - fkhe .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ 2 .* gtheta ./ dpsidrho;
        utheta_he(:,1) = 2 .* utheta_he(:,2) - utheta_he(:,3);
        
        gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
        gtheta2 = pdederive(x,phys.e .* profli.tip .* nhe3p,0,2,2,1) ./ max(1,nhe3p);
        stheta  = min(sign(gtheta1),sign(gtheta2));
        gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
        utheta_he3  = - fkhe3 .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ 2 .* gtheta ./ dpsidrho;
        utheta_he3(:,1) = 2 .* utheta_he3(:,2) - utheta_he3(:,3);
        
        gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
        gtheta2 = pdederive(x,phys.e .* profli.tip .* nz1p,0,2,2,1) ./ nz1p;
        stheta  = min(sign(gtheta1),sign(gtheta2));
        gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
        utheta_z1  = - fkz1 .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ zimp .* gtheta ./ dpsidrho;
        utheta_z1(:,1) = 2 .* utheta_z1(:,2) - utheta_z1(:,3);
        
        gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
        gtheta2 = pdederive(x,phys.e .* profli.tip .* nz2p,0,2,2,1) ./ nz2p;
        stheta  = min(sign(gtheta1),sign(gtheta2));
        gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
        utheta_z2  = - fkz2 .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ zmax .* gtheta ./ dpsidrho;
        utheta_z2(:,1) = 2 .* utheta_z2(:,2) - utheta_z2(:,3);
        
        gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
        gtheta2 = pdederive(x,phys.e .* profli.tip .* nbp,0,2,2,1) ./ max(1,nbp);
        stheta  = min(sign(gtheta1),sign(gtheta2));
        gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
        utheta_b  = - fkb .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ 2 .* gtheta ./ dpsidrho;
        utheta_b(:,1) = 2 .* utheta_b(:,2) - utheta_b(:,3);
        
        %  	  figure(23);clf
        %  	  plot(x,profli.utheta,'r',x,utheta_z1,'b',x,utheta_z2,'c');
        %  	  drawnow
        
        %  	  figure(24);clf
        %  	  plot(x,utheta_main,'r',x,utheta_h,'b',x,utheta_d,'c',x,utheta_t,'m',x,utheta_he,'k');
        %  	  drawnow
        if option.Sn_fraction > 0
            gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
            gtheta2 = pdederive(x,phys.e .* profli.tip .* nwp,0,2,2,1) ./ nwp;
            stheta  = min(sign(gtheta1),sign(gtheta2));
            gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
            utheta_w  = - fkw .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ max(1,z0wavez(profli.tep)) .* gtheta ./ dpsidrho;
            utheta_w(:,1) = 2 .* utheta_w(:,2) - utheta_w(:,3);
            utheta_sn  = - fksn .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ max(1,z0snavez(profli.tep)) .* gtheta ./ dpsidrho;
            utheta_sn(:,1) = 2 .* utheta_sn(:,2) - utheta_sn(:,3);
        else
            gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
            gtheta2 = pdederive(x,phys.e .* profli.tip .* nwp,0,2,2,1) ./ nwp;
            stheta  = min(sign(gtheta1),sign(gtheta2));
            gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
            utheta_w  = - fkw .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ max(1,z0wavez(profli.tep)) .* gtheta ./ dpsidrho;
            utheta_w(:,1) = 2 .* utheta_w(:,2) - utheta_w(:,3);
        end
        
        % calcul du moment poloidal injecte par IdN
        if option.solid_rotation <  -1
            % dans le sens de la vistesse diamagnetique electronique (couplage par les collisions sur les electrons)
            fpnbi1      = (real(zs.pnbi) ./ max(eps,real(zs.pnbi) + imag(zs.pnbi))) * ve;
            fjnbi1      = (real(zs.inbicd) ./ max(eps,real(zs.inbicd) + imag(zs.inbicd))) * ve;
            pfast1      = fpnbi1 .* profli.pnbi .* ((real(zs.esup_nbi) ./ max(1,trapz(x,profli.vpr .* profli.pnbi,2))) * ve);
            nfast1      = max(1,pfast1 .* 2 ./ (phys.e .* option.einj));
            utheta_nbi1 = fjnbi1 .* profli.jnbicd ./ nfast1 ./ btot ./ phys.e;
            
            % calcul du moment poloidal injecte par IdN
            % dans le sens de la vistesse diamagnetique electronique (couplage par les collisions sur les electrons)
            % must be account for 2 NBI
            fpnbi2      = (imag(zs.pnbi) ./ max(eps,real(zs.pnbi) + imag(zs.pnbi))) * ve;
            fjnbi2      = (imag(zs.inbicd) ./ max(eps,real(zs.inbicd) + imag(zs.inbicd))) * ve;
            pfast2      = fpnbi2 .* profli.pnbi .* ((imag(zs.esup_nbi) ./ max(1,trapz(x,profli.vpr .* profli.pnbi,2))) * ve);
            nfast2      = max(1,pfast2 .* 2 ./ (phys.e .* option.einj));
            utheta_nbi2 = fjnbi2 .* profli.jnbicd ./ nfast2 ./ btot ./ phys.e;
            if gas_nbi == -1
                Unbi   = phys.mp .*  (ones(size(cons.ftnbi)) * ve) .* nfast1 .* utheta_nbi1 + ...
                    phys.mp .*  (ones(size(cons.ftnbi)) * ve) .* nfast2 .* utheta_nbi2;
            elseif gas_nbi == 3
                Unbi   = phys.mp .*  ((2 .* (1-real(cons.ftnbi)) + 3 .* real(cons.ftnbi)) * ve) .* nfast1 .* utheta_nbi1 + ...
                    phys.mp .*  ((2 .* (1-imag(cons.ftnbi)) + 3 .* imag(cons.ftnbi)) * ve) .* nfast2 .* utheta_nbi2;
            else
                % still OK for 5 and 11; no boron or He3 NBI injector planned
                Unbi   = phys.mp .*  ((2 .* (1-real(cons.ftnbi)) + 1 .* real(cons.ftnbi)) * ve) .* nfast1 .* utheta_nbi1 + ...
                    phys.mp .*  ((2 .* (1-imag(cons.ftnbi)) + 1 .* imag(cons.ftnbi)) * ve) .* nfast2 .* utheta_nbi2;
            end
            %disp('ut = all + nbi');
            
        else
            Unbi = zeros(size(utheta_main));
            %disp('ut = all');
        end
        
        if option.Sn_fraction > 0
            Utor    = phys.mp .*  (      nHp .* utheta_h + 2 .* nDp  .* utheta_d  +  ...
                3 .* nTp .* utheta_t + 4 .* nhep .* utheta_he +  ...
                aimp .* nz1p .* utheta_z1  + amax .* nz2p .* utheta_z2 +  ...
                (1 - option.Sn_fraction) .* 183.84 .* nwp  .* utheta_w + ...
                option.Sn_fraction .* 118.71 .* nwp  .* utheta_sn  + ...
                5 .* nbp .* utheta_b + 2 .* nhe3p .* utheta_he3);
        else
            Utor    = phys.mp .*  (      nHp .* utheta_h + 2 .* nDp  .* utheta_d  +  ...
                3 .* nTp .* utheta_t + 4 .* nhep .* utheta_he +  ...
                aimp .* nz1p .* utheta_z1  + amax .* nz2p .* utheta_z2 + 183.84 .* nwp  .* utheta_w + ...
                5 .* nbp .* utheta_b + 2 .* nhe3p .* utheta_he3);
        end
        %  	  figure(22);clf
        %  	  plot(x,nHp .* utheta_h + 2 .* nDp  .* utheta_d  + 3 .* nTp .* utheta_t,'r', ...
        %  	      x,4 .* nhep .* utheta_he,'b',x,nz1p .* utheta_z1,'c', x, nz2p .* utheta_z2,'m',x,183.84 .* nwp  .* utheta_w,'k',x,Unbi ./ phys.mp,'g');
        %  	  xlabel('x')
        %  	  ylabel('Utor');
        %  	  title('HDT = r, He = b , imp1 = c, imp2 = m, W = k, NBI = g');
        %  	  drawnow
        
        % the sign is due to definition of dpsidrho
        Utor    = - option.signe .* profli.fdia .* profli.r2i .* dpsidrho .* (Utor + Unbi);
        
    else
        % the sign is due to definition of dpsidrho
        Utor    = option.signe .* profli.fdia .* profli.r2i .* dpsidrho .* Mtor .* utheta_main;
        %disp('ut = main');
        
    end
else
    
    Utor    = zeros(size(profli.fdia));
    %disp('ut = 0');
end

% calcul du champ electrique radial (Er gradient(rho))
if option.Sn_fraction > 0
     Ptor    = phys.mp .* pdederive(x,profli.tip .* nHp,0,2,2,1) + ...
        2 .* phys.mp .* pdederive(x,profli.tip .* nDp,0,2,2,1) + ...
        3 .* phys.mp .* pdederive(x,profli.tip .* nTp,0,2,2,1) + ...
        2 .* phys.mp .* pdederive(x,profli.tip .* nhep,0,2,2,1) + ...
        3/2 .* phys.mp .* pdederive(x,profli.tip .* nhe3p,0,2,2,1) + ...
        aimp ./ zimp .* phys.mp .* pdederive(x,profli.tip .* nz1p,0,2,2,1) + ...
        amax ./ zmax .* phys.mp .* pdederive(x,profli.tip .* nz2p,0,2,2,1) + ...
        11/5 .* phys.mp .* pdederive(x,profli.tip .* nbp,0,2,2,1) + ...
        ((1 - option.Sn_fraction) .* 183.84 ./ max(1,z0wavez(profli.tep)) + ...
        option.Sn_fraction .* 118.71./ max(1,z0snavez(profli.tep))) .* ...
        phys.mp .* pdederive(x,profli.tip .* nwp,0,2,2,1);
   
else
    Ptor    = phys.mp .* pdederive(x,profli.tip .* nHp,0,2,2,1) + ...
        2 .* phys.mp .* pdederive(x,profli.tip .* nDp,0,2,2,1) + ...
        3 .* phys.mp .* pdederive(x,profli.tip .* nTp,0,2,2,1) + ...
        2 .* phys.mp .* pdederive(x,profli.tip .* nhep,0,2,2,1) + ...
        3/2 .* phys.mp .* pdederive(x,profli.tip .* nhe3p,0,2,2,1) + ...
        aimp ./ zimp .* phys.mp .* pdederive(x,profli.tip .* nz1p,0,2,2,1) + ...
        amax ./ zmax .* phys.mp .* pdederive(x,profli.tip .* nz2p,0,2,2,1) + ...
        11/5 .* phys.mp .* pdederive(x,profli.tip .* nbp,0,2,2,1) + ...
        183.84 ./ max(1,z0wavez(profli.tep)) .* phys.mp .* pdederive(x,profli.tip .* nwp,0,2,2,1);
end
Ptor    = Ptor ./ (profli.rmx(:,end) * ve);			
% utilisation de la formule 8.22 de Helander 	
% ici er = grad(phi_elect) = Er/ grad_rho		
% profli.er  = (Ptor + Rtor .* profli.r2i .* dpsidrho  - Utor)  ./ Mtor; => cette formule est fausse
profli.er  = (profli.omega .* dpsidrho + Ptor./ Mtor - Utor ./ Mtor) .* profli.grho;

%  figure(21);clf
%  plot(profli.rmx',profli.er' ./ profli.grho','r',profli.rmx',profli.omega' .* dpsidrho','b',profli.rmx',Ptor'./ Mtor','c',profli.rmx',-Utor' ./ Mtor','m')
%  xlabel('\rho')
%  ylabel('V/m')
%  title('E_r = r, Rotation = b , Pressure = c  , Poloidal = m');
%  drawnow

switch option.mode_vtheta
case 'same v_tor'

  % calculde la roration toroidale
  profli.vtor    = omega.* rmax;
  inter          = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (zimp .* nz1p);
  warning off
  omega_z1       = (profli.er - inter) ./ dpsidrho; 
  profli.utheta  = (profli.vtor    - omega_z1 .* rmax) ./ profli.fdia .* rmax .* option.signe;
  profli.utheta(:,1) = 0;
  warning on
  profli.vtheta  = profli.utheta .* bpol;

otherwise

  % calculde la roration toroidale
  inter          = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (zimp .* nz1p);
  warning off
  omega_z1       = (profli.er - inter) ./ dpsidrho; 
  profli.vtor    = omega_z1 .* rmax + profli.utheta .* option.signe .* profli.fdia ./ rmax;
  profli.vtor(:,1) = 2 .*  profli.vtor(:,2) - profli.vtor(:,3);
  warning on

end

% limite de vitesse a la vitesse du son dans le plasma
cslim    = sqrt(1.602176462e-19 .* (profli.tip + profli.zeff .* profli.tep)./1.6726485e-27 ./ aj);
indbad           = find(abs(profli.vtor) > cslim);
if ~isempty(indbad)
	profli.vtor(indbad)   = sign(profli.vtor(indbad)) .* cslim(indbad);
	profli.er(indbad)     = inter(indbad) + dpsidrho(indbad) .* (profli.vtor(indbad) -   ...
                                option.signe .* profli.utheta (indbad).* profli.fdia(indbad) ./ rmax(indbad)) ./ rmax(indbad);
end
% calcul de omega_(ExB)
% definition de nclass
inter            = profli.er ./ rmax  ./ max(eps,bpol) .* grho;
inter(:,1)       = 0;
profli.web       = abs(rmax .^ 2 .* bpol .^ 2 ./ btot .* pdederive(x,inter,0,2,2,1) ./ ...
                       max(eps,abs(pdederive(x,profli.psi,0,2,2,1))));
profli.web(:,1)  = 0;

% recopie pour amortissement
profli.omega  = omega;


% calcul du flux et de coefficient de transport
% les sources de moment externes
profli.rot_nbi     =              abs(real(psnbi)) .* ((real(snbi) ./ max(1,trapz(x,profli.vpr .* abs(real(psnbi)),2))) * ve) + ...
                     sqrt(-1) .*  abs(imag(psnbi)) .* ((imag(snbi) ./ max(1,trapz(x,profli.vpr .* abs(imag(psnbi)),2))) * ve);
% cacul du moment total
Rtor    = Mtor    .* profli.omega ./ profli.r2i;
Rtor(~isfinite(Rtor)) = 0;
rot_n0             = abs(profli.s0m .* Rtor);
profli.rot_n0      = rot_n0 .* ((zs.sn0fr ./ max(1,trapz(x,profli.vpr .* rot_n0,2))) * ve);
profli.rot_lh      = (abs(profli.jlh + profli.jeccd)) .* ((slh ./ max(1,trapz(x,profli.vpr .*  ...
                     (abs(profli.jlh + profli.jeccd)),2))) * ve);
% le ripple est inclu dans le transport 
% juste pour etre rigoureux , cette source est negligeable
rot_epar = abs(profli.epar) .* ((s_epar ./ max(1,trapz(x,profli.vpr .* abs(profli.epar),2))) * ve);
% le flux de moment
warning off
grot = cumtrapz(x, profli.vpr .* (rot_epar + real(profli.rot_nbi) + imag(profli.rot_nbi) +  ...
                profli.rot_n0 + profli.rot_lh)  - z0dxdt(Rtor .* profli.vpr,cons.temps) +   ...
               ((z0dxdt(profli.rmx(:,end),cons.temps) ./ profli.rmx(:,end)) * x) .* ...
                pdederive(x,Rtor .* profli.vpr,0,2,2,1),2)./ max(eps,profli.vpr_tor .* profli.grho2);
warning on
grot(:,1) = 0;
profli.frot        = grot;
profli.rtor        = Rtor;

% calcul des coefficient de transport associe
%Drot   = max(1e-3,min(30,max(max(profli.xii,profli.xie),profli.dn) .* (min(10,max(0.1,zs.taue ./ max(1e-6,taufr))) * ve)));
Drot   = profli.xii;
Rtord1 = pdederive(x,Rtor,0,2,2,1);
%  Drot   = max(profli.xie,profli.xii);
gedge  = trapz(x,Rtor .* profli.vpr,2) ./ profli.vpr_tor(:,end) ./ profli.grho2(:,end) ./ max(1e-6 ,taufr);
warning off
Dedge  = abs(gedge(:,end)  ./ (Rtord1(:,end) ./ profli.rmx(:,end)));
%warning on
Dedge(abs(Rtord1(:,end))<= eps) = Drot(abs(Rtord1(:,end))<= eps,end);
Drot(:,end) = Dedge;
Drot(Drot > 30)   = 30;
Drot(Drot < 1e-3) = 1e-3;
warning off
Vrot = - (grot + Drot .* Rtord1 ./ (profli.rmx(:,end) * ve)) ./ Rtor;
warning on
Vrot(abs(Rtor) <= eps) = 0;
Vrot(abs(Rtord1(:,end))<= eps) = 0;
Vrot = min(300,max(-300,Vrot));
profli.drot        = Drot;
profli.vrot        = Vrot;

% recopie pour amortissement
profli.omega  = omega;

%  %  
%      figure(21);
%    plot(cons.temps,snbi,'r',cons.temps,sturb,'b',cons.temps,slh,'g',cons.temps,sripicrh,'m',cons.temps,sriplh,'c',cons.temps,sripth,'k',cons.temps,s_epar,'r-.')
%  %  plot(x,alpha,'r',x,beta,'b')
%      drawnow
%  %    keyboard
