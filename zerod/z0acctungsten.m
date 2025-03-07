% cacul de l'accumulation de tungstene
% [nwp,nwm,zu1w,zu2w] = z0acctungsten(post.z0dinput.option,post.z0dinput.geo,post.z0dinput.cons,post.zerod,post.profil0d)
function [nwp,nwm,zu1w,zu2w,tleak,fw,fwth_fit,ya0,cw_edge,fraction_prompt,zwavep] = z0acctungsten(option,geo,cons,zerod,profil0d)

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



% compatibilite mode fast pour plot
if isfield(profil0d,'temps') && (length(cons.temps) ~= length(profil0d.temps))
        times = cons.temps;
        temps = profil0d.temps;
	% modification des donnees
	noms = fieldnames(cons);
	for l=1:length(noms)
		nomc = noms{l};
		val  = getfield(cons,nomc);
		valn  = interp1(times,val,temps,'linear');
		indbad      = find(any(~isfinite(valn),2));
		if ~isempty(indbad)
			valn(indbad,:) = ones(length(indbad),1) * val(end,:);
		end
		cons = setfield(cons,nomc,valn);
	end
	cons.temps = temps;
	
	noms = fieldnames(geo);
	for l=1:length(noms)
		nomc = noms{l};
		val  = getfield(geo,nomc);
		if ~isempty(val)
			valn  = interp1(times,val,temps,'linear');
			indbad      = find(any(~isfinite(valn),2));
			if ~isempty(indbad)
				valn(indbad,:) = ones(length(indbad),1) * val(end,:);
			end
			geo = setfield(geo,nomc,valn);
		end
	end
	
	noms = fieldnames(zerod);
	for l=1:length(noms)
		nomc = noms{l};
		val  = getfield(zerod,nomc);
		if length(val) == length(times)
			val  = interp1(times,val,temps,'nearest');
			zerod = setfield(zerod,nomc,val);
		end
	end
	
end

% adding sources of tungsten due to heating systems
cw_offset = option.cw_offset + option.cw_ecrh .* abs(cons.pecrh ./ 1e6) + option.cw_icrh .* abs(cons.picrh ./ 1e6) +  ...
            option.cw_lhcd .* abs(cons.plh  ./ 1e6) + option.cw_nbi1 .* real(cons.pnbi)  ./ 1e6 + option.cw_nbi1 .* imag(cons.pnbi)  ./ 1e6;

x      = profil0d.xli;
vt     = ones(size(profil0d.tep,1),1);
ve     = ones(size(x));
rext         = profil0d.Raxe + geo.a * x;
btor         = profil0d.fdia ./ rext;
grho         = abs((profil0d.rmx(:,end) * ve)./ max(eps,pdederive(x,rext,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,profil0d.psi,0,2,2,1) ./ rext .* grho ./ (profil0d.rmx(:,end) * ve);
ut           = atan(abs(bpol(:,end) ./  btor(:,end)));

% longueur de ionisation pres des plaques
[svi1s,svi2s,svcx,svrec,sii,sss,Ass] = z0sectionh(zerod.telim,zerod.telim .* zerod.tibord ./ zerod.tebord);
%% equilibre entre 1s et 2s pour le neutres de centre 
alphas = zerod.nelim .* sss ./ ( zerod.nelim .* sss + Ass);
% etat d'equilibre 1s/2s
sviss  = svi2s .* alphas + svi1s .* (1 - alphas);
s3  = sviss  + sii + svcx;
delta_s = sqrt(2 .* phys.e .* zerod.telim  ./ zerod.meff ./ phys.mp) ./ zerod.nelim ./  s3 ./ ut;

% temperature de fuite
% from ref :DIVIMP simulation ..., A. Jarvinen et al, Physica Sripta T145 (2011) p 014013-
% utilisation de (z0wavez(te_bord) +  z0wavez(te_lim)) ./ 2
if option.Sn_fraction > 0
    zwsn_u     = (1 - option.Sn_fraction) .* z0wavez(zerod.tebord) + option.Sn_fraction .* z0snavez(zerod.tebord) ;
    zwsn_t     = (1 - option.Sn_fraction) .* z0wavez(zerod.telim) + option.Sn_fraction .* z0snavez(zerod.telim) ;   
    tleak      = 3.8e-9 .* (abs(option.ftwleak) .* zwsn_u + (1 - abs(option.ftwleak)) .* zwsn_t) .* sqrt(zerod.nelim .* delta_s);
else
    tleak      = 3.8e-9 .* (abs(option.ftwleak) .* z0wavez(zerod.tebord) + (1 - abs(option.ftwleak)) .* z0wavez(zerod.telim)) .* sqrt(zerod.nelim .* delta_s);
end
fw         = max(1e-308,exp(- (tleak ./ zerod.telim) .^ 2));
if option.ftwleak < 0
	% from ref :DIVIMP simulation ..., A. Jarvinen et al, Physica Sripta T145 (2011) p 014013-
	% etalonage de la fuite
	fw_fit    = [10 1 0.3253 1.5101e-16 1e-45];
	fw_leak   = [1 1 0.035  4.0000e-06 1e-308];
	fwth_fit  = min(1,exp(interp1(log(fw_fit),log(fw_leak),log(fw),'pchip','extrap')));
	fwth_fit(~isfinite(fwth_fit)) = 0;
	fwth_fit(fw < 1e-45) = 0;
else
	fwth_fit = fw;
end

switch option.gaz
    case 1
        zg = 1;
    case 2
        zg = 1;
    case 3
        zg = 1;
    case 4
        zg = 2;
    case 5
        zg = mean((1  + 4 .* cons.iso)   ./  (1 + 2 .* cons.iso));
    case 11
        zg = mean((1  + 25 .* cons.iso)   ./ (1 + 5.* cons.iso));
end

switch option.yield_model
case 'fit'
    % if Sn
    if option.Sn_fraction > 0
        warning('No fit for Sn evaporation'),
    end
    
	% ne donne pas la bonne valeur !
    switch option.gaz
        case   1
            ya0 = z0yield('H','C',zerod.telim);
        case 2
            ya0 = z0yield('D','C',zerod.telim);
        case 3
            ya0 = z0yield('D','C',zerod.telim) ./ (1 + cons.iso) + z0yield('T','C',zerod.telim) .* cons.iso ./ (1 + cons.iso);
        case 5
            ya0 = z0yield('D','C',zerod.telim) ./ (1 + cons.iso) + z0yield('He','C',zerod.telim) .* cons.iso ./ (1 + cons.iso);
        case 11
            error('This option is not valid for boron');
        case 4
            ya0 = z0yield('He','C',zerod.telim);
    end
	% from ref :DIVIMP simulation ..., A. Jarvinen et al, Physica Sripta T145 (2011) p 014013-
	% etalonage de la concentration
	te_lim 	= [106.9 	7.9 		3.8];
	Cw_div  = [1.8e-4 	5.6e-8		8.5e-10];
	ya0_w = exp(interp1(log(cat(2,1e6,te_lim,300 .* phys.k .* phys.e)),log(cat(2,1,Cw_div,1e-308)),log(zerod.telim),'pchip','extrap'));
	ya_0  = ya0 ./ max(1e-308,z0yield('D','W',zerod.telim)) .* ya0_w;
	ya_0(~isfinite(ya_0)) = 0;
    ya0 = ya_0;

case 'Javev'

    % initialisation a 0
    ya0 = zeros(size(zerod.nelim));
    % concentration de l'impurete zimp
    czimp = profil0d.nzp(:,end) ./ profil0d.nep(:,end);
    % concentration de l'imurete zmax
    czmax = 0.01 .* abs(option.fzmax_div) + option.rimp .* czimp;
    % la concentration pour le self sputtering : prise en compte de ce qui n'est pas redepose sur place.
    cws =  profil0d.nwp(:,end) ./ profil0d.nep(:,end);
    % concentration en He
    switch option.gaz
        case 5
            che  =  option.frhe0 .* profil0d.nep(:,end);
            che3 =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
        otherwise
            che  =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
            che3 =  zeros(size(profil0d.nep(:,end)));
    end
    % le reste pour HDT
    if option.Sn_fraction > 0
        zwsn_u     = (1 - option.Sn_fraction) .* z0wavez(zerod.tebord) + option.Sn_fraction .* z0snavez(zerod.tebord) ;
        creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* zwsn_u - 2 .* che - 2 .* che3);
    else
        creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zerod.tebord) - 2 .* che - 2 .* che3);
    end
    % H D T
    switch option.gaz
        case 11
            cB  = zerod.nTm ./zerod.nem;
            creste = max(0,crest - 5 .* cB);
            cT  = zeros(size(zerod.n1m));
            cD  = creste .* zerod.nDm ./zerod.n1m;
            cH  = creste .* max(0,zerod.n1m - zerod.nDm) ./zerod.n1m;
        otherwise
            cT  = creste .* zerod.nTm ./zerod.n1m;
            cD  = creste .* zerod.nDm ./zerod.n1m;
            cB  = zeros(size(zerod.n1m));
            cH  = creste .* max(0,zerod.n1m - zerod.nDm - zerod.nTm) ./zerod.n1m;
    end
           
    % le yield effectif est la somme des yields ponderee :
    % il y une seule temperature sur le limiteur
    ya0  = czimp .* z0yield2(zerod.telim,zerod.telim,option.zimp,'W',1) + czmax .* z0yield2(zerod.telim,zerod.telim,option.zmax,'W',1) + ...
        cws   .* z0yield2(zerod.telim,zerod.telim,'W','W',1) + che .* z0yield2(zerod.telim,zerod.telim,'He','W',1) + ...
        cH    .* z0yield2(zerod.telim,zerod.telim,'H','W',1) + cD  .* z0yield2(zerod.telim,zerod.telim,'D','W',1)  + ...
        cT    .* z0yield2(zerod.telim,zerod.telim,'T','W',1) + che3 .* z0yield2(zerod.telim,zerod.telim,'He3','W',1) + ...
        cB    .* z0yield2(zerod.telim,zerod.telim,'B','W',1);
    
    % if Sn
    if option.Sn_fraction > 0
       ya0  = (1 - option.Sn_fraction) .* ya0  + option.Sn_fraction .* ( ...
           czimp .* z0yield2(zerod.telim,zerod.telim,option.zimp,'Sn',1) + czmax .* z0yield2(zerod.telim,zerod.telim,option.zmax,'Sn',1) + ...
           cws   .* z0yield2(zerod.telim,zerod.telim,'Sn','Sn',1) + che .* z0yield2(zerod.telim,zerod.telim,'He','Sn',1) + ...
           cH    .* z0yield2(zerod.telim,zerod.telim,'H','Sn',1) + cD  .* z0yield2(zerod.telim,zerod.telim,'D','Sn',1)  + ...
           cT    .* z0yield2(zerod.telim,zerod.telim,'T','Sn',1)  + che3 .* z0yield2(zerod.telim,zerod.telim,'He3','Sn',1) + ...
           cB    .* z0yield2(zerod.telim,zerod.telim,'B','Sn',1) +  ...
           (1 - option.Sn_fraction) .* cws .* (z0yield2(zerod.telim,zerod.telim,'W','Sn',1) + z0yield2(zerod.telim,zerod.telim,'Sn','W',1)));
    end
case 'Matsunami'

    % initialisation a 0
    ya0 = zeros(size(zerod.nelim));
    % concentration de l'impurete zimp
    czimp = profil0d.nzp(:,end) ./ profil0d.nep(:,end);
    % concentration de l'imurete zmax
    czmax = 0.01 .* abs(option.fzmax_div) + option.rimp .* czimp;
    % la concentration pour le self sputtering : prise en compte de ce qui n'est pas redepose sur place.
    cws =  profil0d.nwp(:,end) ./ profil0d.nep(:,end);
    % concentration en He
    switch option.gaz
        case 5
            che  =  option.frhe0 .* profil0d.nep(:,end);
            che3 =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
        otherwise
            che  =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
            che3 =  zeros(size(profil0d.nep(:,end)));
    end
    % le reste pour HDT
    if option.Sn_fraction > 0
        zwsn_u     = (1 - option.Sn_fraction) .* z0wavez(zerod.tebord) + option.Sn_fraction .* z0snavez(zerod.tebord) ;
        creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* zwsn_u - 2 .* che);
    else
        creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zerod.tebord) - 2 .* che);
    end
    % H D T
    switch option.gaz
        case 11
            cB  = zerod.nTm ./zerod.nem;
            creste = max(0,crest - 5 .* cB);
            cT  = zeros(size(zerod.n1m));
            cD  = creste .* zerod.nDm ./zerod.n1m;
            cH  = creste .* max(0,zerod.n1m - zerod.nDm) ./zerod.n1m;
        otherwise
            cT  = creste .* zerod.nTm ./zerod.n1m;
            cD  = creste .* zerod.nDm ./zerod.n1m;
            cB  = zeros(size(zerod.n1m));
            cH  = creste .* max(0,zerod.n1m - zerod.nDm - zerod.nTm) ./zerod.n1m;
    end
    
    % le yield effectif est la somme des yields ponderee :
    % il y une seule temperature sur le limiteur
    ya0  = czimp .* z0yield2(zerod.telim,zerod.telim,option.zimp,'W',0) + czmax .* z0yield2(zerod.telim,zerod.telim,option.zmax,'W',0) + ...
        cws   .* z0yield2(zerod.telim,zerod.telim,'W','W',0) + che .* z0yield2(zerod.telim,zerod.telim,'He','W',0) + ...
        cH    .* z0yield2(zerod.telim,zerod.telim,'H','W',0) + cD  .* z0yield2(zerod.telim,zerod.telim,'D','W',0)  + ...
        cT    .* z0yield2(zerod.telim,zerod.telim,'T','W',0) + che3 .* z0yield2(zerod.telim,zerod.telim,'He3','W',0) + ...
        cB    .* z0yield2(zerod.telim,zerod.telim,'B','W',0);
    % if Sn
    if option.Sn_fraction > 0
       ya0  = (1 - option.Sn_fraction) .* ya0  + option.Sn_fraction .* ( ...
           czimp .* z0yield2(zerod.telim,zerod.telim,option.zimp,'Sn',0) + czmax .* z0yield2(zerod.telim,zerod.telim,option.zmax,'Sn',0) + ...
           cws   .* z0yield2(zerod.telim,zerod.telim,'Sn','Sn',0) + che .* z0yield2(zerod.telim,zerod.telim,'He','Sn',0) + ...
           cH    .* z0yield2(zerod.telim,zerod.telim,'H','Sn',0) + cD  .* z0yield2(zerod.telim,zerod.telim,'D','Sn',0)  + ...
           cT    .* z0yield2(zerod.telim,zerod.telim,'T','Sn',0) + che3 .* z0yield2(zerod.telim,zerod.telim,'He3','W',0) + ...
           cB    .* z0yield2(zerod.telim,zerod.telim,'B','W',0) +  ...
           (1 - option.Sn_fraction) .* cws .* (z0yield2(zerod.telim,zerod.telim,'W','Sn',0) + z0yield2(zerod.telim,zerod.telim,'Sn','W',0)));
    end

otherwise
	error('unknown sputtering yield model')
end
% concentration au bord
% utilisation de la formule 6.81 p 2196 du livre de Stangeby : "the plasma boundary ..."
% uniquement en divertor, le limiteur est suppose ne pas provoquer d'accumulation de tungstene.
cw_edge = ya0 .* (min(1,fwth_fit) .* zerod.xpoint);

% model de redepostion 
% ref : Stangeby NF 51 (2011) 063001
% ref D. Naujoks et al, NF 36 (1996)  p671-687
% la section efficace de ionisation de W est complexe  a calculer avec ions + electrons + effet cinetique
tp   = [0 	1 	10 	100 	1000 	1e5];
svp  = [1e-10	1e-8	1.5e-7  5e-7	2e-7	1e-7];
svwi = pchip(tp,svp,zerod.telim) .* 1e-6; % m^-3  * s^-1
% self sputtering dominant ou par impurete legere : E0 = telim
lioniz   =  sqrt(2 .* zerod.telim .* phys.e ./ (183.84 .* phys.ua)) ./ zerod.nelim ./  svwi;
% 1 fois ionise
larmor_w =  4.57e-3 .* sqrt(183.84) .* sqrt(zerod.telim ./ 1e3) ./ geo.b0;
p        = lioniz ./ larmor_w;
% from ERO with all effect
pp       = [0 	0.2 	0.6 	2.2 	3.2 	6	1000];
fp       = [1	0.85	0.65	0.3	0.25	0.15	0.1];   
fraction_prompt =  (1 - option.Sn_fraction) .* min(1-eps,max(0.1,pchip(pp,fp,max(0,min(p,1000)))));
if option.cw_factor < 0
	cw_edge = cw_edge .* (1 - fraction_prompt);
end

% Calcul du profil n_w avec la formule 5.9 et 5.10 du livre de Helander : "Collisional Transport ..."
[zwavep,zwavep2] = z0wavez(profil0d.tep);
if option.Sn_fraction > 0
   [zsnavep,zsnavep2] = z0snavez(profil0d.tep);
   zwavep = (1 - option.Sn_fraction) .* zwavep + option.Sn_fraction .* zsnavep;
   zwavep2 = (1 - option.Sn_fraction) .* zwavep2 + option.Sn_fraction .* zsnavep2;
end

% effet du piedestal sur le derivees (comme pour le bootstrap)
indh    = find(zerod.modeh);
nipd1   = pdederive(x,profil0d.nip,0,2,2,1);
tipd1   = pdederive(x,profil0d.tip,0,2,2,1);
omegad1 = pdederive(x,profil0d.omega ,0,2,2,1) ;
if ~isempty(indh) && (option.gradient_Wacc ~= 0)
	dx            = x(end) - x(end-1);
	nipd1(indh,end-1)   = option.gradient_Wacc .* (profil0d.nip(indh,end) - profil0d.nip(indh,end-1)) ./ dx;
	tipd1(indh,end-1)   = option.gradient_Wacc .* (profil0d.tip(indh,end) - profil0d.tip(indh,end-1)) ./ dx;
	omegad1(indh,end-1) = option.gradient_Wacc .* (profil0d.omega(indh,end) - profil0d.omega(indh,end-1)) ./ dx;
end


% terme du au flux d echaleur electronique (decontamination)
%lnldi   = 17.3 - 0.5 .* log(profil0d.nep./1e20) + log(profil0d.tep ./1e3);
%Aw = 184;
%tau_iz  = 6.6e17 .* sqrt(Aw) .* (profil0d.tip ./1e3) .^(3/2) ./ profil0d.nip ./ lnldi ./ zwavep .^ 4;
%btot    = sqrt(profil0d.bpol .^ 2 + (profil0d.fdia ./ profil0d.Raxe) .^ 2);
%omega_i =  zg .* phys.e .* btot ./ phys.mp ./ (zerod.meff * ve);
% le modele complet depend du rapport de la convection ??? la diffusion.
% il est trop complexe pour etre implemente dans METIS.
% ref : R. Dux et al, PPCF 45 (2003) p 1815-1825
% ref : M. . Tokar et al, Nucl. Fux 37 (1997) p 1691- 
% Nous utilisons l'expression classique (en cylindrique) + un terme decontaminant provenant du chauffage des electrons
% le flux depend du transfert de chaleur de la source electronique vers les ions, 
%  mais il est translater en equivalent d'un effet de flux d'ion (modification de formule 5.9).
%Ke = profil0d.qei .* phys.mp .* (zerod.meff * ve) .* omega_i .^ 2 .* tau_iz ./ (zg .* phys.e .^ 2 .* profil0d.nip .* profil0d.tip  .* profil0d.tep); 
%fei = ((profil0d.qe -  profil0d.qei)  - (profil0d.qi  + profil0d.qei)) ./ max(1,profil0d.qe + profil0d.qi);
%fei = abs(profil0d.qei) ./ max(1,profil0d.qe + profil0d.qi);
fei = profil0d.qei ./ max(1,profil0d.qe + profil0d.qi);
fei(:,1) = fei(:,2);
fei = tanh(fei);
%figure(21);clf;plot(cons.temps,fei);drawnow
% le terme de rotation vient de :
% C. angioni et al, Nuclear Fusion 54 (2014) 083028
% Y. Camenen et al, Physics of Plasmas 16, 012503 (5009)
vthi = sqrt(2 .* phys.e .* profil0d.tip ./ (183.84 .* phys.ua));
%Cu   = option.rot_acc .* ( 2 .* 183.84 ./ (zerod.meff * ve) - zwavep);
% version normalisee (pour etre compatible avec le terme en Ti et ni, de l'ordre de 1)
Cu   = option.rot_acc .* ( 2 .* 183.84 ./ (zerod.meff * ve) - zwavep) ./  (2 .* 183.84);
%up   = - pdederive(x,profil0d.omega ,0,2,2,1) .* profil0d.Raxe ./ vthi ./ profil0d.epsi;
% probleme de normalisation (on cherche grad(X)/x au lieu R/L_X)
% il est decomtaminent avec un profil de rotation decroissant du centre au bord (typique avec NBI).
%up   = -pdederive(x,profil0d.omega ,0,2,2,1) ./ vthi ./ profil0d.Raxe;
up   = - omegad1 ./ vthi ./ profil0d.Raxe;
up(:,1) = 0;
%  figure(21);
%  clf;
%  plot(cons.temps,Cu(:,1:19) .* up(:,1:19),'r', ...
%       cons.temps,pdederive(x(:,1:19),profil0d.nip(:,1:19),0,2,2,1) ./ profil0d.nip(:,1:19),'g', ...
%       cons.temps,(1 - 1 ./ zwavep(:,1:19) - (3 ./ 2) - option.heat_acc .* fei(:,1:19)) .* pdederive(x(:,1:19),profil0d.tip(:,1:19),0,2,2,1) ./ profil0d.tip(:,1:19),'b');
%  drawnow

switch option.acc_col
case 'on'
    if isfield(profil0d,'nwp')
        %ref: S. Breton PoP 25 (2018) 012303
        %figure(21);clf;
        meffp         = phys.ua .* zerod.meff * ones(size(profil0d.xli));
        lnl_i         = 17.3 - 0.5 .* log(profil0d.nip ./ 1e20) + 1.5 .* log(profil0d.tip ./ 1e3);
        tau_main_main = 6.6e17 .* sqrt(meffp ./ phys.ua) .* (profil0d.tip ./ 1e3) .^ (3/2) ./ profil0d.nip ./ lnl_i ./ zg .^ 4;
        %subplot(2,2,1);semilogy(profil0d.xli,tau_main_main);
        nu_main_main  = 1 ./ tau_main_main;
        vth_main      = sqrt(2 .* profil0d.tip .* phys.e ./ meffp);
        nustar_main   = nu_main_main .* profil0d.qjli .* profil0d.Raxe ./ vth_main ./ profil0d.epsi .^ (3/2);
        nustar_main(:,1) = nustar_main(:,2);
        %subplot(2,2,2);semilogy(profil0d.xli,nustar_main);
        alpha_z       = profil0d.nwp .* zwavep2 ./ profil0d.nip .* zg .^ 2;
        %subplot(2,2,3);plot(profil0d.xli,alpha_z);
        g             = nustar_main .* profil0d.epsi .^ (3/2);
        K             = 1 - 0.52 .* alpha_z ./ (0.59 + alpha_z + 1.34 ./ g .^ 2);
        H             = -0.5 + (0.29  + 0.68 .* alpha_z) ./ (0.59 + alpha_z + 1.34 ./ g .^ 2);
        %subplot(2,2,4);plot(profil0d.xli,K,'r',profil0d.xli,H,'b');drawnow
        f             = Cu .* up  + K .* nipd1 ./ profil0d.nip + (H - 1 ./ zwavep - option.heat_acc .* fei) .* tipd1 ./ profil0d.tip;
        % effect of anormal transport
        faccu    = 1;
        if option.faccu < 0
            D_W      = profil0d.xie;
        else
            D_W      = profil0d.xii;
        end
        wci_main = zg .* phys.e .* sqrt((profil0d.fdia ./ profil0d.Raxe) .^ 2 +profil0d.bpol .^ 2) ./ meffp;
        rho_i    = vth_main ./ wci_main;
        Dc       = rho_i .^ 2 ./ tau_main_main .* profil0d.qjli .^ 2;
        KmD      = K + abs(option.faccu) .* D_W ./ Dc;
        %      figure(21);clf;
        %      subplot(2,1,1)
        %      plot(profil0d.xli,D_W,'b',profil0d.xli,Dc,'r');
        %      subplot(2,1,2)
        %      plot(profil0d.xli,K,'b',profil0d.xli,KmD,'r');
        %      drawnow
        KmD      = max(eps,KmD);
    else
        % formule 5.9 modifiee + terme de rotation
        %f    = Cu .* up  + pdederive(x,profil0d.nip,0,2,2,1) ./ profil0d.nip + (1 - 1 ./ zwavep - (3 ./ 2) - option.heat_acc .* fei) .* pdederive(x,profil0d.tip,0,2,2,1) ./ profil0d.tip;
        f    = Cu .* up  + nipd1 ./ profil0d.nip + (1 - 1 ./ zwavep - (3 ./ 2) - option.heat_acc .* fei) .* tipd1 ./ profil0d.tip;
        KmD    = 1;
        faccu = option.faccu;
    end
    otherwise
        % formule 5.9 modifiee + terme de rotation
        %f    = Cu .* up  + pdederive(x,profil0d.nip,0,2,2,1) ./ profil0d.nip + (1 - 1 ./ zwavep - (3 ./ 2) - option.heat_acc .* fei) .* pdederive(x,profil0d.tip,0,2,2,1) ./ profil0d.tip;
        f    = Cu .* up  + nipd1 ./ profil0d.nip + (1 - 1 ./ zwavep - (3 ./ 2) - option.heat_acc .* fei) .* tipd1 ./ profil0d.tip;
        KmD    = 1;
        faccu = option.faccu;
end
%  figure(21);clf
%  %  %plot3((cons.temps*ve)',(vt*x)',fei');
%  plot(cons.temps,cw_edge);
%    drawnow


fact     = zwavep ./ KmD .* f;
inte_f   = cumtrapz(x(end:-1:1),fact(:,end:-1:1),2);
inte_f   = inte_f(:,end:-1:1);

%figure(22);clf;plot(cons.temps,inte_f(:,1:19),'b',cons.temps,f(:,1:19),'r')
% valeur finale
% cas de donnees experimentales
mode_exp_profile_provided = false;
if isappdata(0,'NW_SHAPE_EXP') & isfield(profil0d,'nep');
	nwp_exp   = getappdata(0,'NW_SHAPE_EXP');
	if all(nwp_exp.nwp(:)< 0)
	    nwp_shape = max(1,interp1_ex(nwp_exp.temps,abs(nwp_exp.nwp),cons.temps,'nearest','extrap'));
	    neshape   = pchip(nwp_exp.x,abs(nwp_shape),profil0d.xli);
	    mode_exp_profile_provided = true;
	elseif all(nwp_exp.nwp(:)< 1e13)
	    nwp_shape = max(eps,interp1_ex(nwp_exp.temps,nwp_exp.nwp,cons.temps,'nearest','extrap'));
	    neshape   = pchip(nwp_exp.x,nwp_shape,profil0d.xli);
	    neshape   = (neshape ./ max(eps,neshape(:,end) * ve));
	else
	    nwp_shape = max(1e13,interp1_ex(nwp_exp.temps,nwp_exp.nwp,cons.temps,'nearest','extrap'));
	    neshape   = pchip(nwp_exp.x,nwp_shape,profil0d.xli);
	    neshape   = (neshape ./ max(1,neshape(:,end) * ve));	
	end
else
	neshape  = (profil0d.nep ./ max(1,profil0d.nep(:,end) * ve)) .^ option.fne_acc;
end
if mode_exp_profile_provided == true;
	nwp      = neshape .* exp(inte_f .* faccu);
elseif length(option.cw_factor) > 1
	nwp      = (zerod.xpoint * ve) .* abs(option.cw_factor*ve) .* neshape .* ((zerod.nebord .* cw_edge) * ve) .* exp(inte_f .* faccu)  + ...
   		   (cw_offset * ve) .* neshape .* (zerod.nebord  * ve) .* exp(inte_f .* faccu);
else
	nwp      = (zerod.xpoint * ve) .* abs(option.cw_factor) .* neshape .* ((zerod.nebord .* cw_edge) * ve) .* exp(inte_f .* faccu)  + ...
   		   (cw_offset * ve) .* neshape .* (zerod.nebord  * ve) .* exp(inte_f .* faccu);
end

%  switch option.delay_wacc
%  case 1
%    for k = 1:21
%         nwp(:,k) = sgolayfilt(nwp(:,k),1,3);
%    end
%  end
%  
%  figure(32);clf
%  subplot(2,2,1)
%  plot(x,neshape);
%  subplot(2,2,2)
%  plot(cons.temps,cw_edge);
%  subplot(2,2,3);
%  plot(x,inte_f);
%  subplot(2,2,4);
%  plot(x,nwp)
%  drawnow

% securite
nwp = max(1,min(profil0d.nep ./ 74,nwp));

% valeur moyenne
nwm  = trapz(x,nwp .* profil0d.vpr,2) ./ trapz(x,profil0d.vpr,2);

zu1w = trapz(x,zwavep .* nwp .* profil0d.vpr,2) ./ trapz(x,profil0d.nzp .* profil0d.vpr,2);
zu2w = trapz(x,zwavep .^ 2 .* nwp .* profil0d.vpr,2) ./ trapz(x,profil0d.nzp .* profil0d.vpr,2);


%  figure(21)
%  plot(cons.temps,nwm,cons.temps,zerod.nwm);
%  hold on
%  drawnow

%figure(21);clf
%plot(cons.temps,trapz(x,zwavep .^ 2 .* nwp .* profil0d.vpr,2) ./ trapz(x,profil0d.nep .* profil0d.vpr,2));
%hold on
%drawnow

% verification du flux nul 
return
pion = profil0d.nip .* profil0d.tip .* phys.e;
pw   = nwp .* profil0d.tip .* phys.e;
flux = pdederive(x,pion,0,2,2,1) ./ pion  - pdederive(x,pw,0,2,2,1) ./ pw ./ zwavep - 3./2 .* pdederive(x,profil0d.tip,0,2,2,1) ./ profil0d.tip;
figure(21)
semilogy(x,abs(flux)./(zerod.n0a*ve));
drawnow
%keyboard


