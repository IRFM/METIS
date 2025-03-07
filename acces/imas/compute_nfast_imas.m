function nfast  = compute_nfast_imas(data_zerod,profil0d,option,cons,model)

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


% 0 template
nfast = zeros(size(model));

% resample data zerod
temps = data_zerod.temps;
if length(temps ) > 1
    noms = fieldnames(data_zerod);
    for k=1:length(noms)
        var = data_zerod.(noms{k});
	if length(temps) == length(var)
		data_zerod.(noms{k}) = interp1(temps,var,profil0d.temps,'linear','extrap');
        end 
    end
    data_zerod.temps = profil0d.temps;
end
% resample cons
temps = cons.temps;
if length(temps ) > 1
    noms = fieldnames(cons);
    for k=1:length(noms)
        var = cons.(noms{k});
	if length(temps) == length(var)
		cons.(noms{k}) = interp1(temps,var,profil0d.temps,'linear','extrap');
        end 
    end
    cons.temps = profil0d.temps;
end

% useful vector
ve  = ones(size(profil0d.xli));
vt  = ones(size(profil0d.temps));

% ATTENTION : normalisation des mass sur le gaz principal 
% density of fast ions
psup_alpha  = profil0d.pfus;
psup_alpha  = psup_alpha ./ (max(1, trapz(profil0d.xli,psup_alpha .* profil0d.vpr,2)) * ve) .*  (data_zerod.esup_fus * ve)  ./ (3/2);
nfast_alpha = zerod_fast_ions_density(profil0d.nep ,profil0d.tep ,profil0d.zeff , ...
	      data_zerod.meff * ve,4,2,3.56e6,psup_alpha);
talpha      = psup_alpha ./ max(1,nfast_alpha) ./  phys.e;

switch gas_nbi
    case -1
        minj = ones(size(cons.ftnb));
    case 11
        minj = 1 .* (1-cons.ftnbi) + 11 .* cons.ftnbi;
    case {3,5}
        minj = 2 .* (1-cons.ftnbi ) + 3 .* cons.ftnbi ;
    otherwise
        minj = 2 .* (1-cons.ftnbi ) + 1 .* cons.ftnbi ;
end
psup_nbi1  = real(profil0d.nbinesource );
psup_nbi1  = psup_nbi1 ./ (max(1, trapz(profil0d.xli,psup_nbi1 .* profil0d.vpr ,2)) * ve) .*  (real(data_zerod.esup_nbi )* ve) ./ (3/2);
nfast_nbi1 = zerod_fast_ions_density(profil0d.nep ,profil0d.tep ,profil0d.zeff , ...
	      data_zerod.meff * ve,minj*ve,1,option.einj,psup_nbi1);

psup_nbi2  = imag(profil0d.nbinesource );
psup_nbi2  = psup_nbi2 ./ (max(1, trapz(profil0d.xli,psup_nbi2 .* profil0d.vpr ,2)) * ve) .*  (imag(data_zerod.esup_nbi ) * ve) ./ (3/2);
nfast_nbi2 = zerod_fast_ions_density(profil0d.nep ,profil0d.tep ,profil0d.zeff , ...
	     data_zerod.meff * ve ,minj*ve,1,option.einj2,psup_nbi2);

tnbi      = (psup_nbi1 + psup_nbi2) ./ max(1,nfast_nbi1 +nfast_nbi2) ./  phys.e;

% choix du minoritaire
switch option.mino
    case 'He3'
        ag = 3;
        zg = 2;
        lg = 7.92e-3;
    case 'T'
        ag = 3;
        zg = 1;
        lg = 7.92e-3;
    case 'He4'
        ag = 4;
        zg = 2;
        lg = 4.55e-3;
    case 'D'
        ag = 2;
        zg = 1;
        lg = 6.46e-3;
    case 'B'
        ag = 11;
        zg = 5;
        lg = 4.576e-3 .* sqrt(11) / 5;
    otherwise
        ag = 1;
        zg = 1;
        lg = 4.576e-3;
end
psup_icrh   = profil0d.picrh ;
psup_icrh   = psup_icrh ./ (max(1, trapz(profil0d.xli,psup_icrh .* profil0d.vpr ,2)) * ve) .*  (data_zerod.esup_icrh  *ve) ./ (3/2);
nfast_icrh = zerod_fast_ions_density(profil0d.nep ,profil0d.tep ,profil0d.zeff , ...
	     data_zerod.meff * ve,ag,zg,data_zerod.einj_icrh*ve ,psup_icrh);
ticrh       = psup_icrh ./ max(1,nfast_icrh) ./  phys.e;
warning off
% pas de meilleurs formulation pour le moment
nfast_lh  = max(0,data_zerod.esup_lh ) ./ (data_zerod.einj_lh  .* phys.e) .* 2 ./ data_zerod.vp ;
warning on
psup_lh   = profil0d.plh ;
nfast_lh  = psup_lh ./ (max(1, trapz(profil0d.xli,psup_lh .* profil0d.vpr ,2)) * ve) .*  (nfast_lh * ve);
psup_lh   = psup_lh ./ (max(1, trapz(profil0d.xli,psup_lh .* profil0d.vpr ,2)) * ve) .*  (data_zerod.esup_lh * ve)  ./ (3/2);
tlh       = psup_lh ./ max(1,nfast_lh) ./  phys.e;


% hydrogenoid density
nDp   = zeros(size(profil0d.n1p));
nTp   = zeros(size(profil0d.n1p));
nHp   = zeros(size(profil0d.n1p));

% helium
nHep4  = zeros(size(profil0d.nhep));
nHep3  = zeros(size(profil0d.nhep));

% Boron
nBp   = zeros(size(profil0d.n1p));


% compute corrected density depending on minority scheme
switch option.mino
    case 'H'
        nHp = nfast_icrh;
    case 'T'
        nTp = nfast_icrh;
    case 'He3'
        nHep3 = nfast_icrh;
    case 'He4'
        nHep4 = nfast_icrh;
    case 'B'
        nBp = nfast_icrh;
    otherwise
        error(sprintf('minority species %s not yet implemanted',option.mino));
end

% correction of fast NBI ions
switch gas_nbi
    case -1
        nHp = nHp + nfast_nbi;
        nHp = nHp + nfast_nbi2;
    case 3
        nTp = nTp + real(cons.ftnbi * ve)  .* nfast_nbi1;
        nDp = nDp + (1 - real(cons.ftnbi * ve)) .* nfast_nbi1;
        nTp = nTp + imag(cons.ftnbi * ve)  .* nfast_nbi2;
        nDp = nDp + (1 - imag(cons.ftnbi * ve)) .* nfast_nbi2;
    case 5
        nDp = nDp + (1 - real(cons.ftnbi * ve)) .* nfast_nbi1;
        nDp = nDp + (1 - imag(cons.ftnbi * ve)) .* nfast_nbi2;
        nHep3 = nHep3 + real(cons.ftnbi * ve)  .* nfast_nbi1;
        nHep3 = nHep3 + imag(cons.ftnbi * ve)  .* nfast_nbi2;
    case 11
        nHp = nHp + (1 - real(cons.ftnbi * ve)) .* nfast_nbi;
        nHp = nHp + (1 - imag(cons.ftnbi * ve)) .* nfast_nbi2;
        nBp = nBp + real(cons.ftnbi * ve)  .* nfast_nbi1;
        nBp = nBp + imag(cons.ftnbi * ve)  .* nfast_nbi2;
        
    otherwise
        nHp = nHp + real(cons.ftnbi * ve)  .* nfast_nbi1;
        nDp = nDp + (1 - real(cons.ftnbi * ve)) .* nfast_nbi1;
        nHp = nHp + real(cons.ftnbi * ve)  .* nfast_nbi2;
        nDp = nDp + (1 - real(cons.ftnbi * ve)) .* nfast_nbi2;
end

% contribution from fast alpha
nHep4 = nHep4 + nfast_alpha;


% fill template 
nfast(:,:,1) = nHp;
nfast(:,:,2) = nDp;
nfast(:,:,3) = nTp;
nfast(:,:,4) = nHep3;
nfast(:,:,5) = nHep4;
nfast(:,:,6) = nBp;
