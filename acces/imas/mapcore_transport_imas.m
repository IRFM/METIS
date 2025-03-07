% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coretransp = mapcore_transport_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coretransp,sigma_B0_eff)
 
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

% improve precision
[A_el,Z_el,name_el] = chargemasse;
dd   = abs(Z_el - z0dstruct.z0dinput.option.zimp);
mask = (dd == min(dd));
aimp = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(aimp)
    aimp = 7/3 .* z0dstruct.z0dinput.option.zimp;
end
dd   = abs(Z_el - z0dstruct.z0dinput.option.zmax);
mask = (dd == min(dd));
amax = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(amax)
    amax = 7/3 .* z0dstruct.z0dinput.option.zmax;
end
% % impurete principale
zimp = z0dstruct.z0dinput.option.zimp;
% aimp = ceil(zimp .* (7/3));
% 
% % 2ieme impurete
zmax = z0dstruct.z0dinput.option.zmax;
% amax = ceil(zmax .* (7/3));



%% IMAS part start here
% precomputation
ve = ones(size(profil0d.xli));
nDp   = max(1e13,profil0d.n1p .* ((interp1_imas(data_zerod.temps,data_zerod.nDm,profil0d.temps,'pchip','extrap')./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
nTp   = max(1e13,profil0d.n1p .* ((interp1_imas(data_zerod.temps,data_zerod.nDm,profil0d.temps,'pchip','extrap')./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
switch z0dstruct.z0dinput.option.gaz
    case 11
        nbp   = nTp;
        nTp   = zeros(size(nTp));
        nHp   = max(1e13,profil0d.n1p - nDp);
    otherwise
        nbp   = zeros(size(nTp));
        nHp   = max(1e13,profil0d.n1p - nTp - nDp);
end
switch z0dstruct.z0dinput.option.gaz
    case 5
        nhep3  = max(1e13,profil0d.nhep);
        nhep   = max(0,z0dstruct.z0dinput.option.frhe0 .* profil0d.nep);
    otherwise
        nhep  = max(1e13,profil0d.nhep);
        % at this stage for this case
        nhep3 = zeros(size(nhep));
end
%
switch z0dstruct.z0dinput.option.mino
    case 'He3'
        switch z0dstruct.z0dinput.option.gaz
            case 4
                nHe3m = z0dstruct.z0dinput.option.cmin .* data_zerod.nhem;
                nHem  = max(0,data_zerod.nhem - nHe3m);
            case 5
                nHe3m = data_zerod.nhem;
                nHem  = z0dstruct.z0dinput.option.frhe0 .* data_zerod.nem;
                
            otherwise
                nHe3m = z0dstruct.z0dinput.option.cmin .* data_zerod.n1m;
                nHem  = max(0,data_zerod.nhem - nHe3m);
        end
    otherwise
        nHem  = data_zerod.nhem;
        nHe3m = 0 .* nHem;
end
frhe3  = nHe3m ./ max(1e11,nHe3m + nHem);
frhe3  = interp1_imas(data_zerod.temps,frhe3,profil0d.temps,'pchip','extrap') * ve;
switch z0dstruct.z0dinput.option.gaz
   case 5
       % nothing
   otherwise
       nhep3        = nhep .* (1 - frhe3);       
       nhep         = nhep .* frhe3;
end
% nHp   = max(1e13,profil0d.n1p - nTp - nDp);
% nhep  = max(1e13,profil0d.nhep);
nz1p  = max(1e13,profil0d.nzp);
nz2p  = max(1e11,profil0d.nzp .* z0dstruct.z0dinput.option.rimp);   
nwp  = max(1,profil0d.nwp);
% une seule vitesse d'ensemble
% Mtor    = phys.mp .*  max(1e13,nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + ...
%           ceil(7/3 .* z0dstruct.z0dinput.option.zimp) .* nz1p + ceil(7/3 .* z0dstruct.z0dinput.option.zmax) .* nz2p + 183.84 .* nwp);
if z0dstruct.z0dinput.option.Sn_fraction > 0
    Mtor    = phys.mp .*  max(1e13,nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + aimp .* nz1p + amax .* nz2p +  ...
              (1 - option.Sn_fraction) .*183.84 .* nwp) + option.Sn_fraction .*118.71 .* nwp + 3.02 .* nhep3 + 11 .* nbp;
else
    Mtor    = phys.mp .*  max(1e13,nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + aimp .* nz1p + amax .* nz2p + 183.84 .* nwp + 3.02 .* nhep3 + 11 .* nbp);
end

% generic fields
coretransp.ids_properties.comment = 'METIS transport coefficients and flux';
coretransp.time		    = profil0d.temps;
rb0 =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
coretransp.vacuum_toroidal_field.r0          = mean(z0dstruct.z0dinput.geo.R);
coretransp.vacuum_toroidal_field.b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ coretransp.vacuum_toroidal_field.r0;

% model METIS
% initialisation substructures
model = coretransp.model{1};
model.identifier.name        = 'METIS';
model.identifier.index       = 1;
model.identifier.description = 'transport coefficients and flux identified from METIS kinetic profiles and METIS sources';
model.flux_multiplier        = 0;


% loop on time slices
% null profile
zz              = zeros(size(profil0d.xli));
%% SOURCE terms
for k=1:length(profil0d.temps)
    % initialisation substructures
    profiles_1d = model.profiles_1d{1};
    profiles_1d.grid_d                         = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid_d);
    profiles_1d.grid_v                         = profiles_1d.grid_d;
    profiles_1d.grid_flux                      = profiles_1d.grid_d;
    profiles_1d.conductivity_parallel          = 1./max(1e-307,profil0d.eta(k,:));
    profiles_1d.electrons.particles.d          = profil0d.dn(k,:);
    profiles_1d.electrons.particles.v          = - profil0d.vn(k,:);
    profiles_1d.electrons.particles.flux       = profil0d.ge(k,:);
    profiles_1d.electrons.energy.d             = profil0d.xie(k,:);
    profiles_1d.electrons.energy.v             = zz;
    profiles_1d.electrons.energy.flux          = profil0d.qe(k,:);
    profiles_1d.total_ion_energy.d             = profil0d.xii(k,:);
    profiles_1d.total_ion_energy.v             = zz;
    profiles_1d.total_ion_energy.flux          = profil0d.qi(k,:);
    profiles_1d.momentum_tor.d                 = profil0d.drot(k,:);
    profiles_1d.momentum_tor.v                 = - profil0d.vrot(k,:);
    profiles_1d.momentum_tor.flux              = profil0d.frot(k,:) ./ Mtor(k,:) .* profil0d.r2i(k,:) .* profil0d.Raxe(k,:);
    profiles_1d.time                           = profil0d.temps(k);
    model.profiles_1d{k} = profiles_1d;
end
% final assignement
coretransp.model{1} = model;

