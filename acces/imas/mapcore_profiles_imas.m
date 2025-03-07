%% FILL CORE_PROFILE IDS WITH ONLY INTERESTING VARIABLES
function coreprof = mapcore_profiles_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coreprof,sigma_B0_eff)

%% IMAS part start here
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


% precomputation
amin  =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.a,profil0d.temps,'pchip','extrap');

% Sn
if ~isfield(z0dstruct.z0dinput.option,'Sn_fraction')
    z0dstruct.z0dinput.option.Sn_fraction = 0;
end

% isotopic composition for option.gaz == 5
if z0dstruct.z0dinput.option.gaz == 5
    nHe3onD = real(z0dstruct.z0dinput.cons.iso);
    nTonD   = imag(z0dstruct.z0dinput.cons.iso);
    warning('nHe3onD & nTonD not yet used !');
else
    nHe3onD = zeros(size(z0dstruct.z0dinput.cons.iso));
    nTonD   = real(z0dstruct.z0dinput.cons.iso);
end
z0dstruct.z0dinput.cons.iso = real(z0dstruct.z0dinput.cons.iso);

%
% densite ioniques
%
ve    = ones(size(profil0d.xli));
nDm   = interp1_imas(data_zerod.temps,data_zerod.nDm,profil0d.temps,'pchip','extrap');
nTm   = interp1_imas(data_zerod.temps,data_zerod.nTm,profil0d.temps,'pchip','extrap');
nDp   = max(1,profil0d.n1p .* ((nDm./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
nTp   = max(1,profil0d.n1p .* ((nTm./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
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
%nHp   = max(1,profil0d.n1p - nTp - nDp);
%nhep  = max(1,profil0d.nhep);
nz1p  = max(1,profil0d.nzp);
nz2p  = max(1,profil0d.nzp .* z0dstruct.z0dinput.option.rimp);   
nwp   = max(1,profil0d.nwp);
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
nions  = NaN * ones(length(profil0d.temps),length(profil0d.xli),10);
nions(:,:,1) = nHp;
nions(:,:,2) = nDp;
nions(:,:,3) = nTp;
switch z0dstruct.z0dinput.option.gaz
   case 5
        nions(:,:,4) = nhep3;
        nions(:,:,5) = nhep;
   otherwise
       nions(:,:,4) = nhep .* frhe3;
       nions(:,:,5) = nhep .* (1 - frhe3);
       %nhep3        = nhep .* (1 - frhe3);       
       %nhep         = nhep .* frhe3;
end
nions(:,:,6) = nz1p;
nions(:,:,7) = nz2p;
nions(:,:,8) = nwp;
if z0dstruct.z0dinput.option.Sn_fraction > 0
   nions(:,:,8) = (1-z0dstruct.z0dinput.option.Sn_fraction) .* nwp;
   nions(:,:,9) = z0dstruct.z0dinput.option.Sn_fraction .* nwp; 
else
   nions(:,:,9) = zeros(size(nwp));
end
nions(:,:,10) = nbp;

% ion masse and charge table
% H
A(1) = 1; 
Z(1) = 1;
label{1} = 'H';
% D
A(2) = 2; 
Z(2) = 1;
label{2} = 'D';
% T
A(3) = 3; 
Z(3) = 1;
label{3} = 'T';
% He3
A(4) = 3; 
Z(4) = 2;
label{4} = 'He3';
% He4
A(5) = 4; 
Z(5) = 2;
label{5} = 'He4';
% imp 1
Z(6) = z0dstruct.z0dinput.option.zimp;
[A(6),label{6}] = getafromz_imas(Z(6));
% imp 2
Z(7) = z0dstruct.z0dinput.option.zmax;
[A(7),label{7}] = getafromz_imas(Z(7));
% W
Z(8) = 74;
[A(8),label{8}] = getafromz_imas(Z(8));
% Sn
Z(9) = 50;
[A(9),label{9}] = getafromz_imas(Z(9));
% boron
Z(10) = 5;
A(11) = 11;
label{10} = 'B';

% enforced Zeff et ne -> improved precision
% due to resampling, small error on electroneutrality is introduce and tiny difference in Zeff
nizi2 = zeros(size(profil0d.nep));
nizi  = zeros(size(profil0d.nep));
for k= 1:length(Z)
    if Z(k) == 74
        nizi	= nizi  + z0wavez(profil0d.tep)      .* squeeze(nions(:,:,k));
        nizi2	= nizi2 + z0wavez(profil0d.tep) .^ 2 .* squeeze(nions(:,:,k));
    elseif Z(k) == 50
        nizi	= nizi  + z0snavez(profil0d.tep)      .* squeeze(nions(:,:,k));
        nizi2	= nizi2 + z0snavez(profil0d.tep) .^ 2 .* squeeze(nions(:,:,k));
    else
        nizi	= nizi  + Z(k)     .* squeeze(nions(:,:,k));
        nizi2	= nizi2 + Z(k) .^2 .* squeeze(nions(:,:,k));
    end
end
fprintf('Relative difference on electoneutrality = %g\n',sqrt(mean((profil0d.nep(:) - nizi(:)) .^ 2 ./ profil0d.nep(:) .^ 2 ./ size(profil0d.nep,1))));
fprintf('Relative difference on Zeff = %g\n',sqrt(mean((profil0d.zeff(:) - nizi2(:) ./ nizi(:)) .^ 2./ size(profil0d.nep,1))));
profil0d.nep   = nizi;
profil0d.zeff  = nizi2 ./ nizi;
%  figure;
%  subplot(2,1,1)
%  zplotprof(gca,profil0d.temps,profil0d.xli,profil0d.nep,'color','r');
%  zplotprof(gca,profil0d.temps,profil0d.xli,nizi,'color','b','marker','o','linestyle',':');
%  subplot(2,1,2)
%  zplotprof(gca,profil0d.temps,profil0d.xli,profil0d.zeff,'color','r');
%  zplotprof(gca,profil0d.temps,profil0d.xli,nizi2 ./ nizi,'color','b','marker','o','linestyle',':');
%  keyboard


% rotation 
[rtor,vtor,vpol,omega,mtor] = z0rot_imas(data_zerod,profil0d,z0dstruct.z0dinput.option,frhe3,z0dstruct.z0dinput.geo,z0dstruct.z0dinput.cons);

% fast particles densities (do not take into account addition of Sn !
[psupra,psupra_para,psupra_perp,nfast] = nfast4imas(data_zerod,profil0d,z0dstruct.z0dinput.option,frhe3,z0dstruct.z0dinput.geo,z0dstruct.z0dinput.cons);

% neutrals
% cold first and hot second
switch z0dstruct.z0dinput.option.gaz
    case 1
        amn = [1,1];
        zn  = [1,1];
        ion_index = [1,1];
        neutlab   = {'H2','H2'};
    case 2
        amn = [2,2];
        zn  = [1,1];
        ion_index = [2,2];
        neutlab   = {'D2','D2'};
    case 3
        amn = [2,3,2,3];
        zn  = [1,1,1,1];
        ion_index = [2,3,2,3];
        neutlab   = {'D2','T2','D2','T2'};
    case 4
        amn = [2,2];
        zn  = [4,4];
        ion_index = [5,5];
        neutlab   = {'He4','He4'};
    case 5 
        amn = [2,3,2,3];
        zn  = [1,2,1,2];
        ion_index = [2,4,2,4];
        neutlab   = {'D2','He3','D2','He3'};
        
     case 11 
        amn = [1,11,1,11];
        zn  = [1,5,1,5];
        ion_index = [1,10,1,10];
        neutlab   = {'H2','B','H2','B'};
        
   otherwise
        error('plasma compostion not yet implemented in METIS4IMAS');
end
%
n_neutrals  = NaN * ones(length(profil0d.temps),length(profil0d.xli),length(amn));
t_neutrals  = NaN * ones(length(profil0d.temps),length(profil0d.xli),length(amn));
% utilisation de la vitesse du son (plus stable numeriquement, generalement peu differente)
t_hot  = (profil0d.tip + profil0d.tep .* profil0d.zeff)./ 2;
t_cold = interp1_imas(data_zerod.temps,data_zerod.telim,profil0d.temps,'pchip','extrap') * ones(size(profil0d.xli));

%
switch z0dstruct.z0dinput.option.gaz
    case {3,5}
        iso = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.cons.iso,profil0d.temps,'pchip','extrap');
        iso = iso * ones(size(profil0d.xli));
        n_neutrals(:,:,1)   = profil0d.n0m ./ (1+iso);
        t_neutrals(:,:,1)   = t_cold;
        n_neutrals(:,:,3)   = profil0d.n0  ./ (1+iso);
        t_neutrals(:,:,3)   = t_hot;
        n_neutrals(:,:,2)   = profil0d.n0m ./ (1+iso) .* iso;
        t_neutrals(:,:,2)   = t_cold;
        n_neutrals(:,:,4)   = profil0d.n0  ./ (1+iso) .* iso;
        t_neutrals(:,:,4)   = t_hot;
        % take into account molecule (for H and D) and number of charge for He4.
        n_neutrals = n_neutrals ./ 2;
     case 11
        iso = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.cons.iso,profil0d.temps,'pchip','extrap');
        iso = iso * ones(size(profil0d.xli));
        n_neutrals(:,:,1)   = profil0d.n0m ./ (1+iso) / 2;
        t_neutrals(:,:,1)   = t_cold;
        n_neutrals(:,:,3)   = profil0d.n0  ./ (1+iso) / 2;
        t_neutrals(:,:,3)   = t_hot;
        n_neutrals(:,:,2)   = profil0d.n0m ./ (1+iso) .* iso / 5;
        t_neutrals(:,:,2)   = t_cold;
        n_neutrals(:,:,4)   = profil0d.n0  ./ (1+iso) .* iso / 5;
        t_neutrals(:,:,4)   = t_hot;
    otherwise
        n_neutrals(:,:,1)   = profil0d.n0m;
        t_neutrals(:,:,1)   = t_cold;
        n_neutrals(:,:,2)   = profil0d.n0;
        t_neutrals(:,:,2)   = t_hot;
        % take into account molecule (for H and D) and number of charge for He4.
        n_neutrals = n_neutrals ./ 2;
end

% generic fields
coreprof.ids_properties.comment = 'METIS profiles 1D';
coreprof.time		    = profil0d.temps;
rb0 =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
coreprof.vacuum_toroidal_field.r0          = mean(z0dstruct.z0dinput.geo.R);
coreprof.vacuum_toroidal_field.b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ coreprof.vacuum_toroidal_field.r0;

%% global_quantities
coreprof.global_quantities.ip                     = interp1_imas(data_zerod.temps,data_zerod.ip,profil0d.temps,'pchip','extrap');
coreprof.global_quantities.current_non_inductive  = interp1_imas(data_zerod.temps,data_zerod.ini,profil0d.temps,'pchip','extrap');
coreprof.global_quantities.current_bootstrap      = interp1_imas(data_zerod.temps,data_zerod.iboot + data_zerod.ifus,profil0d.temps,'pchip','extrap');
coreprof.global_quantities.v_loop                 = interp1_imas(data_zerod.temps,data_zerod.vloop,profil0d.temps,'pchip','extrap');
coreprof.global_quantities.li                     = interp1_imas(data_zerod.temps,data_zerod.li,profil0d.temps,'pchip','extrap');
coreprof.global_quantities.li_3                   = interp1_imas(data_zerod.temps,data_zerod.li,profil0d.temps,'pchip','extrap');
coreprof.global_quantities.beta_tor               = interp1_imas(data_zerod.temps,data_zerod.betan,profil0d.temps,'pchip','extrap') .*  ...
                                                                 coreprof.global_quantities.ip ./ coreprof.vacuum_toroidal_field.b0 ./ amin ./ 1e6;
coreprof.global_quantities.beta_tor_norm          = interp1_imas(data_zerod.temps,100 .* data_zerod.betan,profil0d.temps,'pchip','extrap');
coreprof.global_quantities.beta_pol               = interp1_imas(data_zerod.temps,data_zerod.betaptot,profil0d.temps,'pchip','extrap');
coreprof.global_quantities.energy_diamagnetic     = interp1_imas(data_zerod.temps,data_zerod.wdia,profil0d.temps,'pchip','extrap');


% loop on time slices
% null profile
zz              = zeros(size(profil0d.xli));
%model_ion   = coreprof.profiles_1d{1}.ion{1};
model_ion   = ids_allocate('core_profiles','profiles_1d/ion',1);
model_ion   = model_ion{1};
%model_state = model_ion.state;
model_state = ids_allocate('core_profiles','profiles_1d/ion/state',1);
%model_state = model_state{1};
%model_element = model_ion.element{1};
model_element = ids_allocate('core_profiles','profiles_1d/ion/element',1);
model_element = model_element{1};

%% SOURCE terms
for k=1:length(profil0d.temps)
    profiles_1d = coreprof.profiles_1d{1};
    if isfield(coreprof,'profiles_2d')
        profiles_2d = coreprof.profiles_2d{1};
    else
        try
            profiles_2d = ids_allocate('core_profiles','profiles_2d',1);
        catch
            profiles_2d = [];
        end
    end
    if isfield(coreprof,'statistics')
        statistics  = coreprof.statistics{1};
    else
        try
            statistics = ids_allocate('core_profiles','statistics',1);
        catch
            statistics = [];
        end
    end
    profiles_1d.grid                                       = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
    %
    profiles_1d.electrons.temperature                      = profil0d.tep(k,:);
    profiles_1d.electrons.density_thermal                  = profil0d.nep(k,:);
    profiles_1d.electrons.density_fast                     = nfast(k,:,9);
    profiles_1d.electrons.density                          = profiles_1d.electrons.density_thermal + profiles_1d.electrons.density_fast;
    profiles_1d.electrons.pressure                         = phys.e .* profil0d.tep(k,:) .* profil0d.nep(k,:) + psupra(k,:,9);
    profiles_1d.electrons.pressure_fast_perpendicular      = psupra_perp(k,:,9);
    profiles_1d.electrons.pressure_fast_parallel           = psupra_para(k,:,9);
    %profiles_1d.electrons.velocity_tor                     = not computed inside METIS
    %profiles_1d.electrons.velocity_pol                     = not computed inside METIS
    %profiles_1d.electrons.velocity.radial                  = not computed inside METIS
    %profiles_1d.electrons.velocity.diamagnetic             = not computed inside METIS
    %profiles_1d.electrons.velocity.parallel                = not computed inside METIS
    %profiles_1d.electrons.velocity.poloidal                = not computed inside METIS
    %profiles_1d.electrons.velocity.toroidal                = not computed inside METIS
    %
    profiles_1d.t_i_average            = profil0d.tip(k,:);
    profiles_1d.n_i_total_over_n_e     = profil0d.nip(k,:) ./ max(1,profil0d.nep(k,:));
    profiles_1d.momentum_tor           = rtor(k,:);
    profiles_1d.zeff                   = profil0d.zeff(k,:);
    profiles_1d.pressure_ion_total     = phys.e .* sum(nions(k,:,[1:8,10]),3) .* profil0d.tip(k,:) + sum(psupra(k,:,[1:8,10]),3);
    profiles_1d.pressure_thermal       = phys.e .* profil0d.tep(k,:) .* profil0d.nep(k,:) + phys.e .* profil0d.nip(k,:) .* profil0d.tip(k,:);
    profiles_1d.pressure_perpendicular = sum(psupra_perp(k,:,:),3) + profiles_1d.pressure_thermal;
    profiles_1d.pressure_parallel      = sum(psupra_para(k,:,:),3) + profiles_1d.pressure_thermal;
    profiles_1d.j_total                = profil0d.jeff(k,:);
    profiles_1d.j_tor                  = profil0d.jli(k,:);
    profiles_1d.j_ohmic                = profil0d.jeff(k,:) - profil0d.jni(k,:);
    profiles_1d.j_non_inductive        = profil0d.jni(k,:);
    profiles_1d.current_parallel_inside = cumtrapz(profil0d.xli,profil0d.jni(k,:) .* profil0d.spr(k,:),2);
    profiles_1d.j_bootstrap            = profil0d.jboot(k,:);
    profiles_1d.conductivity_parallel  = 1./max(1e-307,profil0d.eta(k,:));
    profiles_1d.e_field_parallel       = profil0d.epar(k,:);
    %
    profiles_1d.e_field.radial         = profil0d.er(k,:);
    %profiles_1d.e_field.diamagnetic    = not computed inside METIS
    profiles_1d.e_field.parallel       = profil0d.epar(k,:);
    %profiles_1d.e_field.poloidal       = not computed inside METIS
    %profiles_1d.e_field.toroidal       = not computed inside METIS
    %
    profiles_1d.q                      = profil0d.qjli(k,:);
    profiles_1d.magnetic_shear         = pdederive(profil0d.rmx(k,:),profil0d.qjli(k,:),0,2,2,1) ./ profil0d.qjli(k,:) .* profil0d.rmx(k,:);
    profiles_1d.time                   = profil0d.temps(k);
    
%     % check q phi and psi
%     qcheck = abs(pdederive(profiles_1d.grid.psi,profil0d.phi(k,:),1,2,2,1)) ./ 2 ./ pi;
%     phi_check = profiles_1d.grid.rho_tor .^ 2 .* pi .* coreprof.vacuum_toroidal_field.b0(k);
%     qcheck2 = abs(pdederive(profiles_1d.grid.psi,phi_check,1,2,2,1)) ./ 2 ./ pi;
%     mean(sqrt((qcheck - profiles_1d.q) .^ 2)) 
%     mean(sqrt((qcheck2 - profiles_1d.q) .^ 2)) 
%     mean(sqrt((qcheck2 -qcheck) .^ 2)) 
% 
    
    % filling ion structure for each species
    % initialise substructure
    ion = model_ion;
    element = model_element;
    for l = 1:10
        % filling ion  structure
        element.a            = A(l);
        element.z_n          = Z(l);
        ion.element{1}       = element;
        ion.z_ion            = Z(l);
        ion.label            = label{l};
        ion.name             = label{l};
        %ion.neutral_index    = 2 kinds of neutral can gives one ion, must be a vector
        ion.temperature      = profil0d.tip(k,:); %dépend de la définition de la température <<<<<<<<<<<<<<<<<<<<<<<
        ion.density_thermal  = nions(k,:,l);
        ion.density_fast     = nfast(k,:,l);
        ion.density          = ion.density_thermal + ion.density_fast;
        ion.pressure         = phys.e .* nions(k,:,l) .* profil0d.tip(k,:) + psupra(k,:,l); % problème de definition !!!
        ion.pressure_fast_perpendicular = psupra_perp(k,:,l);
        ion.pressure_fast_parallel      = psupra_para(k,:,l);
        ion.velocity_tor                = vtor(k,:,l);
        ion.velocity_pol                = vpol(k,:,l);
        %ion.velocity.radial             = not computed inside METIS
        %ion.velocity.diamagnetic        = not computed inside METIS
        %ion.velocity.parallel           = not computed inside METIS
        ion.velocity.poloidal           = vpol(k,:,l);
        ion.velocity.toroidal           = vtor(k,:,l);
        ion.state                       = mapionstate(Z(l),A(l),label{l},profil0d.tep(k,:),profil0d.tip(k,:),squeeze(nions(k,:,l)),squeeze(nfast(k,:,l)), ...
            squeeze(psupra(k,:,l)),squeeze(psupra_perp(k,:,l)),squeeze(psupra_para(k,:,l)),squeeze(vtor(k,:,l)),squeeze(vpol(k,:,l)),model_state);
        ion.temperature_fit.source      = ' ';
        ion.density_fit.source          = ' ';
        ion.state{1}.density_fit.source = ' ';
        profiles_1d.ion{l}              = ion;
        
        if l == 9
            % correction infilling ion  structure for Sn separately
            ion.density_fast     = zeros(size(nfast(k,:,end))); % no fast ion for Sn
            ion.density          = ion.density_thermal + ion.density_fast;
            ion.pressure         = phys.e .* nions(k,:,end) .* profil0d.tip(k,:); %  no supra for Sn as for W + psupra(k,:,l); % problème de definition !!!
            ion.pressure_fast_perpendicular = zeros(size(psupra_perp(k,:,end))); % no fast ion for Sn
            ion.pressure_fast_parallel      = zeros(size(psupra_para(k,:,end)));% no fast ion for Sn
            ion.state                       = mapionstate(Z(end),A(end),label{end},profil0d.tep(k,:),profil0d.tip(k,:),squeeze(nions(k,:,end)),zeros(size(squeeze(nfast(k,:,end)))), ...
                                              zeros(size(squeeze(psupra(k,:,end)))),zeros(size(squeeze(psupra_perp(k,:,end)))),zeros(size(squeeze(psupra_para(k,:,end)))),squeeze(vtor(k,:,end)),squeeze(vpol(k,:,end)),model_state);
        end
        
    end
    %
    profiles_1d.ion{end+1}          = ion;
    
    % filling neural  structure for each species
    neutral = profiles_1d.neutral{1};
    element = neutral.element{1};
    for l = 1:length(amn)
        % filling ion  structure
        element.a            = amn(l);
        element.z_n          = zn(l);
        neutral.element{1}       = element;
        neutral.label            = neutlab{l};
        neutral.name             = neutlab{l};
        neutral.ion_index        = ion_index(l);
        neutral.temperature      = t_neutrals(k,:,l);
        neutral.density_thermal  = n_neutrals(k,:,l);
        neutral.density          = n_neutrals(k,:,l);
        neutral.density_fast     = zz;
        neutral.pressure         = phys.e .* t_neutrals(k,:,l) .* n_neutrals(k,:,l);
        neutral.pressure_fast_perpendicular = zz;
        neutral.pressure_fast_parallel      = zz;
        %neutral.velocity.radial             = not computed inside METIS
        %neutral.velocity.diamagnetic        = not computed inside METIS
        %neutral.velocity.parallel           = not computed inside METIS
        %neutral.velocity.poloidal           = not computed inside METIS
        %neutral.velocity.toroidal           = not computed inside METIS
        profiles_1d.neutral{l}              = neutral;
    end
    
    coreprof.profiles_1d{k} = profiles_1d;
    coreprof.profiles_2d{k} = profiles_2d;
    coreprof.statistics{k} = statistics;
end

