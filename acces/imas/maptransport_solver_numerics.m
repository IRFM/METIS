function transport_solver_numerics = maptransport_solver_numerics(z0dstruct,data_zerod,profil0d,texte_diary ...
    ,error_flag,transport_solver_numerics,run,occurrence,sigma_B0_eff,sigma_bvac_r)


%% FOR METIS TIME-INDEPENDENT VARIABLES WHICH BECOME TIME-DEPENDENT IN THE IDS
ntime = length(profil0d.temps);
vt = ones(size(profil0d.temps));

% compatibility between various call mode
if ~isfield(z0dstruct,'profil0d')
    z0dstruct.profil0d = z0dstruct.profil;
end
% Sn
if ~isfield(z0dstruct.z0dinput.option,'Sn_fraction')
    z0dstruct.z0dinput.option.Sn_fraction = 0;
end


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
%amin                   = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.a,profil0d.temps,'pchip','extrap');
%magnetic_shear         = pdederive(profil0d.xli,profil0d.qjli,0,2,2,1) ./ profil0d.qjli .* (ones(size(profil0d.qjli,1),1) * profil0d.xli);
asser                  = interp1_imas(data_zerod.temps,data_zerod.asser,profil0d.temps,'nearest','extrap');
ip                  = interp1_imas(data_zerod.temps,data_zerod.ip,profil0d.temps,'pchip','extrap');
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
nhep  = max(1,profil0d.nhep);
nz1p  = max(1,profil0d.nzp);
nz2p  = max(1,profil0d.nzp .* z0dstruct.z0dinput.option.rimp);
if  z0dstruct.z0dinput.option.Sn_fraction > 0
    nwp   = (1 - z0dstruct.z0dinput.option.Sn_fraction) .* max(1,profil0d.nwp);
    nsnp  = z0dstruct.z0dinput.option.Sn_fraction .* max(1,profil0d.nwp);    
else
    nwp   = max(1,profil0d.nwp);
    nsnp  = zeros(size(nwp));
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
nions(:,:,9) = nsnp;
nions(:,:,10) = nsnp;

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
%Sn
Z(9) = 50;
[A(9),label{9}] = getafromz_imas(Z(9));
% boron
Z(10) = 5;
A(11) = 11;
label{10} = 'B';


% rotation
[rtor,vtor,vpol,omega,mtor] = z0rot_imas(data_zerod,profil0d,z0dstruct.z0dinput.option,frhe3,z0dstruct.z0dinput.geo,z0dstruct.z0dinput.cons);
% time derivative
if length(profil0d.temps) > 1
   drtordt = z0dxdt(rtor,profil0d.temps);
else
   drtordt = z0dxdt(z0rot_imas(z0dstruct.zerod,z0dstruct.profil0d,z0dstruct.z0dinput.option,ones(size(z0dstruct.profil0d.temps)) *frhe3,z0dstruct.z0dinput.geo,z0dstruct.z0dinput.cons),z0dstruct.profil0d.temps);
   drtordt = interp1_imas(z0dstruct.profil0d.temps,drtordt,profil0d.temps,'pchip','extrap');
end

% fast particles densities
[psupra,psupra_para,psupra_perp,nfast] = nfast4imas(data_zerod,profil0d,z0dstruct.z0dinput.option,frhe3,z0dstruct.z0dinput.geo,z0dstruct.z0dinput.cons);

% computation of time derivative
profil0d_dt = [];
if length(profil0d.temps) > 1
	tloc = profil0d.temps;
        noms = fieldnames(profil0d);
	for k = 1:length(noms)
		var = profil0d.(noms{k});
		if size(var,1) == length(tloc) && all(isfinite(var(:)))
 			profil0d_dt.(noms{k}) = z0dxdt(var,tloc);
		else
			profil0d_dt.(noms{k}) = var;
		end
	end
	% in time vector store time step
	profil0d_dt.temps = interp1_imas(tloc(1:end-1),diff(tloc),tloc,'pchip','extrap');

else
	tloc = z0dstruct.profil0d.temps;
	tnow = profil0d.temps;
        noms = fieldnames(z0dstruct.profil0d);
	for k = 1:length(noms)
		var = z0dstruct.profil0d.(noms{k});
		if (size(var,1) == length(tloc)) && all(isfinite(var(:)))
 			profil0d_dt.(noms{k}) = interp1_imas(tloc,z0dxdt(var,tloc) ,tnow,'pchip','extrap');
		else
			profil0d_dt.(noms{k}) = var;
		end
	end	
	% in time vector store time step
	profil0d_dt.temps = interp1_imas(tloc(1:end-1),diff(tloc),tnow,'pchip','extrap');
end

%% DEFINE EMPTY CODE AND IDS_PROPERTIES STRUCTURES
if ~isfield(transport_solver_numerics,'code')
    transport_solver_numerics.code = code_empty_imas;
end
if ~isfield(transport_solver_numerics,'ids_properties')
    transport_solver_numerics.ids_properties = ids_properties_empty_imas;
end

%% COMMENT
if ~defined_imas(transport_solver_numerics.ids_properties.comment)
    transport_solver_numerics.ids_properties.comment = 'IMAS implementation of METIS';
end

%% HOMOGENEOUS_TIME
if ~defined_imas(transport_solver_numerics.ids_properties.homogeneous_time)
    transport_solver_numerics.ids_properties.homogeneous_time = 1;
end

% time
transport_solver_numerics.time = profil0d.temps;
% simple
rb0 =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
transport_solver_numerics.vacuum_toroidal_field.r0 = mean(z0dstruct.z0dinput.geo.R);
transport_solver_numerics.vacuum_toroidal_field.b0 = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ transport_solver_numerics.vacuum_toroidal_field.r0;

% factor on psi
factor_psi = (profil0d.psi(end,1) - profil0d.psi(end,end)) ./ ...
    (z0dstruct.profil0d.psi(end,1) - z0dstruct.profil0d.psi(end,end));


% ion data structure (derivatives_1d)
for l=1:length(label)
    % substructure initialisation
    ion{l} = transport_solver_numerics.derivatives_1d{1}.ion{1};
    ion{l}.a 					= A(l);
    ion{l}.z_ion				= Z(l);
    ion{l}.z_n                  = Z(l);
    ion{l}.label				= label{l};
    ion{l}.name					= label{l};
    ion{l}.multiple_states_flag               = 0;
    ion{l}.particles.identifier.name 	        = label{l};
    ion{l}.particles.identifier.index         = l;
    ion{l}.particles.identifier.description	= label{l};
    ion{l}.energy.identifier.name	        = label{l};
    ion{l}.energy.identifier.index            = l;
    ion{l}.particles.iterations_n             = data_zerod.nb;
    ion{l}.energy.iterations_n                = data_zerod.nb;
end
% ion data structure (boundary_conditions_1d)
for l=1:length(label)
    % substructure initialisation
    ion2{l} = transport_solver_numerics.boundary_conditions_1d{1}.ion{1};
    ion2{l}.a 					= A(l);
    ion2{l}.z_ion				= Z(l);
    ion2{l}.z_n                  = Z(l);
    ion2{l}.label				= label{l};
    ion2{l}.name				= label{l};
    ion2{l}.multiple_states_flag               = 0;
    ion2{l}.particles.identifier.name 	        = label{l};
    ion2{l}.particles.identifier.index         = l;
    ion2{l}.particles.identifier.description	= label{l};
    ion2{l}.energy.identifier.name	        = label{l};
    ion2{l}.energy.identifier.index            = l;
    ion2{l}.particles.iterations_n             = data_zerod.nb;
    ion2{l}.energy.iterations_n                = data_zerod.nb;
end
% ion data structure (convergence.equations{1})
for l=1:length(label)
    % substructure initialisation
    ion3{l} = transport_solver_numerics.convergence.equations{1}.ion{1};
    ion3{l}.a 					= A(l);
    ion3{l}.z_ion				= Z(l);
    ion3{l}.z_n                  = Z(l);
    ion3{l}.label				= label{l};
    ion3{l}.name				= label{l};
    ion3{l}.multiple_states_flag               = 0;
    ion3{l}.particles.identifier.name 	        = label{l};
    ion3{l}.particles.identifier.index         = l;
    ion3{l}.particles.identifier.description	= label{l};
    ion3{l}.energy.identifier.name	        = label{l};
    ion3{l}.energy.identifier.index            = l;
    ion3{l}.particles.iterations_n             = data_zerod.nb;
    ion3{l}.energy.iterations_n                = data_zerod.nb;
end
% ni / ne
nione = profil0d.nip ./ max(1,profil0d.nep);

% loop on time
for k=1:length(transport_solver_numerics.time)
    if length(profil0d) == 1
        index = length(z0dstruct.profil0d.temps);
    else
        index = k;
    end
    % initialisation unused structure 
    transport_solver_numerics.boundary_conditions_ggd{k} = transport_solver_numerics.boundary_conditions_ggd{1};
    % initialisation of substructures of convergence
    transport_solver_numerics.convergence.equations{k} = transport_solver_numerics.convergence.equations{1};
    
    % substructure initialisation
    transport_solver_numerics.derivatives_1d{k} = transport_solver_numerics.derivatives_1d{1};
    %time
    transport_solver_numerics.derivatives_1d{k}.time         = profil0d.temps(k);
    transport_solver_numerics.convergence.equations{k}.time  = profil0d.temps(k);
    
    % grid
    transport_solver_numerics.derivatives_1d{k}.grid  = mapprofiles1d_grid_imas(profil0d,k,transport_solver_numerics.derivatives_1d{1}.grid);
    transport_solver_numerics.derivatives_1d{k}.grid.psi = factor_psi .* transport_solver_numerics.derivatives_1d{k}.grid.psi;
    % time derivative
    transport_solver_numerics.derivatives_1d{k}.drho_tor_dt = profil0d_dt.rmx(k,:);
    transport_solver_numerics.derivatives_1d{k}.d_dvolume_drho_tor_dt  = profil0d_dt.vpr_tor(k,:);
    transport_solver_numerics.derivatives_1d{k}.dpsi_dt = factor_psi .*  profil0d_dt.psi(k,:);
 
    
    % only for full simulation; data are not available in time to time run
    if length(profil0d.temps) > 1
        transport_solver_numerics.derivatives_1d{k}.d_dt.n_i_total_over_n_e = z0dxdt_imas(nione,profil0d.temps,k);
        transport_solver_numerics.derivatives_1d{k}.d_dt.pressure_ion_total = z0dxdt_imas(phys.e .* profil0d.nip(:,:) .* profil0d.tip(:,:) + sum(psupra(:,:,:),3),profil0d.temps,k);
        transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.density_fast  = z0dxdt_imas(nfast(:,:,end),profil0d.temps,k);
        transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.pressure = z0dxdt_imas(phys.e .* profil0d.tep(:,:) .* profil0d.nep(:,:) + psupra(:,:,9),profil0d.temps,k);
        transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.pressure_fast_perpendicular = z0dxdt_imas(psupra_perp(:,:,9),profil0d.temps,k);
        transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.pressure_fast_parallel = z0dxdt_imas(psupra_perp(:,:,9),profil0d.temps,k);
    end
    transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.temperature = profil0d_dt.tep(k,:);
    transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.density_thermal  = profil0d_dt.nep(k,:);
    if length(profil0d.temps) > 1
	transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.density = transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.density_thermal +  ...   
	                                                                     transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.density_fast;
    else
	transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.density = transport_solver_numerics.derivatives_1d{k}.electrons.d_dt.density_thermal;
    end
    %
    transport_solver_numerics.derivatives_1d{k}.dpsi_drho_tor   = factor_psi .* pdederive(profil0d.rmx(k,:),profil0d.psi(k,:),0,2,2,1);
    transport_solver_numerics.derivatives_1d{k}.dpsi_drho_tor(end) = -(2*pi) .* phys.mu0 .* factor_psi .* profil0d.rmx(k,end) .* interp1_imas(data_zerod.temps,data_zerod.ip,profil0d.temps(k),'pchip','extrap') ./ profil0d.C2(k,end);

    %datak.prof.dpsidt_phi  = datak.prof.dpsidt3p - gene.x ./ datak.equi.rhomax .* datak.equi.drhomaxdt .*datak.prof.psid1;
    inter = (ones(size(z0dstruct.profil0d.temps)) * z0dstruct.profil0d.xli) ./ ...
        (z0dstruct.profil0d.rmx(:,end) * ones(size(z0dstruct.profil0d.xli))) .* ...
        (z0dxdt(z0dstruct.profil0d.rmx(:,end),z0dstruct.profil0d.temps) * ones(size(z0dstruct.profil0d.xli))) + ...
        (ones(size(z0dstruct.profil0d.temps)) * z0dstruct.profil0d.xli) ./  ...
        ((z0dstruct.profil0d.fdia(:,end) ./ z0dstruct.profil0d.Raxe(:,end)) * ones(size(z0dstruct.profil0d.xli))) .* ...
        (z0dxdt(z0dstruct.profil0d.fdia(:,end) ./ z0dstruct.profil0d.Raxe(:,end),z0dstruct.profil0d.temps) * ones(size(z0dstruct.profil0d.xli)));
    %
    transport_solver_numerics.derivatives_1d{k}.dpsi_dt_cphi = transport_solver_numerics.derivatives_1d{k}.dpsi_dt - ...
        factor_psi .* interp1_imas(z0dstruct.profil0d.temps,inter,profil0d.temps(k),'pchip','extrap')  .*  ...
        transport_solver_numerics.derivatives_1d{k}.dpsi_drho_tor .* z0dstruct.profil0d.rmx(k,end);

    transport_solver_numerics.derivatives_1d{k}.dpsi_dt_crho_tor_norm = transport_solver_numerics.derivatives_1d{k}.dpsi_dt - transport_solver_numerics.derivatives_1d{k}.drho_tor_dt .* transport_solver_numerics.derivatives_1d{k}.dpsi_drho_tor;

    transport_solver_numerics.derivatives_1d{k}.d2psi_drho_tor2 = pdederive(profil0d.rmx(k,:),profil0d.psi(k,:),1,0,2,2);
    if z0dstruct.z0dinput.option.cronos_regul == 4
        % near magnetic axis
        q_not = profil0d.fdia(k,end) ./  profil0d.Raxe(k,end) ./ (profil0d.psi(k,1) - profil0d.psi(k,2)) .* ...
            factor_psi .* profil0d.rmx(k,end) .^ 2 .* profil0d.xli(1,2) .^ 2 ./ 2;
        transport_solver_numerics.derivatives_1d{k}.d2psi_drho_tor2(1) = - factor_psi .* profil0d.fdia(k,end) ./  profil0d.Raxe(k,end) ./ ...
            q_not .*  profil0d.rmx(k,end) .^ 2;
    end
    
    transport_solver_numerics.derivatives_1d{k}.ion = ion;
    transport_solver_numerics.derivatives_1d{k}.d_drho_tor_norm.n_i_total_over_n_e  = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),nione(k,:),0,2,2,1);
    transport_solver_numerics.derivatives_1d{k}.d2_drho_tor_norm2.n_i_total_over_n_e = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),nione(k,:),1,0,2,2);
    %
    transport_solver_numerics.derivatives_1d{k}.d_drho_tor_norm.pressure_ion_total   = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),phys.e .* profil0d.nip(k,:) .* profil0d.tip(k,:) + sum(psupra(k,:,:),3),0,2,2,1);
    transport_solver_numerics.derivatives_1d{k}.d2_drho_tor_norm2.pressure_ion_total = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),phys.e .* profil0d.nip(k,:) .* profil0d.tip(k,:) + sum(psupra(k,:,:),3),1,0,2,2);
    
    transport_solver_numerics.derivatives_1d{k}.electrons.d_drho_tor_norm.temperature = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),profil0d.tep(k,:),0,2,2,1);
    transport_solver_numerics.derivatives_1d{k}.electrons.d_drho_tor_norm.density_thermal = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),profil0d.nep(k,:),0,2,2,1);
    transport_solver_numerics.derivatives_1d{k}.electrons.d_drho_tor_norm.density_fast = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),nfast(k,:,9),0,2,2,1);
    transport_solver_numerics.derivatives_1d{k}.electrons.d_drho_tor_norm.density = transport_solver_numerics.derivatives_1d{k}.electrons.d_drho_tor_norm.density_thermal + ...
                                                                                    transport_solver_numerics.derivatives_1d{k}.electrons.d_drho_tor_norm.density_fast;
                                                                                    
    transport_solver_numerics.derivatives_1d{k}.electrons.d_drho_tor_norm.pressure = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),phys.e .* profil0d.tep(k,:) .* profil0d.nep(k,:) + psupra(k,:,9),0,2,2,1);
    transport_solver_numerics.derivatives_1d{k}.electrons.d_drho_tor_norm.pressure_fast_perpendicular = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),psupra_perp(k,:,9),0,2,2,1);
    transport_solver_numerics.derivatives_1d{k}.electrons.d_drho_tor_norm.pressure_fast_parallel = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),psupra_para(k,:,9),0,2,2,1);
    transport_solver_numerics.derivatives_1d{k}.electrons.d2_drho_tor_norm2.temperature = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),profil0d.tep(k,:),1,0,2,2);
    transport_solver_numerics.derivatives_1d{k}.electrons.d2_drho_tor_norm2.density_thermal = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),profil0d.nep(k,:),1,0,2,2);
    transport_solver_numerics.derivatives_1d{k}.electrons.d2_drho_tor_norm2.density_fast = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),nfast(k,:,9),1,0,2,2);
    transport_solver_numerics.derivatives_1d{k}.electrons.d2_drho_tor_norm2.density = transport_solver_numerics.derivatives_1d{k}.electrons.d2_drho_tor_norm2.density_thermal + ...
                                                                                      transport_solver_numerics.derivatives_1d{k}.electrons.d2_drho_tor_norm2.density_fast;
    transport_solver_numerics.derivatives_1d{k}.electrons.d2_drho_tor_norm2.pressure = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),phys.e .* profil0d.tep(k,:) .* profil0d.nep(k,:) + psupra(k,:,9),1,0,2,2);
    transport_solver_numerics.derivatives_1d{k}.electrons.d2_drho_tor_norm2.pressure_fast_perpendicular = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),psupra_perp(k,:,9),1,0,2,2);
    transport_solver_numerics.derivatives_1d{k}.electrons.d2_drho_tor_norm2.pressure_fast_parallel = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),psupra_para(k,:,9),1,0,2,2);
    
    % ion sub structure
    for l=1:length(label)
        %
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.temperature = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),profil0d.tip(k,:),0,2,2,1);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.density_thermal = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),nions(k,:,l),0,2,2,1);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.density_fast = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),nfast(k,:,l),0,2,2,1);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.density = transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.density_thermal + ...
                                                                                     transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.density_fast;
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.pressure = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),phys.e .* nions(k,:,l) .* profil0d.tip(k,:) + psupra(k,:,l),0,2,2,1);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.pressure_fast_perpendicular = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),psupra_perp(k,:,l),0,2,2,1);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.pressure_fast_parallel = pdederive(profil0d.rmx(k,:),psupra_para(k,:,l),0,2,2,1);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.velocity_tor = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),vtor(k,:,l),0,2,2,1);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_drho_tor_norm.velocity_pol = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),vpol(k,:,l),0,2,2,1);
        
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.temperature = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),profil0d.tip(k,:),1,0,2,2);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.density_thermal = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),nions(k,:,l),1,0,2,2);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.density_fast = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),nfast(k,:,l),1,0,2,2);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.density = transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.density_thermal + ...
                                                                                       transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.density_fast;
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.pressure = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),phys.e .* nions(k,:,l) .* profil0d.tip(k,:) + psupra(k,:,l),1,0,2,2);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.pressure_fast_perpendicular = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),psupra_perp(k,:,l),1,0,2,2);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.pressure_fast_parallel = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),psupra_para(k,:,l),1,0,2,2);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.velocity_tor = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),vtor(k,:,l),1,0,2,2);
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d2_drho_tor_norm2.velocity_pol = pdederive(profil0d.rmx(k,:) ./ profil0d.rmx(k,end),vpol(k,:,l),1,0,2,2);
        
        % only for full simulation; data are not available in time to time run
        transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.temperature = profil0d_dt.tip(k,:);
        if length(profil0d.temps) > 1
            transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.density_thermal = z0dxdt_imas(nions(:,:,l),profil0d.temps,k);
            transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.density_fast = z0dxdt_imas(nfast(:,:,l),profil0d.temps,k);
            transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.density      = transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.density_thermal + ...
                                                                                   transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.density_fast;
            transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.pressure = z0dxdt_imas(phys.e .* nions(:,:,l) .* profil0d.tip(:,:) + psupra(:,:,l),profil0d.temps,k);
            transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.pressure_fast_perpendicular = z0dxdt_imas(psupra_perp(:,:,l),profil0d.temps,k);
            transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.pressure_fast_parallel = z0dxdt_imas(psupra_para(:,:,l),profil0d.temps,k);
            transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.velocity_tor = z0dxdt_imas(vtor(:,:,l),profil0d.temps,k);
            transport_solver_numerics.derivatives_1d{k}.ion{l}.d_dt.velocity_pol = z0dxdt_imas(vpol(:,:,l),profil0d.temps,k);
        end
    end
    
    % initialisation substructure 
    transport_solver_numerics.boundary_conditions_1d{k} = transport_solver_numerics.boundary_conditions_1d{1};
    %
    transport_solver_numerics.boundary_conditions_1d{k}.time = profil0d.temps(k);
    if asser(k) == 0
        transport_solver_numerics.boundary_conditions_1d{k}.current.identifier.name         = 'ip';     
        transport_solver_numerics.boundary_conditions_1d{k}.current.identifier.index        = 2;
        transport_solver_numerics.boundary_conditions_1d{k}.current.identifier.description  = 'plasma current';
        transport_solver_numerics.boundary_conditions_1d{k}.current.value                   = ip(k);
    else
        transport_solver_numerics.boundary_conditions_1d{k}.current.identifier.name         = 'Psi_edge';
        transport_solver_numerics.boundary_conditions_1d{k}.current.identifier.index        =  1;
        transport_solver_numerics.boundary_conditions_1d{k}.current.identifier.description  = 'edge poloidal flux';
        transport_solver_numerics.boundary_conditions_1d{k}.current.value                   =  profil0d.psi(k,end);
    end
    transport_solver_numerics.boundary_conditions_1d{k}.current.rho_tor_norm = 1;
    
    
    transport_solver_numerics.boundary_conditions_1d{k}.ion = ion2;
    % ion sub structure
    for l=1:length(label)
	transport_solver_numerics.boundary_conditions_1d{k}.ion{l}.particles.value = nions(k,end,l);
	transport_solver_numerics.boundary_conditions_1d{k}.ion{l}.particles.rho_tor_norm = 1;
	transport_solver_numerics.boundary_conditions_1d{k}.ion{l}.energy.value = profil0d.tip(k,end);
	transport_solver_numerics.boundary_conditions_1d{k}.ion{l}.energy.rho_tor_norm = 1;
    end
    
    transport_solver_numerics.boundary_conditions_1d{k}.energy_ion_total.identifier.name = 'temperature';
    transport_solver_numerics.boundary_conditions_1d{k}.energy_ion_total.identifier.index = 1;
    transport_solver_numerics.boundary_conditions_1d{k}.energy_ion_total.identifier.description = 'temperature (ev)';
    transport_solver_numerics.boundary_conditions_1d{k}.energy_ion_total.value = profil0d.tip(k,end);
    transport_solver_numerics.boundary_conditions_1d{k}.energy_ion_total.rho_tor_norm = 1;
    
    transport_solver_numerics.boundary_conditions_1d{k}.momentum_tor.identifier.name = 'momentum';
    transport_solver_numerics.boundary_conditions_1d{k}.momentum_tor.identifier.index = 1;
    transport_solver_numerics.boundary_conditions_1d{k}.momentum_tor.identifier.description = 'momentum (kg.m.s^-1)';
    transport_solver_numerics.boundary_conditions_1d{k}.momentum_tor.value = profil0d.rtor(k,end);
    transport_solver_numerics.boundary_conditions_1d{k}.momentum_tor.rho_tor_norm = 1;
    
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.particles.identifier.name = 'density';
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.particles.identifier.index = 1;
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.particles.identifier.description = 'density (m^-3)';
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.particles.value= profil0d.nep(k,end);
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.particles.rho_tor_norm = 1;
    
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.energy.identifier.name = 'temperature';
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.energy.identifier.index = 1;
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.energy.identifier.description = 'temperature (ev)';
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.energy.value = profil0d.tep(k,end);
    transport_solver_numerics.boundary_conditions_1d{k}.electrons.energy.rho_tor_norm = 1;
    
    transport_solver_numerics.convergence.equations{k}.current.iterations_n = fix(interp1_imas(data_zerod.temps,data_zerod.difcurconv,profil0d.temps(k),'pchip','extrap'));
    transport_solver_numerics.convergence.equations{k}.ion = ion3;
    % ion sub structure
    for l=1:length(label)
        transport_solver_numerics.convergence.equations{k}.ion{l}.particles.iterations_n = data_zerod.nb;
        transport_solver_numerics.convergence.equations{k}.ion{l}.energy.iterations_n    = data_zerod.nb;
    end
    
    transport_solver_numerics.convergence.equations{k}.energy_ion_total.iterations_n = data_zerod.nb;
    transport_solver_numerics.convergence.equations{k}.electrons.particles.iterations_n = data_zerod.nb;
    transport_solver_numerics.convergence.equations{k}.electrons.energy.iterations_n = data_zerod.nb;
    
end

% data in convergences
transport_solver_numerics.convergence.time_step.data = profil0d_dt.temps;
transport_solver_numerics.convergence.time_step.time = profil0d.temps;

% new version
% time step
transport_solver_numerics.time_step.data = cat(1,diff(profil0d.temps),-9.99999e99);
transport_solver_numerics.time_step.time = profil0d.temps;

% solver
switch z0dstruct.z0dinput.option.mode_expo_inte
case 1
	transport_solver_numerics.solver.name         = 'EXPO_INTE';
	transport_solver_numerics.solver.description  = 'exponential integrateur';
	transport_solver_numerics.solver.index        = 1;
otherwise
	transport_solver_numerics.solver.name         = 'IFD';
	transport_solver_numerics.solver.description  = 'implicit finite difference';
	transport_solver_numerics.solver.index        = 0;
end

% coordinate
transport_solver_numerics.primary_coordinate.name   = 'rho_tor';
transport_solver_numerics.primary_coordinate.index  = 2;
transport_solver_numerics.primary_coordinate.description = 'toroidal flux coordinate';


% solver_1d 
%  initialisation substructure
solver_1d = transport_solver_numerics.solver_1d{1};
equation = solver_1d.equation{1};
primary_quantity = equation.primary_quantity;
boundary_condition = equation.boundary_condition{1};

for k=1:length(transport_solver_numerics.time)
   % case eviolution for time derivative
   if length(profil0d) == 1
        index = length(z0dstruct.profil0d.temps);
   else
        index = k;
   end
   % start filling
   solver_1d.time                = profil0d.temps(k);
   solver_1d.grid                = mapprofiles1d_grid_imas(profil0d,k,solver_1d.grid);
   solver_1d.grid.psi            = factor_psi .* solver_1d.grid.psi;

   % others data
   solver_1d.drho_tor_dt                    = profil0d_dt.rmx(k,:);
   solver_1d.d_dvolume_drho_tor_dt          = profil0d_dt.vpr_tor(k,:);


   %%%%%%%
   % flux
   %%%%%%%
   primary_quantity.identifier.name         = 'flux';
   primary_quantity.identifier.index        = 1; % convention ?
   primary_quantity.identifier.description  = 'core_profiles/grid/psi';
   primary_quantity.profile                 = factor_psi .*  profil0d.psi(k,:);
   primary_quantity.d_dr                    = transport_solver_numerics.derivatives_1d{k}.dpsi_drho_tor;
   primary_quantity.d2_dr2                  = pdederive(profil0d.rmx(k,:),profil0d.psi(k,:),1,0,2,2);
   primary_quantity.d_dt                    = profil0d_dt.psi(k,:);
   primary_quantity.d_dt_cphi               = transport_solver_numerics.derivatives_1d{k}.dpsi_dt_cphi;
   primary_quantity.d_dt_cr                 = primary_quantity.d_dt - solver_1d.d_dvolume_drho_tor_dt .* transport_solver_numerics.derivatives_1d{k}.dpsi_drho_tor; 
   equation.primary_quantity                = primary_quantity;

   % centre
   boundary_condition.identifier.name       = 'centre_symetry';
   boundary_condition.identifier.index      = 2;
   boundary_condition.dentifier.description = 'null derivative at the centre';
   boundary_condition.type.name             = 'centre_symetry';
   boundary_condition.type.index            = 2;
   boundary_condition.type.description      = 'null derivative at the centre';
   boundary_condition.value                 = 0;
   boundary_condition.position              = 0;
   equation.boundary_condition{1}           = boundary_condition;

   % bord
    if asser(k) == 0
        boundary_condition.identifier.name         = 'ip';     
        boundary_condition.identifier.index        =  2;
        boundary_condition.identifier.description   =  'plasma current';
        boundary_condition.type.name               = 'ip';     
        boundary_condition.type.index              =  2;
        boundary_condition.type.description        =  'plasma current';
        boundary_condition.value                   =  ip(k);
        boundary_condition.position                =  profil0d.rmx(k,end);
    else
        boundary_condition.identifier.name         = 'Psi_edge';     
        boundary_condition.identifier.index        =  1;
        boundary_condition.identifier.description   =  'edge poloidal flux';
        boundary_condition.type.name               = 'Psi_edge';     
        boundary_condition.type.index              =  1;
        boundary_condition.type.description        =  'edge poloidal flux';
        boundary_condition.value                   =  profil0d.psi(k,end);
        boundary_condition.position                =  profil0d.rmx(k,end);
     end
     equation.boundary_condition{2}           = boundary_condition;

   % tranport coefficient
   % D
   equation.coefficient{1}.profile          = 1./max(1e-307,profil0d.eta(k,:));
   equation.coefficient = equation.coefficient(1);
   % V
   %equation.coefficient{2}.profile          =

   % not used by METIS
   % convergence
   equation.convergence.iterations_n       = fix(interp1_imas(data_zerod.temps,data_zerod.difcurconv,profil0d.temps(k),'pchip','extrap'));
   solver_1d.equation{1} = equation;
   % end of flux

   %%%%%%%%%%%
   % density
   %%%%%%%%%%%
   primary_quantity.identifier.name         = 'density';
   primary_quantity.identifier.index        = 2; % convention ?
   primary_quantity.identifier.description  = 'core_profiles/profiles_1d/electrons/density';
   primary_quantity.profile                 = profil0d.nep(k,:);
   primary_quantity.d_dr                    = pdederive(profil0d.rmx(k,:),profil0d.nep(k,:),0,2,2,1);
   primary_quantity.d2_dr2                  = pdederive(profil0d.rmx(k,:),profil0d.nep(k,:),1,0,2,2);
   primary_quantity.d_dt                    = profil0d_dt.nep(k,:);
   %primary_quantity.d_dt_cphi               = 
   %primary_quantity.d_dt_cr                 = 
   equation.primary_quantity                = primary_quantity;

   % centre
   boundary_condition.identifier.name       = 'centre_symetry';
   boundary_condition.identifier.index      = 2;
   boundary_condition.identifier.description = 'null derivative at the centre';
   boundary_condition.type.name             = 'centre_symetry';
   boundary_condition.type.index            = 2;
   boundary_condition.type.description      = 'null derivative at the centre';
   boundary_condition.value                 = 0;
   boundary_condition.position              = 0;
   equation.boundary_condition{1}           = boundary_condition;

   % bord
   boundary_condition.identifier.name         = 'ne_bord';     
   boundary_condition.identifier.index        =  2;
   boundary_condition.identifier.description   =  'edge electron density';
   boundary_condition.type.name               = 'ne_bord';     
   boundary_condition.type.index              =  2;
   boundary_condition.type.description        =  'edge electron density';
   boundary_condition.value                   =  profil0d.nep(k,end);
   boundary_condition.position                =  profil0d.rmx(k,end);
   equation.boundary_condition{2}             = boundary_condition;

   % tranport coefficient
   % D
   equation.coefficient{1}.profile          = profil0d.dn(k,:);
   % V
   equation.coefficient{2} = equation.coefficient{1};
   equation.coefficient{2}.profile          = -profil0d.vn(k,:);

   % not used by METIS
   %equation.control_parameters.integer{1}.value  = -999999999;
   %equation.control_parameters.real{1}.value     = -9.99999e99;

   % convergence
   equation.convergence.iterations_n       = 1;
   solver_1d.equation{2} = equation;
   % end of density 


   %%%%%%%%%%%%%%%%%%%%%%
   % electron tempreature
   %%%%%%%%%%%%%%%%%%%%%%
   primary_quantity.identifier.name         = 'electron_temperature';
   primary_quantity.identifier.index        = 3; % convention ?
   primary_quantity.identifier.description  = 'core_profiles/profiles_1d/electrons/temperature';
   primary_quantity.profile                 = profil0d.tep(k,:);
   primary_quantity.d_dr                    = pdederive(profil0d.rmx(k,:),profil0d.tep(k,:),0,2,2,1);
   primary_quantity.d2_dr2                  = pdederive(profil0d.rmx(k,:),profil0d.tep(k,:),1,0,2,2);
   primary_quantity.d_dt                    = profil0d_dt.tep(k,:);
   %primary_quantity.d_dt_cphi               = 
   %primary_quantity.d_dt_cr                 = 
   equation.primary_quantity                = primary_quantity;

   % centre
   boundary_condition.identifier.name       = 'centre_symetry';
   boundary_condition.identifier.index      = 2;
   boundary_condition.identifier.description = 'null derivative at the centre';
   boundary_condition.type.name             = 'centre_symetry';
   boundary_condition.type.index            = 2;
   boundary_condition.type.description      = 'null derivative at the centre';
   boundary_condition.value                 = 0;
   boundary_condition.position              = 0;
   equation.boundary_condition{1}           = boundary_condition;

   % bord
   boundary_condition.identifier.name         = 'te_bord';     
   boundary_condition.identifier.index        =  2;
   boundary_condition.identifier.description   =  'edge electron temperature';
   boundary_condition.type.name               = 'te_bord';     
   boundary_condition.type.index              =  2;
   boundary_condition.type.description        =  'edge electron dtemperature';
   boundary_condition.value                   =  profil0d.nep(k,end);
   boundary_condition.position                =  profil0d.rmx(k,end);
   equation.boundary_condition{2}             = boundary_condition;

   % tranport coefficient
   % D
   equation.coefficient{1}.profile          = profil0d.xie(k,:);
   equation.coefficient = equation.coefficient(1);
   % V
   %equation.coefficient{2}.profile          = 

   % not used by METIS
   %equation.control_parameters.integer{1}.value  = -999999999;
   %equation.control_parameters.real{1}.value     = -9.99999e99;

   % convergence
   equation.convergence.iterations_n       = 1;
   solver_1d.equation{3} = equation;
   % end of electron temperature

   %%%%%%%%%%%%%%%%%%%%%%
   % ion tempreature
   %%%%%%%%%%%%%%%%%%%%%%
   primary_quantity.identifier.name         = 'electron_temperature';
   primary_quantity.identifier.index        = 4; % convention ?
   primary_quantity.identifier.description  = 'core_profiles/profiles_1d/electrons/temperature';
   primary_quantity.profile                 = profil0d.tip(k,:);
   primary_quantity.d_dr                    = pdederive(profil0d.rmx(k,:),profil0d.tip(k,:),0,2,2,1);
   primary_quantity.d2_dr2                  = pdederive(profil0d.rmx(k,:),profil0d.tip(k,:),1,0,2,2);
   primary_quantity.d_dt                    = profil0d_dt.tip(k,:);
   %primary_quantity.d_dt_cphi               = 
   %primary_quantity.d_dt_cr                 = 
   equation.primary_quantity                = primary_quantity;

   % centre
   boundary_condition.identifier.name       = 'centre_symetry';
   boundary_condition.identifier.index      = 2;
   boundary_condition.identifier.description = 'null derivative at the centre';
   boundary_condition.type.name             = 'centre_symetry';
   boundary_condition.type.index            = 2;
   boundary_condition.type.description      = 'null derivative at the centre';
   boundary_condition.value                 = 0;
   boundary_condition.position              = 0;
   equation.boundary_condition{1}           = boundary_condition;

   % bord
   boundary_condition.identifier.name         = 'te_bord';     
   boundary_condition.identifier.index        =  2;
   boundary_condition.identifier.description   =  'edge electron temperature';
   boundary_condition.type.name               = 'te_bord';     
   boundary_condition.type.index              =  2;
   boundary_condition.type.description        =  'edge electron dtemperature';
   boundary_condition.value                   =  profil0d.nip(k,end);
   boundary_condition.position                =  profil0d.rmx(k,end);
   equation.boundary_condition{2}             = boundary_condition;

   % tranport coefficient
   % D
   equation.coefficient{1}.profile          = profil0d.xii(k,:);
   equation.coefficient = equation.coefficient(1);
   % V
   %equation.coefficient{2}.profile          = 

   % not used by METIS
   %equation.control_parameters.integer{1}.value  = -999999999;
   %equation.control_parameters.real{1}.value     = -9.99999e99;

   % convergence
   equation.convergence.iterations_n       = 1;
   solver_1d.equation{4} = equation;
   % end of ion temperature 

   %%%%%%%%%%%%%%%%%%%%%%
   % rotation
   %%%%%%%%%%%%%%%%%%%%%%
   primary_quantity.identifier.name         = 'momentum_tor';
   primary_quantity.identifier.index        = 5; % convention ?
   primary_quantity.identifier.description  = 'core_profiles/profiles_1d/momentum_tor';
   primary_quantity.profile                 = rtor(k,:);
   primary_quantity.d_dr                    = pdederive(profil0d.rmx(k,:),rtor(k,:),0,2,2,1);
   primary_quantity.d2_dr2                  = pdederive(profil0d.rmx(k,:),rtor(k,:),1,0,2,2);
   primary_quantity.d_dt                    = drtordt(k,:);
   %primary_quantity.d_dt_cphi               = 
   %primary_quantity.d_dt_cr                 = 
   equation.primary_quantity                = primary_quantity;

   % centre
   boundary_condition.identifier.name       = 'centre_symetry';
   boundary_condition.identifier.index      = 2;
   boundary_condition.identifier.description = 'null derivative at the centre';
   boundary_condition.type.name             = 'centre_symetry';
   boundary_condition.type.index            = 2;
   boundary_condition.type.description      = 'null derivative at the centre';
   boundary_condition.value                 = 0;
   boundary_condition.position              = 0;
   equation.boundary_condition{1}           = boundary_condition;

   % bord
   boundary_condition.identifier.name         = 'te_bord';     
   boundary_condition.identifier.index        =  2;
   boundary_condition.identifier.description   =  'edge electron temperature';
   boundary_condition.type.name               = 'te_bord';     
   boundary_condition.type.index              =  2;
   boundary_condition.type.description        =  'edge electron dtemperature';
   boundary_condition.value                   =  rtor(k,end);
   boundary_condition.position                =  profil0d.rmx(k,end);
   equation.boundary_condition{2}             = boundary_condition;

   % tranport coefficient
   % D
   equation.coefficient{1}.profile          = profil0d.drot(k,:);
   % V
   equation.coefficient{2} = equation.coefficient{1};
   equation.coefficient{2}.profile          = - profil0d.vrot(k,:);

   % not used by METIS
   %equation.control_parameters.integer{1}.value  = -999999999;
   %equation.control_parameters.real{1}.value     = -9.99999e99;

   % convergence
   equation.convergence.iterations_n       = 1;
   solver_1d.equation{5} = equation;
   % end of rotation 

   % recopy for each time step
   transport_solver_numerics.solver_1d{k} = solver_1d;  %Numerics related to 1D radial solver, for various time slices. {dynamic}  struct_array [max_size=unbounded]

   %
   transport_solver_numerics.restart_files{k}.names{1}            = '/dev/null';
   transport_solver_numerics.restart_files{k}.descriptions{1}     = 'Why code does nott used IDSs for restarting ? METIS needs only IDSs for restarting'; 
   transport_solver_numerics.restart_files{k}.names{2}            = '/dev/random';
   transport_solver_numerics.restart_files{k}.descriptions{2}     = 'sometime, this file provide restart data'; 
   transport_solver_numerics.restart_files{k}.time                 = profil0d.temps(k);




end


function out = z0dxdt_imas(x,t,k)

out = z0dxdt(x,t);
out = out(k,:);