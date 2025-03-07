% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcore_sources_full_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,coresource)

% preprocessing
warning off
pei = pdederive(profil0d.xli,profil0d.qei,2,2,2,1) ./ profil0d.vpr;
warning on
pei(:,1) = 2 .* pei(:,2) - pei(:,3);
% improve precision on conservated quantities
qei = cumtrapz(profil0d.xli,pei .* profil0d.vpr,2);
pei = pei ./ ((qei(:,end) ./ profil0d.qei(:,end)) * ones(1,length(profil0d.xli)));


%% IDENTIFIER of the source
coresource.identifier.name  = 'total';
coresource.identifier.index = 1;
coresource.identifier.description='Total source: combines all sources; i.e. sum of sources and wells of power, current, particles and momentum from METIS models';
% null profile
zz              = zeros(size(profil0d.xli));
% initialisation substructure
profiles_1d = coresource.profiles_1d{1};
global_quantities = coresource.global_quantities{1};
ion = profiles_1d.ion{1};
element = ion.element{1};
%% SOURCE terms
for k=1:length(profil0d.temps)
    
    profiles_1d.grid                         = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
    profiles_1d.electrons.particles          = real(profil0d.nbinesource(k,:)) + imag(profil0d.nbinesource(k,:)) +  ...
        profil0d.s0(k,:) + profil0d.s0m(k,:) + profil0d.spellet(k,:);
    % there is a problem with profil0d.source_el: it can be modified
    % internaly in METIS to ensure the convergence. It is better to
    % recompute the sum localy.
    source_tot = profil0d.pohm(k,:)   - profil0d.prad(k,:)  - profil0d.pbrem(k,:) - profil0d.pcyclo(k,:) - ...
                 profil0d.pioniz(k,:) + profil0d.picrh(k,:) + profil0d.pfweh(k,:) + profil0d.pnbi(k,:)   +  ...
                 profil0d.pfus(k,:)   + profil0d.pecrh(k,:) + profil0d.plh(k,:);
    
    source_ion = profil0d.pnbi_ion(k,:) + profil0d.picrh_ion(k,:) + profil0d.pfus_ion(k,:) - profil0d.pioniz_i(k,:);
    
    source_el = source_tot - source_ion;
    
    profiles_1d.electrons.energy             = source_el  - pei(k,:);
    profiles_1d.total_ion_energy             = source_ion + pei(k,:);
    
%     profiles_1d.electrons.energy             = profil0d.source_el(k,:)  - pei(k,:);
%     profiles_1d.total_ion_energy             = profil0d.source_ion(k,:) + pei(k,:);
    %
    switch z0dstruct.z0dinput.option.gaz
        case 5
            % assuming main reaction D-He3 only
            %
            element.a                                = 4;
            element.z_n                              = 2;
            %
            ion.z_ion                                = 2;
            ion.label                                = 'He4+2';
            ion.name                                 = 'He4+2';
            ion.particles                            = profil0d.salf(k,:);
            ion.energy                               = profil0d.pfus_ion(k,:);
            ion.element{1}                           = element;
            profiles_1d.ion{1}                       = ion;
            %
            element.a                                = 2;
            element.z_n                              = 1;
            %
            ion.z_ion                                = 1;
            ion.label                                = 'D+';
            ion.name                                 = 'D+';
            ion.particles                            = - profil0d.salf(k,:);
            ion.energy                               = -1.602176462e-19 .* profil0d.tip(k,:) .* profil0d.salf(k,:);
            ion.element{1}                           = element;
            profiles_1d.ion{2}                       = ion;
            %
            element.a                                = 11;
            element.z_n                              = 5;
            %
            ion.z_ion                                = 5;
            ion.label                                = 'He3+2';
            ion.name                                 = 'He3+2';
            ion.particles                            = - profil0d.salf(k,:);
            ion.energy                               = -1.602176462e-19 .* profil0d.tip(k,:) .* profil0d.salf(k,:);
            ion.element{1}                           = element;
            profiles_1d.ion{3}                       = ion;
            
        case 11
            %
            element.a                                = 4;
            element.z_n                              = 2;
            %
            ion.z_ion                                = 2;
            ion.label                                = 'He+2';
            ion.name                                 = 'He+2';
            ion.particles                            = profil0d.salf(k,:);
            ion.energy                               = profil0d.pfus_ion(k,:);
            ion.element{1}                           = element;
            profiles_1d.ion{1}                       = ion;
            %
            element.a                                = 1;
            element.z_n                              = 1;
            %
            ion.z_ion                                = 1;
            ion.label                                = 'H+';
            ion.name                                 = 'H+';
            ion.particles                            = - profil0d.salf(k,:) / 3;
            ion.energy                               = -1.602176462e-19 .* profil0d.tip(k,:) .* profil0d.salf(k,:) / 3;
            ion.element{1}                           = element;
            profiles_1d.ion{2}                       = ion;
            %
            element.a                                = 11;
            element.z_n                              = 5;
            %
            ion.z_ion                                = 5;
            ion.label                                = 'B+5';
            ion.name                                 = 'B+5';
            ion.particles                            = - profil0d.salf(k,:) / 3;
            ion.energy                               = -1.602176462e-19 .* profil0d.tip(k,:) .* profil0d.salf(k,:) / 3;
            ion.element{1}                           = element;
            profiles_1d.ion{3}                       = ion;
           
        otherwise
            %
            element.a                                = 4;
            element.z_n                              = 2;
            %
            ion.z_ion                                = 2;
            ion.label                                = 'He+2';
            ion.name                                 = 'He+2';
            ion.particles                            = profil0d.salf(k,:);
            ion.energy                               = profil0d.pfus_ion(k,:);
            ion.element{1}                           = element;
            profiles_1d.ion{1}                       = ion;
            %
            element.a                                = 2;
            element.z_n                              = 1;
            %
            ion.z_ion                                = 1;
            ion.label                                = 'D+';
            ion.name                                 = 'D+';
            ion.particles                            = - profil0d.salf(k,:);
            ion.energy                               = -1.602176462e-19 .* profil0d.tip(k,:) .* profil0d.salf(k,:);
            ion.element{1}                           = element;
            profiles_1d.ion{2}                       = ion;
            %
            element.a                                = 3;
            element.z_n                              = 1;
            %
            ion.z_ion                                = 1;
            ion.label                                = 'T+';
            ion.name                                = 'T+';
            ion.particles                            = - profil0d.salf(k,:);
            ion.energy                               = -1.602176462e-19 .* profil0d.tip(k,:) .* profil0d.salf(k,:);
            ion.element{1}                           = element;
            profiles_1d.ion{3}                       = ion;
    end
    %
    profiles_1d.momentum_tor                 = profil0d.rot_nbi(k,:) + profil0d.rot_n0(k,:) + profil0d.rot_lh(k,:);
    profiles_1d.j_parallel                   = profil0d.jni(k,:);
    profiles_1d.conductivity_parallel        = 1./max(1e-307,profil0d.eta(k,:));
    profiles_1d.time                         = profil0d.temps(k);
    
    % donnees complementaires
    profiles_1d.current_parallel_inside       = cumtrapz(profil0d.xli,profil0d.spr(k,:) .* profiles_1d.j_parallel,2);
    profiles_1d.torque_tor_inside             = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.momentum_tor,2);
    profiles_1d.total_ion_power_inside        = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.total_ion_energy,2);
    profiles_1d.electrons.particles_inside    = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.particles,2);
    profiles_1d.electrons.power_inside        = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.energy,2);
    
    coresource.profiles_1d{k} = profiles_1d;
end

% adding global data
for k=1:length(profil0d.temps)
    global_quantities.total_ion_power = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource.profiles_1d{k}.total_ion_energy,2);
    global_quantities.torque_tor = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource.profiles_1d{k}.momentum_tor,2);
    global_quantities.current_parallel = trapz(profil0d.xli,profil0d.spr(k,:) .* coresource.profiles_1d{k}.j_parallel,2);
    global_quantities.time = profil0d.temps(k);
    global_quantities.electrons.particles = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource.profiles_1d{k}.electrons.particles,2);
    global_quantities.electrons.power = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource.profiles_1d{k}.electrons.energy,2);
    global_quantities.power = global_quantities.electrons.power + global_quantities.total_ion_power;
    
    coresource.global_quantities{k} = global_quantities;
    
end
