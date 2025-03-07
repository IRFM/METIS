% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcore_sources_fusion_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,coresource)

%% IDENTIFIER of the source
coresource.identifier.name  = 'fusion';
coresource.identifier.index = 6;
coresource.identifier.description='FUSION source of power, current and momentum from METIS model';
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
    profiles_1d.electrons.particles          = zz;
    profiles_1d.electrons.energy             = max(0,profil0d.pfus(k,:) - profil0d.pfus_ion(k,:));
    profiles_1d.total_ion_energy             =  profil0d.pfus_ion(k,:);
    switch z0dstruct.z0dinput.option.gaz
        case 5
            % assuming main reaction D-He3 only
            %
            element.a                                = 4;
            element.z_n                              = 2;
            %
            ion.z_ion                                = 2;
            ion.label                                = 'He4+2';
            ion.name                                = 'He4+2';
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
            ion.name                                = 'He3+2';
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
            ion.name                                = 'He+2';
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
            ion.name                                = 'H+';
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
            ion.name                                = 'B+5';
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
            ion.name                                = 'He+2';
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
            ion.name                                = 'D+';
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
    profiles_1d.momentum_tor                 = zz;
    profiles_1d.j_parallel                   = profil0d.jfus(k,:);
    profiles_1d.conductivity_parallel        = zz;
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


