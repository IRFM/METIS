% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function [coresource_ohm,coresource_neo] = mapcore_sources_ohm_boot_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,coresource)

%% output initilalisation
coresource_ohm = coresource;
coresource_neo = coresource;

%% 1 - OHM
%% IDENTIFIER of the source
coresource_ohm.identifier.name  = 'ohmic';
coresource_ohm.identifier.index = 7;
coresource_ohm.identifier.description='ohmic power from METIS model';
% null profile
zz              = zeros(size(profil0d.xli));
% initialisation substructure
profiles_1d = coresource.profiles_1d{1};
global_quantities = coresource.global_quantities{1};
%% SOURCE terms
for k=1:length(profil0d.temps)

	profiles_1d.grid                         = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
	profiles_1d.electrons.particles          = zz;
	profiles_1d.electrons.energy             = profil0d.pohm(k,:);
        profiles_1d.total_ion_energy             = zz;
        profiles_1d.momentum_tor                 = zz;
        profiles_1d.j_parallel                   = zz;
        profiles_1d.conductivity_parallel        = zz;
        profiles_1d.time                         = profil0d.temps(k);
        
        % donnees complementaires
	profiles_1d.current_parallel_inside       = cumtrapz(profil0d.xli,profil0d.spr(k,:) .* profiles_1d.j_parallel,2);
	profiles_1d.torque_tor_inside             = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.momentum_tor,2);
	profiles_1d.total_ion_power_inside        = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.total_ion_energy,2);
	profiles_1d.electrons.particles_inside    = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.particles,2);
	profiles_1d.electrons.power_inside        = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.energy,2);

        
        
	coresource_ohm.profiles_1d{k} = profiles_1d;
end

% adding global data
for k=1:length(profil0d.temps)
      global_quantities.total_ion_power = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource_ohm.profiles_1d{k}.total_ion_energy,2);
      global_quantities.torque_tor = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource_ohm.profiles_1d{k}.momentum_tor,2);
      global_quantities.current_parallel = trapz(profil0d.xli,profil0d.spr(k,:) .* coresource_ohm.profiles_1d{k}.j_parallel,2);
      global_quantities.time = profil0d.temps(k);
      global_quantities.electrons.particles = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource_ohm.profiles_1d{k}.electrons.particles,2);
      global_quantities.electrons.power = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource_ohm.profiles_1d{k}.electrons.energy,2);
      global_quantities.power = global_quantities.electrons.power + global_quantities.total_ion_power;

      coresource_ohm.global_quantities{k} = global_quantities;
      
end


%% 2 - neoclassical
%% IDENTIFIER of the source
coresource_neo.identifier.name  = 'bootstrap_current';
coresource_neo.identifier.index = 13;
coresource_neo.identifier.description='bootstrap source of current and parallel electric conductivity from METIS model';
% null profile
zz              = zeros(size(profil0d.xli));
%% SOURCE terms
for k=1:length(profil0d.temps)

	profiles_1d.grid                         = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
	profiles_1d.electrons.particles          = zz;
	profiles_1d.electrons.energy             = zz;
        profiles_1d.total_ion_energy             = zz;
        profiles_1d.momentum_tor                 = zz;
        profiles_1d.j_parallel                   = profil0d.jboot(k,:);
        profiles_1d.conductivity_parallel        = 1./max(1e-307,profil0d.eta(k,:));
        profiles_1d.time                         = profil0d.temps(k);
        
        % donnees complementaires
	profiles_1d.current_parallel_inside       = cumtrapz(profil0d.xli,profil0d.spr(k,:) .* profiles_1d.j_parallel,2);
	profiles_1d.torque_tor_inside             = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.momentum_tor,2);
	profiles_1d.total_ion_power_inside        = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.total_ion_energy,2);
	profiles_1d.electrons.particles_inside    = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.particles,2);
	profiles_1d.electrons.power_inside        = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.energy,2);


	coresource_neo.profiles_1d{k} = profiles_1d;
end


% adding global data
for k=1:length(profil0d.temps)
      global_quantities.total_ion_power = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource_neo.profiles_1d{k}.total_ion_energy,2);
      global_quantities.torque_tor = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource_neo.profiles_1d{k}.momentum_tor,2);
      global_quantities.current_parallel = trapz(profil0d.xli,profil0d.spr(k,:) .* coresource_neo.profiles_1d{k}.j_parallel,2);
      global_quantities.time = profil0d.temps(k);
      global_quantities.electrons.particles = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource_neo.profiles_1d{k}.electrons.particles,2);
      global_quantities.electrons.power = trapz(profil0d.xli,profil0d.vpr(k,:) .* coresource_neo.profiles_1d{k}.electrons.energy,2);
      global_quantities.power = global_quantities.electrons.power + global_quantities.total_ion_power;

      coresource_neo.global_quantities{k} = global_quantities;
      
end
