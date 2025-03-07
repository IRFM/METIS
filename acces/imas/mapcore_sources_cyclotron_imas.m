% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcore_sources_cyclotron_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,coresource)

%% IDENTIFIER of the source
coresource.identifier.name  = 'synchrotron_radiation';
coresource.identifier.index = 9;
coresource.identifier.description='Cyclotron radiation source of power from METIS model (current source is neglected)';
% null profile
zz              = zeros(size(profil0d.xli));
% initialisation substructure
profiles_1d = coresource.profiles_1d{1};
global_quantities = coresource.global_quantities{1};
%% SOURCE terms
for k=1:length(profil0d.temps)

	profiles_1d.grid                         = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
	profiles_1d.electrons.particles          = zz;
	profiles_1d.electrons.energy             = - profil0d.pcyclo(k,:);

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
