% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcore_sources_equipartition_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,coresource)


% preprocessing
warning off
pei = pdederive(profil0d.xli,profil0d.qei,2,2,2,1) ./ profil0d.vpr;
warning on
pei(:,1) = 2 .* pei(:,2) - pei(:,3);
% improve precision on conservated quantities
qei = cumtrapz(profil0d.xli,pei .* profil0d.vpr,2);
pei = pei ./ ((qei(:,end) ./ profil0d.qei(:,end)) * ones(1,length(profil0d.xli)));


%% IDENTIFIER of the source
coresource.identifier.name  = 'collisional_equipartition';
coresource.identifier.index = 11;
coresource.identifier.description='Equipartition energy exchange between electron and ions';
% null profile
zz              = zeros(size(profil0d.xli));
% initialisation substructure
profiles_1d = coresource.profiles_1d{1};
global_quantities = coresource.global_quantities{1};
%% SOURCE terms
for k=1:length(profil0d.temps)
    
    profiles_1d.grid                         = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
    profiles_1d.electrons.particles          = zz;
    profiles_1d.electrons.energy             = - pei(k,:);
    profiles_1d.total_ion_energy             =   pei(k,:);
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
