% mapping du IDS radfiation
function radiation = mapradiation_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,radiation,sigma_B0_eff)

% genric part of IDS
radiation.ids_properties.comment = 'METIS radition';
radiation.time		    = profil0d.temps;
rb0 =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
radiation.vacuum_toroidal_field.r0          = mean(z0dstruct.z0dinput.geo.R);
radiation.vacuum_toroidal_field.b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ radiation.vacuum_toroidal_field.r0;

% Sn
if ~isfield(z0dstruct.z0dinput.option,'Sn_fraction')
    z0dstruct.z0dinput.option.Sn_fraction = 0;
end
if ~isfield(z0dstruct.z0dinput.option,'te_max')
    z0dstruct.z0dinput.option.te_max = 1e5;
end

% initialise empty gdd
for k=1:length(profil0d.temps)
  radiation.grid_gg{k} = radiation.grid_ggd{1};
end

%  identifier:
%  Name 	Index 	Description
%  unspecified 	0 	Unspecified emission process
%  nuclear_decay 	6 	Emission from nuclear decay
%  bremsstrahlung 	8 	Emission from bremsstrahlung
%  synchrotron_radiation 	9 	Emission from synchrotron radiation
%  line_radiation 	10 	Emission from line radiation
%  recombination 	11 	Emission from recombination
%  process_model = radiation.process{1};
%  identifier_model = process_model.identifier;
%  global_quantities_model = process_model.global_quantities{1};
%  profiles_1d_model       = process_model.profiles_1d{1};
% null profile
zz              = zeros(size(profil0d.xli));

% compute recombinaison
if length(data_zerod.temps) == 1
    switch z0dstruct.z0dinput.option.gaz
        case 11
            nboron = z0dstruct.zs.nTm;
        otherwise
            nboron = zeros(size(profil0d.temps));
    end
    
    taue = interp1(z0dstruct.zs.temps,z0dstruct.zs.taue,z0dstruct.profil.temps,'nearest','extrap');
    [pradhep,pradimpp,pradmaxp,pradwp,pradhdt,pradb] = zrad0improved(z0dstruct.profil,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax, ...
                                                     z0dstruct.z0dinput.option.rimp,0,taue,z0dstruct.z0dinput.option.Sn_fraction, ...
                                                     max(105,z0dstruct.z0dinput.option.te_max/1e3 + 5),nboron);
    pradhdt = interp1(z0dstruct.profil.temps,pradhdt,profil0d.temps,'nearest','extrap');
    prad    = data_zerod.prad;
    pradsol = data_zerod.pradsol;                                                 
else
    switch z0dstruct.z0dinput.option.gaz
        case 11
            nboron = data_zerod.nTm;
            if length(nboron) ~= length(profil0d.temps)
                nboron =   interp1_imas(z0dstruct.z0dinput.cons.temps,nboron,profil0d.temps,'pchip','extrap');
            end
        otherwise
            nboron = zeros(size(profil0d.temps));
    end
    
    taue = interp1(data_zerod.temps,data_zerod.taue,profil0d.temps,'nearest','extrap');
    [pradhep,pradimpp,pradmaxp,pradwp,pradhdt,pradb] = zrad0improved(profil0d,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax, ...
                                                     z0dstruct.z0dinput.option.rimp,0,taue,z0dstruct.z0dinput.option.Sn_fraction, ...
                                                     max(105,z0dstruct.z0dinput.option.te_max/1e3 + 5),nboron);
    prad    = interp1(data_zerod.temps,data_zerod.prad,profil0d.temps,'nearest','extrap');
    pradsol = interp1(data_zerod.temps,data_zerod.pradsol,profil0d.temps,'nearest','extrap');
end

% initialisatin sub structures
process = radiation.process{1};
profiles_1d = process.profiles_1d{1};
global_quantities = process.global_quantities{1};
%  line_radiation 	10 	Emission from line radiation
identifier.name         = 'line_radiation';
identifier.index        = 10;
identifier.description  = 'Emission from line radiation';
process.identifier      = identifier;
% loop on time
for k=1:length(profil0d.temps)
        % re-compute it  (to separate from line)
        %profiles_1d                           = profiles_1d_model
        profiles_1d.grid                       = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
        profiles_1d.electrons.emissivity       = min(0,- profil0d.prad(k,:) + pradhdt(k,:));
        profiles_1d.electrons.power_inside     = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.emissivity,2);
        profiles_1d.emissivity_ion_total       = zz;
        profiles_1d.power_inside_ion_total     = zz;
        profiles_1d.emissivity_neutral_total   = zz;
        profiles_1d.power_inside_neutral_total = zz;
        %ion: {[1x1 struct]}
        %neutral: {[1x1 struct]}
        profiles_1d.time = profil0d.temps(k);
        process.profiles_1d{k} = profiles_1d;

        global_quantities.inside_lcfs.power_ion_total     = 0;
        global_quantities.inside_lcfs.power_neutral_total = 0;
        global_quantities.inside_lcfs.power_electrons     = profiles_1d.electrons.power_inside(end);
        global_quantities.inside_lcfs.power               = global_quantities.inside_lcfs.power_ion_total + ...
        global_quantities.inside_lcfs.power_neutral_total + global_quantities.inside_lcfs.power_electrons;

        global_quantities.inside_vessel.power_ion_total     = 0;
        global_quantities.inside_vessel.power_neutral_total = 0;
        global_quantities.inside_vessel.power_electrons     = min(0,- prad(k) - pradsol(k)  + trapz(profil0d.xli,profil0d.vpr(k,:) .* pradhdt(k,:),2));
        global_quantities.inside_vessel.power               = global_quantities.inside_vessel.power_ion_total  + ...
        global_quantities.inside_vessel.power_neutral_total + global_quantities.inside_vessel.power_electrons;

        global_quantities.time       = profil0d.temps(k);
	process.global_quantities{k} = global_quantities;             
end
% fill structure
radiation.process{1}      = process;

%  bremsstrahlung 	8 	Emission from bremsstrahlung
identifier.name         = 'bremsstrahlung';
identifier.index        = 8;
identifier.description  = 'Emission from bremsstrahlung';
process.identifier      = identifier;
% loop on time
for k=1:length(profil0d.temps)
        %profiles_1d                           = profiles_1d_model
        profiles_1d.grid                       = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
        profiles_1d.electrons.emissivity       = - profil0d.pbrem(k,:);
        profiles_1d.electrons.power_inside     = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.emissivity,2);
        profiles_1d.emissivity_ion_total       = zz;
        profiles_1d.power_inside_ion_total     = zz;
        profiles_1d.emissivity_neutral_total   = zz;
        profiles_1d.power_inside_neutral_total = zz;
        %ion: {[1x1 struct]}
        %neutral: {[1x1 struct]}
        profiles_1d.time = profil0d.temps(k);
        process.profiles_1d{k} = profiles_1d;

        global_quantities.inside_lcfs.power_ion_total     = 0;
        global_quantities.inside_lcfs.power_neutral_total = 0;
        global_quantities.inside_lcfs.power_electrons     = profiles_1d.electrons.power_inside(end);
        global_quantities.inside_lcfs.power               = global_quantities.inside_lcfs.power_ion_total + ...
        global_quantities.inside_lcfs.power_neutral_total + global_quantities.inside_lcfs.power_electrons;

        global_quantities.inside_vessel.power_ion_total     = 0;
        global_quantities.inside_vessel.power_neutral_total = 0;
        global_quantities.inside_vessel.power_electrons     = profiles_1d.electrons.power_inside(end);
        global_quantities.inside_vessel.power               = global_quantities.inside_vessel.power_ion_total  + ...
        global_quantities.inside_vessel.power_neutral_total + global_quantities.inside_vessel.power_electrons;

        global_quantities.time       = profil0d.temps(k);
	process.global_quantities{k} = global_quantities;             
end
% fill structure
radiation.process{2}      = process;



%  synchrotron_radiation 	9 	Emission from synchrotron radiation
identifier.name         = 'synchrotron_radiation';
identifier.index        = 9;
identifier.description  = 'Emission from synchrotron radiation';
process.identifier      = identifier;
% loop on time
for k=1:length(profil0d.temps)
        %profiles_1d                           = profiles_1d_model
        profiles_1d.grid                       = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
        profiles_1d.electrons.emissivity       = - profil0d.pcyclo(k,:);
        profiles_1d.electrons.power_inside     = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.emissivity,2);
        profiles_1d.emissivity_ion_total       = zz;
        profiles_1d.power_inside_ion_total     = zz;
        profiles_1d.emissivity_neutral_total   = zz;
        profiles_1d.power_inside_neutral_total = zz;
        %ion: {[1x1 struct]}
        %neutral: {[1x1 struct]}
        profiles_1d.time = profil0d.temps(k);
        process.profiles_1d{k} = profiles_1d;

        global_quantities.inside_lcfs.power_ion_total     = 0;
        global_quantities.inside_lcfs.power_neutral_total = 0;
        global_quantities.inside_lcfs.power_electrons     = profiles_1d.electrons.power_inside(end);
        global_quantities.inside_lcfs.power               = global_quantities.inside_lcfs.power_ion_total + ...
        global_quantities.inside_lcfs.power_neutral_total + global_quantities.inside_lcfs.power_electrons;

        global_quantities.inside_vessel.power_ion_total     = 0;
        global_quantities.inside_vessel.power_neutral_total = 0;
        global_quantities.inside_vessel.power_electrons     = profiles_1d.electrons.power_inside(end);
        global_quantities.inside_vessel.power               = global_quantities.inside_vessel.power_ion_total  + ...
        global_quantities.inside_vessel.power_neutral_total + global_quantities.inside_vessel.power_electrons;

        global_quantities.time       = profil0d.temps(k);
	process.global_quantities{k} = global_quantities;             
end
% fill structure
radiation.process{3}      = process;


%  recombination 	11 	Emission from recombination
identifier.name         = 'recombination';
identifier.index        = 11;
identifier.description  = 'Emission from recombination';
process.identifier      = identifier;
% loop on time
for k=1:length(profil0d.temps)
        %profiles_1d                           = profiles_1d_model
        profiles_1d.grid                       = mapprofiles1d_grid_imas(profil0d,k,profiles_1d.grid);
        profiles_1d.electrons.emissivity       = - pradhdt(k,:);
        profiles_1d.electrons.power_inside     = cumtrapz(profil0d.xli,profil0d.vpr(k,:) .* profiles_1d.electrons.emissivity,2);
        profiles_1d.emissivity_ion_total       = zz;
        profiles_1d.power_inside_ion_total     = zz;
        profiles_1d.emissivity_neutral_total   = zz;
        profiles_1d.power_inside_neutral_total = zz;
        %ion: {[1x1 struct]}
        %neutral: {[1x1 struct]}
        profiles_1d.time = profil0d.temps(k);
        process.profiles_1d{k} = profiles_1d;

        global_quantities.inside_lcfs.power_ion_total     = 0;
        global_quantities.inside_lcfs.power_neutral_total = 0;
        global_quantities.inside_lcfs.power_electrons     = profiles_1d.electrons.power_inside(end);
        global_quantities.inside_lcfs.power               = global_quantities.inside_lcfs.power_ion_total + ...
        global_quantities.inside_lcfs.power_neutral_total + global_quantities.inside_lcfs.power_electrons;

        global_quantities.inside_vessel.power_ion_total     = 0;
        global_quantities.inside_vessel.power_neutral_total = 0;
        global_quantities.inside_vessel.power_electrons     = profiles_1d.electrons.power_inside(end);
        global_quantities.inside_vessel.power               = global_quantities.inside_vessel.power_ion_total  + ...
        global_quantities.inside_vessel.power_neutral_total + global_quantities.inside_vessel.power_electrons;

        global_quantities.time       = profil0d.temps(k);
	process.global_quantities{k} = global_quantities;             
end
% fill structure
radiation.process{4}      = process;
