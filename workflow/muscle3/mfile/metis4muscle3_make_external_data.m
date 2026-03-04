function [data,mode_change,LCFS_nbpts] = metis4muscle3_make_external_data(internal_state)

persistent count_invalid_equi_ids
if isempty(count_invalid_equi_ids)
    count_invalid_equi_ids = 0;
    %dbstop if all error
end

% output initialisation
data = internal_state.external_data.data;
if ~isfield(data,'coordinates_change')
    data.coordinates_change.temps        = [];
    data.coordinates_change.rho_tor_norm = [];
    data.coordinates_change.x            = [];
end

% to get more compact code
used_ports = internal_state.workflow.used_ports;

% search for data use in coordinate change
mode_change.equilibrium = false;
mode_change.z0dstruct = false;
mode_change.metis = false;
mode_change.approximation = false;

% if external equilibrium is available
[timestamp_equi,equi_ids,equi_used_port] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,'equilibrium',used_ports);
% get equilibrium status
if ~isempty(equi_ids)
    % for testing
    % if rand(1) >= 0.5
    %     equi_ids.code.output_flag = 3;
    % else
    %     equi_ids.code.output_flag = 0;
    % end
    if isempty(equi_ids.code.output_flag)
        equi_is_valid = true;
        count_invalid_equi_ids = 0;
    elseif equi_ids.code.output_flag == 0
        equi_is_valid = true;
        count_invalid_equi_ids = 0;
    else
        fprintf('Equilibrium IDS is flag %d and will be not used\n',equi_ids.code.output_flag);
        equi_is_valid = false;
        count_invalid_equi_ids = count_invalid_equi_ids + 1;
        if count_invalid_equi_ids > 1
            fprintf('No valid equilibrium IDS has been received since %d call of METIS actor\n',count_invalid_equi_ids);
        end
    end
else
    equi_is_valid = true;
    count_invalid_equi_ids = 0;
end
if ~isempty(equi_ids) && equi_is_valid
    mode_change.equilibrium = true; 
    % add new time to coordinate change
    % select time
    time_loc = get_time_external(internal_state.actor.settings.value.metis,equi_ids,timestamp_equi);
    % data at one time
    coordinates_change_1t.temps        = time_loc;
    rho_tor_norm = equi_ids.time_slice{end}.profiles_1d.rho_tor_norm(:)';
    x_loc        = equi_ids.time_slice{end}.profiles_1d.r_outboard(:)' -  ...
                   equi_ids.time_slice{end}.profiles_1d.r_inboard(:)';
    x_loc        = (x_loc - min(x_loc)) ./ (max(x_loc) - min(x_loc));
    % resample on nbpts1d
    coordinates_change_1t.x            = linspace(0,1,internal_state.actor.settings.value.metis.nbpts1d);
    coordinates_change_1t.rho_tor_norm = interp1_imas(x_loc,rho_tor_norm,coordinates_change_1t.x,'pchip','extrap');                                    
    % manage timed data accumulation  
    data.coordinates_change = manage_external_data_accumulation(data.coordinates_change,coordinates_change_1t);
    disp('METIS4MUSCLE3: using input equilibrium IDS for coordinate change applied to external data');
    
elseif isfield(internal_state.metis.output,'z0dstruct') &&  ...
        ~isempty(internal_state.metis.output.z0dstruct)
    % if z0dstruct data available
    mode_change.z0dstruct = true;
    data.coordinates_change.temps         = internal_state.metis.output.z0dstruct.z0dinput.cons.temps;
    if isfield(internal_state.metis.output.z0dstruct,'profil')
        rmx   = internal_state.metis.output.z0dstruct.profil.rmx;
        x_loc = internal_state.metis.output.z0dstruct.profil.xli;
    else
        rmx   = internal_state.metis.output.z0dstruct.profil0d.rmx;
        x_loc = internal_state.metis.output.z0dstruct.profil0d.xli;
    end
    rho_tor_norm  = (rmx - rmx(:,1) * ones(1,size(rmx,2))) ./ ...
                                            ((rmx(:,end) - rmx(:,1)) * ones(1,size(rmx,2)));
    % resample on nbpts1d
    data.coordinates_change.x            = linspace(0,1,internal_state.actor.settings.value.metis.nbpts1d);
    data.coordinates_change.rho_tor_norm = interp1_imas(x_loc.',rho_tor_norm.',data.coordinates_change.x.','pchip','extrap').';                                
    disp('METIS4MUSCLE3: using actor METIS internal data for coordinate change applied to external data');

elseif isfield(internal_state.metis.loaded_data,'post') && ...
        ~isempty(internal_state.metis.loaded_data.post)
    % if METIS data available
    mode_change.metis = true;
    if isfield(internal_state.metis.loaded_data.post.profil0d,'temps')
        data.coordinates_change.temps = internal_state.metis.loaded_data.post.profil0d.temps;
    else
        data.coordinates_change.temps = internal_state.metis.loaded_data.post.z0dinput.cons.temps;       
    end
    rmx = internal_state.metis.loaded_data.post.profil0d.rmx;
    rho_tor_norm  = (rmx - rmx(:,1) * ones(1,size(rmx,2))) ./ ...
                                            ((rmx(:,end) - rmx(:,1)) * ones(1,size(rmx,2)));
    x_loc =  internal_state.metis.loaded_data.post.profil0d.xli;   
    % resample on nbpts1d
    data.coordinates_change.x            = linspace(0,1,internal_state.actor.settings.value.metis.nbpts1d);
    data.coordinates_change.rho_tor_norm = interp1_imas(x_loc.',rho_tor_norm.',data.coordinates_change.x.','pchip','extrap').';
    disp('METIS4MUSCLE3: using METIS loaded data for coordinate change applied to external data');
    
else
    % otherwise make some approximation
    mode_change.approximation = true;
    data.coordinates_change.temps = cat(1,internal_state.workflow.clock.t_cur - internal_state.workflow.clock.dt, ...
                                          internal_state.workflow.clock.t_cur,internal_state.workflow.clock.t_next);
    data.coordinates_change.rho_tor_norm  = ones(length(data.coordinates_change.temps),1) * linspace(0,1,internal_state.actor.settings.value.metis.nbpts1d);
    data.coordinates_change.x             = linspace(0,1,internal_state.actor.settings.value.metis.nbpts1d);
    disp('METIS4MUSCLE3: linear coordinate change applied to external data');
    
end

% get sign for current and momentum sources
if isempty(equi_ids) || ~equi_is_valid
    [timestamp_cp,cp_ids,cp_used_port] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,'plasma_profiles',used_ports);
    if isempty(cp_ids)
        [timestamp_cp,cp_ids,cp_used_port] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,'core_profiles',used_ports);
    end
end

if  ~isempty(equi_ids) && equi_is_valid
    % get b0 sign
    sign_b0 = sign(equi_ids.vacuum_toroidal_field.b0);
    sign_ip = sign(equi_ids.time_slice{end}.global_quantities.ip);
    sign_current_source  = sign_ip;
    % in METIS the rotation direction is count positive in plamsa current
    % direction. 
    sign_momentum_source = sign_ip;
    % test against METIS setting
    if isfield(internal_state.metis.output,'z0dstruct') &&  ...
            ~isempty(internal_state.metis.output.z0dstruct)
        if internal_state.metis.output.z0dstruct.z0dinput.option.signe ~= ...
                (sign_b0*sign_ip)
            fprintf('METIS4MUSCLE3: incoherent setting between METIS and inputed IDSs;\n option.sign should be %d\n',(sign_b0*sign_ip));
        end
    end    
elseif ~isempty(cp_ids)
    % get b0 sign
    sign_b0 = sign(cp_ids.vacuum_toroidal_field.b0);
    sign_ip = sign(cp_ids.global_quantities.ip);
    sign_current_source  = sign_ip;
    % in METIS the rotation direction is count positive in plamsa current
    % direction. 
    sign_momentum_source = sign_ip;
    % test against METIS setting
    if isfield(internal_state.metis.output,'z0dstruct') &&  ...
            ~isempty(internal_state.metis.output.z0dstruct)
        if internal_state.metis.output.z0dstruct.z0dinput.option.signe ~= ...
                (sign_b0*sign_ip)
            fprintf('METIS4MUSCLE3: incoherent setting between METIS and inputed IDSs;\n option.sign should be %d\n',(sign_b0*sign_ip));
        end
    end        
else
    % use METIS setting
    if isfield(internal_state.metis.output,'z0dstruct') &&  ...
            ~isempty(internal_state.metis.output.z0dstruct)
        sign_metis = internal_state.metis.output.z0dstruct.z0dinput.option.signe;
        orientation = internal_state.metis.output.z0dstruct.z0dinput.option.orientation;
        if isfield(internal_state.metis.output.z0dstruct.z0dinput.option,'COCOS')
            cocos = internal_state.metis.output.z0dstruct.z0dinput.option.COCOS;
        else
            %assume native COCOS of IMAS DDv4
            cocos = 17;
        end
    elseif isfield(internal_state.metis.loaded_data,'post') &&  ...
            ~isempty(internal_state.metis.loaded_data.post)
        sign_metis = internal_state.metis.loaded_data.post.z0dinput.option.signe;
        orientation = internal_state.metis.loaded_data.post.z0dinput.option.orientation;
        if isfield(internal_state.metis.loaded_data.post.z0dinput.option,'COCOS')
            cocos = internal_state.metis.loaded_data.post.z0dinput.option.COCOS;
        else
            % assume native COCOS of IMAS DDv4
            cocos = 17;
        end  
    else
        % assume native COCOS of IMAS DDv4 and ITER configuration
        sign_metis = 1;
        orientation = -1;
        cocos = 17;
    end
    %
    switch cocos
        case {1,11,2,12,7,17,8,18}
            srphiz = 1;
        otherwise
            srphiz = -1;
    end
    %
    sign_current_source  = orientation .*  sign_metis .* srphiz;
    sign_momentum_source = orientation .*  sign_metis .* srphiz;
    % sign_b0 = orientation .* srphiz;

end

% process each case separately as each one has different importation code
% 1- external equilibrium and/or current diffusion
if (internal_state.external_data.state.equilibrium) || ...
        (internal_state.external_data.state.current_diffusion)
    %
    ids_name = internal_state.actor.settings.value.metis_external_data.equilibrium.ids_name;
    [timestamp,data_ids,equi_used_port] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if isempty(data_ids)
        ids_name = internal_state.actor.settings.value.metis_external_data.current_diffusion.ids_name;
        [timestamp,data_ids,equi_used_port] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    end
  
    if ~isempty(data_ids) && equi_is_valid
        %
        % get boundary length
        LCFS_nbpts = length(data_ids.time_slice{end}.boundary.outline.r);
        %
        if (internal_state.external_data.state.equilibrium) && isfield(data,'equilibrium') && ~isempty(data.equilibrium)
            previous_data = data.equilibrium;
        elseif (internal_state.external_data.state.current_diffusion) && isfield(data,'equilibrium') && ~isempty(data.current_diffusion)
            previous_data = data.current_diffusion;
        else
            previous_data = [];
        end
        %
        time_equi = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        %
        new_data = make_metis4muscle3_external_equilibrium(previous_data,data_ids,time_equi, ...
            internal_state.workflow.clock.dt, ...
            internal_state.actor.settings.value.metis.psioffset, ...
            internal_state.actor.settings.value.metis.numberofequilibriumtimeslices);
        if internal_state.external_data.state.equilibrium
            data.equilibrium = new_data;
        end
        if internal_state.external_data.state.current_diffusion
            data.current_diffusion = new_data;
        end
    else
        LCFS_nbpts = NaN;
    end
else
    LCFS_nbpts = NaN;
end

% 3 - electron temperature
if internal_state.external_data.state.electron_temperature
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.electron_temperature.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % select time
        data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        % coordinates change
        if ~isempty(data_ids.profiles_1d{end}.electrons.temperature)
            [data_1t.x,data_1t.te] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.profiles_1d{end}.electrons.temperature(:)', ...
                data.coordinates_change);
            % on initialisation
            if ~isfield(data,'electron_temperature') || isempty(data.electron_temperature)
                data.electron_temperature = manage_external_data_init(data_1t,internal_state.workflow.clock);
            else
                data.electron_temperature = manage_external_data_accumulation(data.electron_temperature,data_1t);
            end
        else
            if isfield(data,'electron_temperature')
                data = rmfield(data,'electron_temperature');
            end
            warning('No external data for Te !');
        end
    end
end

% 4 - electron density
if internal_state.external_data.state.electron_density
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.electron_density.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % select time
        data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        % coordinates change
        if ~isempty(data_ids.profiles_1d{end}.electrons.density)
            [data_1t.x,data_1t.ne] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.profiles_1d{end}.electrons.density(:)', ...
                data.coordinates_change);
            % on initialisation
            if ~isfield(data,'electron_density') || isempty(data.electron_density)
                data.electron_density = manage_external_data_init(data_1t,internal_state.workflow.clock);
            else
                data.electron_density = manage_external_data_accumulation(data.electron_density,data_1t);
            end
        else
            if isfield(data,'electron_density')
                data = rmfield(data,'electron_density');
            end
            warning('No external data for ne !');
        end
    end
end

% 5 - ion temperature (t_i_average)
if internal_state.external_data.state.ion_temperature
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.ion_temperature.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % select time
        data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        % coordinates change
        if ~isempty(data_ids.profiles_1d{end}.t_i_average)
            [data_1t.x,data_1t.ti] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.profiles_1d{end}.t_i_average(:)', ...
                data.coordinates_change);
            % on initialisation
            if ~isfield(data,'ion_temperature') || isempty(data.ion_temperature)
                data.ion_temperature = manage_external_data_init(data_1t,internal_state.workflow.clock);
            else
                data.ion_temperature = manage_external_data_accumulation(data.ion_temperature,data_1t);
            end
        else
            if isfield(data,'ion_temperature')
                data = rmfield(data,'ion_temperature');
            end
            warning('No external data for Ti !');
        end
    end
end

% 6 - Zeff
if internal_state.external_data.state.charge_effective
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.charge_effective.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % select time
        data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        % coordinates change
        if ~isempty(data_ids.profiles_1d{end}.zeff)
            [data_1t.x,data_1t.zeff] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.profiles_1d{end}.zeff(:)', ...
                data.coordinates_change);
            % on initialisation
            if ~isfield(data,'charge_effective') || isempty(data.charge_effective)
                data.charge_effective = manage_external_data_init(data_1t,internal_state.workflow.clock);
            else
                data.charge_effective = manage_external_data_accumulation(data.charge_effective,data_1t);
            end
        else
            if isfield(data,'charge_effective')
                data = rmfield(data,'charge_effective');
            end
            warning('No external data for Zeff !');
        end
    end
end

% 7 - LHCD
if internal_state.external_data.state.LHCD
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.LHCD.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % select time
        data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        % index of the source
        index_source = get_source_index(data_ids,{'LH','LHCD'});
        if ~isempty(index_source)
            % coordinates change
            [data_1t.x,data_1t.jlh] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.j_parallel(:)', ...
                data.coordinates_change);
            %
            [data_1t.x,data_1t.plh] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.electrons.energy(:)', ...
                data.coordinates_change);
            % set correct sign
            data_1t.jlh = sign_current_source .* data_1t.jlh;
            % on initialisation
            if ~isfield(data,'LHCD') || isempty(data.LHCD)
                data.LHCD = manage_external_data_init(data_1t,internal_state.workflow.clock);
            else
                data.LHCD = manage_external_data_accumulation(data.LHCD,data_1t);
            end
        else
            if isfield(data,'LHCD')
                data = rmfield(data,'LHCD');
            end
            warning('No external data for LHCD !');
        end
    end
end

% 8 - ECRH
if internal_state.external_data.state.ECRH
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.ECRH.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % select time
        data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        % index of the source
        index_source = get_source_index(data_ids,{'EC','ECRH','ECCD'});
        if ~isempty(index_source)
            % coordinates change
            [data_1t.x,data_1t.jeccd] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.j_parallel(:)', ...
                data.coordinates_change);
            %
            [data_1t.x,data_1t.peccd] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.electrons.energy(:)', ...
                data.coordinates_change);
            % set correct sign
            data_1t.jeccd = sign_current_source .* data_1t.jeccd;
            % on initialisation
            if ~isfield(data,'ECRH') || isempty(data.ECRH)
                data.ECRH = manage_external_data_init(data_1t,internal_state.workflow.clock);
            else
                data.ECRH = manage_external_data_accumulation(data.ECRH,data_1t);
            end
        else
            if isfield(data,'ECRH')
                data = rmfield(data,'ECRH');
            end
            warning('No external data for ECRH !');
        end
    end
end

% 9 - ICRH
if internal_state.external_data.state.ICRH
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.ICRH.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % select time
        data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        % index of the source
        index_source = get_source_index(data_ids,{'IC','ICRH','FW','FWCD'});
        if ~isempty(index_source)
            % coordinates change
            [data_1t.x,data_1t.jfwcd] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.j_parallel(:)', ...
                data.coordinates_change);
            %
            [data_1t.x,pel] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.electrons.energy(:)', ...
                data.coordinates_change);
            %
            [data_1t.x,data_1t.pion] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.total_ion_energy(:)', ...
                data.coordinates_change);
            % manage heating mode for METIS internal model
            if all(data_1t.pion==0) && any(pel~=0)
                data_1t.pfw = pel;
                data_1t.pel = zeros(size(pel));
            else
                data_1t.pfw = zeros(size(pel));
                data_1t.pel = pel;
            end
            % set correct sign
            data_1t.jfwcd = sign_current_source .* data_1t.jfwcd;
            % on initialisation
            if ~isfield(data,'ICRH') || isempty(data.ICRH)
                data.ICRH = manage_external_data_init(data_1t,internal_state.workflow.clock);
            else
                data.ICRH = manage_external_data_accumulation(data.ICRH,data_1t);
            end
        else
            if isfield(data,'ICRH')
                data = rmfield(data,'ICRH');
            end
            warning('No external data for ICRH !');
        end
    end
end

% 10 - NBI
if internal_state.external_data.state.NBI
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.NBI.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % select time
        data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        % index of the source
        index_source = get_source_index(data_ids,{'NBI','NBICD'});
        if ~isempty(index_source)
            % coordinates change
            [data_1t.x,data_1t.jnbicd] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.j_parallel(:)', ...
                data.coordinates_change);
            %
            [data_1t.x,data_1t.pel] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.electrons.energy(:)', ...
                data.coordinates_change);
            %
            [data_1t.x,data_1t.pion] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.total_ion_energy(:)', ...
                data.coordinates_change);
            % set correct sign
            data_1t.jnbicd = sign_current_source .* data_1t.jnbicd;
            % on initialisation
            if ~isfield(data,'NBI') || isempty(data.NBI)
                data.NBI = manage_external_data_init(data_1t,internal_state.workflow.clock);
            else
                data.NBI = manage_external_data_accumulation(data.NBI,data_1t);
            end
        else
            if isfield(data,'NBI')
                data = rmfield(data,'NBI');
            end
            warning('No external data for NBI !');
        end       
    end
end

% 11 - PLINE
if internal_state.external_data.state.PLINE
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.NBI.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % select time
        data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
        % index of the source
        index_source = get_source_index(data_ids,{'LINE','LINERADIATION','LINE_RADIATION','PLINE','PRAD'});
        if ~isempty(index_source)
            
            % coordinates change
            [data_1t.x,data_1t.prad] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                data_ids.source{index_source}.profiles_1d{end}.electrons.energy(:)', ...
                data.coordinates_change);
            % on initialisation
            if ~isfield(data,'PLINE') || isempty(data.PLINE)
                data.PLINE = manage_external_data_init(data_1t,internal_state.workflow.clock);
            else
                data.PLINE = manage_external_data_accumulation(data.PLINE,data_1t);
            end
        else
            if isfield(data,'PLINE')
                data = rmfield(data,'PLINE');
            end
            warning('No external data for PLINE !');
        end        
    end
end

% 12 - RUNAWAY electrons
if internal_state.external_data.state.runaway
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.runaway.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % get other IDSs
        % test if external equilibrium is available
        [~,data_ids_eq] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,'equilibrium',used_ports);
        % test if external core_profiles/plasma_profiles is available
        [~,data_ids_cp] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,'plasma_profiles',used_ports);
        if ~isempty(data_ids_eq) &&  ~isempty(data_ids_cp)
            % select time
            data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
            % index of the source
            index_source = get_source_index(data_ids,{'RUNAWAY','RUNAWAYS'});
            if ~isempty(index_source)
                % coordinates change
                [data_1t.x,data_1t.jrun] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                    data_ids.source{index_source}.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                    data_ids.source{index_source}.profiles_1d{end}.j_parallel(:)', ...
                    data.coordinates_change);
                
                %
                [~,data_1t.johm] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                    data_ids_cp.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                    data_ids_cp.profiles_1d{end}.j_ohmic(:)', ...
                    data.coordinates_change);
                %
                [~,data_1t.rho_tor] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                    data_ids_eq.time_slice{end}.profiles_1d.rho_tor_norm(:)', ...
                    data_ids_eq.time_slice{end}.profiles_1d.rho_tor(:)', ...
                    data.coordinates_change);
                %
                [~,data_1t.spr] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                    data_ids_eq.time_slice{end}.profiles_1d.rho_tor_norm(:)', ...
                    data_ids_eq.time_slice{end}.profiles_1d.darea_drho_tor(:)', ...
                    data.coordinates_change);
                %
                % set correct sign
                data_1t.jrun = sign_current_source .* data_1t.jrun;
                data_1t.johm = sign_current_source .* data_1t.johm;
                %
                data_1t.irun =  trapz(data_1t.rho_tor(:)',data_1t.spr(:)' .* data_1t.jrun(:)',2);
                data_1t.iohm =  trapz(data_1t.rho_tor(:)',data_1t.spr(:)' .* data_1t.johm(:)',2);
                data_1t.ip   =  sign_current_source .* data_ids_cp.global_quantities.ip(end);
                % on initialisation
                if ~isfield(data,'runaway') || isempty(data.runaway)
                    data.runaway = manage_external_data_init(data_1t,internal_state.workflow.clock);
                else
                    data.runaway = manage_external_data_accumulation(data.runaway,data_1t);
                end
            else
                if isfield(data,'runaway')
                    data = rmfield(data,'runaway');
                end
                warning('No external data for runaway !');
            end
        end
    end
end

% 13 - momentum
if internal_state.external_data.state.momentum
    % test if data are realy available
    ids_name = internal_state.actor.settings.value.metis_external_data.momentum.ids_name;
    [timestamp,data_ids] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,ids_name,used_ports);
    % if there is some data on this port
    if ~isempty(data_ids)
        data_1t = [];
        % get other IDSs
        % test if external equilibrium is available
        [~,data_ids_eq] = get_ids_data_from_ports(internal_state.actor.io.data,internal_state.actor.io.list,'equilibrium',used_ports);
        if ~isempty(data_ids_eq)
            % select time
            data_1t.temps = get_time_external(internal_state.actor.settings.value.metis,data_ids,timestamp);
            % coordinates change
            if (isfield(data_ids.profiles_1d{end},'momentum_tor') && ~isempty(data_ids.profiles_1d{end}.momentum_tor)) || ...
                    (isfield(data_ids.profiles_1d{end},'momentum_phi') && ~isempty(data_ids.profiles_1d{end}.momentum_phi))
                if isfield(data_ids.profiles_1d{end},'momentum_tor')
                    [data_1t.x,data_1t.momentum_tor] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                        data_ids.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                        data_ids.profiles_1d{end}.momentum_tor(:)', ...
                        data.coordinates_change);
                else
                    [data_1t.x,data_1t.momentum_tor] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                        data_ids.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                        data_ids.profiles_1d{end}.momentum_phi(:)', ...
                        data.coordinates_change);
                end
                %
                [~,data_1t.ne] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                    data_ids.profiles_1d{end}.grid.rho_tor_norm(:)', ...
                    data_ids.profiles_1d{end}.electrons.density(:)', ...
                    data.coordinates_change);
                %
                [~,data_1t.r2i] = metis4muscle3_rhotornorm2xli(data_1t.temps, ...
                    data_ids_eq.time_slice{end}.profiles_1d.rho_tor_norm(:)', ...
                    data_ids_eq.time_slice{end}.profiles_1d.gm1(:)', ...
                    data.coordinates_change);
                % set correct sign
                data_1t.momentum_tor = sign_momentum_source .* data_1t.momentum_tor;
                %
                omega = data_1t.momentum_tor ./ (data_1t.ne ./ max(data_1t.ne(:))) .* data_1t.r2i;
                data_1t.omega = omega ./ max(eps,max(abs(omega(:))));
                
                % on initialisation
                if ~isfield(data,'momentum') || isempty(data.momentum)
                    data.momentum = manage_external_data_init(data_1t,internal_state.workflow.clock);
                else
                    data.momentum = manage_external_data_accumulation(data.momentum,data_1t);
                end
            else
                if isfield(data,'momentum')
                    data = rmfield(data,'momentum');
                end
                warning('No external data for momentum !');
            end
        end
    end
end



