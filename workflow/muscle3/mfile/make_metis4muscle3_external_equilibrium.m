function data_out = make_metis4muscle3_external_equilibrium(data_in,equi_ids,time_equi,dt,psioffset,nbslices)

% get cocos depending on DD version
imas_version = getappdata(0,'IMAS_VERSION');
% only official IMAS COCOS are taken into account
switch imas_version
    case 4
        cocos = 17;
    case 3
        cocos = 11;
    otherwise
        error('unhandled IMAS version');
end

% correction of poloidal flux offset
switch cocos
    case 11
        sign_cocos = -1;
    otherwise
        sign_cocos =  1;
end
% create_external_equilibrium_from_imas is in COCOS 11 !
equi_ids.time_slice{end}.profiles_1d.psi = -sign_cocos .* (equi_ids.time_slice{end}.profiles_1d.psi - sign_cocos .* 2 .* pi .* psioffset);
if ~isempty(equi_ids.time_slice{end}.profiles_2d)
    for k=1:length(equi_ids.time_slice{end}.profiles_2d)
        if ~isempty(equi_ids.time_slice{end}.profiles_2d{k}.psi)
            equi_ids.time_slice{end}.profiles_2d{k}.psi = -sign_cocos .* (equi_ids.time_slice{end}.profiles_2d{k}.psi - sign_cocos .* 2 .* pi .* psioffset);
        end
    end
end
equi_ids.time_slice{end}.global_quantities.psi_axis = -sign_cocos .* (equi_ids.time_slice{end}.global_quantities.psi_axis - sign_cocos .* 2 .* pi .* psioffset);
equi_ids.time_slice{end}.global_quantities.psi_boundary = -sign_cocos .* (equi_ids.time_slice{end}.global_quantities.psi_boundary - sign_cocos .* 2 .* pi .* psioffset);
if ~isempty(equi_ids.time_slice{end}.global_quantities.psi_external_average)
    equi_ids.time_slice{end}.global_quantities.psi_external_average = -sign_cocos .* (equi_ids.time_slice{end}.global_quantities.psi_external_average - sign_cocos .* 2 .* pi .* psioffset);
end
if ~isempty(equi_ids.time_slice{end}.boundary.psi)
    equi_ids.time_slice{end}.boundary.psi = -sign_cocos .* (equi_ids.time_slice{end}.boundary.psi - sign_cocos .* 2 .* pi .* psioffset);
end
% there is other substructure with psi data, but in any case it will no be
% used inside METIS.

% backward compatibiliy with DDv 3:
equi_ids.time_slice{end}.profiles_1d.j_tor = equi_ids.time_slice{end}.profiles_1d.j_phi;
if ~isempty(equi_ids.time_slice{end}.profiles_2d)
    for k=1:length(equi_ids.time_slice{end}.profiles_2d)
        if ~isempty(equi_ids.time_slice{end}.profiles_2d{k}.j_phi)
            equi_ids.time_slice{end}.profiles_2d{k}.j_tor = equi_ids.time_slice{end}.profiles_2d{k}.j_phi;
        end
    end
end
% problem with boundary and x point description
equi_ids.time_slice{end}.boundary.x_point = {};
if ~isempty(equi_ids.time_slice{end}.contour_tree.node)
    for k=1:length(equi_ids.time_slice{end}.contour_tree.node)
        switch equi_ids.time_slice{end}.contour_tree.node{k}.critical_type
            case 1
                % this a X_point
                equi_ids.time_slice{end}.boundary.x_point{end+1}.r = equi_ids.time_slice{end}.contour_tree.node{k}.r;
                equi_ids.time_slice{end}.boundary.x_point{end+1}.z = equi_ids.time_slice{end}.contour_tree.node{k}.z;
            case 0
                % this is the magnetic axis
            case 2
                % this is the LCFS
            otherwise
                error('something new has been created; snow flake ?')
        end
    end
end

% update time in equi ids
equi_ids.time(end) = time_equi;
equi_ids.time_slice{end}.time = time_equi;

% here maybe some quantities should be changed of sign due for
% compatibility with JT-60SA simulator: this is purely internal to METIS
% case initialisation
if isempty(data_in)
    data_out = equi_ids;
    % set 3 pseudo time slices on initialisation
    data_out.time = cat(1,time_equi - 2 .* dt,time_equi - dt,time_equi);
    data_out.time_slice   = equi_ids.time_slice(end); % enforce only one time slice !
    data_out.time_slice(2) = equi_ids.time_slice(end);
    data_out.time_slice(3) = equi_ids.time_slice(end);
    data_out.time_slice{1}.time = time_equi - 2 .* dt;
    data_out.time_slice{2}.time = time_equi- dt;
    data_out.time_slice{3}.time = time_equi;   
    % B0
    data_out.vacuum_toroidal_field.b0 = data_out.vacuum_toroidal_field.b0(end) * ones(3,1);
    %
elseif data_in.time(end) == time_equi
    % convergence
    data_out = data_in;
    data_out.time(end) = time_equi;
    data_out.time_slice(end) = equi_ids.time_slice(end);
    data_out.time_slice{end}.time = time_equi;   
    % B0
    data_out.vacuum_toroidal_field.b0(end) = equi_ids.vacuum_toroidal_field.b0(end) .* ...
         equi_ids.vacuum_toroidal_field.r0(end) ./ data_out.vacuum_toroidal_field.r0(end);
    %
elseif data_in.time(end) < time_equi
    % evolution
    data_out = data_in;
    data_out.time(end+1) = time_equi;
    data_out.time_slice(end+1) = equi_ids.time_slice(end);
    data_out.time_slice{end+1}.time = time_equi;   
    % B0
    data_out.vacuum_toroidal_field.b0(end+1) = equi_ids.vacuum_toroidal_field.b0(end) .* ...
              equi_ids.vacuum_toroidal_field.r0(end) ./ data_out.vacuum_toroidal_field.r0(end);
    % limite number of time slices
    if length(data_out.time) > nbslices
        indkeep = (length(data_out.time) - nbslices +1 ):length(data_out.time);
        data_out.time = data_out.time(indkeep);
        data_out.vacuum_toroidal_field.b0 = data_out.vacuum_toroidal_field.b0(indkeep);
        data_out.time_slice = data_out.time_slice(indkeep);
    end
else
    % rewind
    if any(data_in.time) == time_equi
        % replace
        indr = find(data_in.time == time_equi,1);
        data_out = data_in;
        %data_out.time(indr) = time_equi;
        %data_out.time_slice{end+1}.time = time_equi;   
        data_out.time_slice(indr) = equi_ids.time_slice(end);
        data_out.vacuum_toroidal_field.b0(indr) = equi_ids.vacuum_toroidal_field.b0(end) .* ...
                  equi_ids.vacuum_toroidal_field.r0(end) ./ data_out.vacuum_toroidal_field.r0(end);
    else
        % insert
        indm = find(data_in.time < time_equi);
        indp = find(data_in.time > time_equi);
        %
        data_out.time = cat(1,data_in.time(indm),equi_time,data_in.time(indp));
        data_out.time_slice = data_in.time_slice(indm);
        data_out.time_slice(end+1) = equi_ids.time_slice(end);
        indadd = length(data_out.time_slice) + indp -indp(1) + 1;
        data_out.time_slice(indadd) = data_in.time_slice(indp);
        data_out.vacuum_toroidal_field = data_in.vacuum_toroidal_field;
        data_out.vacuum_toroidal_field.b0 = cat(1,data_in.vacuum_toroidal_field.b0(indm), ...
            equi_ids.vacuum_toroidal_field.b0(end) .* ...
            equi_ids.vacuum_toroidal_field.r0(end) ./  ...
            data_out.vacuum_toroidal_field.r0(end), ...
            data_in.vacuum_toroidal_field.b0(indp));
        % adding missing field
        noms = fieldnames(data_in);
        for k=1:length(noms)
            if ~isfield(data_out,noms{k})
                data_out.(noms{k}) = data_in.(noms{k});
            end
        end
        % limite number of time slices
        if length(data_out.time) > nbslices
            indkeep = (length(data_out.time) - nbslices +1 ):length(data_out.time);
            data_out.time = data_out.time(indkeep);
            data_out.vacuum_toroidal_field.b0 = data_out.vacuum_toroidal_field.b0(indkeep);
            data_out.time_slice = data_out.time_slice(indkeep);
        end
    end
end




