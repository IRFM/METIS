function  [input,LCFS_nbpts] = get_metis_input_from_muscle3(metis_data,actor,clock_in,used_ports)

% initialisation
input = [];
% select input data source
is_ps = false;
for k=1:length(used_ports)
    if actor.io.status.(sprintf('pulse_schedule_%s',used_ports{k}))
        is_ps = true;
    end
end
if is_ps
    % get ids data
    [timestamp_ps,ps_ids,ps_used_port] = get_ids_data_from_ports(actor.io.data,actor.io.list,'pulse_schedule',used_ports);
    fprintf('METIS4MUSCLE3: using pulse_schedule IDS data as METIS input data, reading from port %s\n',ps_used_port);
    % input data will be get from pulse_schedule IDS
    [input.option,input.time,input.cons1t,input.geo1t,input.sepa1t] = ids2metis1t(ps_ids);
    % when data are coming from pulse_schedule, the two following fields
    % are not valid or not available
    input.cons1t.flux = 0;
    input.cons1t.vloop = NaN; % this is no a real field of cons1t
    % idem for geo1t
    input.geo1t.vp   = NaN;
    input.geo1t.sp   = NaN;
    input.geo1t.sext = NaN;
    % update datation using clock information
    if input.time ~= clock_in.t_cur
        fprintf('Information: missmatch between MUSCLE3 message time_stamp (%g s) and pulse_schedule internal time(%g s): forced to use MUSCLE3 message time_stamp.\n', ...
            clock_in.t_cur,input.time);
    end
    input.time = clock_in.t_cur;
    input.cons1t.temps =  clock_in.t_cur;
elseif ~isempty(metis_data.loaded_data.post)
    disp('METIS4MUSCLE3: using data from loaded METIS simulation file as METIS input data');
    % z0dinput from METIS simulation file
    z0dinput = metis_data.loaded_data.post.z0dinput;
    % securite
    if isfield(z0dinput.exp0d,'Rsepa') && isempty(z0dinput.exp0d.Rsepa)
        z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
        z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
    end
    if isfield(z0dinput.exp0d,'XDURx') && isempty(z0dinput.exp0d.XDURx)
        z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURx');
        z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURt');
        z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURv');
    end
    % extract one time data
    if clock_in.t_cur <= max(z0dinput.cons.temps)
        input.cons1t  = zerod_get_one_time(z0dinput.cons,z0dinput.cons.temps,clock_in.t_cur,'linear','extrap');
        input.geo1t   = zerod_get_one_time(z0dinput.geo,z0dinput.cons.temps,clock_in.t_cur,'linear','extrap');
        input.exp0d1t = zerod_get_one_time(z0dinput.exp0d,z0dinput.cons.temps,clock_in.t_cur,'linear','extrap');
    else
        input.cons1t  = zerod_get1t(z0dinput.cons,length(z0dinput.cons.temps));
        input.geo1t   = zerod_get1t(z0dinput.geo,length(z0dinput.cons.temps));    
        input.exp0d1t = zerod_get1t(z0dinput.exp0d,length(z0dinput.cons.temps));
    end
    % LCFS 
    input.sepa1t = [];
    % securite
    if isfield(input.exp0d1t,'Rsepa') && isfield(input.exp0d1t,'Zsepa') && ...
            ~isempty(input.exp0d1t.Rsepa) && ~isempty(input.exp0d1t.Zsepa) && ...
            all(size(input.exp0d1t.Rsepa) == size(input.exp0d1t.Zsepa))
        input.sepa1t.Rsepa = input.exp0d1t.Rsepa;
        input.sepa1t.Zsepa = input.exp0d1t.Zsepa;
    else
        input.sepa1t = [];
    end
    % update datation using clock information
    input.time = clock_in.t_cur;
    input.cons1t.temps =  clock_in.t_cur;
    % option
    input.option = z0dinput.option;

else
    error('No input data available for METIS');
end

% overload option with parameters defined in settings
model = metis4imas(1);
if isfield(actor.settings.value,'metis_option') && ~isempty(actor.settings.value.metis_option)
    noms = fieldnames(actor.settings.value.metis_option);
    for k=1:length(noms)
        metis_name = noms{k};
        if ~isfield(input.option,noms{k})
            % search for label
            nomsl = fieldnames(model.label);
            for l=1:length(nomsl)
                switch  model.label.(nomsl{l})
                    case noms{k}
                        metis_name = nomsl{l};
                        break
                end
            end
        end
        if isfield(input.option,metis_name) && ~isempty(actor.settings.value.metis_option.(noms{k}))
            input.option.(metis_name) = actor.settings.value.metis_option.(noms{k});
        elseif isempty(actor.settings.value.metis_option.(noms{k}))
            fprintf('METIS option parameter %s defined in settings is empty (maybe the type in settings not match the type of METIS option).\n',noms{k});
        else
            fprintf('METIS option parameter %s defined in settings do not match any existing METIS parameter.\n',noms{k});
        end
    end
end

% LCFS_nbpts
if isfield(input,'sepa1t')
    if isfield(input.sepa1t,'Rsepa')
        LCFS_nbpts = size(input.sepa1t.Rsepa,2);
    else
        LCFS_nbpts = Nan;
    end
else
    LCFS_nbpts = Nan;
end


