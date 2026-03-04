function data = metis4muscle3_read(actor,direction,overlay_on_off,data)

% default state
if (nargin < 2) || isempty(direction)
    direction = 'input';
end
if nargin < 4
    data = [];
end

% mdoel
msg_model.timestamp = [];
msg_model.next_timestamp = [];
msg_model.data_raw = [];
msg_model.data     = [];
msg_model.settings = [];
%
% look if ports are connected
name = fieldnames(actor.io.list);
for k=1:length(name)
    if ~isfield(data,name{k})
        data.(name{k}) = msg_model;
        % generate empty data structure
        for l=1:length(actor.io.list.(name{k}).list_data_name)
            data.(name{k}).data.(actor.io.list.(name{k}).list_data_name{l}) = [];
        end
    end
    if actor.io.status.(name{k})
        switch actor.io.list.(name{k}).direction
            case direction
                fprintf('reading data on port %s\n',name{k});
                warning off
                if ~overlay_on_off
                    msg =  struct(actor.instance.receive_with_settings(name{k}));
                else
                    msg =  struct(actor.instance.receive(name{k}));
                end
                warning on
                if isfield(msg,'timestamp')
                    data.(name{k}).timestamp = msg.timestamp;
                end
                if isfield(msg,'next_timestamp')
                    data.(name{k}).next_timestamp = msg.next_timestamp;
                end
                if isfield(msg,'data')
                    data.(name{k}).data_raw = msg.data;
                    % process data
                    switch class(data.(name{k}).data_raw)
                        case 'py.NoneType'
                            % nothing to do in Matlab
                        otherwise
                            for l=1:length(actor.io.list.(name{k}).list_data_name)
                                if iscell(data.(name{k}).data_raw)
                                    value =  data.(name{k}).data_raw{l};
                                elseif (length(actor.io.list.(name{k}).list_data_name) == length(data.(name{k}).data_raw)) && ...
                                        (length(actor.io.list.(name{k}).list_data_name) >1)
                                    value =  data.(name{k}).data_raw(l);
                                else
                                    value =  data.(name{k}).data_raw;
                                end
                                switch actor.io.list.(name{k}).list_type{l}
                                    case {'empty'}
                                        data.(name{k}).data.(actor.io.list.(name{k}).list_data_name{l}) = [];
                                    case {'int','float'}
                                        data.(name{k}).data.(actor.io.list.(name{k}).list_data_name{l}) = double(value);
                                    case {'str','char','string'}
                                        data.(name{k}).data.(actor.io.list.(name{k}).list_data_name{l}) = char(value);
                                    case 'ids'
                                        % strange test due to Python object
                                        if ~isempty(value) && (length(value) > 1)
                                            data.(name{k}).data.(actor.io.list.(name{k}).list_data_name{l}) =  imas_deserialize(uint8(value), ...
                                                strtrim(actor.io.list.(name{k}).list_data_name{l}));
                                        else
                                            fprintf('Empty data for IDS %s\n',actor.io.list.(name{k}).list_data_name{l});
                                        end
                                    case 'matlab_serialized_structure'
                                        data.(name{k}).data.(actor.io.list.(name{k}).list_data_name{l}) = get_matlab_serialized_structure(value);
                                    otherwise
                                        warning(sprintf('not yet implemented type %s !',actor.io.list.(name{k}).list_type{l}));
                                        data.(name{k}).data.(actor.io.list.(name{k}).list_data_name{l}) = value;
                                end
                            end
                    end
                end
                if isfield(msg,'settings')
                    warning off
                    data.(name{k}).settings = struct(msg.settings);
                    warning on
                    % need to preprocesse settings
                    if isempty(fieldnames(data.(name{k}).settings))
                        % set to empty is simple to manage in Matlab
                        data.(name{k}).settings =[];
                    end
                end
        end
    end
end
