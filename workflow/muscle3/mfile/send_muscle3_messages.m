function send_muscle3_messages(internal_state,synchronization_flush)

% test input 
if (nargin <= 1) || isempty(synchronization_flush)
    synchronization_flush = 1;
end

% time data
time_stamp      = internal_state.workflow.clock.t_cur; 
time_stamp_next = internal_state.workflow.clock.t_next; 

if ~isfinite(time_stamp_next)
   % sen message to muscle3 to indicate the end of the run
   time_stamp_next = py.None();
end

% look if ports are connected
noms = fieldnames(internal_state.actor.io.list);
for k=1:length(noms)
    % loop on connected port
    if internal_state.actor.io.status.(noms{k})
        % only ouputed data
        switch internal_state.actor.io.list.(noms{k}).direction
            case 'output'
                list = internal_state.actor.io.list.(noms{k}).list_type;
                data_message = {};
                for l=1:length(list)
                    % send message
                    fprintf('Making message for %s, ',internal_state.actor.io.list.(noms{k}).list_data_name{l});
                    switch list{l}
                        case 'float'
                            var = create_float_message(internal_state,internal_state.actor.io.list.(noms{k}).list_data_name{l});
                        case 'int'
                            var = create_int_message(internal_state,internal_state.actor.io.list.(noms{k}).list_data_name{l});
                        case 'str'
                            var = create_str_message(internal_state,internal_state.actor.io.list.(noms{k}).list_data_name{l});
                        case 'ids'
                            if isfield(internal_state.imas.ids_data,internal_state.actor.io.list.(noms{k}).list_data_name{l})
                                var = create_ids_message(internal_state,internal_state.actor.io.list.(noms{k}).list_data_name{l},l);
                            else
                                var = make_python_bytes([]);
                            end
                        case 'matlab_serialized_structure'
                            var = create_matlab_message(internal_state,internal_state.actor.io.list.(noms{k}).list_data_name{l});
                        case 'empty'
                            var = create_matlab_message(internal_state,internal_state.actor.io.list.(noms{k}).list_data_name{l});                            
                        otherwise
                            error('Data type not yet implemented for output port');
                    end
                    data_message{end+1} = var;
                end
                if length(data_message) > 1
                    msg = py.libmuscle.Message(time_stamp,time_stamp_next,data_message);
                else
                    msg = py.libmuscle.Message(time_stamp,time_stamp_next,var);
                end
                fprintf('and sending on port %s\n',noms{k});
                internal_state.actor.instance.send(noms{k}, msg);
                % problem with MUSCLE3 synchronisation !!!!
                switch upper(internal_state.actor.io.list.(noms{k}).operator)
                    case 'O_I'
                        if synchronization_flush == 1
                            disp('force MUSCLE3 workflow to flush the data preventing dead lock with output port operator O_I');
                            drawnow;      % flush matlab buffer (including graphics)
                            py.os.sync(); % flush systemp buffer
                        end
                end
        end
    end
end

function var = create_float_message(internal_state,data_name)

switch data_name
    case {'time'}
        var =  double(internal_state.workflow.clock.t_cur); 
    otherwise
        var = double(NaN);
end


function var = create_empty_message(internal_state,data_name)

var =  [];

function var = create_int_message(internal_state,data_name)

disp('not yet implemented')
keyboard

switch data_name
    case {'status'}
        var =  NaN; 
    otherwise
        var = NaN;
end

function var = create_str_message(internal_state,data_name)

disp('not yet implemented')
keyboard


function var = create_ids_message(internal_state,data_name,l)

% get IDS Matlab data structure
var_ids = internal_state.imas.ids_data.(data_name);
if iscell(var_ids)
        var_ids = var_ids{l};
end
% serialized data
var_serial = imas_serialize(var_ids,data_name);
% make a Python bytes stream
%var =  py.bytes(uint8(var_serial));
var = make_python_bytes(var_serial);

function var = create_matlab_message(internal_state,data_name)

switch data_name
    % case not METIS not yet implemented
    otherwise
        % METIS data structure
        eval(sprintf('%s = internal_state.metis.onetime.(data_name);',data_name));
end

% serialize Matlab sructure
tmpf = sprintf('%s.mat',tempname);
save(tmpf,data_name);
fid = fopen(tmpf,'r');
var_serial = fread(fid,Inf,'char*1');
fclose(fid);
delete(tmpf);

% make a Python bytes stream
%var =  py.bytes(uint8(var_serial));
var = make_python_bytes(var_serial);





