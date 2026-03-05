function [timestamp,data,port_name,real_used_port] = get_ids_data_from_ports(data_input,ports_list,ids_name,used_ports)


% initailisation
timestamp = [];
data = [];
port_name = '';
real_used_port = '';

% correspondance table (with used_ports)
operator.f_init    = 'F_INIT';
operator.input_f   = 'F_INIT';
operator.input_s   = 'S';
operator.output_f  = 'O_F';
operator.output_i  = 'O_I';
operator.in_f      = 'F_INIT';
operator.in_s      = 'S';
operator.out_f     = 'O_F';
operator.out_i     = 'O_I';

% loop on ports (new implementation to scan all ports if multiple
% connexion)
nomp = fieldnames(ports_list);
% loop on used ports (last read overwrite previous read data)
for k=1:length(used_ports)
    for l=1:length(nomp)
        switch ports_list.(nomp{l}).direction
            % not output
            case 'input'
                % selection on used ports in actor
                switch used_ports{k}
                    case {'input_f','input_s','in_f','in_s'}
                        if any(strcmp(ids_name, ports_list.(nomp{l}).list_data_name)) && ...
                                strcmp(operator.(used_ports{k}), ports_list.(nomp{l}).operator)
                            port_name = nomp{l};
                            real_used_port = used_ports{k};
                            %
                            break
                        end
                end
            case 'init'
                switch used_ports{k}
                    case 'f_init'
                        if any(strcmp(ids_name, ports_list.(nomp{l}).list_data_name))
                            port_name = nomp{l};
                            real_used_port = used_ports{k};
                            %
                            break
                        end
                end
        end
    end
    if ~isempty(port_name) && ~isempty(data_input.(port_name).data.(ids_name))
        if isempty(data)
            fprintf('copying data pointer from port %s for IDS %s (selected port = %s)\n',nomp{l},ids_name,real_used_port);
        else
            fprintf('copying data pointer from port %s for IDS %s and overwrinting previous data (selected port  = %s)\n',nomp{l},ids_name,real_used_port);
        end
        %
        data      = data_input.(port_name).data.(ids_name);
        timestamp = data_input.(port_name).timestamp;
    end
    
end

% if data isempty and ids name start with plasma_ then try with core_
if isempty(data) && startsWith(ids_name,'plasma_')
    ids_name = strrep(ids_name,'plasma_','core_');
    [timestamp,data,port_name,real_used_port] = get_ids_data_from_ports(data_input,ports_list,ids_name,used_ports);
end




    
