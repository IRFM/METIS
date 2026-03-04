function port_list = change_port_operator(port_list,operator_init,operator_input,operator_output)

if nargin < 2
    operator_init  = '';
    operator_input  = '';
    operator_output = '';
elseif nargin < 3
    operator_input  = '';
    operator_output = '';
elseif nargin < 4
    operator_output = '';
end

% update operator if needed
if nargin >= 2
    noms = fieldnames(port_list);
    for k=1:length(noms)
        switch port_list.(noms{k}).direction
            case 'init'
                if ~isempty(operator_init)
                    if isstruct(operator_init)
                        port_list.(noms{k}).operator = operator_init.(noms{k});
                    else
                        port_list.(noms{k}).operator = operator_init;
                    end
                end
            case 'input'
                if ~isempty(operator_input)
                    if isstruct(operator_input)
                        port_list.(noms{k}).operator = operator_input.(noms{k});
                   else
                        port_list.(noms{k}).operator = operator_input;
                    end
                end
            case 'output'
                if ~isempty(operator_output)
                    if isstruct(operator_output)
                        port_list.(noms{k}).operator = operator_output.(noms{k});
                    else
                        port_list.(noms{k}).operator = operator_output;
                    end
                end
        end
    end
end
 