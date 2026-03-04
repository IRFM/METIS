function ports = metis4muscle3_make_ports(io)


% Ports need to match the MUSCLE3 yaml description on the other side !
if ~isfield(io,'ports') || isempty(io.ports)
    ports = py.dict();
    operator = [];
    % make list by operator
    name = fieldnames(io.list);
    for k=1:length(name)
        switch upper(io.list.(name{k}).operator)
            case 'NONE'
                % not a port on real actor
            otherwise
                if isfield(operator,io.list.(name{k}).operator)
                    operator.(io.list.(name{k}).operator){end+1} = name{k};
                else
                    operator.(io.list.(name{k}).operator){1} = name{k};
                end
        end
    end
    % loop on ports
    name = fieldnames(operator);
    for k=1:length(name)
        ports{py.getattr(@py.ymmsl.Operator,name{k})} = py.list(operator.(name{k}));
    end
end
