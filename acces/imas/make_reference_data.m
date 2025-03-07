function ps_out = make_reference_data(ps_in)

if isfield(ps_in.flux_control.i_plasma.reference,'data')
    % nothing to do
    ps_out = ps_in;
    return
end

noms = fieldnames(ps_in);
for k=1:length(noms)
    switch noms{k}
        case {'ids_properties','code','time'}
            % no action on these field
            ps_out.(noms{k}) = ps_in.(noms{k});
        otherwise
            ps_out.(noms{k}) = add_reference_data_field(ps_in.(noms{k}));
    end
end


function out = add_reference_data_field(in)

if iscell(in)
    if isempty(in)
        out = in;
    else
        for k = 1:length(in)
            out{k} = add_reference_data_field(in{k});
        end
    end
elseif isstruct(in)
    if isempty(in)
        out = in;
        return
    elseif length(in) > 1
        for k=1:length(in)
            out(k) = add_reference_data_field(in(k));
        end
    else
        noms = fieldnames(in);
        for k = 1:length(noms)
            switch noms{k}
                case 'reference'
                    if isstruct(in.(noms{k}))
                        if isfield(in.(noms{k}),'data')
                            % it is OK
                            out.(noms{k}) = in.(noms{k});
                        else
                            error('this case should never happen');
                        end
                    else
                            out.(noms{k}).data = in.(noms{k});                        
                    end
                otherwise
                    out.(noms{k}) = add_reference_data_field(in.(noms{k}));
            end
        end
    end
else
        out = in;    
end
        