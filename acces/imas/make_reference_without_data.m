function ps_out = make_reference_without_data(ps_in)

if ~isfield(ps_in.flux_control.i_plasma.reference,'data')
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
            ps_out.(noms{k}) = remove_reference_data_field(ps_in.(noms{k}));
    end
end


function out = remove_reference_data_field(in)

if iscell(in)
    if isempty(in)
        out = in;
    else
        for k = 1:length(in)
            out{k} = remove_reference_data_field(in{k});
        end
    end
elseif isstruct(in)
    if isempty(in)
        out = in;
        return
    elseif length(in) > 1
        for k=1:length(in)
            out(k) = remove_reference_data_field(in(k));
        end
    else
        noms = fieldnames(in);
        for k = 1:length(noms)
            switch noms{k}
                case 'reference'
                    if isstruct(in.(noms{k}))
                        if isfield(in.(noms{k}),'data')
                             out.(noms{k}) = in.(noms{k}).data;                        
                       else
                            error('this case should never happen');
                        end
                    else
                            % it is OK
                            out.(noms{k}) = in.(noms{k});
                    end
                otherwise
                    out.(noms{k}) = remove_reference_data_field(in.(noms{k}));
            end
        end
    end
else
        out = in;    
end
        