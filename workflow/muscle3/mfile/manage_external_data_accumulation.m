function  data = manage_external_data_accumulation(data,data_1t)


% METIS convention = temps for time
if isempty(data)
     data = data_1t;
elseif isempty(data.temps)
     data = data_1t;
elseif data_1t.temps > data.temps(end)
    % just add it at this stage, maybe latter we will limit the length is
    % there is som memory size problem
    % evolution
    noms = fieldnames(data);
    for k=1:length(noms)
        switch noms{k}
            case 'x'
                data.(noms{k}) = data_1t.(noms{k})(1,:);
            otherwise
                data.(noms{k})(end+1,:) = data_1t.(noms{k});
        end
    end
elseif data_1t.temps == data.temps(end)
    % convergence
    noms = fieldnames(data);
    for k=1:length(noms)
        switch noms{k}
            case 'x'
                data.(noms{k}) = data_1t.(noms{k})(1,:);
            otherwise
                data.(noms{k})(end,:) = data_1t.(noms{k});
        end
    end
elseif  data_1t.temps < data.temps(end) 
    % rewind
    mask = (data_1t.temps > data.temps);
    noms = fieldnames(data);
    for k=1:length(noms)
        switch noms{k}
            case 'x'
                data.(noms{k}) = data_1t.(noms{k})(1,:);
            otherwise
                if any(mask)
                    data.(noms{k}) = cat(1,data.(noms{k})(mask,:),data_1t.(noms{k}));
                else
                    data.(noms{k}) = data_1t.(noms{k});
                end
        end
    end
    
end
