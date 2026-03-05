function codeparam = encode_codeparam(data)

% headear
codeparam=sprintf('%s\n%s\n','<?xml version="1.0" encoding="utf-8"?>','<data>');

% loop on fields
noms = fieldnames(data);
for k=1:length(noms)
    if isempty(data.(noms{k}))
        codeparam=sprintf('%s<%s/>\n',codeparam,noms{k});
    else
        if isnumeric(data.(noms{k}))
            val = strtrim(sprintf('%16.16g',data.(noms{k})));
        else
            val = data.(noms{k});
        end
        codeparam=sprintf('%s<%s>%s</%s>\n',codeparam,noms{k},val,noms{k});
    end
end

% closing the nested structure
codeparam=sprintf('%s</data>\n',codeparam);