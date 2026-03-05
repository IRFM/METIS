function xml = simple_struct2xml(info,xml)

noms = fieldnames(info);
if nargin < 2
    xml = '<?xml version="1.0" encoding="utf-8"?>';
    xml = sprintf('%s\n%s',xml,'<ROOT>');
end
for k=1:length(noms)
    if ischar(info.(noms{k}))
        xml = sprintf('%s\n<%s>%s</%s>',xml,noms{k},deblank(info.(noms{k})),noms{k});
    elseif length(info.(noms{k})) > 1
        data_i = info.(noms{k});
        for l=1:length(info.(noms{k}))
            xml = sprintf('%s\n<table_%s>',xml,noms{k});
            if iscell(data_i)
                val = data_i{l};
            else
                val = data_i(l);                
            end
            if ischar(val)
                xml = sprintf('%s\n<%s>%s</%s>',xml,noms{k},deblank(val),noms{k});
            else
                xml = sprintf('%s\n<%s>%g</%s>',xml,noms{k},val,noms{k});
            end
            xml = sprintf('%s\n</table_%s>',xml,noms{k});
        end
    elseif iscell(info.(noms{k}))
        val = info.(noms{k});
        xml = sprintf('%s\n<%s>%g</%s>',xml,noms{k},val{:},noms{k});
    elseif isstruct(info.(noms{k}))
        xml = sprintf('%s\n<%s>',xml,noms{k});
        xml = simple_struct2xml(info.(noms{k}),xml);
        xml = sprintf('%s\n</%s>\n',xml,noms{k});        
    else
        xml = sprintf('%s\n<%s>%g</%s>',xml,noms{k},info.(noms{k}),noms{k});
    end
end
if nargin < 2
    xml = sprintf('%s\n%s',xml,'</ROOT>');
end

