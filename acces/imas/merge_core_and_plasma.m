function ids = merge_core_and_plasma(root_name,second)

if nargin < 2
    ids  = ids_gen(sprintf('core_%s',root_name));
    pids = ids_gen(sprintf('plasma_%s',root_name));
else
    ids = root_name;
    pids = second;
end
% merge both
ids = merge_ids(ids,pids);


function ids = merge_ids(ids,pids)

if ~isstruct(pids)
    if iscell(pids)
        for l=1:min(length(ids),length(pids))
            ids{l} = merge_ids(ids{l},pids{l});
        end
    elseif isempty(ids) && ~isempty(pids)
        ids = pids;
    end
    return
elseif isstruct(pids) && ~isstruct(ids)
    ids = pids;
    return
end

noms = fieldnames(pids);
for k=1:length(noms)
    if isfield(ids,noms{k})
        if iscell(ids.(noms{k}))
            for l=1:length(ids.(noms{k}))
                ids.(noms{k}){l} = merge_ids(ids.(noms{k}){l},pids.(noms{k}){l});
            end
        else
            ids.(noms{k}) = merge_ids(ids.(noms{k}),pids.(noms{k}));
        end
    else
        ids.(noms{k}) = pids.(noms{k});
    end
end

