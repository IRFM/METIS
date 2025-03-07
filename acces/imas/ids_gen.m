% ids_gen(IDSname)
% 
% Generate an IDS of specified type with default field values.
% Matlab version to overcome obsolescence of IMAS mexfile
% 
% IDSname : name of the IDS to generate.
%
% NOTE: To initialise IDSs ids_gen has been deprecated in favor of
% ids_init.
%
% NOTE: The array of structures in the resulting IDS will be
% of size 1.
%
% See also: ids_init
function ids = ids_gen(ids_name)

ids = ids_init(ids_name);
ids = ids_gen_allocate(ids,ids_name,'');

function out = ids_gen_allocate(in,ids_name,pathinids)

% test 
if ~isempty(in)
    out = in;
    noms = fieldnames(in);
    for k=1:length(noms)
        if (isempty(in.(noms{k})) && iscell(in.(noms{k}))) || isstruct(in.(noms{k}))
            if ~isempty(pathinids)
                pids = sprintf('%s/%s',pathinids,noms{k});
            else
                pids = noms{k};
            end
            out.(noms{k}) = ids_gen_allocate(in.(noms{k}),ids_name,pids);
        end
    end
elseif ~ischar(in)
    try
        inter = ids_allocate(ids_name,pathinids,1);
        % next level
        out{1} = ids_gen_allocate(inter{1},ids_name,pathinids);
    catch
        out = [];
    end   
end
    