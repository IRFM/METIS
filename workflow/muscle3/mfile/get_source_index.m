function index_source = get_source_index(ids_in,source_names)

% for core_sources selection of the source
index_source = NaN;
for k=1:length(ids_in.source)
    name = deblank(upper(ids_in.source{k}.identifier.name));
    if ~isempty(name)
        if ~isempty(strmatch(name,source_names,'exact'))
            index_source = k;
            break
        end
    end
end
