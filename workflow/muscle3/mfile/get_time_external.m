function time_loc = get_time_external(settings,data_ids,timestamp)

% select time
if settings.idstime
    if isfield(data_ids,'time') && ~isempty(data_ids.time)
        time_loc = data_ids.time - settings.timeoffset;
    elseif isfield(data_ids,'profiles_1d') &&  ~isempty(data_ids.profiles_1d{end}.time)
        time_loc = data_ids.profiles_1d{end}.time - settings.timeoffset;
    elseif isfield(data_ids,'time_slice') &&  ~isempty(data_ids.time_slice{end}.time)
        time_loc = data_ids.time_slice{end}.time - settings.timeoffset;
    else
        time_loc = timestamp - settings.timeoffset;
    end
else
    time_loc = timestamp - settings.timeoffset;
end
