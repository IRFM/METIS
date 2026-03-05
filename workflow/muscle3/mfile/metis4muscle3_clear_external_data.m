function metis4muscle3_clear_external_data

% removing app data
noms = fieldnames(getappdata(0));
lcount = 0;
for k = 1:length(noms)
    if ~isempty(findstr(noms{k},'_EXP'))
        rmappdata(0,noms{k});
        lcount = lcount + 1;
    end
end
if lcount > 0
    fprintf('METIS4MUSCLE3: %d external METIS data have been removed\n',lcount);
end

% clear external data for mechanism mimicking ST
if isappdata(0,'LAST_TIME_ST')
    rmappdata(0,'LAST_TIME_ST');
end
if isappdata(0,'INDICE_INV_ST')
    rmappdata(0,'INDICE_INV_ST');
end
if isappdata(0,'MINIMAL_TAU_ST')
    rmappdata(0,'MINIMAL_TAU_ST');
end
