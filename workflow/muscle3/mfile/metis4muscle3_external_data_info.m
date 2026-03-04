function metis4muscle3_external_data_info(logger)

% information sur les donnees externes
noms = fieldnames(getappdata(0));
first = true;
for k = 1:length(noms)
	if ~isempty(strfind(noms{k},'_EXP'))
        if first
            first = false;
            logger.info('METIS4MUSCLE3: external data in used !')
        end
		fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));		
	end
end
