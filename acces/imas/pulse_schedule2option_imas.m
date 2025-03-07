% cette fonction extrait les options de metis du cpo pulse_schedule
function options = pulse_schedule2option_imas(pulse_schedule)

% add codeparam to options
info = metis4imas(1);
options = info.valeur;
codeparam = pulse_schedule.code.parameters;
if ~isempty(codeparam)
	tpn = tempname;
	fid = fopen(tpn,'w');
	if iscell(codeparam)
        ind = max(find(codeparam{end} == '>'));
		fprintf(fid,'%s\n',codeparam{end}(1:ind));
	else
        ind = max(find(codeparam  == '>'));
		fprintf(fid,'%s\n',codeparam(1:ind));
	end
	fclose(fid);
	info  = xml_read(tpn);
	delete(tpn);
	noms = fieldnames(info);
	for k=1:length(noms)
		if isfield(options,noms{k})
			options.(noms{k}) = info.(noms{k});
		end
	end
end 
