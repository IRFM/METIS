% cette fonction extrait les options de metis du cpo pulse_schedule
function options = pulse_schedule2option_imas(pulse_schedule)

% add codeparam to options
info_all = metis4imas(1);
options = info_all.valeur;
info = [];
codeparam = pulse_schedule.code.parameters;
if ~isempty(codeparam)
    try
        info = decode_codeparam(codeparam);
    catch
        disp('switch to xml_read to decode the code parameters');
    end
    if isempty(info)
        disp('error in decode_codeparam, use alternative methode');
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
    end
	noms = fieldnames(info);
	for k=1:length(noms)
		if isfield(options,noms{k})
			options.(noms{k}) = info.(noms{k});
		end
	end
end 


