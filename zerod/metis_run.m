function metis_run

close_0plot;
% recopie machine dans option
evalin('base','z0dinput.option.machine = z0dinput.machine;');
evalin('base','z0dinput.option.shot = z0dinput.shot;');

% information sur les donnees externes
noms = fieldnames(getappdata(0));
text_warn = '';
for k = 1:length(noms)
	if findstr(noms{k},'_EXP')
		fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
		text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
	end
end
if isdeployed && ~isempty(text_warn)
  warndlg(text_warn,'Pay attention: external data used');
end

% pour eviter les incoherence sur le mode restart
evalin('base','clear  z0dstruct');
  
% overwrite parameters whith reference parameters if defined
evalin('base','z0dinput.option = z0doverwriteparam(z0dinput.option);');

% sort du mode evolution
evalin('base','z0dinput.option.evolution = 0;');

% recopie machine dans option
evalin('base','zerod_machine_name;');


evalin('base','[post.zerod,void,post.profil0d] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);post.z0dinput = z0dinput;z0plotsc;', ...
	'error(lasterr,''Error during the run of Metis simulator !'');');
