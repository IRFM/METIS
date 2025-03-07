% call metis with option to fit energy content and LH current drive efficiency.
function metis_fitlh

drawnow
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



z0dinput = evalin('base','z0dinput');
plh = z0dinput.cons.plh;
indplh  = find(plh >0.1e5);
if isempty(indplh)
        z0dinput.thyb = [min(z0dinput.cons.temps),max(z0dinput.cons.temps)];
else
	plhm    = mean(plh(indplh));
	indok   = find(plh > (0.3 .* plhm));
	tmin    = min(z0dinput.cons.temps(indok));
	tmax    = max(z0dinput.cons.temps(indok));
        z0dinput.thyb = [tmin,tmax];
end
z0dinput.option.lhmode = 1;
z0dinput.option.scaling = 4;
z0dinput.option.evolution = 0;
z0dinput.option.machine = z0dinput.machine;
z0dinput.option.shot = z0dinput.shot;
% recopie machine dans option
zassignin('base','z0dinput',z0dinput);
fprintf('METIS_fitLH: call metis with option to fit energy content and LH current drive efficiency (between  %g s and %g s)\n',z0dinput.thyb(1),z0dinput.thyb(2));
evalin('base','[post.zerod,void,post.profil0d] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d,z0dinput.thyb);post.z0dinput = z0dinput;z0plotsc;', ...
	'errordlg(lasterr,''Error during the run of Metis simulator !'');');
