% script d'appel du zerod pour le fit de l'efficacite de LH
% recherche des temps avec LH
plh = z0dinput.cons.plh;
indplh  = find(plh >0.1e5);
if isempty(indplh)
	defaultanswer={num2str(min(z0dinput.cons.temps)),num2str(max(z0dinput.cons.temps))};
else
	plhm    = mean(plh(indplh));
	indok   = find(plh > (0.3 .* plhm));
	tmin    = min(z0dinput.cons.temps(indok));
	tmax    = max(z0dinput.cons.temps(indok));
	defaultanswer={num2str(tmin),num2str(tmax)};

end
prompt={'begining time (s)','final time (s)'};
name='LH time interval for efficiency fitting';
numlines=1;
answer=inputdlg(prompt,name,numlines,defaultanswer);
pause(1);

if ~isempty(answer)
	z0dinput.thyb = [str2num(answer{1}),str2num(answer{2})];
	option = z0dinput.option;
	option.lhmode = 1;
	option.scaling = 4;
	[post.zerod,void,post.profil0d] = zerodfast(option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d,z0dinput.thyb);
	post.z0dinput = z0dinput;
	z0plotsc;
end
