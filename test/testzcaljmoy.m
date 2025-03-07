% script de test de zcaljmoy
function testzcaljmoy

root = getappdata(0,'root');
if isempty(root)
	zineb_path;
	root = getappdata(0,'root');
end

% ajout des path pour chargement automatique du fichier
if isdir(fullfile(root,'certification','fullruns'))
	addpath(fullfile(root,'certification','fullruns'),'-begin');
end
if isdir(fullfile(root,'certification','cronosfunctions'))
	addpath(fullfile(root,'certification','cronosfunctions'),'-begin');
end
if isdir(fullfile(root,'certification','modules'))
	addpath(fullfile(root,'certification','modules'),'-begin');
end

t=dir(fullfile(root,'certification','modules'));

for k=1:size(t,1)
 	data = struct([]);
	if (t(k).name(1) ~= '.') & ~isdir(t(k).name) 
		try
			data= load(fullfile(root,'certification','modules',t(k).name));
		catch
			data = struct([]);
		end	
	end
	
	if ~isempty(data)
		test = data.test;
		noms = fieldnames(test);
		for l = 1:length(noms)
			ud = test.(noms{l});
			datak = zget1t(ud.data,length(ud.data.gene.temps));
			datan=zprofequi(ud.param.phys,ud.param.gene,datak);
			x    = ud.param.gene.x;
			figure
			plot(x,datak.equi.jmoy,'b',x,datak.prof.jmoy,'r',x,datan.prof.jmoy,'k');
			xlabel('x')
			ylabel('Jmoy A/m^2');
			title(noms{l});
			legend('equi ref','prof ref','test')
			drawnow
		end		
	end
end