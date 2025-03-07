% script de lecture des donnees experuimatales de densite pour METIS
switch post.z0dinput.machine
case 'TS'
	% recherche des occurence tprof
	shot = fix(post.z0dinput.shot);
	occtab =[];
	tmintab = [];
	tmaxtab = [];
	occok = NaN;
	tminok = -inf;
	tmaxok = inf;
	for occ=0:9			
		[ip,t] = tsbase(shot+occ/10,'sprofip');
		if ~isempty(t);
			occtab(end+1) = occ;
			tmintab(end+1) = min(t);
			tmaxtab(end+1) = max(t);
		end	
	end
	if isempty(occtab)
		disp('utilisation de gne');
		% lecture de gne
		[ne,tne,xne,cne] = tsbase(shot,'gne');
		ne = zautofilterprof(tne,xne,ne,'wfilt');
		nes.temps = tne;
		nes.ne    = ne;
		nes.x     = xne;
		setappdata(0,'NE_EXP',nes);
	else
		%
		for k =1:length(occtab)
			fprintf('choc #%d.%d de %g a %g\n',shot,occtab(k),tmintab(k),tmaxtab(k));
		end
		
		rep = input('choix de l''occurence ?');
		if ~isempty(rep)
			occok = rep;
		else
			occok = min(occtab);
		end
		% lecture de tprof
		[ne,tne,xne,cne] = tsbase(shot+occok/10,'gprofnefit');
		ne = zautofilterprof(tne,xne,ne,'wfilt');
		nes.temps = tne;
		nes.ne    = ne;
		nes.x     = xne;
		setappdata(0,'NE_EXP',nes);
	end	

otherwise 
	error('unplug device');
end