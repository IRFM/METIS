% calcul rapide avec extraction de point pertinent
%[post.zerod,void,post.profil0d] = zerodslow(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d,[],1001);post.z0dinput = z0dinput;z0plotsc;
function [zs,info,profli]= zerodslow(option,cons,geo,exp0d,thyb,nbte)

if nargin < 5
	thyb = [];
end

%MATLAB:interp1:NaNinY
warning off
% compatibilite
if ~isfield(cons,'xece') & isfield(option,'xece')
 	cons.xece = option.xece .* ones(size(cons.temps));
end

times   = cons.temps;
if length(nbte) == 1
	temps =  linspace(min(times),max(times),nbte)';
else
	temps = nbte(:);
end

% indication de reduction
fprintf('Metis sample in slow mode : %d -> %d -> %d\n',length(times),length(temps),length(times));

% modification des donnees
noms = fieldnames(cons);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(cons,nomc);
	valn  = interp1(times,val,temps,'linear');
	indbad      = find(any(~isfinite(valn),2));
	if ~isempty(indbad)
		valn(indbad,:) = ones(length(indbad),1) * val(end,:);
	end
	cons = setfield(cons,nomc,valn);
end
cons.temps = temps;

noms = fieldnames(geo);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(geo,nomc);
	if ~isempty(val)
		valn  = interp1(times,val,temps,'linear');
		indbad      = find(any(~isfinite(valn),2));
		if ~isempty(indbad)
			valn(indbad,:) = ones(length(indbad),1) * val(end,:);
		end
		geo = setfield(geo,nomc,valn);
	end
end

noms = fieldnames(exp0d);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(exp0d,nomc);
	if ~isempty(val)
		if (size(val,1) == length(times))
			valn  = interp1(times,val,temps,'linear');
			indbad      = find(any(~isfinite(valn),2));
			if ~isempty(indbad)
				valn(indbad,:) = ones(length(indbad),1) * val(end,:);
			end
			exp0d = setfield(exp0d,nomc,valn);
		end
	end
end

% appel du 0d
warning on
[zs,info,profli]= zerod(option,cons,geo,exp0d,1e-2,thyb);
warning off

% modification des donnees de sortie
noms = fieldnames(zs);
for l=1:length(noms)
	nomc = noms{l};
	val  = getfield(zs,nomc);
	if length(val) == length(temps)
		val  = interp1(temps,val,times,'linear');
		zs = setfield(zs,nomc,val);
	end
end
profli.temps = temps;

warning on

