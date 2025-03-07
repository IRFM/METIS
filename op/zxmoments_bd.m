% fonction d'extraction des moments pour les champs de la base de donnees cronos
function param = zxmoments_bd(param,temps,val,kmin,kmax,mask)

% securite
val = real(val);

% extraction des temps pertinents
indnin = max(1,min(kmin,length(val))):max(1,min(kmax,length(val)));
val  = val(indnin);
t    = temps(indnin);
if ~isempty(mask)
    mask = mask(indnin);
end
% calcul des moments
% maximum tout temps
r = max(val);
if isfinite(r)
	param.simMax = java.lang.Double(r);
end
% minimum tout temps
r = min(val);
if isfinite(r)
	param.simMin = java.lang.Double(r);
end
% moyenne tout temps
r = mean(val);
if isfinite(r)
	param.simAvg = java.lang.Double(r);
end
% median tout temps
r = median(val);
if isfinite(r)
	param.simMed = java.lang.Double(r);
end
% valeur la plus probable
r = zlikehood(val);
if isfinite(r)
	param.simProb = java.lang.Double(r);
end
% ecart type tout temps
r = std(val);
if isfinite(r)
	param.simStd = java.lang.Double(r);
end

% donnees sur le plateau
if ~isempty(mask)
	val     = val(find(mask));
else
	val     = [];
end
if isempty(val)
	return
end

% calcul des moments
% maximum tout temps
r = max(val);
if isfinite(r)
	param.plateauMax = java.lang.Double(r);
end
% minimum tout temps
r = min(val);
if isfinite(r)
	param.plateauMin = java.lang.Double(r);
end
% moyenne tout temps
r = mean(val);
if isfinite(r)
	param.plateauAvg = java.lang.Double(r);
end
% median tout temps
r = median(val);
if isfinite(r)
	param.plateauMed = java.lang.Double(r);
end
% valeur la plus probable
r = zlikehood(val);
if isfinite(r)
	param.plateauProb = java.lang.Double(r);
end
% ecart type tout temps
r = std(val);
if isfinite(r)
	param.plateauStd = java.lang.Double(r);
end



% recherche de la valeur la plus probable
function  sout = zlikehood(sin);

sin = sin(isfinite(sin));
if isempty(sin)
	sout = NaN;
	return
end

stdin = std(sin);
if stdin == 0
	[n,x]= hist(sin);
else
	[n,x]= hist(sin,fix(pi./std(sin) .*( max(sin) - min(sin)))+1);
end

ind  = find(n == max(n));
if ~isempty(ind)
	sout = mean(x(ind));
else
	sout = NaN;
end
