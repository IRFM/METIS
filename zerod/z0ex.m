% extraction des donnees scalaires
function zsmat = z0ex(zsmat,zs,tmin,tmax)

indok = find((zs.temps >= tmin) & (zs.temps <= tmax));

noms =  fieldnames(zs);
for k = 1:length(noms);
    nomc  = noms{k};
    if ~isfield(zsmat,nomc)
           var   = getfield(zs,nomc);
           zsmat = setfield(zsmat,nomc,msvar(var,indok));
    else
           var    = getfield(zs,nomc);
           varmat = getfield(zsmat,nomc);
           if isempty(varmat)
               varmat(end + 1) = NaN;
           else
               varmat(end + 1) = msvar(var,indok);
           end
           zsmat  = setfield(zsmat,nomc,varmat);
    end
end


function x=msvar(x,ind)

if all(size(x)==1)
       return;
end
x = x(ind);
indnok = find(~isfinite(x));
x(indnok) = [];
if isempty(x)
    x = NaN;
else 
    x =mean(x);
end

