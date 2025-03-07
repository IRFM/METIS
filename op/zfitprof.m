% reechantillionne spacialement les profils
% rend des profils monotones positif 
function vout=zfitprof(xin,vin,xout)

if size(xin,1) == 1
	xin = ones(size(vin,1),1) * xin;
end
xoutm = ones(size(vin,1),1) * xout;

ind   = find(~isfinite(vin));
if ~isempty(ind)
	vin(ind) =zeros(1,length(ind));
end

% filtrage ?
vv =zeros(size(vin));
for k = 1:size(vin,2)
	vv(:,k) = medfilt1(vin(:,k),3);
end

% prolongation de l'intervalle
comp    = ones(size(vin,1),1);
xin     = cat(2,-comp,xin,1.5*comp);
vv      = cat(2,0.*comp,vv,0.*comp);
vvp     = tsplinet(xin,vv,xoutm);
dvvdxx  = pdederive(xoutm,vvp,0,2,2,1);

% distribution des derivees >0
count  = 100;
while (count >0)&(~isempty(find(dvvdxx > 0)))
	comp    = zeros(size(dvvdxx,1),1);
	dvvdxxg = cat(2,dvvdxx(:,2:end),comp);
	dvvdxxd = cat(2,comp,dvvdxx(:,1:(end-1)));
	dvvdxx  = dvvdxx .* (dvvdxx <= 0) + 0.5 .* dvvdxxg .* (dvvdxxg > 0) + 0.5 .* dvvdxxd .* (dvvdxxd > 0);
	count   = count - 1;
end
dvvdxx  = dvvdxx .* (dvvdxx <= 0) ;

% reconstitution du profil
vvc     = cumtrapz(xout',dvvdxx')';
vvc     = vvc  + ((vvp(:,end).* (vvp(:,end)>0) - vvc(:,end)) * ones(1,size(xoutm,2)));
vout    = vvc  .* ((trapz(xout',vvp')' ./ trapz(xout',vvc')') * ones(1,size(xoutm,2)));

ind   = find(~isfinite(vout));
if ~isempty(ind)
	vout(ind) =zeros(1,length(ind));
end
