function vout=zmonotone(xx,vv)

if size(xx,1) ~= 1
	xx  = xx(1,:);
end

ind   = find(~isfinite(vv));
if ~isempty(ind)
	vv(ind) =zeros(1,length(ind));
end

% derivation
dvvdxx  = pdederive(xx,vv,0,2,2,1);

% distribution des derivees >0
count  = 100;
while (count >0)&(~isempty(find(dvvdxx > 0)))
	comp    = zeros(size(dvvdxx,1),1);
	dvvdxxg = cat(2,dvvdxx(:,2:end),comp);
	dvvdxxd = cat(2,comp,dvvdxx(:,1:(end-1)));
	dvvdxx  = dvvdxx .* (dvvdxx <= 0) + 0.5 .* dvvdxxg .* (dvvdxxg > 0) + 0.5 .* dvvdxxd .* (dvvdxxd > 0);
	% derivee nulle au centre
	dvvdxx(:,1:3)  = zeros(size(dvvdxx(:,1:3)));
	count   = count - 1;
end
dvvdxx  = dvvdxx .* (dvvdxx <= 0) ;

% reconstitution du profil
vvc     = cumtrapz(xx',dvvdxx')';
vvc     = vvc  + ((vv(:,end).* (vv(:,end)>0) - vvc(:,end)) * ones(1,length(xx)));
vout    = vvc  .* ((trapz(xx',vv')' ./ trapz(xx',vvc')') * ones(1,length(xx)));

ind   = find(~isfinite(vout));
if ~isempty(ind)
	vout(ind) =zeros(1,length(ind));
end
