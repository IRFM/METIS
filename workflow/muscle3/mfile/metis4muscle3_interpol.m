function [zs,profil,z0dstruct] = metis4muscle3_interpol(post,time)

% nombre de temps interne 
nbte = 3;
%
t0 = time - max(1e-6,min(1,0.01 .* mean(post.z0dinput.geo.a(1)) .* mean(post.z0dinput.geo.R(1))));
temps    	= linspace(t0,time,nbte)';
temps           = cat(1,temps,2 .* temps(end) - temps(end-1));
vt              = ones(size(temps));
z0dinput 	= zerod_interp1t(post.z0dinput,temps);
z0dinput.option.init_geo  = 1;
z0dinput.option.evolution = 0;
z0dinput.option.vloop     = 0;
z0dstruct.zs       = zerod_interp1t(post.zerod,temps);
z0dstruct.info     = post.z0dinput.zsinfo;
z0dstruct.profil   = zerod_interp1t(post.profil0d,temps);
% seul les end-1 temps sont valides
noms = fieldnames(z0dstruct.zs);
for k = 1:length(noms)
    if length(z0dstruct.zs.(noms{k})) > 1
        z0dstruct.zs.(noms{k}) = z0dstruct.zs.(noms{k})(1:(end-1));
    end
end
noms = fieldnames(z0dstruct.profil);
for k = 1:length(noms)
    if size(z0dstruct.profil.(noms{k}),1) > 1
        z0dstruct.profil.(noms{k}) = z0dstruct.profil.(noms{k})(1:(end-1),:);
    end
end
noms = fieldnames(z0dinput.cons);
for k = 1:length(noms)
    if length(z0dinput.cons.(noms{k})) > 1
        z0dinput.cons.(noms{k}) = z0dinput.cons.(noms{k})(1:(end-1));
    end
end
noms = fieldnames(z0dinput.geo);
for k = 1:length(noms)
    if length(z0dinput.geo.(noms{k})) > 1
        z0dinput.geo.(noms{k}) = z0dinput.geo.(noms{k})(1:(end-1));
    end
end
noms = fieldnames(z0dinput.exp0d);
for k = 1:length(noms)
    if size(z0dinput.exp0d.(noms{k}),1) > 1
        z0dinput.exp0d.(noms{k}) = z0dinput.exp0d.(noms{k})(1:(end-1),:);
    end
end
z0dstruct.z0dinput = z0dinput;


% sorties
noms = fieldnames(z0dstruct.zs);
zs = [];
for k=1:length(noms)
    zs.(noms{k}) = z0dstruct.zs.(noms{k})(end);
end
profil = [];
noms = fieldnames(z0dstruct.profil);
for k=1:length(noms)
    profil.(noms{k}) = z0dstruct.profil.(noms{k})(end,:);
end



% interpolation 1 temps d'une structure
function struct_out = zerod_interp1t(struct_in,temps)

if ~isstruct(struct_in)
    if isempty(struct_in)
        struct_out = struct_in;        
    elseif isnumeric(struct_in)
	if size(struct_in,1) == 1
	    struct_out = struct_in;
	else
	    struct_out = struct_in(end)  .* ones(size(temps));
	end
    elseif iscell(struct_in)
         struct_out = struct_in{end};
    else
        struct_out = struct_in;    
    end
else
    noms = fieldnames(struct_in);
    if isfield(struct_in,'temps')
        temps_loc = struct_in.temps;
    else
        temps_loc = [];
    end
    for k = 1:length(noms)
        if isstruct(struct_in.(noms{k}))
            struct_out.(noms{k})  = zerod_interp1t(struct_in.(noms{k}),temps);
        elseif isempty(temps_loc)
            struct_out.(noms{k})  = zerod_interp1t(struct_in.(noms{k}),temps);
        elseif length(struct_in.temps) == size(struct_in.(noms{k}),1)
            struct_out.(noms{k}) = interp1(temps_loc,struct_in.(noms{k}),temps,'linear',NaN);
            near = interp1(temps_loc,struct_in.(noms{k}),temps,'nearest','extrap');
            struct_out.(noms{k})(~isfinite(struct_out.(noms{k}))) = near(~isfinite(struct_out.(noms{k})));
        else
            struct_out.(noms{k})  = zerod_interp1t(struct_in.(noms{k}),temps);
        end
    end
    if ~isempty(temps_loc)
        struct_out.temps = temps;
    end
end
 
 
 
function [Rsepa,Zsepa] = resample_sepa(R,Z,nbp)

% reservation memoire
Rsepa = NaN * ones(size(R,1),nbp);
Zsepa = NaN * ones(size(Z,1),nbp);

% mise sur les memes angles
Rc   = (max(R,[],2) + min(R,[],2)) ./ 2;
Zc   = (max(Z,[],2) + min(Z,[],2)) ./ 2; 
th   = linspace(0,2.*pi,nbp);
vh   = ones(1,size(R,2));
vt   = ones(size(R,1),1);
cx   = (R - Rc*vh) + sqrt(-1) .* (Z-Zc*vh);
thx  = unwrap(angle(cx),[],2);
thx = angle(cx);
thx = thx .* (thx>=0) + (2 .* pi + thx) .* (thx <0);
rhox = abs(cx);
thx(thx<0) = thx(thx<0) + 2 .* pi;
for k = 1:size(thx,1)
	[thk,indice] = unique(thx(k,:));
	rhok         = rhox(k,indice);
	[thk,indice] = sort(thk); 
	rhok         = rhok(indice);
	rhok = cat(2,rhok,rhok(:,2:end-1),rhok);
	thk = cat(2,thk -2.*pi,thk(:,2:end-1),thk+2.*pi);
	ind_bad = find(diff(thk) <= 0);
	if ~isempty(ind_bad)
	    thk(ind_bad) = [];
	    rhok(ind_bad) = [];
	    if length(thk) < 3
		error('no enougth distinct points to resample LCFS');
	    end
	end
	rhoo   = pchip(thk,rhok,th);
	R      = Rc(k)  + rhoo .* cos(th);
	Z      = Zc(k)  + rhoo .* sin(th);
	Rsepa(k,:) = R;
	Zsepa(k,:) = Z;
end
