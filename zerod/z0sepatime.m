function z0dinput = z0sepatime(input,z0dinput)


if nargin == 0
	input.temps 	= 'vector of time (s)';
	input.rxup      = 'upper triangularity (minor radius unit), X point distance from magnetic axis';
	input.zxup      = 'upper altitude X point (minor radius unit)';
	input.apup      = 'upper separatrix angle (R,X)  (LFS, degrees)';
	input.amup      = 'upper separatrix angle (-R,X) (HFS, degrees)';
	input.ra        = 'major radius (m) [6.2]';
	input.za        = 'altitude of the plasma (m) [6.2]';
	input.a         = 'minor radius (m) [2]';
	input.rxdo      = 'lower triangularity (minor radius unit), X point distance from magnetic axis';
	input.zxdo      = 'lower altitude X point (minor radius unit)';
	input.apdo      = 'lower separatrix angle (R,X)  (LFS, degrees)';
	input.amdo      = 'lower separatrix angle (-R,X) (HFS, degrees)';
	input.b0        = 'maximal toroidal magnetic field on the conductor (T)';
	input.delta     = 'minimal distance between conductor an plasma (m)';
	input.xpoint    = 'first time for X point separatrix';
	z0dinput	= input;
	return
end

% temps
temps = z0dinput.cons.temps;

% reservation
z0dinput.exp0d.Rsepa = NaN * ones(length(temps),201);
z0dinput.exp0d.Zsepa = NaN * ones(length(temps),201);


% boucle sur le champs
% interpolation
noms = fieldnames(input);
moments = [];
for k=1:length(noms)
	switch noms{k}
	case 'temps'
		% rien 
	case 'xpoint'
		% rien 
	otherwise
		moments.(noms{k}) = interp1(input.temps,input.(noms{k}),temps,'linear',input.(noms{k})(end)); 
	end
end


% boucle sur les temps
for k=1:length(temps)
	ind = max(find(input.temps <= temps(k))); 
	if input.xpoint(ind) == 1
		% separatrice x point
		noms = fieldnames(moments);
		option = [];
		for l=1:length(noms)
			option.(noms{l}) = moments.(noms{l})(k); 
		end
		option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
		option.mode      = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
		sepa1  = z0dsepanew2(temps(k),option);
		% calcul des moments
		% centre pour angle d'integration
		rc = mean(sepa1.R,2);
		zc = mean(sepa1.Z,2);
		vc = ones(1,size(sepa1.R,2));
		uc = unwrap(angle((sepa1.R-rc*vc) + sqrt(-1) .* (sepa1.Z  -zc*vc)));
		uc    = uc .* (uc >0) + (uc + 2*pi) .* (uc<= 0);
		uc(:,1)   = uc(:,end) + 2 .* pi;
		xu    = linspace(0,1,length(vc));
		%dudx  = pdederive(xu,uc,2,2,2,1);
		%dudx(:,1) = (dudx(:,1) +dudx(:,end)) ./ 2;
		%dudx(:,end) = dudx(:,1);
		dRdx  = pdederive(xu,sepa1.R,2,2,2,1);
		dZdx  = pdederive(xu,sepa1.Z,2,2,2,1);
		% calcul de R0 et Z0
		maskrmax  = (sepa1.R == max(sepa1.R));
		% recalcul des parametres sur le vecteur final
		rmin  = min(sepa1.R);
		rmax  = max(sepa1.R);
		geo.a = 0.5 .* (rmax - rmin);
		geo.R = 0.5 .* (rmax + rmin);
		zmin  = min(sepa1.Z);
		zmax  = max(sepa1.Z);
		%geo.z0   = (zmax + zmin) ./ 2 + moments.za(k);
		%geo.K    = (zmax - zmin) ./ 2 ./ geo.a;
		geo.z0   = option.za;
		%geo.z0   = (zmax + zmin + sum(sepa1.Z .* maskrmax,2) ./ sum(maskrmax,2)) ./ 3 + option.za;
		geo.K    = (abs(trapz(xu,sepa1.Z .*  dRdx,2) ./ pi ./ geo.a) + (zmax - zmin)) ./ 3 ./ geo.a;
		
		rzmax = sepa1.R(min(find(sepa1.Z == zmax)));
		rzmin = sepa1.R(min(find(sepa1.Z == zmin)));
		uu   =  angle(rzmax - geo.R + sqrt(-1) .* (zmax - geo.z0));
		ul   =  angle(rzmin - geo.R + sqrt(-1) .* (zmin - geo.z0));
		tu   =  abs((acos((rzmax - geo.R) ./ geo.a) - acos(cos(uu))) ./ sin(uu));
		tl   =  abs((acos((rzmin - geo.R) ./ geo.a) - acos(cos(ul))) ./ sin(ul));
		tm   =  (tl + tu) ./ 2;
		d    =   abs(rzmax + rzmin -  2 .* geo.R) ./ 2 ./ geo.a;
		geo.d =  0.6 .* d + 0.4  .* sin(tm);
		R = sepa1.R;
		Z = sepa1.Z;
	else
		% geometrie simple
		geo.a  =  moments.a(k);
		geo.R  =  moments.ra(k);
		geo.z0=  moments.za(k);
		geo.d  =  (moments.rxup(k) + moments.rxdo(k)) ./ 2;
		geo.K  =  (moments.zxup(k) + moments.zxdo(k)) ./ 2;

		% 
		t  = asin(max(0,min(1,geo.d)));
		u  = linspace(0,2.*pi,201);
		vu = ones(size(u));
		R  = geo.R * vu + (geo.a * vu) .* cos(u + t * sin(u));
		Z  = (geo.a .* geo.K) * sin(u) + geo.z0;
	end 

	% correction du b0, calcul sur l'axe geometrique
	geo.b0 = moments.b0(k) .* (1 - geo.a ./ geo.R - moments.delta(k) ./ geo.R);
	
	z0dinput.exp0d.Rsepa(k,:) = R;       % vecteur R des points de la separatrice (m)
	z0dinput.exp0d.Zsepa(k,:) = Z - geo.z0;       % vecteur Z des points de la separtrice (m)
	% recopie dans z0dinput
	noms = fieldnames(geo);
	for l = 1:length(noms)
		z0dinput.geo.(noms{l})(k) = geo.(noms{l});
	end
end

% mise sur les memes angles
R    = z0dinput.exp0d.Rsepa;
Z    = z0dinput.exp0d.Zsepa;
Rc   = z0dinput.geo.R;
Zc   = z0dinput.geo.z0; 
th   = linspace(0,2.*pi,201);
vh   = ones(size(th));
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
	rhoo   = pchip(thk,rhok,th);
	R      = Rc(k)  + rhoo .* cos(th);
	Z      = Zc(k)  + rhoo .* sin(th);
	z0dinput.exp0d.Rsepa(k,:) = R;
	z0dinput.exp0d.Zsepa(k,:) = Z;
end
