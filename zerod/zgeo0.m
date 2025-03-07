% cette fonction calcule le volume du plasma sa section et la surface externe
function [vps,sps,sexts,peris,geo,xpoint,R,Z] = zgeo0(geo,Rsepa,Zsepa,evolution)


% mode de fonctionnement
sepamode = 0;
if nargin  > 2
	if ~isempty(Rsepa) && ~isempty(Zsepa)
		if all(isfinite(Rsepa(:))) && all(isfinite(Zsepa(:)))
			sepamode = 1;
		else
			fprintf('S');
		end
	end
end

% calcul "a peu pres" si non geo.dsponible
vp0   = 2 .* pi .^ 2 .* geo.R .* geo.a .^ 2 .* geo.K;
sp0   = pi .* geo.a .^ 2 .* geo.K;
% cette formule est fausse, mais c'est celle utiliser pour le scaling de rayonnement
sext0 =  2 .* pi .^ 2 .* geo.R .* geo.a .* sqrt((1 + geo.K .^ 2));
peri0 =  sqrt(2) .* pi .* geo.a .* sqrt((1 + geo.K .^ 2));

if sepamode == 1
	% la separatrice doit etre donnee deja centree	en Z ; il ne faut pas soustraire z0.
	% on preserve R*b0
	rb0  = geo.b0 .* geo.R;
	% calcul des moments pour assurer la coherence dans tous les cas avec les moments
	% centre pour angle d'integration
	rc = mean(Rsepa,2);
	zc = mean(Zsepa,2);
	vc = ones(1,size(Rsepa,2));
	%uc = unwrap(angle((Rsepa-rc*vc) + sqrt(-1) .* (Zsepa  -zc*vc)));
	uc = unwrap(angle((Rsepa-rc*vc) + sqrt(-1) .* Zsepa ));
	uc    = uc .* (uc >0) + (uc + 2*pi) .* (uc<= 0);
	uc(:,1)   = uc(:,end) + 2 .* pi;
	xu    = linspace(0,1,length(vc));
	dRdx  = pdederive(xu,Rsepa,2,2,2,1);
	dZdx  = pdederive(xu,Zsepa,2,2,2,1);
	% calcul de R0 et Z0
	maskrmax  = (Rsepa == (max(Rsepa,[],2) * vc));
	% recalcul des parametres sur le vecteur final
	rmin  = min(Rsepa,[],2);
	rmax  = max(Rsepa,[],2);
	geo.a = 0.5 .* (rmax - rmin);
	geo.R = 0.5 .* (rmax + rmin);
	zmin  = min(Zsepa,[],2);
	zmax  = max(Zsepa,[],2);
	% geo.z0   = (zmax + zmin + sum(Zsepa .* maskrmax,2) ./ sum(maskrmax,2)) ./ 3;
	geo.K    = (abs(trapz(xu,Zsepa .*  dRdx,2) ./ pi ./ geo.a) + (zmax - zmin)) ./ 3 ./ geo.a;
	
	rzmax = geo.R;
	rzmin = geo.R;
	for k = 1:size(Zsepa,1)
		rzmax(k) = Rsepa(k,min(find(Zsepa(k,:) == zmax(k))));
		rzmin(k) = Rsepa(k,min(find(Zsepa(k,:) == zmin(k))));
	end
	uu   =  angle(rzmax - geo.R + sqrt(-1) .* zmax);
	ul   =  angle(rzmin - geo.R + sqrt(-1) .* zmin);
	tu   =  abs((acos((rzmax - geo.R) ./ geo.a) - acos(cos(uu))) ./ sin(uu));
	tl   =  abs((acos((rzmin - geo.R) ./ geo.a) - acos(cos(ul))) ./ sin(ul));
	tm   =  (tl + tu) ./ 2;
	d    =   abs(rzmax + rzmin -  2 .* geo.R) ./ 2 ./ geo.a;
	geo.d =  0.6 .* d + 0.4  .* sin(tm);

	% recalcul du champs avec rb0
	geo.b0  = rb0 ./ geo.R;
	
	% on integre su theta directement
	vu = ones(1,size(Rsepa,2));
	vt = ones(size(Rsepa,1),1);
	rho = abs((Rsepa - geo.R * vu) + i .* (Zsepa));
	uu  = unwrap(angle((Rsepa - geo.R * vu) + i .* (Zsepa)),[],2);
	u   = mean(uu,1);
	if (max(std(uu,1)) > 0.001) && (evolution == 1)
				u = uu(end,:);
	end
	R  = Rsepa;
	Z  = Zsepa;

    
    	if (u(1) == u(2))
          	u = u(2:end);
          	R = R(:,2:end);
          	Z = Z(:,2:end);
    	end
    	if (u(end-1) == u(end))
        	u = u(1:(end-1));
        	R = R(:,1:(end-1));
        	Z = Z(:,1:(end-1));
    	end    

	dldu  = sqrt(pdederive(u,R,2,2,2,1) .^2 + pdederive(u,Z,2,2,2,1) .^2 );
	dS = -Z .* pdederive(u,R,2,2,2,1);
	vp = abs(trapz(u,2.*pi.*R.*dS,2));
	sp = abs(trapz(u,dS,2));
	sext = abs(trapz(u,2.*pi.* R .* dldu,2));

	% calcul du premietre
	peri = abs(trapz(u,dldu,2));


	if (max(std(uu,1)) > 0.001) && (evolution == 0)
        R  = Rsepa;
        Z  = Zsepa;
        uu  = unwrap(angle((Rsepa - geo.R * vu) + i .* (Zsepa)),[],2);
		for k=1:size(uu,1)
			u = uu(k,:);
			Rk = R(k,:);
			Zk = Z(k,:);
			if (u(1) == u(2))
				u  = u(2:end);
				Rk = Rk(2:end);
				Zk = Zk(2:end);
			end
			if (u(end-1) == u(end))
				u  = u(1:(end-1));
				Rk = Rk(1:(end-1));
				Zk = Zk(1:(end-1));
			end    
		
			dldu  = sqrt(pdederive(u,Rk,2,2,2,1) .^2 + pdederive(u,Zk,2,2,2,1) .^2 );
			dS = -Zk .* pdederive(u,Rk,2,2,2,1);
			vp(k) = abs(trapz(u,2.*pi.*Rk.*dS,2));
			sp(k) = abs(trapz(u,dS,2));
			sext(k) = abs(trapz(u,2.*pi.* Rk .* dldu,2));
		
			% calcul du premietre
			peri(k) = abs(trapz(u,dldu,2));

		end
	end

	
	% un seul affichage par run
	if evolution  ~= 1
		fprintf('Metis: using separatrix given by points (R,Z) :')
	end
	
	
%  	% detection point x
%  	vu = ones(1,size(R,2));
%  	vt = ones(size(R,1),1);
%  	t  = asin(max(0,min(1,geo.d)));
%  	Rm  = geo.R *vu + (geo.a * vu) .* cos(vt * u + t * sin(u));
%  	Zm  = (geo.a .* geo.K) * sin(u);
	for k = 1:size(Zsepa,1)
		indh     = find(Zsepa(k,:) == zmax(k),1);
		indh     = cat(2,indh - 3, indh - 2,indh - 1,indh,indh + 1,indh + 2,indh + 3); 
		indh     = mod(indh-1,size(Zsepa,2))+1;
		rh       = Rsepa(k,indh);
		zh       = Zsepa(k,indh);
		ph       = polyfit(rh,zh,2);
		eh       = sqrt(mean((zh - polyval(ph,rh)).^ 2)) ./ (max(rh) - min(rh));
		indl     = find(Zsepa(k,:) == zmin(k),1);
		indl     = cat(2,indl - 3, indl - 2,indl - 1,indl,indl + 1,indl + 2,indl + 3);
		indl     = mod(indl-1,size(Zsepa,2))+1;
		rl       = Rsepa(k,indl);
		zl       = Zsepa(k,indl);
		pl       = polyfit(rl,zl,2);
		el       = sqrt(mean((zl - polyval(pl,rl)).^ 2)) ./ (max(rl) - min(rl));
		xpoint(k,1) = max(el > 2e-2 , eh > 2e-2);
	end
else
	% on integre su theta directement
	t  = asin(max(0,min(1,geo.d)));
	u  = linspace(0,2.*pi,201);
	vu = ones(size(u));
	vt = ones(size(t));
	R  = geo.R *vu + (geo.a * vu) .* cos(vt * u + t * sin(u));
	Z  = (geo.a .* geo.K) * sin(u);
	dldu  = sqrt(pdederive(u,R,2,2,2,1) .^2 + pdederive(u,Z,2,2,2,1) .^2 );
	dS = -Z .* pdederive(u,R,2,2,2,1);
	vp = trapz(u,2.*pi.*R.*dS,2);
	sp = trapz(u,dS,2);
	sext = trapz(u,2.*pi.* R .* dldu,2);
	
	% calcul du premietre
	peri = trapz(u,dldu,2);
	% to have fixed number of output
	xpoint = zeros(size(geo.R));
end

% securite
vps = vp;
vps(~isfinite(vps))  = vp0(~isfinite(vps));
vps(vps == 0)        = vp0(vps ==0);
sps = sp;
sps(~isfinite(sps)) = sp0(~isfinite(sps));
sps(sps == 0)       = sp0(sps ==0);
sexts = sext;
sexts(~isfinite(sexts)) = sext0(~isfinite(sexts));
sexts(sexts==0)         = sext0(sexts==0);
peris = peri;
peris(~isfinite(peris)) = peri0(~isfinite(peris));
peris(peris==0)         = peri0(peris==0);
vps    = max(vps,1e-3);
sps    = max(sps,1e-3);
sexts  = max(sexts,1e-3);
peris  = max(peris,1e-3);
