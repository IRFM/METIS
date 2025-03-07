% calcul un separatrice pour un equilique en tenant compte 
% des points de controle. 
% Utilise une description avec 2 courbes de Bezier.
function [R,Z] = fitsepa(Rc,Zc,nbp,uniform)

if nargin < 3	
	nbp = length(Rc);
elseif isempty(nbp)
	nbp = length(Rc);
end
if nargin < 4
	uniform = 1;
end

% point externe
indext = find(Rc == max(Rc),1);
z0     = Zc(indext);
r0     = (min(Rc) + max(Rc)) ./ 2;

% decoupage en partie HFS et LFS
indmax = find(Zc  == max(Zc),1);
indmin = find(Zc  == min(Zc),1);

% points HFS 
ind_hfs = find(((Zc>= z0) & (Rc <=  Rc(indmax))) | ((Zc<= z0) & (Rc <=  Rc(indmin))));
if ~any(indmin == ind_hfs)
	indhfs(end+1) = indmin;
end
if ~any(indmax == ind_hfs)
	indhfs(end+1) = indmax;
end
bph     = fitbezier(Zc(ind_hfs),Rc(ind_hfs));

%  Zch = Zc(ind_hfs);
%  Rch = Rc(ind_hfs);

% points LFS 
ind_lfs = find(((Zc>= z0) & (Rc >=  Rc(indmax))) | ((Zc<= z0) & (Rc >=  Rc(indmin))));
if ~any(indmin == ind_lfs)
	indlfs(end+1) = indmin;
end
if ~any(indmax == ind_lfs)
	indlfs(end+1) = indmax;
end
bpl     = fitbezier(Zc(ind_lfs),Rc(ind_lfs));

%  Zcl = Zc(ind_lfs);
%  Rcl = Rc(ind_lfs);

% clacul de la separatrice
Zh = Zc(ind_hfs)';
Rh = bezierval(bph,Zh);
Zl = Zc(ind_lfs)';
Rl = bezierval(bpl,Zl);

% sortie unifiee
Z = cat(1,Zh,Zl)';
R = cat(1,Rh,Rl)';

% uniformisation
if uniform
	c   = (R - r0) + sqrt(-1) .* (Z - z0);
	th  = unwrap(angle(c),[],2);
	rho = abs(c);
	th(th<0) = th(th<0) + 2 .* pi;
	[th,ind] = sort(th,2);
	rho       = rho(ind);
	rhol = cat(2,rho,rho,rho);
	thl = cat(2,th -2.*pi,th,th+2.*pi);
	indnok = find(any(diff(thl,1,2)<=0,1));
	while(~isempty(indnok))
		thl(indnok) =[];
		rhol(indnok)  = [];
		indnok = find(any(diff(thl,1,2)<=0,1));
	end
	tu   = linspace(0,2.*pi,nbp);
	ru   = pchip(thl,rhol,tu);
	R    = r0 + ru .* cos(tu);
	Z    = z0 + ru .* sin(tu);	
else	
	
	c   = (R - r0) + sqrt(-1) .* (Z - z0);
	th  = unwrap(angle(c),[],2);
	rho = abs(c);
	th(th<0) = th(th<0) + 2 .* pi;
	[th,ind] = sort(th,2);
	R  = R(ind);
	Z  = Z(ind);
end








