% calcul de la fonction de distribution de NBI
% S0 = source en particules/s
% taus = temps de ralentissement (s)
% v0   = vitesse d'injection
% mu0  = pitch angle initial
% vc   = vitesse critique
% vg   = vitesse gamma de Stix
% v    = vitesses ou il faut evaluer la fonction de distribution
% mu   = valeur de v// / v ou il faut evaluer la fonction de distribution
function [v,mu,fvmu,fn,pn] = z0nbi_fvmu(s0,taus,v0,mu0,vc,vg,vmax)

% polynome de Legendre (100)
%npl = 5;           % nombre de polynome de legendre (5 -> jusqu'en mu ^4)
%nmu = 10 .* (npl - 1) + 1;          % nombre de points en mu
%nv  = 201;        % nombre de points en vitesse 
npl = 3;           % nombre de polynome de legendre (5 -> jusqu'en mu ^4)
nmu = 10 .* (npl - 1) + 1;          % nombre de points en mu
nv  = 51;        % nombre de points en vitesse 
mu = linspace(-1,1,nmu);
if nargin <= 6
  vmax = v0;
end
v  = linspace(0,vmax,nv)';
pn = ones(npl,nmu);
fn = zeros(nv,npl);
pn(2,:) = mu;
%
% Calcul par reccurence des polynomes d'ordre superieur
%
for k=2:(npl - 1)
	pn(k+1,:)=((2*k-1)*mu.*pn(k,:)-(k-1)*pn(k-1,:))/k;
end
% normalisation de la fonction
an = ones(npl,1);
for k =1:npl
	an(k) = (2 .* (k-1)  + 1) ./ 2 .* pchip(mu,pn(k,:),mu0); 
end

% vecteur utiles
vv = ones(size(v));
vmu = ones(size(mu));
van  = ones(size(an));

% contribution Fn
n  = 0:(npl-1);
fn = s0 .* taus .* (vv * an') ./ 2 ./ pi ./ ((v * van') .^ 3 + vc .^ 3) .* ...
		     ((v * van')  .^ 3 ./ v0 .^ 3 .* (v0 .^ 3 + vc .^ 3) ./ ((v * van') .^ 3 + vc .^ 3)) .^ ...
		     (vv * (n .* (n + 1) .* vg .^ 3 ./ 3 ./ vc .^ 3)) .* ((v * van') <= v0); 

% fonction de distribution
fvmu   = fn * pn;


