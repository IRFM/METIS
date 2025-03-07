% fonction bootstrap Sauter corrigee (j*b)
function jboot = zsauter0d(x,te,ti,ne,ni,q,zeff,zion,Raxe,ft,epsi,dpsidx,F,modeh,fgg)


% le facteur est un fit avec NClass (l'integrale n'est pas precise au bord)
if nargin < 15
    fgg = 1;
end


ve   = ones(size(x));
vt   = ones(size(te,1),1);
ee           =   1.602176462e-19;

warning off

Zib            = zeff;

lni            = 30.0-log(zion.^3 .* sqrt(ni) ./ (ti.^(1.5)));
lne            = 31.3-log(sqrt(ne)./te);
fact           = Raxe ./ (epsi.^(1.5));
nustare        = fact.*6.921e-18.*ne.*lne.*q./(te.^2).*Zib;
nustari        = fact.*4.90e-18.*ni.*lni.*q./(ti.^2).*(zion.^4);
X              = ft ./(1 + (1-0.1*ft).*sqrt(nustare) + 0.5*(1-ft).*nustare./Zib);
L31            = (1+(1.4 ./ (Zib+1))).*X - (1.9 ./ (Zib+1)).*X.^2 + (0.3 ./ (Zib+1)).*X.^3 + (0.2 ./ (Zib+1)).*X.^4;
X2             = ft ./(1 + 0.26*(1-ft).*sqrt(nustare) + 0.18*(1-0.37*ft).*nustare./sqrt(Zib));
X3             = ft ./(1 + (1+0.6*ft).*sqrt(nustare) + 0.85*(1-0.37*ft).*nustare.*(1+Zib));
L32_1          = ((0.05 + 0.62*Zib)./Zib./(1+0.44*Zib)).*(X2-X2.^4) + ...
                 (1 ./ (1+0.22*Zib)).*(X2.^2 - X2.^4 - 1.2*(X2.^3 - X2.^4)) + ...
                 (1.2 ./ (1 + 0.5*Zib)).*X2.^4;
L32_2          = -((0.56 + 1.93*Zib)./Zib./(1+0.44*Zib)).*(X3-X3.^4) + ...
                 (4.95 ./ (1+2.48*Zib)).*(X3.^2 - X3.^4 - 0.55*(X3.^3 - X3.^4)) - ...
                 (1.2 ./ (1 + 0.5*Zib)).*X3.^4;
L32            = L32_1 + L32_2;
X4             = ft ./(1 + (1 - 0.1*ft).*sqrt(nustare) + 0.5*(1-0.5*ft).*nustare./Zib);
L34            = (1+(1.4 ./ (Zib+1))).*X4 - ...
                 (1.9 ./ (Zib+1)).*X4.^2  + ...
                 (0.3 ./ (Zib+1)).*X4.^3  + ...
                 (0.2 ./ (Zib+1)).*X4.^4;


alpha0         = -1.17*(1 - ft)./(1 - 0.22*ft -0.19*ft.^2);
alpha          = ((alpha0 + 0.25*(1 - ft.^2) .* sqrt(nustari))./(1 + 0.5*sqrt(nustari)) + 0.315*nustari.^2 .* ft.^6) ./ ...
                 (1 + 0.15*nustari.^2 .* ft.^6);


xTe1           = te/1e3;
xne1           = ne;
xni1           = ni;
xTi1           = ti/1e3;
xpe1           = te.*ne/1e3;
xpi1           = ti.*ni/1e3;
xp1            = xpe1 + xpi1;
rpe            = xpe1 ./ max(xp1,1);
dpdpsi1        = pdederive(x,xp1,0,2,2,1);
dTedpsi1       = pdederive(x,xTe1,0,2,2,1);
dTidpsi1       = pdederive(x,xTi1,0,2,2,1);
% correction piedestal (derivee 2 points + air du rectangle / aire triangle)
indh = find(modeh);
if ~isempty(indh) && (fgg ~= 0)
	dx            = x(end) - x(end-1);
	dpdpsi1(indh,end-1)  = fgg .* (xp1(indh,end) - xp1(indh,end-1)) ./ dx;
	dTedpsi1(indh,end-1) = fgg .* (xTe1(indh,end) - xTe1(indh,end-1)) ./ dx;
	dTidpsi1(indh,end-1) = fgg .* (xTi1(indh,end) - xTi1(indh,end-1)) ./ dx;
end

Isau2psi       = L31./xpe1.*dpdpsi1 +L32./xTe1.*dTedpsi1+L34.*alpha.*(1-rpe) ./ rpe./xTi1.*dTidpsi1;
Isau2psi       = -xpe1.* F .*Isau2psi.*ee.*1e3;

jboot          = -Isau2psi./dpsidx;
warning on

% les gradients s'annulent au centre
jboot(:,1)     = 0;