function eta = zeta0(te,ne,q,zeff,r,R,ft,epsi,eta_breakdown)



% formule de Sauter
if nargin < 9
  ne             = max(ne,1e13);
  te             = max(te,13.6);
end
lne            = 31.3-log(sqrt(ne)./(te));
fact           = R ./ (epsi.^(1.5));
nustare        = fact.*6.921e-18.*ne.*lne.*q./(te.^2).*zeff;    
NZ             = 0.58 + 0.74 ./ (0.76+zeff);
sigsptz        = 1.9012e4*te.^(1.5)./zeff./lne./NZ;
X5             = ft ./(1 + (0.55-0.1*ft).*sqrt(nustare) + 0.45*(1-ft).*nustare./(zeff.^1.5));

sigma          = sigsptz.*(1-(1+0.36./zeff).*X5 + 0.59 ./ zeff.*X5.^2 - 0.23./zeff.*X5.^3);
eta            = 1./real(sigma);

if nargin >= 9
  eta = eta + eta_breakdown;
end