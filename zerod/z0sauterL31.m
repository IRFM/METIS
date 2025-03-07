% fonction L31 Sauter 
function L31 = z0sauterL31(x,te,ne,q,zeff,Raxe,ft,epsi)


ve             = ones(size(x));
vt             = ones(size(te,1),1);
ee             = 1.602176462e-19;
Zib            = zeff;
lne            = 31.3-log(sqrt(ne)./te);
fact           = Raxe ./ (epsi.^(1.5));
nustare        = fact.*6.921e-18.*ne.*lne.*q./(te.^2).*Zib;
X              = ft ./(1 + (1-0.1*ft).*sqrt(nustare) + 0.5*(1-ft).*nustare./Zib);
L31            = (1+(1.4 ./ (Zib+1))).*X - (1.9 ./ (Zib+1)).*X.^2 + (0.3 ./ (Zib+1)).*X.^3 + (0.2 ./ (Zib+1)).*X.^4;
