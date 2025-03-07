% fonction de calcul d'un profil a 3 parametres
function jf = z0j3(x,j0,xh,jh)


%blindage
xh = max(0.05,xh) - 0.025;

ve  = ones(size(x));
vt  = ones(size(j0));
lo     = 0.05;
fo     = exp(-(vt * x - xh *ve) .^ 2 ./ lo);
oo     = exp(-(xh *ve) .^ 2 ./ lo);
do     = max(exp(0) - oo,exp(0) - exp(-(vt * ve - xh *ve) .^ 2 ./ lo)) ;
jf     = j0 * (1 - x)  + (fo - oo) ./ do .* (jh * ve - (j0 .* (1 - xh)) * ve ); 
jf     = jf .* ((vt * x) <= (xh * ve)) + max(0,jf) .*  ((vt * x) > (xh * ve));
jf(~isfinite(jf)) = eps;
jf                 = jf + eps; 
