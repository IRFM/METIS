% ZCENTRE calcul le point central d'un profil
%-------------------------------------------
% fichier zcentre.m ->  zcentre2
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule le point central d'un profil
% en cas de NaN, inf ou 0 en ce point. Elle suppose 
% une derivee nulle au centre. Elle utilise une interpolation 
% par spline cubic
%  
% syntaxe  :
%  
%     fp = zcentre2(f);
%    
% entree :
%
%     f    =  vecteur profil de au moins 4 points
%
% sorties :
% 
%     fp   =  vecteur profil prolonge en fp(1)
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.0, du 05/03/2003.
% 
% 
% liste des modifications : 
% 05/03/2003 : cas ou f est petit au centre mais non nulle
%--------------------------------------------------------------
%
function fp = zcentre2(f)
%
% test f(1) < .1 * f(2)
%
if ~isfinite(f(1))|(100*f(1)/f(2) <= 10)
   fp = f;
   fp(1) = (61/46) .* f(2) - (9/23) .* f(3) + (3/46) .* f(4); 
else
   fp = f;
end

