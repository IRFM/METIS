% ZCENTRE calcul le point central d'un profil
%-------------------------------------------
% fichier zcentre.m ->  zcentre
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
%     fp = zcentre(f);
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
% version 1.0, du 23/08/2000.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function fp = zcentre(f)

if ~isfinite(f(1))|(f(1) == 0)
   fp = f;
   fp(1) = (61/46) .* f(2) - (9/23) .* f(3) + (3/46) .* f(4); 
else
   fp = f;
end

