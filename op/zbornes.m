% ZBORNES assure qu'une matrice a des valeurs entre deux bornes
%--------------------------------------------------------------
% fichier zbornes.m ->  zbornes
%
%
% fonction Matlab 5 :
%
% Cette fonction bornes la valeur des donnees et 
% retire les Inf et les NaN.
%  
% syntaxe  :
%  
%     s = zbornes(s,smin,smax,repli);
%    
% entree :
%
%     s     = donnees
%     smin  = valeur minimum que peu prendre s (si smin = -inf,
%             s n'a pas de valeur minimum)
%     smax  = valeur maximum que peu prendre s (si smin = inf,
%             s n'a pas de valeur maximum)
%    
%     repli = valeur de remplacement pour les point ou s est infini
%             ou NaN          
%
% sorties :
% 
%     s     =  donnees bornees
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 22/02/2001.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function  s = zbornes(s,smin,smax,repli)

if isempty(s)
	return
end
ind =find(~isfinite(s));
if ~isempty(ind)
	s(ind) = repli .* ones(1,length(ind));
end
if isfinite(smin)
	s = s .* (s >= smin) + smin .* (s < smin);
end
if isfinite(smax)
	s = s .* (s <= smax) + smin .* (s > smax);
end


