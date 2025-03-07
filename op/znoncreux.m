%  ZNONCREUX  courte description  
%------------------------------------------------------------------------------- 
% fichier :  znoncreux.m  ->  znoncreux 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   s =znoncreux(x,e) 
%  
% entrees :  
%  x = 
%  e = 
%  
% sorties :  
%   s  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  1.9  du  11/03/2002  
%  
%*#auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function s =znoncreux(x,e)


dedx   = pdederive(x,e,0,2,2,1);
if any(dedx >= 0)
	ind              = min(find(dedx < 0));
	if ~isempty(ind)
	  dedx(1:ind)    = dedx(ind) .* ones(1,ind);
	end
end
s     = cumtrapz(x,dedx);
s     = abs(s - s(end) + e(end) .* (e(end)>0));
