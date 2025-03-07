%  ZLISSAGEBORD  courte description
%------------------------------------------------------------------------------- 
% fichier :  zlissagebord.m  ->  zlissagebord 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   s = zlissagebord(x,e,x0) 
%  
% entrees :  
%  x  = 
%  e  = 
%  x0 = 
%  
% sorties :  
%   s  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  1.9  du  11/03/2002  
%  
% #auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function s = zlissagebord(x,e,x0)

if x0 < 1
   ind = min(find(x >= x0));
   if x0 > 0.9
        fin = length(x);
   else
	     fin = min(find(x > 0.95)); 
   end	
   me  = mean(e(ind:fin));

   ee  = spline([x(ind-2),x(ind),1],[e(ind-2),e(ind),me],x);
   ee  = ee .* (ee > 0) + (ee <= 0) .* me;

   s   = e;
   s(ind:end)  = ee(ind:end);
else
   s   = e;
   
end
