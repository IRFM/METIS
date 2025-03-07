%  FCTVAL4  courte description  
%------------------------------------------------------------------------------- 
% fichier :  fctval4.m  ->  fctval4 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   ind = fctval4(val) 
%  
% entrees :  
%  val = 
%  
% sorties :  
%   ind  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  4.1  du  08/04/2008  
%  
%*@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function ind = fctval4(val)

if all(~isfinite(val))
        ind = 0 ;
elseif all(val==0)
        ind = 0 ;
elseif all(val==1)
        ind = 1 ;
elseif all(val==2)
        ind = 2 ;
elseif all(val==3)
        ind = 3 ;
else
        ind = 4 ;
end

