%  ZLINLOG  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zlinlog.m  ->  zlinlog 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   zlinlog 
%  
% entrees :  
%  
%  
% sorties :  
%  
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
function zlinlog

   ho  = gco;
	hf  = gcf;
	hal = findobj(hf,'type','axes');
	if get(ho,'value') == 1
	    set(hal,'yscale','log');
	else
	    set(hal,'yscale','lin');
	end
