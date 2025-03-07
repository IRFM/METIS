%  FSHEAR  courte description  
%------------------------------------------------------------------------------- 
% fichier :  fshear.m  ->  fshear 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   fv = fshear(x,x0,finf,fmin) 
%  
% entrees :  
%  x    = 
%  x0   = 
%  finf = 
%  fmin = 
%  
% sorties :  
%   fv  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  18/11/2005  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function fv = fshear(x,x0,finf,fmin)
fv = 2.*x0.*(finf-fmin).*exp(2.^(1./2))./(-2+2.^(1./2)).*(1-1./2.*x0.*2.^(1./2)./x).*exp(-x./x0.*2.^(1./2))./x+finf;
