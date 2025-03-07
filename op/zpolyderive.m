% ZPOLYDERIVE calcul la derivee d'ordre N a partir d'un fit polynomial
%---------------------------------------------------------------------
% fichier zpolyderive.m ->  zpolyderive
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule la derivee d'ordre N a partir d'un fit polynomial.
%  
% syntaxe  :
%  
%     yp=zpolyderive(x,y,o);
%    
% entree :
%
%     x  = variable
%     y  = y(x) [vecteur]
%     o  = ordre de la derivee
%
% sorties :
% 
%     yp =  derivee d'ordre o de y(x)
%     
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 25/10/2000.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function yp=zpolyderive(x,y,o)


n = fix((length(x)/log(length(x))+1)/2)*2;
p  = polyfit(x,y,n);
pp = ((length(p)-o):-1:1).*p(1:(end-o));
yp = polyval(pp,x);