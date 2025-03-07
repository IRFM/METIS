%  Z0INTERP_CONS  courte description  
%------------------------------------------------------------------------------- 
% fichier :  z0interp_cons.m  ->  z0interp_cons 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   xout = z0interp_cons(tin,xin,tout) 
%  
% entrees :  
%  tin  = 
%  xin  = 
%  tout = 
%  
% sorties :  
%   xout  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  4.0  du  08/03/2007  
%  
%*@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function xout = z0interp_cons(tin,xin,tout)

xout = xin(end) .* ones(size(tout));
for k=length(tin):-1:2
	xout(tout <= tin(k)) = xin(k-1);
end