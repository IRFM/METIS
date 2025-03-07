%  ZID3BESSEL  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zid3bessel.m  ->  zid3bessel 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   [ne,chi2] = zid3bessel(corde,T,x) 
%  
% entrees :  
%  corde = 
%  T     = 
%  x     = 
%  
% sorties :  
%  ne   = 
%  chi2 = 
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
function [ne,chi2] = zid3bessel(corde,T,x)


xx = x ./ max(x);
zz = [2.40483,3.83171,5.13562];
b0 = besselj(0,zz(1) .* x)';
b1 = besselj(1,zz(2) .* x)';
b2 = besselj(2,zz(3) .* x)';
cc = ones(size(x))';
a0 = T * b0;
a1 = T * b1;
a2 = T * b2;
a3 = T * cc;

M  = cat(2,a0,a1,a2,a3); 
coef = M \ corde;
ne   = coef(1) .* b0 +  coef(2) .* b1 +  coef(3) .* b2 +coef(4);
chi2 = sum((M*coef - corde) .^ 2);

