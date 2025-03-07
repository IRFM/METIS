%  INTERSECTION_2_DROITES  courte description  
%------------------------------------------------------------------------------- 
% fichier :  intersection_2_droites.m  ->  intersection_2_droites 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   [in,t,u,ri,zi,re,ze] = intersection_2_droites(ra1,ra2,za1,za2,rb1,rb2,zb1,zb2,tol) 
%  
% entrees :  
%  ra1 = 
%  ra2 = 
%  za1 = 
%  za2 = 
%  rb1 = 
%  rb2 = 
%  zb1 = 
%  zb2 = 
%  tol = 
%  
% sorties :  
%  in = 
%  t  = 
%  u  = 
%  ri = 
%  zi = 
%  re = 
%  ze = 
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
function [in,t,u,ri,zi,re,ze] = intersection_2_droites(ra1,ra2,za1,za2,rb1,rb2,zb1,zb2,tol)

warning off
t = -(zb2 .* rb1 - ra1 .* zb2 + ra1 .* zb1 + za1 .* rb2 - za1 .* rb1 - zb1 .* rb2) ./ (za2 .* rb2 - za2 .* rb1 - za1 .* rb2 + za1 .* rb1 + ra2 .* zb1 - ra2 .* zb2 - ra1 .* zb1 + ra1 .* zb2);
u = (za2 .* ra1 - za2 .* rb1 - ra2 .* za1 + za1 .* rb1 + ra2 .* zb1 - ra1 .* zb1) ./ (za2 .* rb2 - za2 .* rb1 - za1 .* rb2 + za1 .* rb1 + ra2 .* zb1 - ra2 .* zb2 - ra1 .* zb1 + ra1 .* zb2);
warning on

ri  = ra1 + (ra2 -ra1) .* t;
zi  = za1 + (za2 -za1) .* t;
re  = rb1 + (rb2 -rb1) .* u;
ze  = zb1 + (zb2 -zb1) .* u;

in  = (t>=-tol) .* (t<=(1+tol)) .* (u>=-tol) .* (u<=(1+tol)) .* isfinite(t) .* isfinite(u);
