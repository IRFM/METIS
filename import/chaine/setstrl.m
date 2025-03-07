%----------------------------------------------------------------
% fonction setstrl.m
% 
%  ecrit par J-F Artaud, poste 46.78   -    derniere mise a jour le 05/10/94 - version 2
%
%  formate une chaine a une longueur l et complete avec des blanc si necessaire
%    (coupe a droite)
%
%   syntaxe : 
%
%                             ss=setstrl(se,l);
%   avec : 
%
%                           se = chaine de caractere
%
%                           l = longuer souhiatee
%
%                           ss = chaine de caractere formatee
%
%------------------------------------------------------------------
%
function ss=setstrl(se,l)
       
       if length(se)<l,
          comp=char(ones(1,l-length(se))*32);
          ss=[se,comp];
       else
          ss=se(1,1:l);
       end
