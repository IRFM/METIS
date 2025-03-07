%  ZUILISTEPOINT  mis a jour de la liste des coordonnees des points dans la liste des intervalles
%------------------------------------------------------------------------------- 
% fichier : zuilistepoint.m
% 
% fonction Matlab 5 : 
%	mis a jour de la liste des coordonnees des points dans la intervalles
%  
% syntaxe :  
%	zuilistepoint(tag)
%  
% entrees :  
%	tag : tag du formulaire
%  
% sorties :  
%   
%  
% fonction écrite par J-F Artaud , poste 61-19
% version  1.7  du  29/09/2001  
%  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function zuilistepoint(tag)

[hform,hui] = zuiformhandle(tag);
multi=getappdata(hform,'multi_mem');
x=get(hui.cons,'xdata');
y=get(hui.cons,'ydata');
ind = 1:length(x);
dd=sprintf('#%-3d  [%-6.5g, %-6.5g\t]|',[ind;x;y]);
set(hui.liste_valeurs,'string',dd,'min',1,'max',length(x));
