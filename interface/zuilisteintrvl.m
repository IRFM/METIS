%  ZUILISTEINTRVL  mis a jour de la liste des intervalles
%------------------------------------------------------------------------------- 
% fichier : zuilisteintrvl.m  -> zuilisteintrvl
% 
% fonction Matlab 5 : 
%	mis a jour de la liste des intervalles
%  
% syntaxe :  
%	zuilisteintrvl(tag)
%  
% entrees :  
%	tag : tag du formulaire
%  
% sorties :  
%  
% fonction écrite par ?? poste XX-XX  
% version  1.7  du  29/09/2001  
%  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function zuilisteintrvl(tag)

[hform,h] = zuiformhandle(tag) ;
name=getappdata(h.list_intervalle,'libel') ;
xd=getappdata(h.list_intervalle,'xd') ;
xf=getappdata(h.list_intervalle,'xf') ;

dd=[];
for i=1:length(xd)
	aa = sprintf('#%-s \t [%-6.2g - %-6.2g]',name(i,:),xd(i),xf(i)) ;
	dd = strvcat(dd,aa) ;
end
%set(h.list_intervalle,'string',dd,'min',1,'max',length(xd)) ;
set(h.list_intervalle,'string',dd,'min',1,'max',2) ;
