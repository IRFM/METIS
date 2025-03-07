%  ZUICLOSE  ferme toutes les fenetres graphiques sauf "direct" et "plot"
%------------------------------------------------------------------------------- 
% fichier : zuiclose.m  
% 
% 
% fonction Matlab 5 : 
% 
% ferme toutes les fenetres graphiques sauf "direct" et "plot"
%  
% syntaxe :  
%   zuiclose
%  
% entrees :  
%  
% sorties :  
%  
% fonction écrite par J-F Artaud , poste 46-78
% version  2.1  du 20/05/2003 
%  
%  
% liste des modifications :  
%  27/09/2001 -> rajout des tests sur las tag pour ne pas fermer la fenetre
%                principale 'zuidirect' et la fenetre 'zdataplot'
%  20/05/2003 -> correction bug clignotement fermeture 
%
%-------------------------------------------------------------------------------  
%  
function zuiclose

fh      =  getappdata(0,'formulaire_handle') ;

if isempty(fh)
   return
end
   
names   =  fieldnames(fh) ;

for k = 1:length(names)
	h= getfield(fh,names{k}) ;
	if ishandle(h)
		tag = get(h,'tag') ;
		if ~strcmp(tag,'direct') & ~strcmp(tag,'zdataplot') ;
   			zuicloseone(h) ;
		else
			zuiformvisible(h) ;
		end
	end
end
   
