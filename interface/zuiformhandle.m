% ZUIFORMHANDLE retourne la structure des handles d'un formulaire
%-----------------------------------------------------------------
% fichier zuiformhandle.m ->  zuiformhandle
%
%
% fonction Matlab 5 :
% 
% Cette fonction retourne la structure des handles d'un formulaire. 
%
% syntaxe  :
%  
%  [hform,hui] = zuiformhandle(tag);
%
% entrees :
% 
%  tag     =  tag du formulaire
%
% sortie : 
% 
%  hform   = handle du formulaire
%  hui     = structure contenant les handles des uicontrol
%  
% format :      
%
%  hui.<tag_du_uicontrol>
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2 du 10/06/2004.
% 
% 
% liste des modifications : 
%
%   * 10/06/2004 ->  gestion des handles caches
%--------------------------------------------------------------
%
function [hform,hui] = zuiformhandle(tag)

hform  = [];
hui    = [];

if nargin == 0
	warning('zuiformhandle: pas d''objet')
   return;
elseif isempty(tag)
	warning('zuiformhandle: pas de tag')
	return;
else
        etatmem = get(0,'ShowHiddenHandles');
        set(0,'ShowHiddenHandles','on');
	    hform = findobj(0,'type','figure','tag',tag);
        hdform = double(hform);
        [void,indmax] = max(hdform);
        hform = hform(indmax); 
        set(0,'ShowHiddenHandles',etatmem);

end
if ishandle(hform)
	hui=getappdata(hform,'zhandle');
end
