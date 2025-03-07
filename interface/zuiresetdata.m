% ZUIRESETDATA remet aux valeurs par defaut un uiconrol
%------------------------------------------------------
% fichier zuiresetdata.m ->  zuiresetdata
%
% fonction Matlab 5 :
%	Cette fonction remet aux valeurs par defaut un uiconrol. 
%
% syntaxe  :
%	zuiresetdata(hc)
%
% entrees :
%	hc       = handle du uicontrol
%
% sortie : aucune
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 08/02/2001.
% 
% liste des modifications : 
%   *  27/07/2001 -> on recupere le style et l'etat par defaut de l'uicontrol
%                   (C.Passeron)
%
%--------------------------------------------------------------
%
function zuiresetdata(hc)

if ~ishandle(hc)
    return
end

% memorisation des valeurs par defaut
defaut  = getappdata(hc,'init_defaut');
if ~isempty(defaut)
	string  = defaut.string;
	value   = defaut.value;
	style   = defaut.style;
	etat    = defaut.etat;
	
	% mise a jour des proprietes
	set(hc,'style',style,'string',string,'value',value, ...
               'enable',etat);
end


