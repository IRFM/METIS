% ZUIFORMRESET		reinitialise un formulaire aux valeurs par defaut
%------------------------------------------------------
% fichier zuiformreset.m ->  zuiformreset
%
% fonction Matlab 5 :
%	Cette fonction remet le formulaire aux valeurs par defaut 
%	(c'est a dire aux valeurs qu'il avait a la creation). 
%
% syntaxe  :
%	zuiformreset(hfig)
%
% entrees :
%	hfig = handle de la fenetre formulaire
%
% sortie : aucune
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 07/02/2001.
% 
% liste des modifications : 
%	*  27/07/2001 -> on recupere le style et l'etat par defaut des uicontrols (C.Passeron)
%
%-----------------------------------------------------------------------------
%
function zuiformreset(hfig)

% recherche des uicontrols
hui = findobj(hfig,'type','uicontrol');
% boucle sur les objets
for k =1:length(hui)
	hc  = hui(k);

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
end

	
	
