% ZUIEDIT_PARAM_MODEXT_CTRL control uicontrols formulaire des modules externes, parametres généraux
%--------------------------------------------------------------
% 	fichier zuiedit_param_modext_ctrl.m  
%
% fonction Matlab 5 :
%	fonction de control des donnes des uicontrols du formulaire
%	des modules externes , parametres  généraux
%	sous le mode edition du formulaire principal
%
%
% syntaxe  :
%	zuiedit_parma_modext(action) 
%
% entrees
%
% sorties :
%  action       =  tag du uicontrol active
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 03/09/2003.
% 
% liste des modifications : 
%   * 03/09/2003 -> ajout gestion module ripple
%
%   * 12/09/2001  -> changement de la liste des options pour les glacons  (J-F Artaud)
%
%--------------------------------------------------------------
function zuiedit_param_modext_ctrl(action)

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_modext');
% information pour l'assistant
zuicr(hfig,action) ;
try
  hoc = getfield(h,action) ;
catch
  hoc = [];
end
  
if ~strcmp(get(hoc,'style'),'popupmenu')
        % pas de control de valeur (aucun risque d'erreur)
else
        data = zuidata(hoc);
	value = zuidata(hoc)
	
end

