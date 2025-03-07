% ZUIEDIT_PARAM_GENE_SAUV_FCT  gestion des uicontrols des parametres relatifs a la sauvegarde des resultats
%--------------------------------------------------------------
% fichier zuiedit_param_gene_sauv_fct.m ->
%
% fonction Matlab 5 :
%	fonction de controle des uicontrols des
%	parametres relatifs a la sauvegarde des resultats
%	des parametres  généraux
%	sous le mode edition du formulaire principal
%  
% syntaxe :  
%   zuiedit_param_gene_sauv_fct(action)
%  
% entrees :  
%	action : tag du uicontrol active
%  
% sorties : 
%  				
% fonction ecrite par C. Passeron, poste 61 19
% version 1.3, du 10/04/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function zuiedit_param_gene_sauv_fct(action)

%  
[hfig,h] = zuiformhandle('ed_param_sauv') ;
if ~ishandle(hfig)
	return
end
hoc = getfield(h,action) ;

switch lower(action)

% Annulation
	case {'annulation','close'}
		zuicloseone(hfig);	
	
% raz
	case {'init','raz'}
		zuiformvisible(hfig) ;
		zuiformreset(hfig) ;
		zuiuploadform(hfig) ;
		zuireset(h.raz) ;
		
% Validation
	case 'validation'
		zuidownloadform(hfig);
		zuiformcache(hfig) ;
		zuisavenonok;
		zuireset(h.validation);

	otherwise
		warning('ation non prise en compte')
	   
end
