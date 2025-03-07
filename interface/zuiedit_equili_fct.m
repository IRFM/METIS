% ZUIEDIT_PARAM_GENE_FCT gestion callback formulaire des parametres g��aux
%--------------------------------------------------------------
% fichier : zuiedit_param_gene_fct
%
% fonction Matlab 5 :
%	fonction de gestion des callback associes a 
%	chaque uicontrol pour le formulaire des parametres g��aux
% 	sous le mode edition du formulaire principal
%
% syntaxe :
% 		zuiedit_param_gene_fct(action)
%
% entrees :
%	action : tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.7, du 03/10/2001. 
% 
% liste des modifications : 
% * 03/10/2001 -> ajout de la configuration deu profiler (J-F Artaud)
%
%--------------------------------------------------------------

function zuiedit_equili_fct(action)

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_equili');
 
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)


case 'radio_calcequi'
        if ishandle(hfig)
                %hout = zuiedit_data_geom ;
                zuiedit_calcequi ;
                zuireset(h.radio_calcequi)
        end

case 'radio_equi'
        if ishandle(hfig)
                zinitequi ;
                zuireset(h.radio_equi)
        end

case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig);
		zuireset(h.btn_quit);
	end

otherwise

end

