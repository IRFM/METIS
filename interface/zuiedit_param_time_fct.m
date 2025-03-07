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

function zuiedit_param_time_fct(action)

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_time');
 
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)

case 'radio_exec'
	if ishandle(hfig)
		% disp('formulaire parametre/generaux/execution')
		zuiedit_param_gene_funf('zexec','param.gene','Execution','','zuiedit_param_gene_exec_ctrl');
		zuireset(h.radio_exec)
	end

case 'radio_split'
        if ishandle(hfig)
                zuiedit_param_split ;
                zuireset(h.radio_split)
        end

case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig);
		zuireset(h.btn_quit);
	end

otherwise

end

