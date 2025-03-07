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

function zuiedit_param_gene_fct(action)

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_gene');
 
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)

case 'radio_conv'
	if ishandle(hfig)
		% disp('formulaire parametre/generaux/convergence')
		zuiedit_param_gene_funf('zconv','param.gene','Convergence');
		zuireset(h.radio_conv)
	end

case 'radio_multi'
	if ishandle(hfig)
		% disp('formulaire parametre/generaux/convergence')
		zuiedit_param_gene_funf('zmulti','param.cons.neomulti','Multiplicator');
		zuireset(h.radio_multi)
	end

case 'radio_sauv'
	if ishandle(hfig)
		% disp('formulaire parametre/generaux/sauvegarde')
%		zuiedit_param_gene_funf('zsauv','Sauvegarde','','zuiedit_param_gene_sauv_ctrl');
		zuiedit_param_gene_sauv;
		zuireset(h.radio_sauv)
	end

case 'radio_info'
	if ishandle(hfig)
		% disp('formulaire parametre/generaux/informations')
		zuiedit_param_gene_info;
		zuireset(h.radio_info)
	end

case 'radio_config'
	if ishandle(hfig)
		% disp('formulaire parametre/generaux/informations')
		zuiedit_param_gene_funf('zequa','param.gene','Equation configuration');
		zuireset(h.radio_config)
	end

case 'radio_exec'
	if ishandle(hfig)
		% disp('formulaire parametre/generaux/execution')
		zuiedit_param_gene_funf('zexec','param.gene','Execution','','zuiedit_param_gene_exec_ctrl');
		zuireset(h.radio_exec)
	end

case 'radio_profile'
	if ishandle(hfig)
		% disp('formulaire parametre/profile')
		zuiedit_param_gene_funf('zprofile','param.profile','profiler configuration');
		zuireset(h.radio_profile)
	end

case 'radio_from'
        if ishandle(hfig)
                zuiedit_param_from ;
                zuireset(h.radio_from)
        end


case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig);
		zuireset(h.btn_quit);
	end

otherwise

end

