% ZUIEDIT_FCT  	fonction de gestion des callback associes a 
%			chaque uicontrol (autre que popup ou edit)
%--------------------------------------------------------------
% fichier zuiedit_fct.m
%
% fonction Matlab 5 :
%	fonction de gestion des callback associes a chaque uicontrol du 
%	formulaire edition 
%
% syntaxe :
%	zuiedit_fct(action)
%
% entrees :
%  	action : tag du uicontrol active
%
% sorties :
% 
%
% fonction ecrite par C. Passeron, poste 6119
% version 3.0 du 18/12/2004
% 
% liste des modifications :
% 
%  * 16/10/2001 -> ajout de la fonction zcestou sur le bouton d'aide 
%  * 10/12/2002 -> interface en anglais 
%  * 11/12/2002 -> ajout de la mise a jour de la composition plasma
%  * 18/12/2004 -> ajout de l'interface simulink
%  * 05/02/2006 -> remaniement de l'interface
%--------------------------------------------------------------
function zuiedit_fct(action)

if nargin ==0
	action = ' ';
end
% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('edit');

% si la fenetre a deja ete appelee, on l'active
if ishandle(hfig)
        zuiformvisible(hfig);
end
	
% information pour l'assistant
zuicr(hfig,action) ;

% selon ation
switch lower(action)

% Commandes de Parametre
% ----------------------

case 'radio_time'
        if ishandle(hfig)
                zuiedit_param_time;
                zuireset(h.radio_time)
        end

case 'radio_gene'
	if ishandle(hfig)
		zuiedit_param_gene;
		zuireset(h.radio_gene)
	end

%case 'radio_fonction'
%	if ishandle(hfig)
%		zuiedit_param_modext ;
%		zuireset(h.radio_fonction)
%	end

case 'radio_compo'
	if ishandle(hfig)
                zuiedit_compo;
		zuireset(h.radio_compo)
	end

%case 'radio_split'
%	if ishandle(hfig)
%		zuiedit_param_split ;
%		zuireset(h.radio_split)
%	end

case 'radio_plot'
	if ishandle(hfig)
	        %zuiedit_param_plot ;
                zuiedit_graphique ;
		zuireset(h.radio_plot)
	end

%case 'radio_from'
%	if ishandle(hfig)
%		zuiedit_param_from ;
%		zuireset(h.radio_from)
%	end

%case 'radio_intv'
%	if ishandle(hfig)
%		zuiedit_param_intervalle ;
%		zuireset(h.radio_intv)
%	end

case 'radio_equili'
        if ishandle(hfig)
                zuiedit_equili;
                zuireset(h.radio_equili)
        end

% Commandes de data
% -----------------
case 'radio_mode' 
        if ishandle(hfig)
                zuiedit_transport ;
                zuireset(h.radio_mode)
        end
%	if ishandle(hfig)
%		hout = zuiedit_config_mode ;
%		zuireset(h.radio_mode)
%	end

%case 'radio_geo'
%	if ishandle(hfig)
%		hout = zuiedit_data_geom ;
%		zuireset(h.radio_geo)
%	end

case 'radio_cons'
	if ishandle(hfig)
		hout = zuiedit_data_cons ;
		zuireset(h.radio_cons)
	end

case 'radio_asserv'
	if ishandle(hfig)
                zuiedit_asser;
		zuireset(h.radio_asserv)
	end

case 'radio_prof'
	if ishandle(hfig)
		hout = zuiedit_data_prof ;
		zuireset(h.radio_prof)
	end

case 'radio_source'
	if ishandle(hfig)
		zuiedit_source ;
		zuireset(h.radio_source)
	end

%case 'radio_coef'
%	if ishandle(hfig)
%		hout = zuiedit_data_coef ;
%		zuireset(h.radio_coef)
%	end

%case 'radio_bord'
%	if ishandle(hfig)
%		zuireset(h.radio_bord)
%	end

%case 'radio_impur'
%	if ishandle(hfig)
%		zuireset(h.radio_impur)
%	end

case 'radio_neo'
	if ishandle(hfig)
                zuiedit_neo;
		zuireset(h.radio_neo)
	end

case 'radio_exp'
	if ishandle(hfig)
		zuireset(h.radio_exp)
	end

%case 'radio_equi'
%	if ishandle(hfig)
%		zinitequi ;
%		zuireset(h.radio_equi)
%	end

%case 'radio_scal'
%	if ishandle(hfig)
%               zupdatecompo(1);
%		zuireset(h.radio_scal)
%	end

%case 'radio_mhd'
%	if ishandle(hfig)
%		zuireset(h.radio_mhd)
%	end

case 'radio_mhdyna'
        if ishandle(hfig)
                zuiedit_mhdyna;                
                zuireset(h.radio_mhdyna)
        end



case 'radio_evx'
	if ishandle(hfig)
		zuireset(h.radio_evx)
	end

%case 'radio_simulink'
%        zuisimulink_fct('init');
%	if ishandle(hfig)
%		zuireset(h.radio_simulink)
%	end

case 'radio_postrait'
     if ishandle(hfig)
           zuiedit_param_postrait;
           zuireset(h.radio_postrait)
     end

case 'radio_device'
     if ishandle(hfig)
           zuiedit_param_device;
           zuireset(h.radio_device)
     end
% Boutons quit
	case {'btn_quit','close'}
		if ishandle(hfig)
			zuiformcache(hfig);
			zuireset(h.btn_quit) ;
			zuiclose ;

			% on revient a la fenetre zuidirect
			[hform,h] = zuiformhandle('direct') ;
			zuiformvisible(hform) ;
		end	

% Boutons Aide
	case {'aide'}
		if ishandle(hfig)
			%msgbox(' Desole, pas d''aide pour le moment','Aide','help')
			zbrowser('https://wikicronos.partenaires.cea.fr/wiki/index.php/Documentation');
			zuireset(h.aide)
		end	
	
otherwise
	warning('action not taken into account')
	   
end

