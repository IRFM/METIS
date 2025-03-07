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
[hfig,h] = zuiformhandle('ed_transport');
 
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)

case 'radio_transmode'
       if ishandle(hfig)
               hout = zuiedit_config_mode ;
               zuireset(h.radio_transmode)
       end

case 'limites'
       if ishandle(hfig)
               zuiedit_data_cons_limites ;
               zuireset(h.limites) ;
       end

%case 'radio_coef'
%        if ishandle(hfig)
%                zuiedit_transport_coeff;
%                %hout = zuiedit_data_coef ;
%                zuireset(h.radio_coef)
%        end

case 'transport_coeff_mode'
        if ishandle(hfig)
            	zuiedit_transport_coeff_mode;
              	zuireset(h.transport_coeff_mode);
         end

case 'transport_coeff_calc'
	if ishandle(hfig)
		zuiedit_transport_coeff_calc;
		zuireset(h.transport_coeff_calc);
	end

case 'transport_coeff_prof'
        if ishandle(hfig)
       		hout = zuiedit_data_coef ;
                zuireset(h.transport_coeff_prof);
        end

case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig);
		zuireset(h.btn_quit);
	end

otherwise

end

