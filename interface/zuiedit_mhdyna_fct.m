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

function zuiedit_mhdyna_fct(action)

% disp('callback : ')
% disp(action)


info=zinfo;
% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_mhdyna');
 
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)


case 'mhdcalc'
        zuiedit_mhd_calc;
        zuireset(h.mhdcalc) ;

case 'mhdyna'
        % zuiedit_data_cons_mhd ;
        nom         = 'data.cons.stab' ;
        x           = evalin('base','data.gene.temps') ;
        y           = evalin('base',nom) ;
        texte_x     = 'temps' ;
        texte_y     = 'stab';
        var_x       = 'void';
        var_y       = nom ;
        canal       = 1 ;
        code_retour = '' ;
        liste_ref   = {} ;
        var_ref     = {} ;
        texte_prop  = '' ;
        var_prop    = '' ;
        hout = zuieditcons(nom,info.data.cons.stab,x,y,texte_x,texte_y,var_x,var_y, ...
                           canal,code_retour,liste_ref,var_ref,texte_prop,var_prop)  ;
        zuireset(h.mhdyna) ;

case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig);
		zuireset(h.btn_quit);
	end

otherwise

end

