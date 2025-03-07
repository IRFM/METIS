%
% ZUIDEDIT_DATA_ASSER_FCT
%      	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de consignes des asservissements
%--------------------------------------------------------------
% fonction Matlab 5 :
%
% fichier zuiedit_data_asser_fct.m  
%
% syntaxe :
% 
% entrees :
% 
%  action       =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 10/12/2002.
% 
% liste des modifications : 
% * 10/12/2002 : interface en anglais
%
%--------------------------------------------------------------
function zuiedit_data_asser_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_asser') ;

% information pour l'assistant
zuicr(hfig,action) ;

% variables d'entr� de l'�iteur de consignes zuieditcons
x           = evalin('base','data.gene.temps') ;
texte_x     = 'temps' ;
var_x       = 'void';
code_retour = 'abs' ;
liste_ref   = {} ;
var_ref     = {} ;
texte_prop  = '' ;
var_prop    = '' ;


switch lower(action)
case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
        % selon action
        [acte,nomc] = strtok(action,'_');
        nomc          = nomc(2:end);
	if any(nomc=='_')
	   [nomc,numc] = strtok(nomc,'_');
           numc        = str2num(numc(2:end));
	else
	   numc =[];
	end

 	if isempty(numc)	
	     switch acte
	     case 'edit'
		nom     = strcat('data.cons.asser.',nomc) ;
		y       = evalin('base',nom) ;
		texte_y = nomc ;
		var_y   = nom;
		canal   = 1 ;
		hout = zuieditcons(nom,getfield(info.data.cons.asser,nomc),x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	  
		zuireset(getfield(h,strcat('edit_',nomc)));

	     case 'import'
		hout = zuiedit_import_mode(strcat('data.cons.asser.',nomc),'consigne') ;
	        zuireset(getfield(h,strcat('import_',nomc)));
	     otherwise
	     	 disp('Action not taking into account');
	     end
        else
	     switch acte
	     case 'edit'
		nom     = strcat('data.cons.asser.',nomc,sprintf('(:,%d)',numc));
		y       = evalin('base',nom) ;
		texte_y = sprintf('%s(:,%d)',nomc,numc);
		var_y   = strcat('data.cons.asser.',nomc);
		canal   = numc ;
		hout = zuieditcons(nom,getfield(info.data.cons.asser,nomc),x,y,texte_x,texte_y,var_x,var_y,canal, ...
		                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
		zuireset(getfield(h,sprintf('edit_%s_%d',nomc,numc)));
	
	     case 'import'
		nom     = strcat('data.cons.asser.',nomc,sprintf('(:,%d)',numc));
		hout = zuiedit_import_mode(nom,'consigne') ;
		zuireset(getfield(h,sprintf('import_%s_%d',nomc,numc)));

	     case 'prec'
		 rep =questdlg('should we copy the last channel ?', ...
			       'Confirmation', ...
			       'Yes','No','Yes');
		 switch rep
		 case 'Yes'
			 nom = strcat('data.cons.asser.',nomc,sprintf('(:,%d)',numc-1));
			 y = evalin('base',nom) ;
			 nom = strcat('data.cons.asser.',nomc,sprintf('(:,%d)',numc));
			 zassignin('base',nom,y) ;
		 case 'Non'
		 end
		 zuireset(getfield(h,sprintf('prec_%s_%d',nomc,numc)));
	
	     otherwise
	    	 disp('Action not taken into account');
	     end
	end	
end

