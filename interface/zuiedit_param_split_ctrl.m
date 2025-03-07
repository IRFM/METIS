% ZUIEDIT_PARAM_SPLIT_CTRL   control des uicontrols du formulaire de split , parametres  généraux
%--------------------------------------------------------------
% fichier zuiedit_param_split_ctrl.m  
%
% fonction Matlab 5 :
%	fonction de control des donnes des uicontrols du formulaire
%	de split , parametres  généraux
% 	sous le mode edition du formulaire principal
%
% syntaxe  :
%	zuiedit_param _split_ctrl(action) ;
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
function zuiedit_param_split_ctrl(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_split');

% information pour l'assistant
zuicr(hfig,action) ;
hoc = getfield(h,action) ;

% bornes
binf_dtmin = 1.e-6 ;
bsup_dtmin = 1.e-2 ;
binf_dtmax = 1.e-5 ;
bsup_dtmax = 0.1   ;

% selon l'action
switch lower(action)

case 'mode_split'
	zdata = zuidata(hoc) ;
	if zdata==0
		zuienable(h.libel_nb) ;
		zuienable(h.nb_split) ;
		zuienable(h.borne_nb) ;
	else		
		zuidisable(h.libel_nb) ;
		zuidisable(h.nb_split) ;
		zuidisable(h.borne_nb) ;
	end
	
case 'dtmax_split'
	% dtmax 
	zdata_dtmax = zuidata(hoc) ;
	zdata_dtmin = zuidata(h.dtmin_split) ;
	if zdata_dtmax < binf_dtmax 
		zdata_dtmax = binf_dtmax ;
	end
	if zdata_dtmax > bsup_dtmax 
		zdata_dtmax = bsup_dtmax ;
	end
	if zdata_dtmax < zdata_dtmin*10 
		zdata_dtmax = zdata_dtmin*10 ;
	end
	zuidata(h.dtmax_split,zdata_dtmax) ;
	
case 'dtmin_split'
	% dtmin 
	zdata_dtmin = zuidata(hoc) ;
	zdata_dtmax = zuidata(h.dtmax_split) ;
	if zdata_dtmin < binf_dtmin 
		zdata_dtmin = binf_dtmin ;
	end
	if zdata_dtmin > bsup_dtmin 
		zdata_dtmin = bsup_dtmin ;
	end
	if zdata_dtmin > zdata_dtmax/10 
		zdata_dtmin = zdata_dtmax/10 ;
	end
	zuidata(h.dtmin_split,zdata_dtmin) ;
	
case 'nb_split'
	% nb 
	zdata = zuidata(hoc) ;
	if zdata <2
		zuidata(h.nb_split,2) ;
	end
	
case 'equi_split'
	zdata       = zuidata(hoc) ;
	zdata_dtmin = zuidata(h.dtmin_split) ;	
	zdata_dtmax = zuidata(h.dtmax_split) ;
	if zdata < zdata_dtmin 
		zuidata(h.equi_split,zdata_dtmin) ;
	elseif zdata > zdata_dtmax
		zuidata(h.equi_split,zdata_dtmax) ;
	end
		
otherwise
	warning('ation non prise en compte')
           
 	
end
