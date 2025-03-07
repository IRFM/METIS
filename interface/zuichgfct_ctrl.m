% ZUICHGFCT_CTRL  controle du formulaire de saisie nom fonction des modules externes
%------------------------------------------------------------------------------- 
% fichier : zuichgfct_ctrl -> 
% 
% 
% fonction Matlab 5 : 
% 	fonction de controle du formulaire de saisie du nom de la fonction
%	des modules externes , parametres  généraux
% 	sous le mode edition du formulaire principal
%  
% syntaxe :  
%   zuichgfct_ctrl(action)
%  
% entrees :  
%	action : tag du uicontrol active
%  
% sorties :  
%  
% fonction écrite par C. Passeron, poste 61-19
% version  1.7  du  29/09/2001  
%  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function zuichgfct_ctrl(action)

% recupere le handle de la fenetre concernee
%[hfig,h] = zuiformhandle('chgfct');
hfig = gcf ;
h    = getappdata(hfig,'zhandle') ;

% disp('control : ') ;
% disp(action) ;

%selon action
switch lower(action)
	
	case 'popup_nom_fct'
		fct = get(h.popup_nom_fct,'string') ;
		ind = get(h.popup_nom_fct,'value') ;
		fctchar = fct{ind} ;
		zuidata(h.edit_nom_fct,fctchar) ;
		if strcmp(fctchar,'other')
			[fctchar,path]=uigetfile('*.m','function name ?') ;
			if (fctchar==0)
				zuidata(h.edit_nom_fct,zuidata(h.popup_nom_fct)) ;	
			else
				% supprimer .m dans le nom de la fct
				ind=findstr(fctchar,'.m') ;
				if ~isempty(ind)
					fctchar=fctchar(1:(ind-1)) ;
				end
				zuidata(h.edit_nom_fct,fctchar) ;
			end
		end
	case 'edit_autre'
		zuidata(h.popup_nom_fct,zuidata(h.edit_nom_fct),'zuimax') ;
		
	otherwise
		warning('action not taking into account')
end

