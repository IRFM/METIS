% ZUICHGFCT_FCT  controle du formulaire de saisie nom fonction des modules externes
%------------------------------------------------------------------------------- 
% fichier : zuichgfct_fct -> 
% 
% 
% fonction Matlab 5 : 
% 	fonction de controle du formulaire de saisie du nom de la fonction
%	des modules externes , parametres  g�n�raux
% 	sous le mode edition du formulaire principal
%  
% syntaxe :  
%   zuichgfct_fct(action)
%  
% entrees :  
%	action : tag du uicontrol active
%  
% sorties :  
%  
% fonction �crite par C. Passeron, poste 61-19
% version  2.1  du  28/07/2003
%  
% liste des modifications :  
% 30/08/2001 -> ajout de la connexion au module (J-F Artaud)
% 24/09/2001 -> ajout de la recopie des parametres communs
% 19/10/2001 -> changement de la validation
% 28/07/2003 -> petit bug sans consequence si le nom du modules est compose d'espaces
%  
%-------------------------------------------------------------------------------  
function zuichgfct_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
%[hfig,h] = zuiformhandle('chgfct') ;
hfig = gcf ;
h    = getappdata(hfig,'zhandle') ;
htag = get(hfig,'tag') ;
% information pour l'assistant
zuicr(hfig,action) ;

% selon ation
switch lower(action)
	
case {'init','raz'}
	ishandle(hfig)
	if ishandle(hfig)
		zuiformvisible(hfig) ;
	else
		zuichfct ;
		%[hfig,h] = zuiformhandle('chgfct');
		hfig = gcf;
		h    = getappdata(hfig,'zhandle');
	end
	
	zuiformreset(hfig) ;
 	zuidata(h.popup_nom_fct,zuidata(h.edit_nom_fct),'zuimax') ;
	zuireset(h.raz) ;
	
case {'annulation','close'}
	if ishandle(hfig)
		zuicloseone(hfig) ;
		
	end
case 'validation'
	if ishandle(hfig)
		zuiformcache(hfig) ;
 		zuidownloadform(hfig) ;
%		get(getappdata(hfig,'hfct'),'String')
%		get(getappdata(hfig,'hmodule'),'String') 
		nomf = zuidata(h.edit_nom_fct)   ; % nom du module externe
                nomf(nomf <= ' ') = [];
		module  = getappdata(hfig,'module'); % nom de la fonction
		try
		   nb    = evalin('base',strcat('param.nombre.',module));
		catch
		   nb    = 1;
		end
                if ~isempty(nomf)
		     zmajmodext(module,nomf) ;
                end
		zuidata(getappdata(hfig,'hmodule'),zuidata(h.edit_nom_fct)) ;
		if ~isempty(nomf)
			%par = feval(nomf,nb);
			%zassignin('base',strcat('param.cons.',module),par.valeur);
			nbok =  zmergeparam(module,nomf,nb);
			if nbok > 0
				fprintf('=> %d parametres communs\n',nbok);
			end
			
		end
		zuisavenonok ;
		zuireset(h.validation) ;
	end
otherwise
	warning('action non prise en compte')
	
end

