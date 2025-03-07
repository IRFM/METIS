% ZUIDEDIT_PARAM_COMP_FCT  gestion des callbacks formulaire composition 
%--------------------------------------------------------------
%	fichier zuiedit_param_comp_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	des composition , parametres  généraux
%	sous le mode edition du formulaire principal
%
% syntaxe :
%	zuiedit_param_comp_fct(action)
%
% entrees :
%  action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 30/08/2001.
% 
% liste des modifications : 
%  * 30/08/2001 -> correction bug action 'init' (J-F Artaud)
%  * 16/10/2001 -> correction bug lors du rappel de la fenetre sur la mise a jour des 3 premiers gaz
%
%--------------------------------------------------------------
function zuiedit_param_comp_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_comp');
% information pour l'assistant
zuicr(hfig,action) ;
if ~strcmp(action,'init')
	hoc = getfield(h,action) ;
end
   
z = [1,1,1,2,2,3] ;
a = [1,2,3,4,3,7] ;
% selon ation
switch lower(action)

% fonctions equilibre
	case {'annulation','close'}
		zuicloseone(hfig);	
	
	case {'init','raz'}
		zuiformvisible(hfig);
		zuiformreset(hfig);
		zuiuploadform(hfig);
		% 
		% relecture des infos pour les 3 premiers gazs
		%
		ind  = evalin('base','param.compo.a(1)') ;
		set(h.pop_majoritaire,'value',ind);
		ind1  = evalin('base','param.compo.z(2)') ;
		ind   = evalin('base','param.compo.a(2)') ;
		if ind1==2 
			if ind==4
				ind = 4 ;
			elseif ind==3
				ind = 5 ;
			elseif ind==7
				ind = 6 ;
			end
		end
		set(h.pop_minoritaire1,'value',ind);
		ind1  = evalin('base','param.compo.z(3)') ;
		ind   = evalin('base','param.compo.a(3)') ;
		if ind1==2 
			if ind==4
				ind = 4 ;
			elseif ind==3
				ind = 5 ;
			elseif ind==7
				ind = 6 ;
			end
		end
		set(h.pop_minoritaire2,'value',ind);

		zuireset(h.raz);
	
	case 'validation'
		ind = get(h.pop_majoritaire,'Value');
		zassignin('base','param.compo.z(1)',z(ind))
		zassignin('base','param.compo.a(1)',a(ind))
		ind = get(h.pop_minoritaire1,'Value');
		zassignin('base','param.compo.z(2)',z(ind))
		zassignin('base','param.compo.a(2)',a(ind))
		ind = get(h.pop_minoritaire2,'Value');
		zassignin('base','param.compo.z(3)',z(ind))
		zassignin('base','param.compo.a(3)',a(ind))
		
		zuiformcache(hfig) ;
		zuidownloadform(hfig);
		zuisavenonok;
		zuireset(h.validation);

	otherwise
		warning('action non prise en compte')
	   
end

