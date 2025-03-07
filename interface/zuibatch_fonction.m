% ZUIBATCH_FONCTION  	gestion callback du formulaire de soumission en batch
%--------------------------------------------------------------
% fichier : zuibatch_fonction 
%
% fonction Matlab 5 :
%	fonction de gestion des callback associes a  chaque uicontrol 
%  (autre que popup ou edit).
%	Cette fonction est un argument de la fonction de creation d'un formulaire zuibatch
%
% syntaxe :
%	zuibatch_fonction(action)
%
% entrees :
%  action :  tag du uicontrol active
%
% sorties :
%
% fonction ecrite par C. Passeron, psote 61-19
% version  3.0  du  07/02/2005  
% 
% 
% liste des modifications : 
% * 28/11/2001 -> modification de la presentation
% * 18/03/2002 -> ajout du mode reprise sans reinitilisation de l'equilibre et des sources
% * 17/09/2002 -> ajout du pc cronos
% * 07/02/05   -> nouveau systeme TS
%
%--------------------------------------------------------------
function zuibatch_fonction(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('batch');
% information pour l'assistant
zuicr(hfig,action);

% selon ation
switch lower(action)

% Boutons de commandes communs a tous les formulaires
	case {'init','raz'}
		if ishandle(hfig)
			zuiformvisible(hfig);
		else
			zuibatch;
			[hfig,h] = zuiformhandle('batch');
		end
		zuiformreset(hfig);
		zuiuploadform(hfig);
		zuireset(h.raz);
		
	case {'annulation','close'}
		if ishandle(hfig)
			zuiformcache(hfig);
			zuireset(h.annulation);
		end
	
	case 'validation'
		if ishandle(hfig)
			% appel de zineb_batch
			filename = evalin('base','param.edit.currentfile','[]');
			%machine  = zuidata(h.popup_nom_machine_batch);
			machine = []; % patch since the machine menu has been removed from the zuibatch window
			reprise  = zuidata(h.radio_reprise_batch);
			queue    = zuidata(h.popup_nom_queue_batch);
			%cmd = strcat('zineb_batch(','''',filename,''',',num2str(reprise),',''',machine,''',''',queue,''')')
			oks  =evalin('base','param.edit.saveok');
			delete(hfig)
			drawnow
			feval('zineb_batch',filename,num2str(reprise),machine,queue,oks);
			%zuisavenonok;
			%zuireset(h.validation);
			%zuiformcache(hfig);

		end
        case 'radio_reprise_batch'
	
	     %rien

	otherwise
		warning('ation non prise en compte')
end

