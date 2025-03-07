% ZUIRUN_FONCTION  	gestion des callback du formulaire de soumission en interactif de zineb
%--------------------------------------------------------------
% fichier   zuirun.m
%
% fonction Matlab 5 :
%	fonction de gestion des callback associes a chaque uicontrol 
%	(autre que popup ou edit)
%	du formulaire de soumission en interactif de zineb%
%
% syntaxe :
%	zuirun_fonction(action)
%
% entrees : 
%	action :  tag du uicontrol active
%
% sorties :
%
% fonction ecrite par J-F Artaud, poste 46-78
%	version 1.9  , du 18/03/2002.
% 
% 
% liste des modifications : 
%	* 12/09/2001 -> activation de la commande zineb_run (J-F Artaud)
%	* 28/11/2001 -> modification de la presentation
%  * 18/03/2002 -> ajout du mode reprise sans reinitilisation de l'equilibre et des sources
%
%--------------------------------------------------------------
function zuirun_fonction(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end
% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('run');
% information pour l'assistant
zuicr(hfig,action);
%
filerep  = evalin('base','param.gene.file','[]');
if isempty(filerep) 
	zuidisable(h.radio_reprise_run)
end
	
% selon ation
switch lower(action)

% Boutons de commandes communs a tous les formulaires
	case {'init','raz'}
		if ishandle(hfig)
			zuiformvisible(hfig);
		else
			zuirun;
			[hfig,h] = zuiformhandle('run');
		end
		zuiformreset(hfig);
		% rajouter le nom de fichier de simulation
		% zuidata(h.text_loadfile,zuidata(h.text_loadfile));
		zuiuploadform(hfig);
		zuireset(h.raz);
		
	case {'annulation','close'}
		if ishandle(hfig)
			zuiformcache(hfig);
			zuireset(h.annulation);
		end
	
	case 'validation'
		if ishandle(hfig)
			% appel de zineb_run
			filename = evalin('base','param.edit.currentfile','[]');
			reprise  = zuidata(h.radio_reprise_run);
                        %cmd = strcat('zineb_run(','''',filename,''',',num2str(reprise),')')
			zuireset(h.validation);
			zuiformcache(hfig) ;
			oks  =evalin('base','param.edit.saveok');
			cr = feval('zineb_run',filename,num2str(reprise),oks);
			% zuisavenonok;
			zuireset(h.validation);

		end
        case 'radio_reprise_run'
	
	     %rien
		  
	otherwise
		warning('ation non prise en compte')
end

