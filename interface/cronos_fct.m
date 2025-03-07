% CRONOS_FCT  	fonction de gestion des callback du formulaire initial de Cronos
%--------------------------------------------------------------
% fichier cronos_fct.m
% 
% fonction Matlab 5 :
% 	fonction de gestion des callback associes a chaque uicontrol (autre que popup ou edit)
%	du formulaire initial cronos
%
% syntaxe
%	zuidirect_fct(action)
%
% entrees :
%	action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 6119
% version  2.0 , du 12/12/2002.
% 
% liste des modifications : 
% * 19/02/2002 -> test du "popup" pour activer/desactiver la 
%                 capture des fenetres -> images, et references croiseees
% * 22/02/2002 -> ajout de la variable globale de langue
%
%--------------------------------------------------------------

function zuidirect_fct(action)

if nargin ==0
	action = ' ';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('cronos');
% information pour l'assistant
zuicr(hfig,action) ;

% selon ation
switch lower(action)
	   
% 
% Chargement du fichier de travail
case 'radio_francais'
	if ishandle(hfig)
	   setappdata(0,'langue_cronos','francais');
		zuidirect ;
		%zuireset(h.radio_francais) ;
		zuicloseone(hfig) ;
	end

case 'radio_anglais'
	if ishandle(hfig)
	   st = sprintf(' english version version 2.0') ;
	   herror = warndlg(st,'Be care ! ','modal') ;
		%
		% chargement des ressources
		%
		ffza   = which('zuidirect');
		ppza   = fileparts(ffza);
		ppzaa  = strcat(ppza,'_anglais');
		rmpath(ppza)
		addpath(ppzaa);
		clear functions
		
		%
		% lancement de l'interface
		%
	   setappdata(0,'langue_cronos','anglais');
	   setappdata(0,'uicrossref',[]) ; % securite
	   rmappdata(0,'uicrossref') ;
		zuidirect ;
		zuicloseone(hfig) ;
	end
	%st = sprintf(' La version en anglais n''est pas disponible pour le moment') ;
	%herror = warndlg(st,'Désolé','modal') ;
	%zuireset(h.radio_anglais) ;

case 'popup_capture'
	uicrossref = zuidata(h.popup_capture) ;
	setappdata(0,'uicrossref',uicrossref) ;

% Boutons quit
case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig) ;
		zuireset(h.btn_quit) ;
		zuicloseone(hfig) ;
	end	

% Boutons Aide
case {'aide'}
	if ishandle(hfig)
		msgbox(' Sorry, no help at this time','Help','help')
		zuireset(h.aide)
	end
		
otherwise
		warning('action not taken into account')
end

