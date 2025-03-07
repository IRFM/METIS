% ZUIBATCH_CONTROL  	control des popup et edit du formulaire de soumission en batch 
%--------------------------------------------------------------
% fichier : zuibatch_control
%
% fonction Matlab 5 :
%	cette fonction est un argument de la fonction de creation d'un formulaire 
%	elle controle les popup et edit du formulaire de soumission en batch 
%
% syntaxe :
%	zuibatch_control(action)
%
% entrees : 
%  action :  tag du uicontrol active
%
% sorties :
% 
%
% fonction ecrite  par C. Passeron, psote 61-19
% version 1.7  du  29/09/2001  
% 
% liste des modifications : 
% 
% * 17/09/2002 -> ajout du pc cronos
%
%--------------------------------------------------------------
function zuibatch_control(action)

if nargin ==0
	action = ' ';
end
% disp('control : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('batch');

%selon action
switch lower(action)
	
	case 'popup_nom_queue_batch'
		zuidata(h.popup_nom_queue_batch);
		
	%case 'popup_nom_machine_batch'
	%	zuidata(h.popup_nom_machine_batch);
		
	otherwise
end

