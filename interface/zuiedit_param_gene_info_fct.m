% ZUIDEDIT_PARAM_GENE_INFO_FONCTION gestion des callbacks formulaire 'informations' sous 'Parametres generaux - Edition'
%--------------------------------------------------------------
% fichier zuiedit_param_gene_info_fct.m
%
% fonction Matlab 5 :
% cette fonction est un argument de la fonction de creation du formulaire 
% d'informations des parametres generaux sous le mode edition du formulaire
% principal
%
%
% syntaxe :
%	zuiedit_param_gene_info_fct(action)
%
% entrees :
%  action       =  tag du uicontrol active
%
% sorties :
% 
%
% fonction ecrite par C. Passeron, poste 61-19
% version 1.8, du 10/04/2001.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function zuiedit_param_gene_info_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_gene_info');
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)
		
	case {'btn_quit','close'}
		if ishandle(hfig)
			zuiformcache(hfig);
			zuireset(h.btn_quit);
		end
	
	otherwise
		warning('ation non prise en compte')
	   
end

