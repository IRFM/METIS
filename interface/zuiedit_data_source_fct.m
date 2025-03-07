% ZUIDEDIT_DATA_SOURCE_FCT gestion callbacks du formulaired'edition de profils des coefficients de transport
%      	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils des coefficients de transport
%--------------------------------------------------------------------
% fichier zuiedit_data_source_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils des coefficients de transport
%
% syntaxe :
%	zuiedit_data_source_fct(action)
% 
% entrees :
%  action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 02/09/2003.
% 
% liste des modifications : 
% * 02/09/2003 -> ajout de la source ripple (rip)
%
%--------------------------------------------------------------
function zuiedit_data_source_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_source') ;

% information pour l'assistant
zuicr(hfig,action) ;

% selon ation
switch lower(action)

case 'radio_fci'
	hout = zuiedit_data_srcfci ;
	zuireset(h.radio_fci) ;

case 'radio_fce'
	hout = zuiedit_data_srcfce;
	zuireset(h.radio_fce) ;

case 'radio_hyb'
	hout = zuiedit_data_srchyb ;
	zuireset(h.radio_hyb) ;

case 'radio_idn'
	hout = zuiedit_data_srcidn ;
	zuireset(h.radio_idn) ;

case 'radio_rip'
	hout = zuiedit_data_srcrip ;
	zuireset(h.radio_rip) ;

case 'radio_ext'
	hout = zuiedit_data_srcext ;
	zuireset(h.radio_ext) ;

case 'radio_n0'
	hout = zuiedit_data_srcn0 ;
	zuireset(h.radio_n0) ;

case 'radio_autres'
	hout = zuiedit_data_srcautres ;
	zuireset(h.radio_autres) ;

case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

