% ZUICREATE_FCT  	gestion des callback du formulaire Creation
%--------------------------------------------------------------
% fichier zuicreate_fct.m
%
% fonction Matlab 5
%	fonction de gestion des callback associes a 
%	chaque uicontrol (autre que popup ou edit) du formulaire
%	Creation	
%
% syntaxe :
%	zuicreate_fct(action)
% 
% entrees :
%	action       =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 6119
% version  3.0 , du 03/05/2005.
%
% liste des modifications : 
% * 11/12/2002 : interface en anglais
% * 03/05/2005 : ajout de la cration pour 0d
%--------------------------------------------------------------
function zuicreate_fct(action)

if nargin ==0
	action = ' ';
end
% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('create') ;

% si la fenetre a deja ete appelee, on l'active
if ishandle(hfig)
        zuiformvisible(hfig);
end
	
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)

% Commandes de Parametre
% ----------------------
case 'radio_ts'
	zuireset(h.radio_ts) ;
	zuicloseone(hfig) ;
	evalin('base','zuitsacces;'); 

case 'radio_jet'
	zuireset(h.radio_jet)
	zuicloseone(hfig) ;
	evalin('base','zuijetacces2(1);');  

case 'radio_iter'
	zuireset(h.radio_iter)
	zuicloseone(hfig) ;
	evalin('base','zuiiteracces ;');  

case 'radio_sst1'
	zuireset(h.radio_sst1)
	zuicloseone(hfig) ;
	evalin('base','zuisst1acces ;');  
	
case 'radio_east'
	zuireset(h.radio_east)
	zuicloseone(hfig) ;
	evalin('base','zuieastacces ;');  

case 'radio_hl2a'
	zuireset(h.radio_hl2a)
	zuicloseone(hfig) ;
	evalin('base','zuihl2aacces ;'); 
	 
case 'radio_compass'
	zuireset(h.radio_compass)
	zuicloseone(hfig) ;
	evalin('base','zuicompassacces ;');  

case 'radio_ftu'
	zuireset(h.radio_ftu)
	zuicloseone(hfig) ;
	evalin('base','zuiftuacces ;');

case 'radio_tcv'
	zuireset(h.radio_tcv)
	zuicloseone(hfig) ;
	evalin('base','zuitcvacces ;');

case 'radio_diiid'
	zuireset(h.radio_diiid)
	zuicloseone(hfig) ;
	evalin('base','zuidiiidacces ;');

case 'radio_pdb'
	zuireset(h.radio_pdb)
	zuicloseone(hfig) ;
	evalin('base','zuiitpadbacces(1) ;');  

case 'radio_autre'
	zuireset(h.radio_autre)
	zuicloseone(hfig) ;
	evalin('base','zuivideacces ;');  

case 'radio_plus'
	zuireset(h.radio_plus)
	zuicloseone(hfig) ;
	evalin('base','zuiplusacces ;');  

case 'radio_tcronos'
	zuireset(h.radio_tcronos)
	zuicloseone(hfig) ;
	evalin('base','zuipenelope ;');  

case 'radio_tsg'
	zuireset(h.radio_tcronos)
	zuicloseone(hfig) ;
	evalin('base','zuigcronos;');  

case 'radio_0d'
	zuireset(h.radio_0d)
	zuicloseone(hfig) ;
	evalin('base','zuiacces0d(1);');

case {'btn_quit','close'}
	zuireset(h.btn_quit) ;
	zuicloseone(hfig) ;

case {'aide'}
	msgbox(' sorry, no help at this time','Help','help')
	zuireset(h.aide)
	
otherwise
	warning('action not taken into account')
	   
end

