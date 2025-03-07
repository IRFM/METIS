% ZUIASSISTANT_FCT   gestion des callback du formulaire des assistants
%--------------------------------------------------------------
%
% fichier zuiassistant_fct.m
%
% fonction Matlab 5 :
%	fonction de gestion des callback associes a 
%	chaque uicontrol (autre que popup ou edit) du formulaire
%	de assistant
% 
% syntaxe
%
% entrees :
%	action =  tag du uicontrol active
%
% sorties :
%
% fonction ecrite par J-F Artaud, poste 62-15
% version  2.1 , du 17/06/2003.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zuivisu_fct(action)

if nargin ==0
	action = ' ';
end
% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('assistant');
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)
	   
% 
% jivaro 
	case 'radio_jivaro'
	if ishandle(hfig)
                evalin('base','zjivaro')
		zuireset(h.radio_jivaro);
	end

% changement de base temps
	case 'radio_basetempsch'
	if ishandle(hfig)
		evalin('base','zuiplusacces')
		zuireset(h.radio_basetempsch)
	end
      
% changement du nombre de point radiaux
	case 'radio_rayon'
	if ishandle(hfig)
		%disp('a ecrire ');
  		evalin('base',['info= zsignac;option_signac=info.valeur;', ...
                                'zuicreefunform(''zsignac'',''option_signac'',1,0,''[data,param] = zsignac(data,param,option_signac);'');']);
		zuireset(h.radio_rayon)
	end
   	
% export ascii
	case 'radio_exportascii'
	if ishandle(hfig)
		% disp('Appel de polarplot ');
 		evalin('base',['info=zexportascii;option_exportascii= info.valeur;', ...
                               'zuicreefunform(''zexportascii'',''option_exportascii'',1,0,''dirroot = zexportascii(param,data,post,option_exportascii);'');']);
		zuireset(h.radio_exportascii)
	end
	
% importascii
	case 'radio_importascii'
	if ishandle(hfig)
		% disp('Appel de mseplot ');
 		evalin('base','zimportascii;') ;
		zuireset(h.radio_importascii)
	end
	
% etalhts
	case 'radio_etalhts'
	if ishandle(hfig)
		% disp('Appel de mseplot ');
 		evalin('base',' zcalcetalh([],1);') ;
		zuireset(h.radio_etalhts)
	end
	
% rema
	case 'radio_rematune'
	if ishandle(hfig)
 		evalin('base',' ztunerema(1,2,3);') ;
		zuireset(h.radio_rematune)
	end
	
% fci
	case 'radio_fcitune'
	if ishandle(hfig)
 		evalin('base',' ztunefci(1,2,3);') ;
		zuireset(h.radio_fcitune)
	end
% kinezero
	case 'radio_kine'
	if ishandle(hfig)
 		evalin('base',' zkine(1,2,3);') ;
		zuireset(h.radio_kine)
	end

% neident
	case 'radio_neident'
	if ishandle(hfig)
 		evalin('base','zinversionnl;') ;
		zuireset(h.radio_neident)
	end

% eceident
	case 'radio_eceident'
	if ishandle(hfig)
 		evalin('base','zeceverif(1,2,3);') ;
		zuireset(h.radio_eceident)
	end

	case 'radio_0d'
	if ishandle(hfig)
 		evalin('base','zcall0d;') ;
		zuireset(h.radio_0d)
	end
% entropie maximum + contrainte
	case 'radio_nemaxent'
	if ishandle(hfig)
 		evalin('base','zinversionnl2;') ;
		zuireset(h.radio_nemaxent)
	end

% Luke
	case 'radio_luketune'
	if ishandle(hfig)
 		evalin('base',' ztuneluke(1,2,3);') ;
		zuireset(h.radio_luketune)
	end
	
% neo
	case 'radio_neotune'
	if ishandle(hfig)
 		evalin('base',' ztuneneo(1,2,3);') ;
		zuireset(h.radio_neotune)
	end

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
			msgbox(' not available for the moment','Aide','help')
			zuireset(h.aide)
		end	
	otherwise
		warning('action not taking into account')
end

