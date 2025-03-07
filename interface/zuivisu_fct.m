% ZUIVISU_FCT   gestion des callback du formulaire de visualisation
%--------------------------------------------------------------
%
% fichier zuivisu_fct.m
%
% fonction Matlab 5 :
%	fonction de gestion des callback associes a 
%	chaque uicontrol (autre que popup ou edit) du formulaire
%	de visualisation
% 
% syntaxe
%
% entrees :
%	action =  tag du uicontrol active
%
% sorties :
%
% fonction ecrite par C. Passeron, poste 6119
% version  2.2 , du 19/01/2005.
% 
% liste des modifications : 
% * 11/12/2002 : interface en anglais
%--------------------------------------------------------------

function zuivisu_fct(action)

if nargin ==0
	action = ' ';
end
% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('visu');
% information pour l'assistant
zuicr(hfig,action)

% selon ation
switch lower(action)
	   
% 
% zdataplot
	case 'radio_zdataplot'
	if ishandle(hfig)
		% disp('appel de zdataplot');
		hdp = findobj(0,'type','figure','tag','zdataplot');
		if ishandle(hdp)
			figure(hdp);
		else
			zdataplot ;
		end
		zuireset(h.radio_zdataplot);
	end

%  fast zdataplot
	case 'fast_zdataplot'
	if ishandle(hfig)
		% disp('appel de zdataplot');
		hdp = findobj(0,'type','figure','tag','zdataplot');
		if ishandle(hdp)
			figure(hdp);
		elseif exist(fullfile(getenv('HOME'),'/zineb/','fastzdataplot.fig'))
	                % ouverture des clipboard
	                %zgclipp;
	                %zgclipt;
			hgload(fullfile(getenv('HOME'),'/zineb/','fastzdataplot')) ;
		else
			zdataplot;
		        hdp = findobj(0,'type','figure','tag','zdataplot');
			hgsave(hdp,fullfile(getenv('HOME'),'/zineb/','fastzdataplot')) ;
                        fprintf('%s cree ...\n',fullfile(getenv('HOME'),'/zineb/','fastzdataplot'));
		end
		zuireset(h.fast_zdataplot);
	end

% plotverif
	case 'radio_plotverif'
	if ishandle(hfig)
		% disp('appel de plotverif');
		evalin('base','plotverif')
		zuireset(h.radio_plotverif)
	end
      
% compare
	case 'radio_compare'
	if ishandle(hfig)
		% disp('Appel de compare ');
 		evalin('base','compare') ;
		zuireset(h.radio_compare)
	end
   	
% polarplot
	case 'radio_polarplot'
	if ishandle(hfig)
		% disp('Appel de polarplot ');
 		evalin('base','polarplotnc') ;
		zuireset(h.radio_polarplot)
	end
	
% mseplot
	case 'radio_mseplot'
	if ishandle(hfig)
		% disp('Appel de mseplot ');
 		evalin('base','mseplot') ;
		zuireset(h.radio_mseplot)
	end
	
% coefverif
	case 'radio_coefverif'
	if ishandle(hfig)
		% disp('Appel de coefverif ');
 		evalin('base','coefverif') ;
		zuireset(h.radio_coefverif)
	end
	
% neocohere
	case 'radio_neocohere'
	if ishandle(hfig)
		% disp('Appel de neocohere');
 		evalin('base','neocohere');
		zuireset(h.radio_neocohere )
	end

% neoverif
	case 'radio_neoverif'
	if ishandle(hfig)
		% disp('Appel de neoverif ');
 		evalin('base','neoverif') ;
		zuireset(h.radio_neoverif)
	end

% etatverif
	case 'radio_etatverif'
	if ishandle(hfig)
		% disp('Appel de etatverif ');
 		evalin('base','etaverif') ;
		zuireset(h.radio_etatverif)
	end
	
% comparetemp
	case 'radio_comparetemp'
	if ishandle(hfig)
		% disp('Appel de comparetemp');
 		evalin('base','comparetemp') ;
		zuireset(h.radio_comparetemp)
	end

% qplot
	case 'radio_qplot'
	if ishandle(hfig)
		% disp('Appel de qplot');
 		evalin('base','qplot') ;
		zuireset(h.radio_qplot)
	end
	
	case 'radio_power'
	if ishandle(hfig)
 		evalin('base','bilanpuiss') ;
		zuireset(h.radio_power)
	end
	
	case 'radio_stabmhd'
	if ishandle(hfig)
 		evalin('base','tracemhd(post,param)','drawnow;text(0.1,0.5,''no data'');') ;
		zuireset(h.radio_stabmhd)
	end

	case 'radio_modemhd'
	if ishandle(hfig)
 		evalin('base','plotmhd','drawnow;text(0.1,0.5,''no data'');') ;
		zuireset(h.radio_modemhd)
	end

	case 'radio_couche'
	if ishandle(hfig)
 		evalin('base','zplotcouche','drawnow;text(0.1,0.5,''no data'');') ;
		zuireset(h.radio_couche)
	end

	case 'radio_cohere'
	if ishandle(hfig)
 		evalin('base','zcomptprof','drawnow;figure(''name'',''TS coherence'');text(0.1,0.5,''no data'');') ;
		zuireset(h.radio_couche)
	end
	
	case 'radio_wdiff'
	if ishandle(hfig)
 		evalin('base','zcompjeux','1;') ;
		zuireset(h.radio_wdiff)
	end

	case 'radio_neutron'
	if ishandle(hfig)
 		evalin('base','visuneutron','1;') ;
		zuireset(h.radio_neutron)
	end

	case 'radio_bootstrap'
	if ishandle(hfig)
 		evalin('base','visubootstrap','1;') ;
		zuireset(h.radio_bootstrap)
	end
   
	case 'radio_scenario'
	if ishandle(hfig)
 		evalin('base','tsscenario','1;') ;
		zuireset(h.radio_scenario)
	end
   
	case 'radio_postece'
	if ishandle(hfig)
 		evalin('base','zeceplot','1;') ;
		zuireset(h.radio_postece)
	end
   
	case 'radio_ecrh'
	if ishandle(hfig)
 		evalin('base','fceplot','1;') ;
		zuireset(h.radio_ecrh)
	end
 
  	case 'radio_flux'
	if ishandle(hfig)
 		evalin('base','fluxedge','1;') ;
		zuireset(h.radio_flux)
	end
        
 	case 'radio_qali'
	if ishandle(hfig)
 		qali_stabilite;
		zuireset(h.radio_qali)
	end
        
 	case 'radio_itb'
	if ishandle(hfig)
 		zitbplot;
		zuireset(h.radio_itb)
	end


       	case 'radio_sepa'
	if ishandle(hfig)
 		ztssepaplot;
		zuireset(h.radio_sepa)
	end
	
	case 'radio_catalogue'
	if ishandle(hfig)
 		zbatchcronos;
		zuireset(h.radio_catalogue)
	end

	case 'radio_li'
	if ishandle(hfig)
                evalin('base','verifli');
                zuireset(h.radio_li)
	end

 	case 'radio_nustar'
	if ishandle(hfig)
                evalin('base','zplotnustar');
                zuireset(h.radio_nustar)
	end	

 	case 'radio_st'
	if ishandle(hfig)
                evalin('base','visudds');
                zuireset(h.radio_st)
	end

 	case 'radio_movie'
	if ishandle(hfig)
                evalin('base','zplasma_movie(param,data);');
                zuireset(h.radio_movie)
	end
	
 	case 'radio_bootvalid'
	if ishandle(hfig)
                evalin('base','zbootvalidplot;');
                zuireset(h.radio_bootvalid)
	end
	
 	case 'radio_ticxs'
	if ishandle(hfig)
                evalin('base','ztsticxs;');
                zuireset(h.radio_ticxs)
	end

 	case 'radio_confinement'
	if ishandle(hfig)
                evalin('base','sor=zbatchconf(data,param);');
                zuireset(h.radio_confinement)
	end
% Boutons quit
	case {'btn_quit','close'}
	if ishandle(hfig)
		zuiformcache(hfig) ;
		zuireset(h.btn_quit) ;
		zuicloseone(hfig) ;
	end	

% Boutons jeux1
	case {'btn_jeux1'}
	if ishandle(hfig)
		evalin('base','clear jeux1');
		zuireset(h.btn_jeux1) ;
	end	

% Boutons Aide
	case {'aide'}
		if ishandle(hfig)
			msgbox(' Sorry no help for the present time','Help','help')
			zuireset(h.aide)
		end	
	otherwise
		warning('action not taken into account')
end

