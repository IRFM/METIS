% ZUIEDIT_PARAM_GENE_SAUV_FCT controle formulaire des parametres relatifs a la sauvegarde des resultats
%--------------------------------------------------------------
% fichier zuiedit_param_gene_sauv_fct.m ->
%
% fonction Matlab 5 :
%	fonction de controle des uicontrols du formulaire
% 	des parametres relatifs a la sauvegarde des resultats
%	des parametres  généraux
% 	sous le mode edition du formulaire principal
%
% syntaxe :
%	zuiedit_param_gene_sauv_fct(action)
%
% entrees :
%	action : tag du uicontrol active
%
% sorties
%  				
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 23/03/2004.
% 
% liste des modifications : 
% * 23/03/2004 : accepte le cas rapsauve=[] (signifie pas de sauvegarde rapide)
%--------------------------------------------------------------
function zuiedit_param_gene_sauv_fct(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_sauv') ;
if ~ishandle(hfig)
	return
end
zuicr(hfig,action) ;
hoc = getfield(h,action) ;

%selon l'action
switch lower(action)

% File
	case 'edit_file'
	zdata = zuidata(hoc) ;
	% controle des chemnins, repertoire, espace disque
	eval(['[s,t] = unix(''touch ' zdata ''');'])
	if s~=0
		txt=sprintf('le fichier \n %s \n ne peut pas etre cree',zdata);
		herror = errordlg(txt,'ATTENTION');
		zwaitfor(herror)
		zuiformreset(hfig);
		return
	end
%	user = getenv('LOGNAME')
%	eval(['[s,df] = unix(''df /usr/drfc/' user ''')'])

	eval(['[s,t]=unix(''df -k '  zdata ' | egrep -v Filesystem | awk ''''{print $4}'''' '')']);
	txt = sprintf(' il reste %s blocks(1024) disponibles ',t);
	eval(['[s,disque]= unix(''df -k ' zdata ' | egrep -v Filesystem | awk ''''{print $1}'''' | cut -d''''/'''' -f3'')'])
	title = sprintf('sur le disque %s',disque)
	hwarn = warndlg(txt,title,'warn')
	zuidata(h.edit_file,zdata) ;
	
% Rapsauve
	case 'edit_rapsauve'
	zdata = zuidata(hoc) ;
        if ~isempty(zdata)
	   % controle des chemins, repertoire, espace disque
	   eval(['[s,t] = unix(''touch ' zdata ''');']);
	   if s~=0
		txt=sprintf('le fichier \n %s \n ne peut pas etre cree',zdata);
		herror = errordlg(txt,'ATTENTION');
		zwaitfor(herror)
		zuiformreset(hfig);
		return
	   end
	   eval(['[s,t]=unix(''df -k '  zdata ' | egrep -v Filesystem | awk ''''{print $4}'''' '')']);
	   txt = sprintf(' il reste %s blocks(1024) disponibles ',t);
	   eval(['[s,disque]= unix(''df -k ' zdata ' | egrep -v Filesystem | awk ''''{print $1}'''' | cut -d''''/'''' -f3'')'])
	   title = sprintf('sur le disque %s',disque)
	   hwarn = warndlg(txt,title,'warn')
        end
	zuidata(h.edit_rapsauve,zdata) ;
	
end
