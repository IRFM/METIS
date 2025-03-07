% ZUIDIRECT_FCT  	gestion  callback du formulaire principal
%--------------------------------------------------------------
% fichier zuidirect_fct.m
% 
% fonction Matlab 5 :
% cette fonction est utilisee par zuidirect.m
%
% syntaxe :  
%	zuidirect_fct(action)
%
% entrees :
%	action       =  tag du uicontrol active
%
% sorties :
% 
%
% fonction ecrite par C. Passeron, poste 6119
% version  3.0 , du 15/03/2005.
% 
% 
% liste des modifications : 
%
%   * 14/09/2001 -> changement de l'appel de zineb_vup (fonctionnement dans le workspace) J-F Artaud
%   * 14/09/2001 -> changement de l'appel de zineb_update (mode standard avec zinit) J-F Artaud
%   * 24/09/2001 -> mise a jour des champs de la figure lors du changement de nom (comme pour le load)
%   * 24/10/2001 -> petit changement dans l'affichage du fichier de trace
%   * 06/11/2001 -> changement de l'appel de zuisavedata
%   * 13/02/2002 -> remet la struture post dans le workspace
%   * 25/03/2002 -> ajout du mode rapide pour l'appel du postprocessing
%   * 11/12/2002 -> interface anglais
%   * 15/03/2005 -> suppression de la protection hercule/rigel
%--------------------------------------------------------------

function zuidirect_fct(action)

if nargin ==0
	action = ' ';
end
% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('direct');
% information pour l'assistant
zuicr(hfig,action) ;

% selon ation
switch lower(action)
	   
% 
% Chargement du fichier de travail
case 'radio_loadfile'
if ishandle(hfig)
%	if (existbase('param')==1) & (existbase('data')==1)
%		saveok = evalin('base','param.edit.saveok') ;
%		if saveok == 0
%			rep =questdlg('Do you want to come back to the main menue ?', ...
%	         	     'be careful -> unsaved data !', ...
%	         	     'Yes','No','No');
%			switch rep
%			case 'Yes'
%				zuireset(h.radio_loadfile) ;
%				return
%			case 'No'
%				set(h.radio_savefile,'foregroundcolor',[0. 0. 0])
%			end
%		end
%	end
	zuiload;

	zuiuploadform(hfig) ;

	numshot=evalin('base','param.from.shot.num','[]');
	if ~isempty(numshot)
		numchoc = fix(numshot) ;
 		zuidata(h.text_numchoc,num2str(numchoc));
		
		occurence = round((numshot-numchoc)*10);
 		zuidata(h.text_occurence,num2str(occurence));
		
		tps1=evalin('base','param.gene.tdeb');
		tps2=evalin('base','param.gene.tfin');
		tps=[num2str(tps1,2) ' -> ' num2str(tps2,2)];
 		zuidata(h.text_temps,tps);
	end
	zuireset(h.radio_loadfile);
end

% Rename du fichier de travail
case 'radio_renamefile'
if ishandle(hfig)
	zuirename ;
	zuiuploadform(hfig) ;
	numshot=evalin('base','param.from.shot.num','[]');
	if ~isempty(numshot)
		numchoc = fix(numshot) ;
 		zuidata(h.text_numchoc,num2str(numchoc));
		
		occurence = round((numshot-numchoc)*10);
 		zuidata(h.text_occurence,num2str(occurence));
		
		tps1=evalin('base','param.gene.tdeb');
		tps2=evalin('base','param.gene.tfin');
		tps=[num2str(tps1,2) ' -> ' num2str(tps2,2)];
 		zuidata(h.text_temps,tps);
	end
	zuireset(h.radio_renamefile) ;
end
   
   
% Sauvegarde du fichier de travail
case 'radio_savefile'
if ishandle(hfig)
 	zuisavedata('force') ;
	zuireset(h.radio_savefile) ;
end
   
% Creation du fichier de travail
case 'radio_createfile'
if ishandle(hfig)
	if (existbase('param')==1) & (existbase('data')==1)
	   try
				saveok = evalin('base','param.edit.saveok') ;
		catch
				saveok = 1;
		end 
		if saveok == 0
			rep =questdlg('do you want to come back to the main menu ?', ...
	         	     'be careful -> unsaved data !', ...
	         	     'Yes','No','No');
			switch rep
			case 'Yes'
				zuireset(h.radio_createfile) ;
				return
			case 'No'
				set(h.radio_savefile,'foregroundcolor',[0. 0. 0])
			end
		end
	end
	zuicreate ;
	zuireset(h.radio_createfile) ;
end
	
% Edition du fichier de travail
case 'radio_editfile'
if ishandle(hfig)
	zuiformcache(hfig) ;
	zuiedit ;
	zuireset(h.radio_editfile) ;
end
   
% Execution en interactif 
case 'radio_runinter'
if ishandle(hfig)

	if (existbase('param')==1) & (existbase('data')==1)
		saveok = evalin('base','param.edit.saveok') ;
		if saveok == 0
			rep =questdlg('do you want to come back to the main menu ?', ...
	         	     'be careful -> unsaved data !', ...
	         	     'Yes','No','No');
			switch rep
			case 'Yes'
				zuireset(h.radio_runinter) ;
				return
			case 'No'
				set(h.radio_savefile,'foregroundcolor',[0. 0. 0])
			end
		end
	end
  	zuirun ;
	zuireset(h.radio_runinter) ;
end

% Execution en batch LSF
case 'radio_runbatch'
if ishandle(hfig)
	if (existbase('param')==1) & (existbase('data')==1)
		saveok = evalin('base','param.edit.saveok') ;
		if saveok == 0
			rep =questdlg('do you want to come back to the main menu ?', ...
	         	     'be careful -> unsaved data !', ...
	         	     'Yes','No','No');
			switch rep
			case 'Yes'
				zuireset(h.radio_runbatch) ;
				return
			case 'No'
				set(h.radio_savefile,'foregroundcolor',[0. 0. 0])
			end
		end
	end
	zuibatch ;
	zuireset(h.radio_runbatch) ;
end
   
% Rebuilt
case 'radio_rebuilt'
if ishandle(hfig)
 	filename = evalin('base','param.gene.file') ;
	cmd = ['zineb_rebuilt(''' filename ''')'] ;
	eval(cmd) ;
	zuireset(h.radio_rebuilt)
end

% Visualisation
case 'radio_visu'
if ishandle(hfig)
	zuivisu ;
	zuireset(h.radio_visu) ;
end
     
% Affichage du resume
case 'radio_resume'
if ishandle(hfig)
        evalin('base','zresume(param,data,post);','fprintf(''Error in zresume :\n%s'',lasterr);');
	zuireset(h.radio_resume) ;
end
     
% Appel de l'assistant
case 'radio_assistant'
if ishandle(hfig)
	zuireset(h.radio_assistant) ;
        zuiassistant;
end
   
% Mise a jour du code dynamique
case 'radio_majcode'
if ishandle(hfig)
	button = questdlg('Update of the code ',...
     	'confirmation',' Yes ',' No ','No');
     	if strcmp(button,' Yes ')
 		%filename = evalin('base','param.edit.currentfile');
		%cmd = ['zineb_update(''' filename ''')'] ;
     		%eval(cmd) ;
     		evalin('base','zineb_update;');
     	elseif strcmp(button,' No ')
        	disp('??')
     	end				
	zuireset(h.radio_majcode) ;
end
	
% Mise a jour du fichier de travail
case 'radio_majfile'
if ishandle(hfig)
 	%filename = evalin('base','param.gene.origine');
	%cmd = ['zineb_vup(''' filename ''')'] ;
	%eval(cmd) ;
	% modification du 14/09/2001 J-F Artaud
	evalin('base','[cr,param,data,post] = zineb_vup(param,data,post);');
	zuireset(h.radio_majfile)
end

% Gestion du  path
case 'radio_path'
if ishandle(hfig)
	zineb_path ;
	zuireset(h.radio_path) ;
end

% Post-Traitement
case 'radio_posttrait'
if ishandle(hfig)
	if (existbase('param')==1) & (existbase('data')==1)
		saveok = evalin('base','param.edit.saveok') ;
		if saveok == 0
			rep =questdlg('do you want to come back to the main menu ?', ...
	         	     'be careful -> unsaved data !', ...
	         	     'Yes','No','No');
			switch rep
			case 'Yes'
				zuireset(h.radio_posttrait) ;
				return
			case 'No'
				set(h.radio_savefile,'foregroundcolor',[0. 0. 0])
			end
		end
	else 
	   saveok =0;
	end
 	%filename = evalin('base','param.gene.file');
	%cmd = ['zineb_post(''' filename ''')'] ;
	if saveok == 1
	   evalin('base','[cr,post] = zineb_post(param.gene.file,1);');
	else
	   evalin('base','[cr,post] = zineb_post(param.gene.file,0);');
	end
	%eval(cmd) ;
	
	zuireset(h.radio_posttrait)
end
   

% Affichage du Logfile
case 'radio_tracelog'
if ishandle(hfig)
	filename = evalin('base','param.gene.file') ;
 	logf = strcat(filename,'.zineb_logfile');
	if exist(logf)~=0
		[s,t]=unix([getappdata(0,'editeur'),' ',logf,' &']) ;
    		if s ~= 0
    			disp('Error during edition of the logfile :')
    			disp(t)
    		end
	else
		txt=sprintf('the file %s \n does not exist',logf);
		herror = errordlg(txt,'Be care');
		zwaitfor(herror)
	end
	zuireset(h.radio_tracelog)
end
   
% Suivi du Logfile
case 'radio_suivilog'
if ishandle(hfig)
	filename = evalin('base','param.gene.file') ;
 	logf = strcat(filename,'.zineb_logfile');
	if exist(logf)~=0
		[s,t]=	unix([' xterm -geometry 80x24-69+55 -T "fichier de trace" -sb  -e tail -f ', ...
						logf,' &']);
    		if s ~= 0
    			disp('Error during edition of the logfile :')
    			disp(t)
    		end
		root = getappdata(0,'root');
		[s,t]=	unix([' xterm -geometry 80x24-69+50 -T "fichier de trace" -sb  -e ',root,'/op/ztailgrep ',logf,' &']);
    		if s ~= 0
    			disp('Error during edition of the logfile :')
    			disp(t)
    		end
	else
		txt=sprintf('The file \n %s \n doest not exist',logf);
		herror = errordlg(txt,'ATTENTION');
		zwaitfor(herror)
	end
	zuireset(h.radio_suivilog)
end
	
% Chargement du fichier r�ultat
case 'radio_loadresult'
if ishandle(hfig)
	zuiloadresult ;
	zuiuploadform(hfig) ;
	zuireset(h.radio_loadresult)
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
		msgbox(' Sorry, no Help for the moment','Help','help')
		zuireset(h.aide)
	end
		
otherwise
		warning('action not taken into account')
end

% test l'existance d'une variable dans le workspace
function cr = existbase(nom)

cr = evalin('base',strcat('exist(''',nom,''')'));
