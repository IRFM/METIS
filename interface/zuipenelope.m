% interface graphique pour zpenelope
function zuipenelope

% test de donnees a sauvegardeees
if (existbase('param')==1) & (existbase('data')==1) & (0>1)
   try
			saveok = evalin('base','param.edit.saveok') ;
	catch
	      saveok = 1;
   end
	if saveok == 0
		rep =questdlg('do you want to come back to the main menue ?', ...
	         	     'be careful -> unsaved data !', ...
	         	     'Yes','No','No');
		switch rep
		case 'yes'
			return
		end
	end
end


% appel de l'interface du module zpenelope
info = zpenelope;
zassignin('base','option',info.valeur);
h=zuicreefunform('zpenelope','option',1);
set(h,'name','Lecture des donnees Tconos');
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end


hdlg = msgbox('access to the database (Tcronos, Tprof, ...)','be patient ...','help');
drawnow;
evalin('base','[param,data,post,cr]=zpenelope(option);');
if ishandle(hdlg)
	zuicloseone(hdlg);
end
% test creation donnees
ok = evalin('base','isstruct(param)&isstruct(data)');
if ok == 0
	warndlg('No data avialable','zpenelope error');
	return
end
cr = evalin('base','cr','NaN');
if cr ~= 0
	warndlg('data access problem','zpenelope');
	return
end

% position flag d'edition
zuisavenonok;
evalin('base','param.edit.currentfile=param.gene.file;');

% mise a jour menu principal de cronos
[hfig,h] = zuiformhandle('direct');
if ishandle(hfig)
	zuiuploadform(hfig) ;

	numshot=evalin('base','param.from.shot.num','[]');
	if ~isempty(numshot)
		numchoc = fix(numshot) ;
		zuidata(h.text_numchoc,num2str(numchoc));

		occurence = round((numshot-numchoc)*10);
		zuidata(h.text_occurence,num2str(occurence));

		tps1=evalin('base','param.gene.tdeb');
		tps2=evalin('base','param.gene.tfin');
		tps=[num2str(tps1) ' �' num2str(tps2)];
		zuidata(h.text_temps,tps);
	end
end




% test l'existance d'une variable dans le workspace
function cr = existbase(nom)

cr = evalin('base',strcat('exist(''',nom,''')'));
