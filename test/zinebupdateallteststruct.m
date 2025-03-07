% zinebupdateallteststruct met a niveau la structure de donnee de tous les tests CRONOS
%------------------------------------------------------------------------------- 
% fichier :  zinebupdateallteststruct.m  ->  zinebupdateallteststruct
% 
% 
% fonction Matlab 5 : 
% 
% Cette fonction met a niveau la structure de donnee de tous les tests CRONOS
% La mise a jour est interractive pour les simulations completes cronos qui servent
% pour les tests. Il faut modifier les parametres si necesaire et resauver
%  
% syntaxe :  
%   
%        zinebupdateallteststruct
%  
% fonction ecrite par Artaud , poste 62.15 
% version  4.1 du 08/08/2008  
%  
% liste des modifications :  version CVS
%  
%-------------------------------------------------------------------------------  
%  
function zinebupdateallteststruct

% racine CRONOS
root = getappdata(0,'root');
if isempty(root)
	zineb_path;
	root = getappdata(0,'root');
end

% creation temporaire de la variable globale CRONOS_WORK_DIR dans le cas
% des tests
filename=sprintf('%s/certification',root);
setappdata(0,'CRONOS_WORK_DIR',filename);

disp('Start update data structures for CRONOS tests')
disp('update data structures for CRONOS run test')
% Mise a jour des simulmations completes -> point de depart
r = dir(fullfile(root,'certification','cronostestrun','*.mat.gz'));
% boucle sur la liste
for k =1:length(r)
	fprintf('update data structures of %s\n',r(k).name)
	evalin('base','clear data param post'); 
	evalin('base','post = [];');
	evalin('base','zuidirect;');
	zuiload(fullfile(root,'certification','cronostestrun',r(k).name));
	zuisavenonok
	evalin('base','[cr,param,data,post] = zineb_vup(param,data,post);');
	h = msgbox('You have to anwser to dialog boxes and save the CRONOS simulation to finish the update','Instructions');	
	while(evalin('base','param.edit.saveok','0') == 0)
		pause(0.1)
	end
	try
		close(h);
	end
end
% Mise a jour des simulmations completes -> resultat
r = dir(fullfile(root,'certification','cronostestrun_result','*.mat.gz'));
% boucle sur la liste
for k =1:length(r)
	fprintf('update data structures of %s\n',r(k).name)
	evalin('base','clear data param post'); 
	evalin('base','post = [];');
	evalin('base','zuidirect;');
	zuiload(fullfile(root,'certification','cronostestrun_result',r(k).name));
	zuisavenonok
	evalin('base','[cr,param,data,post] = zineb_vup(param,data,post);');
	h = msgbox('You have to anwser to dialog boxes and save the CRONOS simulation to finish the update','Instructions');
	while(evalin('base','param.edit.saveok','0') == 0)
		pause(0.1)
	end	
        try
		close(h);
	end

end
% Mise a jour des tests complets
disp(' ')
disp(' ')
disp('update data structures for CRONOS complet tests')
r = dir(fullfile(root,'certification','fullruns','*.mat'));
% boucle sur la liste
for k =1:length(r)
	fprintf('update data structures of %s\n',r(k).name)
        zinebupdatestruct(r(k).name)
end
% Mise a jour des tests de fonction CRONOS
disp(' ')
disp(' ')
disp('update data structures for CRONOS function tests')
r = dir(fullfile(root,'certification','cronosfunctions','*.mat'));
% boucle sur la liste
for k =1:length(r)
	fprintf('update data structures of %s\n',r(k).name)
        zinebupdatestruct(r(k).name)
end
% Mise a jour des tests de module
disp(' ')
disp(' ')
disp('update data structures for CRONOS module tests')
r = dir(fullfile(root,'certification','modules','*.mat'));
% boucle sur la liste
for k =1:length(r)
	fprintf('update data structures of %s\n',r(k).name)
        zinebupdatestruct(r(k).name)
end
disp('End update data structures for CRONOS tests')

