% ZUILOAD   charge un fichier de resultats.
%----------------------------------------
% fichier zuiload.m ->  zuiload
%
% fonction Matlab 5 :
%	Cette fonction charge un fichier de resultats. 
%
% syntaxe  :
%	zuiloadresult;
%
% entrees :  
%  
% sorties :  
%
% fonction ecrite par Ch. Passeron , poste 61-19
% version 3.0, du 17/01/2005.
% 
% liste des modifications : 
%	* 13/06/2001 -> ajout de la compression/decompression
%	* 09/07/2001 -> correction bug sur compression
%	* 14/09/2001 -> ajout de la gestion du type source/resultat
%	* 14/09/2001 -> ajout de la gestion du jeux de donnees de reference (jeux1)
%	* 14/09/2001 -> ajout de la gestion de la cohernce des noms
%	* 14/09/2001 -> ajout du test des versions 
%	* 20/11/2001 -> ajout de la reconstruction partielle pour la
%                  consultaion des resultats intermdiaires
%	* 05/02/2002 -> efface le nom du fichier si resultat intermediaire
%       * 17/01/2005 -> remplace rm par rm -f
%
%--------------------------------------------------------------
%
function zuiloadresult

% verification des champs origine et file
try 
   type = evalin('base','param.gene.filetype');
catch
   type = '';   
end
if strcmp(type,'resultat')
%	warndlg('You have already an output file', ...
%	        'No output file to load');
    ButtonName = questdlg('You have already an output file',...
                          'action','load again','cancel','return');
    switch ButtonName,
        case 'load again'
        case 'cancel'
	      return
    end
end
	       

% permutation avec jeux1
rep =questdlg('Do you want to store the present data as the reference dataset ??', ...
              'Reference dataset', ...
              'Yes','No','Cancel','No');
switch rep
case 'Yes'
	evalin('base','jeux1.data = data;','jeux1.data=[];');
	evalin('base','jeux1.param = param;','jeux1.param=[];');
	evalin('base','jeux1.post = post;','jeux1.post=[];');
	evalin('base','t1 = jeux1.data.gene.temps;','t1=[];');
case 'Cancel'
	return
case 'No'
	% rien
end


saveok   = evalin('base','param.edit.saveok');

% securite anti ecrasement des donnees en cours
if saveok == 0
	rep =questdlg('Do you want to save the data ?', ...
	              'Be careful -> unsaved data !', ...
	              'Yes','No','Cancel','Cancel');
	switch rep
	case 'Yes'
		disp('saving in progress')
		zuisavedata;
	case 'Cancel'
		return
	case 'No'
		disp('data removed from workspace')
	end
end

% recupere le nom du fichier de resultat
filename = evalin('base','param.gene.file','[]');
%afile    = filename;
%bfile    = evalin('base','param.edit.currentfile','[]');

% retrait des .mat .gz
[pathf,filef,ext1,ver1] = fileparts(filename);
[vpathf,filef,ext2,ver2] = fileparts(filef);
% creation des noms utiles
filename1 = fullfile(pathf,[filef '.mat']);
filename2 = fullfile(pathf,[filef '.mat.gz']);
filename  = fullfile(pathf,filef);


% extrait le chemin
if isempty(filename) |((exist(filename1) ~= 2) & (exist(filename2) ~= 2))
	txt=sprintf(['the file \n %s \n does not exist\n', ... '
	             'Do you want to load intermediate output ? '],filename);
	%herror = errordlg(txt,'ATTENTION');
	%zwaitfor(herror)
	rep =questdlg(txt, ...
	              'Be careful -> No output ''files'' avialable', ...
	              'Yes','No','No');
	switch rep
	case 'Yes'
		evalin('base','[cr,param,data] = zineb_rebuilt(param.gene.origine);');
		cr = evalin('base','cr','-1');
		if cr == 0
			zassignin('base','param.edit.currentfile','');
                        zassignin('base','param.edit.saveok',0);
		end
	end
	return
end


% chargement
[lfile,rmfile,cr]=zgzip(filename,'uncompress');
evalin('base',strcat('load(''',lfile,''')')) ;
if ~isempty(rmfile)
	[voids,voidt]=unix(['rm -f ',rmfile,' >& /dev/null']);
end
%evalin('base',strcat('load(''',filename,''')'));
evalin('base','data = zreduit(param,data,''uncompact'');');

% modifie les champs
zassignin('base','param.edit.currentfile',filename);
zassignin('base','param.edit.saveok',1);
zassignin('base','void',[]);

% verification de la coherence des noms
cr = zveriffilename(filename);


% verifcation de la coherence des versions
zverifversion;

