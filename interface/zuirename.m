% ZUIRENAME change le nom des fichiers et sauve les donnees avec ce nouveau nom
%--------------------------------------------------------------------------------------
% fichier zuirename.m ->  zuirename
%
% fonction Matlab 5 :
%	Cette fonction  change le nom des fichiers et sauve les donnees avec ce nouveau nom . 
%
% syntaxe  :
%	zuisavedata;
%
% entrees :
%
% sorties :
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.0, du 11/12/2002.
% 
% liste des modifications : 
%  * 25/01/2002 -> ajout de la question sur le type de fichier
%  * 11/03/2002 -> effacement du commenatire lors du changement de nom
%  * 11/12/2002 -> interface en anglais
%
%--------------------------------------------------------------
%
function zuirename

% recupere le nom du fichier de sauvegarde
filename = evalin('base','param.edit.currentfile','[]');
saveok   = evalin('base','param.edit.saveok','0');
%filetype_mem = evalin('base','param.gene.filetype','0');
commem     = evalin('base','param.from.creation.com');
disp('old comments :')
disp(commem)
usermem     = evalin('base','param.from.creation.user');
%  %type du jeux de fichier
%  switch filetype_mem
%  case 'source'
%     stype = 'Input';
%  otherwise
%     stype = 'Output';
%  end
%  
%  
%  rep =questdlg(sprintf('Please select file type (%s):',stype), ...
%  	      'File type ?', ...
%  	      'Input','Output','Cancel',stype);
%  switch rep
%  case 'Input'
%    filetype  = 'source';
%  case 'Output'
%    filetype  = 'resultat';
%  case 'Cancel'
%    return
%  otherwise
%     % rien
%  end

% lorsque le fichier change de nom il devient obligatorirement un fichier source
% cela permet d'eviter les conflits et d'avoir une association propre source -> resultat
filetype = 'source';

% changement de numero de choc pour le machine n'ayant pas de vrai numero
% ceci evite de vori apparaitre toutes les simulations avec le meme numero dans la base de donnees
tokname = evalin('base','param.from.machine','[]');
nummem = [];
switch upper(tokname) 
case {'ITER','DEMO','JT60_SA','FAST'}
	nummem = evalin('base','param.from.shot.num','[]');
	if ~ isempty(nummem)
		numnew = fix(datenum(clock).*10000);
		zassignin('base','param.from.shot.num',numnew);
	end
end

% modifie les champ
zassignin('base','param.edit.currentfile','');
zassignin('base','param.edit.saveok',0);
zassignin('base','param.gene.filetype',filetype);
zassignin('base','param.from.creation.com','');
% sauvegarde
zuisavedata;

% restore si annulation
if isempty(evalin('base','param.edit.currentfile','[]'))
	zassignin('base','param.edit.currentfile',filename);
	zassignin('base','param.edit.saveok',saveok);
	zassignin('base','param.gene.filetype',filetype_mem);
   	zassignin('base','param.from.creation.com',commem);
   	zassignin('base','param.from.creation.user',usermem);
	if ~isempty(nummem)
   		zassignin('base','param.from.shot.num',nummem);
	end
end
