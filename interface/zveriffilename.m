% ZVERIFFILENAME  verifie la coherence entre le nom du fichier et les champs file et origine
%------------------------------------------------------------------------------- 
% fichier : zveriffilename.m  ->  zveriffilename
% 		                            zantiext : supprime les extensions du nom du fichier
% 
% fonction Matlab 5 : 
%	verifie la coherence entre le nom du fichier et les champs file et origine  
%
% syntaxe :  
%   cr = zveriffilename(fichier) 
%  
% entrees :  
%	fichier : nom du fichier a traiter
%  
% sorties :  
%   cr = code retour
%  
% fonction �rite par J-F Artaud , poste 46-78
% version  1.7  du  29/09/2001  
%  
%  
% liste des modifications :  
%  
%-------------------------------------------------------------------------------  
%  

function cr = zveriffilename(fichier)

cr =0;

if nargin < 1
	fichier ='';
end
if isempty(fichier)
	disp('Il faut donner un nom de fichier');
	cr = -1;
	return
end


% verification des champs origine et file
try 
   type = evalin('base','param.gene.filetype');
catch
   type = [];   
end

try
   file    = zantiext(evalin('base','param.gene.file'));
   origine = zantiext(evalin('base','param.gene.origine'));
catch
   disp('pas de donnees dans le workspace')
   cr = -2
   return
end

fichier = zantiext(fichier);
fichier = strrep(fichier,'/tmp_mnt','');

if strcmp(origine,fichier)
   % c'est un fichier source (ok)
	if isempty(type)
		zassignin('base','param.gene.filetype','source');
	end
elseif strcmp(file,fichier)	
   % c'est un fichier resultat (ok)
	if isempty(type)
		zassignin('base','param.gene.filetype','resultat');
	end
elseif ~isempty(type)
	% il faut corriger
	if strcmp(type,'source')
		zassignin('base','param.gene.origine',fichier);
	elseif strcmp(type,'resultat')
		zassignin('base','param.gene.file',fichier);
	end
	zuisavenonok;
else
	% correction auto impossible
	% il faut demander
	kmin = evalin('base','param.gene.kmin');
	k    = evalin('base','param.gene.k');
	if  k  > kmin
		% permutation avec jeux1
		rep =questdlg('Cronos file :', ...
		              'file type ?', ...
		              'Input','Output','Cancel','Output');
	else
		% permutation avec jeux1
		rep =questdlg('Cronos file :', ...
		              'file type ?', ...
		              'Input','Output','Cancel','Output');
	end
	switch rep
	case 'Input'
		zassignin('base','param.gene.filetype','source');
		zassignin('base','param.gene.origine',fichier);
		zuisavenonok;
	case 'Output'
		zassignin('base','param.gene.filetype','resultat');
		zassignin('base','param.gene.file',fichier);
		zuisavenonok;
	case 'Cancel'
		% rien
	end

end


% supprime les extention du nom du fichier
function file = zantiext(file)

e='?';
while (~isempty(e))
	[p,f,e] = fileparts(file);
	file    = fullfile(p,f);
end
