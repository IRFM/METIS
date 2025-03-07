% ZGZIP comprime ou decomprime avec gzip les fichier de donnees cronos/zineb
%-----------------------------------------------------------------------------
% fichier zgzip.m ->  zgzip
%
%
% fonction Matlab 5 :
%
% Cette fonction sert a comprimer ou decomprimer avec gzip
% les fichier de donnees cronos/zineb. Elle est securis???e.
%
% syntaxe  :
%
%     [file_out,file_temp,cr]=zgzip(file,mode);
%
% entree :
%
%     file  = nom du fichier matlab a comprimer ou a decomprimer
%
%     mode  = 'compress'   -> comprime le fichier
%             'uncompress' -> decomprime le fichier
%
%     filesturne = nom du fichier temporaire pour ne pas passer pas les montages %                  NFS
%
% sorties :
%
%    file_out    =  nom du fichier produit par la fonction
%                   (pour le charger apres une decompression)
%
%    file_temp   =  nom du fichier temporaire a supprimer apres un load
%                   si necessaire;
%
%    cr          = compte rendu d'execution (0 = ok)
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 17/01/2005.
%
%
% liste des modifications :
%
%     * 20/06/2001 -> correction du bug sur file_temp (sans incidence)
%     * 24/07/2001 -> acceleration du gzip pour cgc (solution provisoire)
%     * 17/01/2005 -> remplace rm par rm -f
%     * 03/05/2005 -> mode pour saturne
%     * 28/05/2005 -> mise a niveau pour saturne
%
%--------------------------------------------------------------
%
function [file_out,file_temp,cr]=zgzip(file,mode,filesaturne)

cr =0;
file_temp ='';

% test des entrees
if nargin <2
   error('ZGZIP :Il faut donne un nom de fichier et un mode')
end
if isempty(file)
   warning('ZGZIP :Nom de fichier vide')
   return
end
if isempty(mode)
   warning('ZGZIP : "mode" vide')
   return
end

% decomposition du nom de fichier
[pathf,name,ext] = fileparts(file);
[vp,name,ext_int ] = fileparts(name);



if strcmp(mode,'compress')
   % la compression n'est plus necessaire en matlab 7.x
   if ~verLessThan('matlab','7.5')
        file_out = file;
        return
   else
   	% test
	if strcmp(ext,'.gz')
	disp('ZGZIP : Le fichier est deja comprime')
	file_out = file;
	file_temp ='';
	cr = 100;
	return
	end
	if isempty(ext)
	ext ='.mat';
	end
	
	% nom complet
	ffile  =fullfile(pathf,[name ext ver]);
	
	% compression du fichier via gzip
	%[s,t] = unix(['gzip -9qf ',ffile,' >& /dev/null']);
	[s,t] = unix(['mv ',ffile,'.gz ',ffile,'.gz.bck']);
		% j'ai supprimer le -9 qui peut poser des probleme sous tru64
	if nargin == 3
	[s,t] = unix(['gzip -f ',filesaturne,'.mat']);
	if s == 0
		[s,t] = unix(['mv -f ',filesaturne,'.mat.gz ',ffile,'.gz']);
	end
	
	else
	[s,t] = unix(['gzip -f ',ffile]);
	end
	
	% traitement des erreurs
	if s ~= 0
	[s,t] = unix(['mv ',ffile,'.gz.bck ',ffile]);
	cr = s;
	disp('ZGZIP : erreur lors de la compression ->')
	disp(t)
	file_out = file;
			file_temp ='';
	return
	else
	%     [vs,t] = unix(['echo rm -f ',ffile,'.bck >& /dev/null'])
	[vs,t] = unix(['rm -f ',ffile,'.gz.bck >& /dev/null']);
	file_out = fullfile(pathf,[name '.mat.gz' ver]);
			file_temp ='';
	end
   end
elseif strcmp(mode,'uncompress')
   if (strcmp(ext,'.mat') && (exist(file) == 2)) || (isempty(ext) && (exist(strcat(file,'.mat')) == 2))
	% c'est un fichier matlab 
	% pas de decompression
	file_out = file;
	file_temp = '';
	cr = 0;
	return
   else
   
	%  numero propre pour les nom de fichier
	if isappdata(0,'pid')
		pid = getappdata(0,'pid');
	else
		pid = [];
	end
	if isempty(pid)
		[pid,pid_sess,user]=getidprocess;
		setappdata(0,'pid',pid);
	end
	% nom de fichier temporaire
	%tmp = getenv('HOME');
	%if strmatch('saturne',getenv('HOSTNAME'))
	%if isdir('/scratch')
	%    tmp = '/scratch';
	%elseif strcmp(getenv('USER'),'cgc')
	%tmp = strcat(tmp,'/cgc_data/tmp');  %mis en commentaire FI 19/02/03
	%            tmp = '/usr/drfc/zineb/data';
	%end
	tmp=zineb_tempdir(0);
	
	file_temp = fullfile(tmp,sprintf('zineb%d.mat.gz',pid));
	
	% nom du fichier a copier
	file_cp = fullfile(pathf,[name '.mat.gz']);
	
	% copie du fichier comprimer
	[s,t] = unix(['rm -f ',file_temp,' >& /dev/null']);
	[s,t] = unix(['cp ',file_cp,' ',file_temp]);
	[s,t] = unix(['chmod 666 ',file_temp]);
	if s ~= 0
		% le fichier n'est pas comprimer
		cr = 0;
		file_out = file;
		file_temp ='';
		return
	else
		% il faut decomprimer le fichier
		%[s,t] = unix(['gzip -dqf ',file_temp,' >& /dev/null']);
		[s,t] = unix(['gzip -dqf ',file_temp]);
		% traitement des erreurs
		if s  ~= 0
			cr = s;
			disp('ZGZIP : erreur lors de la decompression ->')
			disp(t)
			[s,t] = unix(['rm -f ',file_temp]);
			file_out = file;
			return
		else
			file_out = fullfile(tmp,sprintf('zineb%d.mat',pid));
			file_temp = file_out;
		end
	end
   end
else
   disp('ZGZIP : "mode" inconnu')
   file_out = file;
   cr = 101;
end
