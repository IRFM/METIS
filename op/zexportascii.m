% ZEXPORTASCII export les donnees de Cronos sous forme de fichier ascii
%-----------------------------------------------------------------------
%
% fichier zexportasciit.m
%
% fonction Matlab 5 :
%	Cette fonction exporte les donnees de Cronos sous forme ASCII
%       Chaque donnees est associe a un fichier dont le nom reflette le non de la donnees
%       Les fichiers sont mis dans un repertoire de meme nom que le fichier initial (param.gene.file)
%       Le repertoire est tranforme en archive (tar + gzip).
% 
% syntaxe
%        dirroot = zexportascii(param,data,post,option)
% entrees :
%	param    = structure de parametress de cronos
%       data     = structure de donnees de cronos
%       post     = structure du postprocessing de cronos
%       option   = structure des option de la fonction (cf partie autodeclarante)
%
% sorties :
%       dirroot = nom du fichier de sortie
%
% fonction ecrite par J-F Artaud, poste 62-15
% version  3.0 , du 17/01/2005.
% 
% liste des modifications : 
%  * 17/01/2005 -> remplace rm par rm -f
%
%--------------------------------------------------------------
%
function dirroot = zexportascii(param,data,post,option)
langue                  = getappdata(0,'langue_cronos');

if nargin == 0
	
	valeur.doubleopt   = 0;                        % 
	type.doubleopt     = 'integer';     % type de la donnnee
	borne.doubleopt    = {0,1};         % bornes
	defaut.doubleopt   = 0;           % valeur par defaut
        info.doubleopt   = 'format exportation: 0 -> simple precision, 1-> double precision';

	valeur.liste   = '';                        % 
	type.liste     = 'string';     % type de la donnnee
	borne.liste    = '';         % bornes
	defaut.liste   = '';           % valeur par defaut
        info.liste   = 'liste restrictive des donnees exportees (nom d''une variable type "charcell" du workspace), si vide exporte toute les donnees';

	valeur.tmin    = [];                        % 
	type.tmin     = 'float';     % type de la donnnee
	borne.tmin    =  [-inf,inf];       % bornes
	defaut.tmin   = [];           % valeur par defaut
        info.tmin   = 'temps de debut pour les donnees exporte, si vide commence au premier temps du jeu de donnees';

	valeur.tmax    = [];                        % 
	type.tmax     = 'float';     % type de la donnnee
	borne.tmax    = [-inf,inf];         % bornes
	defaut.tmax   = [];           % valeur par defaut
        info.tmax   = 'temps de fin pour les donnees exporte, si vide finit au dernier temps du jeu de donnees';

	valeur.ope    = '';                        % 
	type.ope     = 'string';     % type de la donnnee
	borne.ope    = '             ';         % bornes
	defaut.ope   = '';           % valeur par defaut
        info.ope   = 'operateur unaire et scalaire applique aux donnees (exemple : max,mean ), si vide  aucune operation applique';


	sortie.description = 'Exporte les donnees de Cronos sous forme de fichiers ASCII place dans une archive (tar + gzip)';   % description (une ligne) de la fonction

         if strcmp(langue,'anglais')	
          info.doubleopt   = 'type of export data: 0 -> single precision, 1-> double precision';
          info.liste   = 'export data list (variable name of  type "charcell" in the workspace), if empty export all datas';
          info.tmin   = 'first time of the export data, if empty, start at lowest time';
          info.tmax   = 'last time of the export data, if empty last at the greatest time';
          info.ope   = 'Operators applied to data (for instance : max,mean ), if empty, no ones';
	  sortie.description = 'Save the whole CRONOS datas in an ASCII file [ASCII + tar + gzip]';   % description (une ligne) de la fonction
        end
	
	
	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	
	sortie.help = 'help zuiexportascii';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	dirroot = sortie;
	
	return

end

% compatibilite
doubleopt = option.doubleopt;
if ischar(option.tmin)
   tmin =[];
else
   tmin      = option.tmin;
end
if ischar(option.tmax)
   tmax =[];
else
   tmax      = option.tmax;
end
option.liste = deblank(option.liste);
if ~isempty(option.liste)
   try
      liste = evalin('base',option.liste');
   catch
    if strcmp(langue,'francais')
      error('La variable qui doit contenir la liste des donnees a exporter n''existe pas');
    end
    if strcmp(langue,'anglais')
      error('unknown list variable');
    end
   end
end
ope = deblank(option.ope);

if isempty(doubleopt)
   doubleopt = 0;
end
if isempty(post)
   post.vide =[];
end
if isempty(tmin)
   tmin = -1e40;
end
if isempty(tmax)
   tmax = 1e40;
end

indtok  = find((tmin <= data.gene.temps) &(tmax >= data.gene.temps));

% operation compression
if strcmp(langue,'francais')
   disp('compactage des donnees')
end    
if strcmp(langue,'anglais')
   disp('data pack operation')
end
data=zreduit(param,data,'compress');

% nom des fichiers
if strcmp(langue,'francais')
  disp('creation de la structure archive')
end
if strcmp(langue,'anglais')
  disp('archive directory creation')
end
if strcmp(param.gene.filetype,'source')
   [pp,dirroot] = fileparts(param.gene.origine);
else
   [pp,dirroot] = fileparts(param.gene.file);
end

dirroot = sprintf('%s_ASCII',dirroot);
% creation du repertoire
if ~isdir(fullfile(pp,dirroot))
   mkdir(pp,dirroot);
end
dirroot = fullfile(pp,dirroot)

% ouverture des fichiers complementaires
fileclass          = fullfile(dirroot,'CLASS');
[fidclass,mess]    = fopen(fileclass,'w');
if fidclass < 3  
   error(mess)
end
filesize          = fullfile(dirroot,'SIZE');
[fidsize,mess]    = fopen(filesize,'w');
if fidsize < 3  
   error(mess)
end

% boucle sur les champ de param
disp('Param');
% liste des champs la structures
champ = fieldnames(param);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('param.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while(~isempty(champ))
    %premier champ de la liste
    champc=champ{1};
    champ(1)=[];
    eval(strcat('test=isstruct(',champc,');'));
    if test
    	% cas d'une sous structure -> ajout de champs
    	eval(strcat('champnew=fieldnames(',champc,');'));   	
    	for k=1:length(champnew)
    		% nom complet pour acceder a la variable
	        champnew{k}=strcat(champc,'.',champnew{k});
	   end
	   % ajout a la liste des champs
	   if isempty(champ)
	      champ =champnew;
	   else
	      champ=cat(1,champ,champnew);
	   end
    else
	   % juste pour les tests
      disp(champc);
      % sauvgarde de la variable
      val =[];
      eval(sprintf('val = %s;',champc),'val =[];');
      zssave(dirroot,val,champc,doubleopt,liste,fidclass,fidsize);
    end
end

% boucle sur les champ de param
disp('Data');
% liste des champs la structures
champ = fieldnames(data);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('data.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while(~isempty(champ))
    %premier champ de la liste
    champc=champ{1};
    champ(1)=[];
    eval(strcat('test=isstruct(',champc,');'));
    if test
    	% cas d'une sous structure -> ajout de champs
    	eval(strcat('champnew=fieldnames(',champc,');'));   	
    	for k=1:length(champnew)
    		% nom complet pour acceder a la variable
	        champnew{k}=strcat(champc,'.',champnew{k});
	   end
	   % ajout a la liste des champs
	   if isempty(champ)
	      champ =champnew;
	   else
	      champ=cat(1,champ,champnew);
	   end
    else
	   % juste pour les tests
      disp(champc);
      % sauvgarde de la variable
      val =[];
      eval(sprintf('val = %s;',champc),'val =[];');
      if ~isempty(val) & isnumeric(val)
         val = double(val);
         ss  = length(size(val));
         switch ss 
         case 2
               val = val(indtok,:);
         case 3
               val = val(indtok,:);
         case 4
               val = val(indtok,:);
         otherwise
               error('dimension > 4');
         end
         if ~isempty(ope)
            switch ope
            case 'mean'
               val = mean(val,1);
            otherwise
               error('operateur non pris en charge');
            end
         end
      end
      zssave(dirroot,val,champc,doubleopt,liste,fidclass,fidsize);
    end
end

% boucle sur les champ de post
disp('Post');
% liste des champs la structures
champ = fieldnames(post);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('post.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while(~isempty(champ))
    %premier champ de la liste
    champc=champ{1};
    champ(1)=[];
    eval(strcat('test=isstruct(',champc,');'));
    if test
    	% cas d'une sous structure -> ajout de champs
    	eval(strcat('champnew=fieldnames(',champc,');'));   	
    	for k=1:length(champnew)
    		% nom complet pour acceder a la variable
	        champnew{k}=strcat(champc,'.',champnew{k});
	   end
	   % ajout a la liste des champs
	   if isempty(champ)
	      champ =champnew;
	   else
	      champ=cat(1,champ,champnew);
	   end
    else
	   % juste pour les tests
      disp(champc);
      % sauvgarde de la variable
      val =[];
      eval(sprintf('val = %s;',champc),'val =[];');
      zssave(dirroot,val,champc,doubleopt,liste,fidclass,fidsize);
    end
end

% ajout des infos relatives au variables
filename = fullfile(dirroot,'INFORMATION');
[fid,mess]  = fopen(filename,'w');
if fid <3 
   error(mess)
end

fprintf(fid,'Informations inside these files are private data and\n');
fprintf(fid,'must not be used in a publication without validation.\n\n\n');

if ~isempty(ope)
   fprintf(fid,'Operateur : %s from %g to %g\n\n\n',ope,tmin,tmax);
end

fclose(fid);
filenamens = fullfile(dirroot,'INFORMATIONNS');
[fid,mess]  = fopen(filenamens,'w');
if fid <3 
   error(mess)
end

% boucle sur les champ de post
disp('Information');
% liste des champs la structures
lm = getappdata(0,'langue_cronos');
setappdata(0,'langue_cronos','anglais');
info = zinfo;
setappdata(0,'langue_cronos',lm);
champ = fieldnames(info);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('info.',champ{k});
end
% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while(~isempty(champ))
    %premier champ de la liste
    champc=champ{1};
    champ(1)=[];
    eval(strcat('test=isstruct(',champc,');'));
    if test
    	% cas d'une sous structure -> ajout de champs
    	eval(strcat('champnew=fieldnames(',champc,');'));   	
    	for k=1:length(champnew)
    		% nom complet pour acceder a la variable
	        champnew{k}=strcat(champc,'.',champnew{k});
	   end
	   % ajout a la liste des champs
	   if isempty(champ)
	      champ =champnew;
	   else
	      champ=cat(1,champ,champnew);
	   end
    else
	   % juste pour les tests
      disp(champc);
      % sauvgarde de la variable
      val =[];
      eval(sprintf('val = %s;',champc),'val =''unavaible'';');
      fprintf(fid,'{ %s }->\t\t%s\n',strrep(strrep(champc,'info.',''),'.','_'),val); 
    end
end

fclose(fid);
[s,t]= unix(sprintf('sort %s >> %s',filenamens,filename));
if s ~= 0
   error(t);
end
delete(filenamens);

% tar of directory
if strcmp(langue,'francais')
  disp('Creation de l''archive ...');
end
if strcmp(langue,'anglais')
  disp('Writing tar file ...');
end

[s,t] = unix(sprintf('tar cvf %s.tar %s',dirroot,dirroot));
if s ~= 0
   error(t);
else
   [s,t] = unix(sprintf('rm -rf %s',dirroot));
end
dirroot = strcat(dirroot,'.tar');
% gzip 
if strcmp(langue,'francais')
  disp('Compression de l''archive ...');
end
if strcmp(langue,'anglais')
  disp('zipping the file ...');
end

[s,t] = unix(sprintf('gzip -qf9 %s',dirroot));
if s ~= 0
   error(t);
end
dirroot = strcat(dirroot,'.gz');

fclose(fidsize);
fclose(fidclass);
if strcmp(langue,'francais')
  fprintf('Archive %s cree.\n',dirroot);
end
if strcmp(langue,'anglais')
  fprintf('file %s created.\n',dirroot);
end


function zssave(dirroot,val,champ,doubleopt,liste,fidclass,fidsize)

if ~isempty(liste)
   if isempty(strmatch(champ,liste))
      return
   end
end
nomvar = strrep(champ,'.','_');

% nom des fichiers
filename   = fullfile(dirroot,nomvar);


% info complementaires
fprintf(fidclass,'class.%s = ''%s'';\n',champ,class(val));
fprintf(fidsize,'sizeof.%s = %s;\n',champ,mat2str(size(val)));



if isempty(val)
   [s,t] = unix(sprintf('touch %s',filename));
   if  s ~= 0
      error(t);
   end
elseif ischar(val)
   [fid,mess]  = fopen(filename,'w');
   if fid < 3 
      error(mess);
   end
   for k=1:size(val,1)
      fprintf(fid,'%s\n',val(k,:));
   end
   fclose(fid);
elseif iscell(val)
   val = val(:);
   for k = 1:length(val)
      zssave(dirroot,val{k},sprintf('%s@%d',champ,k),doubleopt,{},fidclass,fidsize);
   end
elseif iscomplex(val)
   val = double(val);
   ss = size(val);
   if length(ss) > 2
      for k = 1:length(ss)
         filename = sprintf('%s_%d',filename,ss(k));
      end
      val = val(:);
   end

   vali = imag(val);
   valr = real(val);

   ind = find(~isfinite(vali));
   if ~isempty(ind)
      vali(ind) = sign(vali(ind)) * 1e38;
   end
   if  doubleopt
      save(strcat(filename,'_I'),'vali','-ASCII','-DOUBLE');
   else
      save(strcat(filename,'_I'),'vali','-ASCII');
   end

   ind = find(~isfinite(valr));
   if ~isempty(ind)
      valr(ind) = sign(valr(ind)) * 1e38;
   end
   if  doubleopt
      save(strcat(filename,'_R'),'valr','-ASCII','-DOUBLE');
   else
      save(strcat(filename,'_R'),'valr','-ASCII');
   end
elseif isnumeric(val)
   val = double(val);
   ss = size(val);
   if length(ss) > 2
      for k = 1:length(ss)
         filename = sprintf('%s_%d',filename,ss(k));
      end
      val = val(:);
   end
   
   ind = find(~isfinite(val));
   if ~isempty(ind)
      val(ind) = sign(val(ind)) * 1e38;
   end
   if  doubleopt
      save(filename,'val','-ASCII','-DOUBLE');
   else
      save(filename,'val','-ASCII');
   end
else 
  if strcmp(langue,'francais')
    fprintf('la donnee %s ne peut pas etre sauver en ASCII\n',champ);
  end
   if strcmp(langue,'anglais')
    fprintf('bad news, you are unlucky, data can not be saved in ASCII\n',champ);
  end
 
end
