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
% version 3.0 , du 17/01/2005.
% 
% liste des modifications : 
%  * 17/01/2005 -> remplace rm par rm -f
%
%--------------------------------------------------------------
%
function [param,data,post] = zimportascii

% etat des donnees
try
    saveok   = evalin('base','param.edit.saveok','0');
catch
    saveok = 0;
end
langue                  = getappdata(0,'langue_cronos');	

% securite anti ecrasement
% Changement pour eviter d'avoir a repondre au 1er chargement
% verification de l'existance des donnees
if (existbase('param')==1) & (existbase('data')==1)
	if saveok == 0
            if strcmp(langue,'francais')
		rep =questdlg('Voulez vous sauver les donnees ?', ...
		              'Attention -> donnees non sauvegardees !', ...
		              'Oui','Non','Abandon','Abandon');
		switch rep
		case 'Oui'
			disp('sauvegarde en cours')
			zuisavedata;
		case 'Abandon'
			return
		case 'Non'
			disp('Donnees ecrasees')
		end
             elseif strcmp(langue,'anglais')
             
 		rep =questdlg('Do you want save data ?', ...
		              'Warning -> unsaved data !', ...
		              'Yes','No','Cancel','Cancel');
		switch rep
		case 'Yes'
			disp('saving in progress')
			zuisavedata;
		case 'Cancel'
			return
		case 'No'
			disp('data overloaded')
		end
            end
	end
	
	% permutation avec jeux1
        if strcmp(langue,'francais')
	  rep =questdlg('Voulez-vous conserver les donnees en memoire comme reference ?', ...
		              'Permutation -> creation d''un jeux de donnees de reference', ...
		              'Oui','Non','Abandon','Non');
	  switch rep
	  case 'Oui'
		disp('permutation en cours')
		evalin('base','jeux1.data = data;','jeux1.data=[];');
		evalin('base','jeux1.param = param;','jeux1.param=[];');
		evalin('base','jeux1.post = post;','jeux1.post=[];');
		evalin('base','t1 = jeux1.data.gene.temps;','t1=[];');
	  case 'Abandon'
		return
	  case 'Non'
		% rien
	  end
        elseif strcmp(langue,'anglais')
	  rep =questdlg('Do you want to keep last data as a reference ?', ...
		              'Permutation -> saving reference data (name : jeux1)', ...
		              'Yes','No','Cancel','No');
	  switch rep
	  case 'Yes'
		disp('permutation in progress')
		evalin('base','jeux1.data = data;','jeux1.data=[];');
		evalin('base','jeux1.param = param;','jeux1.param=[];');
		evalin('base','jeux1.post = post;','jeux1.post=[];');
		evalin('base','t1 = jeux1.data.gene.temps;','t1=[];');
	  case 'Cancel'
		return
	  case 'No'
		% rien
	  end
	end
end

% compatibilite
param = [];
post  = [];
data  = [];

% nom des fichiers
if strcmp(langue,'francais')
  [file,path]=uigetfile('*_ASCII.tar.gz','Nom du fichier a charger ?');
elseif strcmp(langue,'anglais')
  [file,path]=uigetfile('*_ASCII.tar.gz','file name to load ?');
end
if ~ischar(file)
	return
end
filename = fullfile(path,file);
filename = strrep(filename,'/tmp_mnt','');
filename = strrep(filename,'/usr/deneb/cgc','/usr/drfc/cgc');
filetemp = fullfile(path,strcat('copy_',file));

% cp
if strcmp(langue,'francais')
  disp('copie de l''archive');
elseif strcmp(langue,'anglais')
  disp('loading the file');
end
[s,t] = unix(sprintf('cp %s %s',filename,filetemp));
if s ~= 0
   error(t);
end

% gzip 
if strcmp(langue,'francais')
  disp('decompression de l''archive');
elseif strcmp(langue,'anglais')
  disp('unzipping the file');
end

[s,t] = unix(sprintf('gzip -df %s',filetemp));
if s ~= 0
   error(t);
end
filetemp = strrep(filetemp,'.gz','');
if strcmp(langue,'francais')
  disp('desarchivage des fichiers');
elseif strcmp(langue,'anglais')
  disp('untar the file');
end
[s,t] = unix(sprintf('tar xvf %s',filetemp));
if s ~= 0
   error(t);
end
[s,t] = unix(sprintf('rm -f %s',filetemp));
if s ~= 0
   error(t);
end

% racine de la structure
dirroot = strrep(filename,'.tar.gz','');

% liste des fichiers de donnees
if strcmp(langue,'francais')
  disp('lecture de la liste des fichiers');
elseif strcmp(langue,'anglais')
  disp('reading the file list');
end
[s,t] = unix(sprintf('find %s -name "*" -print',dirroot));
if s ~= 0
   error(t);
end
liste = tseparec(t);
ind   = strmatch(fullfile(dirroot,'INFORMATION'),liste,'exact'); 
liste(ind,:)= [];
ind   = strmatch(fullfile(dirroot,'SIZE'),liste,'exact'); 
if ~isempty(ind)
   liste(ind,:)= [];
end
ind   = strmatch(fullfile(dirroot,'CLASS'),liste,'exact'); 
if ~isempty(ind)
   liste(ind,:)= [];
end
liste(1,:)= [];

% lecture des donnees cles
data.gene.temps       = load(fullfile(dirroot,'data_gene_temps'));
param.gene.nbrho      = load(fullfile(dirroot,'param_gene_nbrho'));
param.nombre.fci      = load(fullfile(dirroot,'param_nombre_fci'));
param.nombre.fce      = load(fullfile(dirroot,'param_nombre_fce'));
param.nombre.idn      = load(fullfile(dirroot,'param_nombre_idn'));
param.nombre.hyb      = load(fullfile(dirroot,'param_nombre_hyb'));
param.nombre.glacon   = load(fullfile(dirroot,'param_nombre_glacon'));
param.gene.nbg        = load(fullfile(dirroot,'param_gene_nbg'));
% creation du fichier vide
[cr,data,param]=zinit('',data.gene.temps,param.gene.nbrho,'',[],[], ...
                                 param.nombre.fci,param.nombre.fce,param.nombre.idn, ...
                                 param.nombre.hyb,param.nombre.glacon,param.gene.nbg);

if cr ~=0
	return
end


% boucle sur les fichiers
for k =1:length(liste)
   nomc = liste(k,:);
   nomc(nomc <= ' ') = '';
   [void,varu] = fileparts(nomc);
    if any(varu == '@')
      [varu,indice] =strtok(varu,'@');
      indice(1) =[];
      indice = str2num(indice);
   else
      indice = 1;
   end
   flag = 0;
   if strcmp(varu((end-1):end),'_I')
            flag = 1;
            varu((end-1):end) =[];
   elseif strcmp(varu((end-1):end),'_R')
           flag  =2 ;
           varu((end-1):end) =[];
   end
   varname   = strrep(varu,'_','.');

   % lecture de la donnee
   try
         donnee    = eval(varname);
   catch
         ind = max(find(varname=='.'));
         varname(ind) ='_';
         try
            donnee    = eval(varname);
         catch
               try
                  % nouvelle donnees
                  varname(ind) ='.';
                  eval(sprintf('%s = [];',varname));
                  donnee    = eval(varname);
                  if strcmp(langue,'francais')                  
                     fprintf('nouvelle variable : %s',varname);
                  elseif strcmp(langue,'anglais')  
                    fprintf('new variable : %s',varname);
                  end                
               catch
                  try
                     % nouvelle donnees
                     varname(ind) ='_';
                     eval(sprintf('%s = [];',varname));
                     donnee    = eval(varname);
                     if strcmp(langue,'francais')                  
                        fprintf('nouvelle variable : %s',varname);
                     elseif strcmp(langue,'anglais')  
                       fprintf('new variable : %s',varname);
                     end                
                  catch
                     ind = max(find(varname=='.'));
                      % nouvelle donnees
                     varname(ind) ='_';
                     eval(sprintf('%s = [];',varname));
                     donnee    = eval(varname);
                     if strcmp(langue,'francais')                  
                        fprintf('nouvelle variable : %s',varname);
                     elseif strcmp(langue,'anglais')  
                       fprintf('new variable : %s',varname);
                     end                
                   end
               end
         end   
   end
   disp(varname)
   size_mem   = size(donnee);
   classmem  = class(donnee);
   
   % lecture du fichier
   try 
         donnee_lu = load(nomc);
         type_lu   = 0;
   catch
         [fid,mess]  = fopen(nomc,'r');
         if fid < 3
            error(mess)
         end
         donnee_lu = fread(fid,inf);
         if isempty(donnee_lu)
            donnee_lu = '';
         else
            donnee_lu = tseparec(donnee_lu);
            if size(donnee_lu,1) >1
               donnee_lu(end,:) =[];
            end
         end
        type_lu   = 1;
   end
   
   switch classmem
   
   case 'char'
         if type_lu == 0
            donnee_lu = mat2str(donnee_lu);
         end
         eval(sprintf('%s = donnee_lu;',varname));
   case 'struct'
         if strcmp(langue,'francais')                  
            fprintf('la variable %s de type struct est ignoree\n',varname);
         elseif strcmp(langue,'anglais')  
           fprintf('the variable %s of type struct is ignored\n',varname);
         end                
         
   case 'cell'
         % on essaye le  numerique en premier 
         eval(sprintf('%s{%d} = donnee_lu;',varname,indice));
   case 'single'
         if type_lu == 1
            % donnee compressee
            eval(sprintf('%s = donnee_lu;',varname));
         else
            donnee_lu = single(donnee_lu);
            if all(size_mem >0)
               if prod(size(donnee_lu)) == prod(size_mem)
                  donnee_lu = reshape(donnee_lu,size_mem);
               else
                 if strcmp(langue,'francais')                  
                  fprintf('la variable %s n''a pas les dimensions attendues\n',varname);
                 elseif strcmp(langue,'anglais')  
                  fprintf('dimension mismatch for the variable %s\n',varname);
                 end                
               end  
            end
            switch flag
            case 0
               eval(sprintf('%s = donnee_lu;',varname));
            case 1
               donnee_r = real(eval(varname));
               eval(sprintf('%s = donnee_r + sqrt(-1) .* donnee_lu;',varname));
            case 2
               donnee_i = imag(eval(varname));
               eval(sprintf('%s = donnee_lu + sqrt(-1) .* donnee_i;',varname));
            end
         end
   case 'double'
         if type_lu == 1
            % donnee compressee
            eval(sprintf('%s = donnee_lu;',varname));
        else
           donnee_lu = double(donnee_lu);
           if all(size_mem >0)
               if prod(size(donnee_lu)) == prod(size_mem)
                  donnee_lu = reshape(donnee_lu,size_mem);
               else
                 if strcmp(langue,'francais')                  
                  fprintf('la variable %s n''a pas les dimensions attendues\n',varname);
                 elseif strcmp(langue,'anglais')  
                  fprintf('dimension mismatch for the variable %s\n',varname);
                 end                
               end  
           end
           switch flag
           case 0
              eval(sprintf('%s = donnee_lu;',varname));
           case 1
              donnee_r = real(eval(varname));
              eval(sprintf('%s = donnee_r + sqrt(-1) .* donnee_lu;',varname));
           case 2
               donnee_i = imag(eval(varname));
              eval(sprintf('%s = donnee_lu + sqrt(-1) .* donnee_i;',varname));
           end
        end
   otherwise
         if strcmp(langue,'francais')                  
           disp('Class non prise en compte');
         elseif strcmp(langue,'anglais')
           disp('Class not taken into account');
         end
   end
end

% decompression
if strcmp(langue,'francais')
  disp('decompactage des donnees');
elseif strcmp(langue,'anglais')
  disp('unpack data');
end
data = zreduit(param,data,'uncompact');

% nouvelles donnees
zuisavenonok;


% efface les fichiers temporaires
if strcmp(langue,'francais')
  disp('effacement des fichiers temporaires');
elseif strcmp(langue,'anglais')
  disp('remove tempory file');
end

[s,t] = unix(sprintf('rm -rf %s',dirroot'));
if s ~= 0
   error(t);
end

if strcmp(langue,'francais')
  disp('importation terminee');
elseif strcmp(langue,'anglais')
  disp('happy end');
end

