%fonction donnant le repertoire de travail 
%option = 0 ou 1
%si 0 repertoire de travail specifie par arch.inc
%si 1 repertoire du fichier input (a priori partage) ou repertoire temporaire partage par tous les noeuds 

function jobdir =zineb_tempdir(option)

jobdir=getappdata(0,'TEMPDIR');

if (~isdir(jobdir)) || (option == 1)
    if isappdata(0,'TEMPDIR_EXCHANGE')
   	 jobdir=getappdata(0,'TEMPDIR_EXCHANGE');
	 if ~isdir(jobdir)
    		jobdir=getappdata(0,'CRONOS_WORK_DIR');
	 end
    else
   	 jobdir=getappdata(0,'CRONOS_WORK_DIR');
    end
end

% securite pour execution hors cronos
if ~isdir(jobdir)
  homedir = getenv('HOME');
  jobdir  = fullfile(homedir,'tmp');
  if ~isdir(jobdir)
      mkdir(homedir,'tmp');
      fprintf('Creating directory %s for shared tempory files\n', jobdir);
  end
end
if ~isdir(jobdir)
    jobdir = tempdir;
end