% calcul les sommes de verification de cronos
function [nom,info] = zineb_cksum(nom)


root = getappdata(0,'root');
if isempty(root)
   zineb_path;
   root = getappdata(0,'root');
end
[s,ho] = unix('hostname');
ho(ho<=32) =[];

if nargin ==0
   nom =[];
end
if isempty(nom)
   nom = fullfile(root,'certification','cksum',sprintf('Cronos_md5sum_%s_%g.%s.txt',date,zinebversion,ho));
end
fprintf('creating signature file : %s\n',nom);
[s,t] = unix(sprintf('date > "%s"',nom));
[s,t] = unix(sprintf('uname -a >> "%s"',nom));
[s,t] = unix(sprintf('whoami >> "%s"',nom));
[s,t] = unix(sprintf('echo "%s" >> "%s"',getappdata(0,'root'),nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
fprintf('creating signature for mfile\n');
[s,t] = unix(sprintf('find %s -name "*.m" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
fprintf('creating signature for f77 files\n');
[s,t] = unix(sprintf('find %s -name "*.f" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
[s,t] = unix(sprintf('find %s -name "*.F" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
fprintf('creating signature for f90 files\n');
[s,t] = unix(sprintf('find %s -name "*.f90" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
[s,t] = unix(sprintf('find %s -name "*.F90" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
fprintf('creating signature for include files\n');
[s,t] = unix(sprintf('find %s -name "*.inc" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
[s,t] = unix(sprintf('find %s -name "*.INC" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
fprintf('creating signature for C files\n');
[s,t] = unix(sprintf('find %s -name "*.c" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
fprintf('creating signature for header files\n');
[s,t] = unix(sprintf('find %s -name "*.h" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
fprintf('creating signature for C++ files\n');
[s,t] = unix(sprintf('find %s -name "*.cc" -exec md5sum {} \\; >> %s',root,nom));
[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
fprintf('creating signature for mexfile files\n');
[s,t] = unix(sprintf('find %s -name "*.mex*" -exec md5sum {} \\; >> %s',root,nom));

[s,t] = unix(sprintf('echo "-------------------------------------------------------" >> %s',nom));
[s,t] = unix(sprintf('md5sum %s >> %s',nom,nom));

[s,info] = unix(sprintf('cat %s',nom));
info = strrep(info,strcat(getappdata(0,'root'),'/'),'');
[fid,mess]  = fopen(nom,'w');
if fid < 3
      error(mess);
else
      fprintf(fid,'%s\n',info);      
      fclose(fid);
end
   
if nargin == 0
   [s,t] = unix(sprintf('gzip -f %s',nom));
   [s,t] = unix(sprintf('chmod a-wx %s.gz',nom));
end
