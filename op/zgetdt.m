function [dt,tmin,tmax,res,restot] = zgetdt(file)
% - =45, e =101, 0 = 47 .. 9 = 57

if nargin < 1
   error('il faut donner le nom du fichier ...')
elseif isempty(file)   
   error('il faut donner le nom du fichier ...')
end

% retrait des extention
[path,file,ext,ver] =fileparts(file);
[vpath,file,vext,vver] =fileparts(file);
%nom avec bonne extention
file = fullfile(path,[file '.zineb_logfile']);

% lecture des infos
[s,t]=unix([ 'grep "dt =" ',file]);
if s ~= 0
  disp('erreur unix :')
  disp(t)
  return
end

% retrait du texte
t =strrep(t,'de','');
at = abs(t);
ind = find(((at>=45)&(at<=57)) | (at <=32)| (at == 101));
t=t(ind);

% separtion des lignes
texte = tseparec(t);

% transformation en matrice
res = str2num(texte);

% supression des retour en arriere
indsup = find(diff(res(:,2)) == 0);
restot = res;
res(indsup,:) = [];
dt   = res(:,1);
tmin = res(:,2);
tmax = res(:,3); 


