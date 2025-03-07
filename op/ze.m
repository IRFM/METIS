% script de localisation des erreurs en mode debug
function ze(nh)
%
% appel a l'historique
%
if nargin == 0
   rep = dbstack;
   last = lasterr;
   lasterr('');
else
   liste = getappdata(0,'error_history');
   if isempty(liste)
      return
   end
   k     = min(1, length(liste) - nh);
   rep   = liste{k}.rep;
   last  = lsite{k}.text; 
end
% 
% cas workspace
%
if length(rep) <2
   disp('Pas de point d''arret en cours ....')
   return
end
%
% recherche de la premiere fonction nom matlab
%
for k =2:length(rep)
   if isempty(findstr(which(rep(k).name),matlabroot))
      l = k;
      break
   end
end

disp('----------------------------');
fprintf('Erreur dans la fonction :\n%s @ %d \n%s \n',rep(l).name,rep(l).line,last);
disp('Pile : ')
for k =2:length(rep)
   fprintf('-> %s @ %d\n',rep(k).name,rep(k).line); 
end
disp(' ')
if any(rep(l).name == '(')
   zne(rep(l+1).name,rep(l).line);
else
   zne(rep(l).name,rep(l).line);
end

if nargin > 0
   return
end

hist.text  = last;
hist.stack = rep;
if isappdata(0,'error_history')
   liste = getappdata(0,'error_history');
   liste{end+1} = hist;
   setappdata(0,'error_history',liste);
else
   liste = {hist};
   setappdata(0,'error_history',liste);
end 
lasterr('');
