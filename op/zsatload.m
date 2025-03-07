% regarde la charge de la ferme
function [node,names,uload] = zsatload

% recherche des noeux
node = [];
names = {};
uload = [];
fprintf('\n@:');
for k = 1:32   
   % ping -> neoud actif
   fprintf('.');
   [s,t] = system(sprintf('ping -q -n -c 1 saturne%d',k));
   if s == 0
      % le noeud est actif
      node(end+1) = k;
      [p,r] = strtok(t);      
      names{end+1 } =strtok(r); 
      % test de charge
%      [s,t] = system(sprintf('rsh saturne%d "cat /proc/loadavg"',k));
       [s,t] = system(sprintf('cat /proc/loadavg',k));
     if s == 0 
         [p1,r] = strtok(t);      
         [p2,r] = strtok(r);      
         uload(end+1)   = max([str2num( strtok(r)),str2num(p2),str2num(p1)]) ;
      else
           disp(t);
           uload(end+1) = NaN;
      end
   end
end
fprintf('|\n');

