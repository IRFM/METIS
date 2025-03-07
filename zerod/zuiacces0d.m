% interface graphique pour z0dacces
function zuiacces0d(mode)

if nargin == 1 
   info = z0dacces;
   assignin('base','op0d',info.valeur);
   zuicreefunform('z0dacces','op0d',1,0,'zuiacces0d');
   return
end

evalin('base','[cr,data,param] = z0dacces('''',op0d,post);');
zuisavenonok
disp('Data ready !');
