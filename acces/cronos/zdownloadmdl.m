% cette fonction charge dans le workspace les donnees pour une simulation simulink dans Cronos
function [cr,simu] = zdownloadmdl(simu)

% sortie
cr = 0;

% test des entrees
if isempty(simu)
   return
elseif simu.onoff == 0
   return
elseif isempty(simu.mdlname)
   return
end

% to make it working without display with matlab2012b
eval('[cr,simu] = zdownloadmdl_real(simu);');

