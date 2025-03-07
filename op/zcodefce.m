% ZCODEFCE cette fonction code la consigne de fce
%------------------------------------------------------------------
% fichier zcodefce.m ->  zcodefce
%
%
% fonction Matlab 5 :
%
% Cette fonction code la consigne de puissance de fce
% 
% syntaxe  :
%  
%      [cons] = zcodefce(puiss,toro,polo)
%
% entree :
%
%     puiss   = puissance (W)
%     toro    = angle toroidal d'injection (degres)
%     polo    = angle poloidal d'injection (degres)
% 
% sortie :
%
%     cons    = consigne de puissance de fce (complexe et codee)
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 15/09/2003.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function cons = zcodefce(puiss,toro,polo)
   
   angle_mul    = (360 + rem(polo,360)) .* 1e-10  + ...
                  fix((360 + rem(toro,360)) .* 1e4) .* 1e-7;
   cons         = puiss .* exp(sqrt(-1) .* angle_mul);


