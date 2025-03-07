% ZDECODEFCE cette fonction decode la consigne de fce
%------------------------------------------------------------------
% fichier zdecodefce.m ->  zdecodefce
%
%
% fonction Matlab 5 :
%
% Cette fonction decode la consigne de puissance de fce
% 
% syntaxe  :
%  
%      [puiss,toro,polo] = zdecodefce(cons);
%
% entree :
%
%      cons    = consigne de puissance de fce (complexe et codee)
% 
% sortie :
%
%     puiss   = puissance (W)
%     toro    = angle toroidal d'injection (degres)
%     polo    = angle poloidal d'injection (degres)
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 15/09/2003.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function [puiss,toro,polo] = zdecodefce(cons)
   
   
   puiss       = abs(cons);
   angle_mul  = angle(cons);
   toro = fix(abs(angle_mul) .* 1e7) .* 1e-4 - 360;
   polo = (abs(angle_mul) .* 1e7 - fix(abs(angle_mul) .* 1e7)) .* 1e3 - 360;


