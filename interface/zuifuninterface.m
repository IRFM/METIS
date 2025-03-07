% ZUIFUNINTERFACE surcouche a zuicreefunform pour les modules extene de cronos
%-------------------------------------------------------------------------------------
% fichier zuifuninterface.m ->  zuifuninterface
%
% fonction Matlab 5 :
%	Cette fonction cree les formulaires associes aux modules extene de cronos 
%	en appelant zuicreefunform. 
%
% syntaxe  :
%   zuifuninterface(clef)
%
% entrees :
%	clef      = clef associe au module (fci,hyb,fce,idn, ...)
%
% sortie : 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.6, du 28/08/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zuifuninterface(clef)

% test des entrees
if nargin == 0
   return
elseif isempty(clef)
   return
end

% recupere les donnees
try 
    fonction = evalin('base',strcat('param.fonction.',clef));
catch
    return
end
if isempty(fonction)
    return
end
try
    nombre = evalin('base',strcat('param.nombre.',clef));
catch
    nombre =1;
end
if isempty(nombre)
    nombre =1
end
    
racine = strcat('param.cons.',clef);

% appel du generateur d'interface dynamique
zuicreefunform(fonction,racine,nombre);
