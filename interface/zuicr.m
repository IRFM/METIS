% ZUICR met dans la structure root les informations destinees a l'assistant
%----------------------------------------------------------------------------
% fichier zuicr.m ->  zuicr
%
%
% fonction Matlab 5 :
% 
% Cette fonction met dans la structure root 
% les informations destinees a l'assistant. 
%
% syntaxe  :
%  
%  zuicr(hform,action)
%
% entrees :
% 
%  hform       = handle du formulaire
%  action      = action associee au callback
%
% sortie : aucune
%
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 08/02/2001.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zuicr(hform,action)
coderetour.handle = hform;
coderetour.action = action;
setappdata(0,'coderetour',coderetour)
