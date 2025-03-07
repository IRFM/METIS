%  ZLISTCOEF  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zlistcoef.m  ->  zlistcoef 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   liste = zlistcoef 
%  
% entrees :  
%  
%  
% sorties :  
%   liste  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  1.9  du  11/03/2002  
%  
%*#auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function liste = zlistcoef

s  = what(fullfile(getappdata(0,'root'),'coef'));
pm     = strrep(s.m,'.m','');
pmex   = strrep(s.mex,strcat('.',mexext),'');

liste  = strvcat(pm{:},pmex{:},'other');
liste = cellstr(liste);
