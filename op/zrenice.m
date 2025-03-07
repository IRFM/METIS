%  ZRENICE  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zrenice.m  ->  zrenice 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   zrenice 
%  
% entrees :  
%  
%  
% sorties :  
%  
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
function zrenice

    [pid,pid_sess,user]=getidprocess;
    unix(sprintf('renice 15 %d',pid));
