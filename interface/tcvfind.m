%  TCVFIND  courte description  
%------------------------------------------------------------------------------- 
% fichier :  tcvfind.m  ->  tcvfind 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   [st,err] = tcvfind(data,index,nom) 
%  
% entrees :  
%  data  = 
%  index = 
%  nom   = 
%  
% sorties :  
%  st  = 
%  err = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  18/11/2005  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function [st,err] = tcvfind(data,index,nom)
   
   ind = strmatch(nom,index,'exact');
   if isempty(ind)
      disp(nom)
      index
      err = -999;
      st =[];
   else
      stc = data{ind};
      err = stc.err; 
      st  = stc.data;
   end     
         
