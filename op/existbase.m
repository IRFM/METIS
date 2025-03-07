%  EXISTBASE  courte description
%------------------------------------------------------------------------------- 
% fichier :  existbase.m  ->  existbase 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   cr = existbase(nom) 
%  
% entrees :  
%  nom = 
%  
% sorties :  
%   cr  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  2.2  du  13/11/2003 
%  
% #auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * 13/11/2003 -> ajou du support des cell et sous structure
%  
%-------------------------------------------------------------------------------  
%  
function cr = existbase(nom)

if any(nom == '.') | any(nom == '{')
  try
     evalin('base',sprintf('size(%s);',nom));
      cr = 1;
  catch
      cr = 0;
  end
else
   cr = evalin('base',strcat('exist(''',nom,''')'));
end
