%  CLOSE_0PLOT  courte description  
%------------------------------------------------------------------------------- 
% fichier :  close_0plot.m  ->  close_0plot 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   close_0plot 
%  
% entrees :  
%  
%  
% sorties :  
%  
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
function close_0plot

h = findobj(0,'type','figure');
for kh = 1:length(h)
   tag = get(h(kh),'tag');
   if strmatch('z0plot',tag)
      close(h(kh))
   end
end
