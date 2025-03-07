%  NOM DE LA FONCTION  courte description  
%------------------------------------------------------------------------------- 
% fichier : nom_du_fichier -> nom_fonction_principale, nom_fonction_locale 1,... 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function .... 
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
errorbar(rho414([1:8]),mean(x414(a414,[1:8])),mean(x414(a414,[1:8]))-min(x414(a414,[1:8])),max(x414(a414,[1:8]))-mean(x414(a414,[1:8])),'rs-')
hold on
errorbar(rho67([1 3:8]),mean(x67(a67,[1 3:8])),mean(x67(a67,[1 3:8]))-min(x67(a67,[1 3:8])),max(x67(a67,[1 3:8]))-mean(x67(a67,[1 3:8])),'bo-')
xlabel('rho')
ylabel('Te (keV)')
errorbar(rho414mhd([1:8]),mean(x414(b414,[1:8])),mean(x414(b414,[1:8]))-min(x414(b414,[1:8])),max(x414(b414,[1:8]))-mean(x414(b414,[1:8])),'g+-')
