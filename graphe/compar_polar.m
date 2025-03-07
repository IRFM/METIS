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
load figure_comparaison_polar

subplot(2,1,1)

plot(rhoprof,qprof,rho,qcron)
ylabel('q')
legend('Tprof','cronos')
title('choc 30555, t=25 s')
subplot(2,1,2)

plot(1:5,post.polar.af(596,:),1:5,afmes(28709,:),'ro',1:5,af1)

legend('cronos','mesure','tprof')
ylabel('af')
