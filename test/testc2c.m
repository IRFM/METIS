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
t=data.gene.temps;
x = param.gene.x;
c2cd1    = rpdederive(x,data.equi.c2c,2,2,2,1);
figure;
subplot(2,1,1);
plot(t,data.equi.c2c(:,95:end),'b',t,jeux1.data.equi.c2c(:,95:end),'r')
subplot(2,1,2)
plot(t,c2cd1(:,end));
