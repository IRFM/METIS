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
% version  1.9  du  11/03/2002  
%  
%*#auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
figure(1)
clear ind
ind(1)    = min(find(data.gene.temps > jetmse.tqEFTMx(1)));
ind(2)    = min(find(data.gene.temps > jetmse.tqEFTMx(2)));
ind(3)    = min(find(data.gene.temps > jetmse.tqEFTMx(3)));
plot(jetmse.xefit,-jetmse.JEFTMx(1,:),'r',x,-jetmse.JEFTMx(2,:),'b--',x,...
     -jetmse.JEFTMx(3,:),'k:',x,jeux1.data.equi.jmoy(ind(1),:),'ro',...
     x,jeux1.data.equi.jmoy(ind(2),:),'b*',x,jeux1.data.equi.jmoy(ind(3),:),'k+')
hold on

plot(jetmse2.xefit,-jetmse2.JEFTMx(1,:),'m',x,-jetmse2.JEFTMx(2,:),'g--',x,...
     -jetmse2.JEFTMx(3,:),'c:',x,data.equi.jmoy(ind(1),:),'mo',...
     x,data.equi.jmoy(ind(2),:),'g*',x,data.equi.jmoy(ind(3),:),'c+')
