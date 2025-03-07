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
addpath /usr/drfc/cgc/matlab5/zineb/v1.9
zineb_path;
zuiload('/usr/drfc/cgc/cgc_data/zineb/data/zineb51644jet0de58880a62260_resultat')
jetmse = load(['/usr/drfc/cgc/cgc_data/jet/data/51644/qMSEx51644.mat']);


x      = param.gene.x;

figure(1)
clear ind
ind(1)    = min(find(data.gene.temps > jetmse.tqEFTMx(1)));
ind(2)    = min(find(data.gene.temps > jetmse.tqEFTMx(2)));
ind(3)    = min(find(data.gene.temps > jetmse.tqEFTMx(3)));
plot(jetmse.xefit,-jetmse.JEFTMx(1,:),'r',x,-jetmse.JEFTMx(2,:),'b--',x,...
     -jetmse.JEFTMx(3,:),'k:',x,data.prof.jmoy(ind(1),:),'ro',...
     x,data.prof.jmoy(ind(2),:),'b*',x,data.prof.jmoy(ind(3),:),'k+')
	  
zuiload('/usr/drfc/cgc/cgc_data/zineb/data/zineb51643jet0de58880a62330_resultat')
jetmse2 = load(['/usr/drfc/cgc/cgc_data/jet/data/51643/qMSEx51643.mat']);
figure(2)
clear ind
ind(1)    = min(find(data.gene.temps > jetmse.tqEFTMx(1)));
ind(2)    = min(find(data.gene.temps > jetmse.tqEFTMx(2)));
ind(3)    = min(find(data.gene.temps > jetmse.tqEFTMx(3)));
plot(jetmse2.xefit,-jetmse2.JEFTMx(1,:),'r',x,-jetmse2.JEFTMx(2,:),'b--',x,...
     -jetmse2.JEFTMx(3,:),'k:',x,data.prof.jmoy(ind(1),:),'ro',...
     x,data.prof.jmoy(ind(2),:),'b*',x,data.prof.jmoy(ind(3),:),'k+')
