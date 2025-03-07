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
load /usr/drfc/artaud/zineb/data/zineb51782_20010503.mat
load /usr/drfc/cgc/matlab5/tcron/JET/51782/prof51782
load /usr/drfc/cgc/matlab5/tcron/JET/51782/efit51782
load /usr/drfc/cgc/matlab5/tcron/JET/51782/temp51782
R1 = squeeze(data.equi.R(1,:,:));
Z1 = squeeze(data.equi.Z(1,:,:));
t  = data.gene.temps;

ind1 = iround(tefit,t(1));


h    = figure;
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])
plot(R1(1:10:end,:)',Z1(1:10:end,:)',rbnd(ind1,:),zbnd(ind1,:),'ro',rmag(ind1),zmag(ind1),'bo')
axis('equal')
xlabel('R(m)')
ylabel('Z(m)')
