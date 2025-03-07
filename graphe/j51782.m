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
load /usr/drfc/cgc/matlab5/tcron/JET/51782/temp51782
load /usr/drfc/cgc/matlab5/tcron/JET/51782/qMSEx51782
load /usr/drfc/cgc/matlab5/tcron/JET/lecture/toto
load /usr/drfc/cgc/matlab5/tcron/JET/lecture/titi
x=param.gene.x;
t=data.gene.temps;

ind = iround(t,tqEFTMx(4));

h    = figure;
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])
subplot(121)
plot(x,data.prof.jmoy(ind,:),'r',x,-JEFTMx(4,:),'ro',x,data.source.idn.j(ind,:),x,data.source.jboot(ind,:));

hold on
plot([0 1],[0 0])

xlabel('sqrt(Phi/pi/B0)');
ylabel('J (A/m^2)');
legend('J_m_o_y','J_M_S_E','J_N_B_I','J_b_o_o_t');
title('current profiles')
subplot(122)
ind2 = iround(trhonth,tqEFTMx(4));
ind3 = iround(tcx,tqEFTMx(4));

plot(x,data.prof.te(ind,:),rhonth(ind2,:),terhonth(ind2,:),'ro',x,data.prof.ti(ind,:),rhoncx(ind3,:),ticx(ind3,:),'bo')
title('temperature profiles')
xlabel('sqrt(Phi/pi/B0)');
ylabel('eV');
legend('Te ','Te Lidar ','Ti ','Ti CX');

