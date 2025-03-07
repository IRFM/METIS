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
inte = varlog.asser_new.data.inte;
eq   = varlog.asser_new.data.eq;
phyb = varlog.asser_new.data.phyb;



parametre = param.cons.asser;

nparam.K = ones(size(varlog.t))*parametre.K;
parametre.gi=0.4;
test= sum( nparam.K .* (parametre.gp .* (eq + parametre.gi .* inte)),2);

figure(6)
clf
k=2;
inc=10;
int=150;
plot(parametre.qpos,parametre.qref,'g*')
hold on
plot(gene.x((indq(k)-inc):(indq(k)+inc)),data.prof.q(int,(indq(k)-inc):(indq(k)+inc)),'ro',gene.x,data.prof.q(int,:))
mq = mean(prof.q(int,(indq(k)-inc):(indq(k)+inc)))
plot(gene.x(indq(k)),mq,'k*')
