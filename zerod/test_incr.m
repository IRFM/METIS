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
xdep = linspace(0,1);
incr = [0.1,0.5,1,2,3,5,10]';
vx=ones(size(xdep));
vi =ones(size(incr));
xdep = vi*xdep;
incr = incr*vx;
hitb =  (1 + incr .* max(0,min(1,abs(xdep))) .^ 3);
aitb = (incr.*xdep+1) ./hitb;
figure(19);plot(xdep',hitb',xdep',aitb')
