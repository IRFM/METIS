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
% version  3.1  du  20/02/2006  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
Ip = data.gene.ip/1e6;
P=(data.gene.paddfci+data.gene.paddohm+data.gene.paddhyb)/1e6;
Bt=data.geo.b0;
R = data.geo.r0;
a=data.geo.a;
k=1;
ne = data.gene.nbar/1e19;
Meff=2;

tau1 = 0.023*(Ip.^0.96).*(Bt.^0.03).*(ne.^-0.4).*(R.^1.83).*((R./a).^0.06) ...
          .*(k.^0.64).*(Meff.^0.2).*(P.^-0.73);      


for k=1:length(data.gene.temps)

  tauL(k) =H97(Ip(k),Bt(k),ne(k),P(k),R(k),1,a(k),2);

end
