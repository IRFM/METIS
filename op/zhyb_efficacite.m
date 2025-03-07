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
  cst    = 0.75*sqrt(2.1+3)/4*7;

  nbar = mean(data.gene.nbar);
  Zeffv = mean(data.gene.zeffm);
  phiLH = mean(angle(data.cons.hyb(:,1)))*360/230;
  B     = mean(data.geo.b0)
  etaLHv = cst/sqrt(nbar/1e19+3)*B./(Zeffv+5)*(1-phiLH);% scaling etaLH
  PLH = 3e6;
  R0  = 2.4;

  ILHv     = etaLHv*1e19.*PLH./nbar./R0/1e6;
