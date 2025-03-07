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
zineb_path;
addpath /usr/drfc/cgc/matlab5/tpion/v4.0
load last_pionmex2
zvar(7)=0;
 [piglob,spoyfl,spdrf1,spdrf2,spdrf3,spdrfe,spdic,spdec,sw1,sw2,sw3, ...
  swz1,swz2,swz3,sfusd,sfd1,sfd2,swdf,swdzf]                     =  ...
  pionmex2(zvar,time,curr,qax,PICRH,FREQ,Aph,rho,s,Vol,Rout,shift0,Rin,Bout,...
          Bin,Acharg,Amass,Tempe,Tempi,Dense,Densi0);
