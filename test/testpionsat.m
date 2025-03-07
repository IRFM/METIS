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
zineb_compile pionmex4.f
load testpion_saturne

   [piglob,spoyfl,spdrf1,spdrf2,spdrf3,spdrfe,spdic,spdec,sw1,sw2,sw3, ...
  swz1,swz2,swz3,sfusd,sfd1,sfd2,swdf,swdzf,flo1,flo2,flo3,flo4,flo5,flo6,...
  sig1,sig2,sig3,sig4,sig5,sig6,va1,va2,va3,va4,va5,va6,fno1,fno2,fno3,fno4,...
 fno5,fno6,nvm,fc,fl,fr,tzna,sd,tzf,nhres,nharmi,is,nres] =  ...
  pionmex4(zvar,time,curr,qax,PICRH,FREQ,Aph,rho,s,Vol,Rout,shift0,Rin,Bout,...
          Bin,Acharg,Amass,Tempe,Tempi,Dense,Densi0,nion1,nion2,nion3,nion4,nion5);
