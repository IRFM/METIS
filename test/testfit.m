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

   nbou = 0;
   FAIL = 1;
   while FAIL ~= 0 & nbou < 50

   try
      [NX,NY,LAMDA,MU,C,FP,DEFRANK,FAIL] =  ...
       zfitdata(nRpla,nZpla,Psi,w./max(max(abs(Psi))),s,[-1,1,0,1]);
      s=s*1.05;
      nbou = nbou+1;

   catch
       error('Possibly s too small for zfitdata, s=s*1.05');

   end
   end