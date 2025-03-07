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
cd /usr/drfc/cgc/matlab5/zineb/v2.2/equi/helena
for k=1:6
      [C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
      r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,kkbig,li,dimerc,drmerc,balcrit,ifail]= ...
      helmex77a(psipin,ptot,jmoy,geo.r0,geo.b0,ip,2,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R(1:k:end),geo.Z(1:k:end),init,b,df2,kkbig);
  erreurvp(k,1) = Vp(end)-Vp(end-1);
  erreurvp(k,2) = Vp(end-1)-Vp(end-2);
end
      [C2,C3,rho2m,R,Z,rho,Vp1,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
      r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,kkbig,li,dimerc,drmerc,balcrit,ifail]= ...
      helmex77(psipin,ptot,jmoy,geo.r0,geo.b0,ip,geo.mode,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R,geo.Z,init,b,df2,kkbig);


      [C2,C3,rho2m,R,Z,rho,Vp2,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp,rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
      r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,kkbig,li,dimerc,drmerc,balcrit,ifail]= ...
      helmex77(psipin,ptot,jmoy,geo.r0,geo.b0,ip,1,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R,geo.Z,init,b,df2,kkbig);
cd /usr/drfc/cgc/matlab5/zineb/vup/equi/helena

for k=1:6

      [C2,C3,rho2m,R,Z,rho,Vp,drhoor,qpsi,psipout,ftra,Fdia,psiT,B2,invB2,Sp, ...
      rav,oor,drhoav,xshift,xrad,xell,xtriapos,xtrianeg,jout,ipout,pout, ...
      r2,oor2,r2tau2,drho2ob2,r3otau3,r3otau,BPR,BPZ,b,df2,kkbig,li,ifail]= ...
      helmex77(psipin,ptot,jmoy,geo.r0,geo.b0,ip, ...
               geo.mode,geo.a,geo.e1,geo.trh1,geo.trb1,geo.R(1:k:end),geo.Z(1:k:end), ...
               init,b,df2,kkbig);
  erreurvp(k,1) = Vp(end)-Vp(end-1);
  erreurvp(k,2) = Vp(end-1)-Vp(end-2);
end
