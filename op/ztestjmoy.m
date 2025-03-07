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
for ind =param.gene.kmin:param.gene.k
datakp1 = zget1t(data,ind);
x = param.gene.x;
[jmoy,jmoync,jmoy1,jmoy2,c2cd1]  = zcaljmoy(param.gene.x,datakp1.prof.psi,datakp1.equi.vpr,datakp1.equi.grho2r2,datakp1.equi.rhomax,datakp1.equi.ri, ...
                 param.phys.mu0,param.gene.creux,datakp1.equi.rhoRZ,datakp1.equi.c2c);
[jmoyold,jmoyncold]  = zcaljmoy(param.gene.x,datakp1.prof.psid1,datakp1.equi.vpr,datakp1.equi.grho2r2,datakp1.equi.rhomax,datakp1.equi.ri, ...
                 param.phys.mu0,param.gene.creux);


figure(16);
clf
subplot(2,1,1)
plot(x,jmoy,'r',x,jmoy1,'m',x,jmoy2,'b',x,datakp1.equi.jmoy,'c',x,datakp1.source.totale.j,'g');
subplot(2,1,2)
plot(x,datakp1.equi.vpr,x,datakp1.equi.c2c,x,c2cd1);
pause
end
