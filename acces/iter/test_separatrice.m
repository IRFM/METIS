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
geom(1) = 1.687;
geom(2) = 0.466;
geom(3) = 0.0;
geom(4) = 0.0;
geom(5) = 2.001;
geom(6) = 0.568;
geom(7) = 22.46;
geom(8) = 67.92;


lim2 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;

if geom(8) > lim2

  geom(8) = lim2-eps;

end

mode = 1;

[x,y]=separatrice(geom,mode);
[xa,ya]=separatrice_matfile(geom,mode);

figure;
subplot(2,2,1)
plot(x,y,'b.',xa,ya,'ro');
subplot(2,2,2)
plot(x,y,'b',xa,ya,'r');
subplot(2,2,3)
plot(x-xa,'.')
subplot(2,2,4)
plot(y - ya,'.');
