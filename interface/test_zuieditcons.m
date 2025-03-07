%  TEST_ZUIEDITCONS   script de test de zuieditcons
%------------------------------------------------------------------------------- 
% fichier : test_zuieditcons.m -> test_zuieditcons
%				                      zuieditcons 
% 
% fonction Matlab 5 : 
%	script de test de zuieditcons
% 
% syntaxe :  
%  
% entrees :  
%  
% sorties :  
%  
% fonction écrite par J-F Artaud , poste 46-78
% version  1.7  du  29/09/2001  
%  
%  
% liste des modifications :  
%  
%-------------------------------------------------------------------------------  
% 
 
% edition de la consigne ip
nom     = 'Ip';
info    = zinfo;
aide    = info.data.cons.ip;
x       = data.gene.temps;
y       = data.cons.ip;
texte_x = 'temps (s)';
texte_y = 'Ip (MA)'; 
var_x   = 'void';
var_y   = 'data.cons.ip';
canal   = 1;
hout=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal);
zwaitfor(hout)
figure
plot(x,y,'b',data.gene.temps,data.cons.ip,'o');

% edition du module de Hyb(2) 
nom     = 'consigne de puissance du coupleur 2 de l''hybride ';
info    = zinfo;
aide    = info.data.cons.hyb;
x       = data.gene.temps;
y       = data.cons.hyb(:,1);
texte_x = 'temps (s)';
texte_y = 'Puissance (W)'; 
var_x   = 'void';
var_y   = 'data.cons.hyb';
canal   = 2;
hout=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,'abs');
zwaitfor(hout)
figure
plot(x,y,'b',data.gene.temps,abs(data.cons.hyb(:,2)),'o');

% edition du n// de Hyb(2) 
nom     = 'N// du coupleur 2 de l''hyvbride';
info    = zinfo;
aide    = info.data.cons.hyb;
x       = data.gene.temps;
y       = data.cons.hyb(:,1);
texte_x = 'temps (s)';
texte_y = 'N// (su)'; 
var_x   = 'void';
var_y   = 'data.cons.hyb';
canal   = 2;
if rand(1) >0.5
	hout=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,'angle');
else
	hout=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,'degres');
end	
zwaitfor(hout)
figure
plot(x,angle(y),'b',data.gene.temps,angle(data.cons.hyb(:,2)),'o');

% edition du profil de depot electronique de l'hybreide (2)
nom     = 'profil de depot de puissance de l''hybride sur les electrons';
info    = zinfo;
aide    = info.data.source.hyb.el;
x       = param.gene.x;
[u,s,v] = svd(data.prof.te,0);
y       = v(:,1).*s(1,1).*max(u(:,1));
texte_x = 'rho';
texte_y = 'source.hyb.el (W)'; 
var_x   = 'void';
var_y   = 'data.source.hyb.el';
canal   = 1;
liste_ref = 'te1   |te2  |te3  |vide';
var_ref   = {{'x','v(:,1)',':'}, ...
	            {'x','v(:,2)',':'}, ...
	            {'x','v(:,3)',':'}, ...
	            {'[]','[]',''}};
texte_prop = 'hybride';
var_prop = 'data.cons.hyb';	            
hout=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,'',liste_ref,var_ref,texte_prop,var_prop);
zwaitfor(hout)
figure
plot(x,y,'b',param.gene.x,data.source.hyb.el,'o');
