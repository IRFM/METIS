% ZCREGULAR regularise un profil au centre
%---------------------------------------------------------------------
% fichier zcregular.m ->  zcregular
%
%
% fonction Matlab 5 :
%
% version fonctionnent aussi avec des complexe de zregular
%
% Cette fonction regularise les profil au centre. Elle modifie les 4
% points les plus centraux. Elle utilise un polynome d'odre 3 de derivee
% nulle au centre et passant par les 3 points adjacants aux points modifies.
% Elle est utiliser pour les densite de courants calculer a partir des formules
% analytique ou au centre le numerateur et le denominateur d'un quotient tendent 
% vers 0. 
%  
% syntaxe  :
%  
%     yy=zcregular(x,y);
%    
% entree :
%
%     x  = variable
%     y  = y(x) [vecteur]
%
% sorties :
% 
%     yy =  profil regularise au centre
%     
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.5, du 12/06/2001.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function yy=zcregular(x,y)

n =5;

x1 = x(n);
x2 = x(n+1);
x3 = x(n+2);

y1 = y(n);
y2 = y(n+1);
y3 = y(n+2);


a = -(x2.^2.*y3-x2.^2.*y1+x3.^2.*y1+y2.*x1.^2-y3.*x1.^2-y2.*x3.^2)./ ...
     (x2-x3)./(x2.^2.*x3.^2-x1.^2.*x2.^2+x2.*x1.^3-x2.*x3.*x1.^2+x3.*x1.^3-x3.^2.*x1.^2);
b = (-x1.^3.*y3+x1.^3.*y2-x3.^3.*y2+y3.*x2.^3-y1.*x2.^3+y1.*x3.^3)./ ...
    (x1.^3.*x2.^2-x1.^3.*x3.^2-x1.^2.*x2.^3+x1.^2.*x3.^3-x3.^3.*x2.^2+x3.^2.*x2.^3);
d = (x2.^2.*x1.^3.*y3-x3.^2.*x1.^3.*y2-y3.*x1.^2.*x2.^3+y2.*x1.^2.*x3.^3+ ...
     x3.^2.*y1.*x2.^3-x2.^2.*y1.*x3.^3)./ ...
    (x1.^3.*x2.^2-x1.^3.*x3.^2-x1.^2.*x2.^3+x1.^2.*x3.^3-x3.^3.*x2.^2+x3.^2.*x2.^3);
c = 0;

p = [a,b,c,d];

yy = y;
yy(1:(n-1)) = polyval(p,x(1:(n-1)));
