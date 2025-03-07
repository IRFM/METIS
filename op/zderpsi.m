% cette fonction calcule la derivee de pis dans le repere de l'equilibre afin de prevenir des effets de bord
function [psid1,psid2,c2cd1]=zderpsi(x,psi,rhorz,rhomax,c2c)

mode = 'spline';

xeq = double(rhorz);

psieq = interp1(x.*rhomax,psi,xeq,mode);

psid1     = rpdederive(xeq,psieq,0,1,2,1);
psid2     = rpdederive(xeq,psid1,2,0,2,1);

psid1     = interp1(xeq,psid1,x.*rhomax,mode).*rhomax;
psid2     = interp1(xeq,psid2,x.*rhomax,mode).*rhomax.^ 2;


if nargin > 4
   c2ceq = interp1(x.*rhomax,c2c,xeq,mode);
   c2cd1 = rpdederive(xeq,c2ceq,2,1,2,1);
   c2cd1 = interp1(xeq,c2cd1,x.*rhomax,mode).*rhomax;
   c2cd1(1:3)  = c2cd1(4);
end
