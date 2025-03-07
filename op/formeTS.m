function [Rparoi,Zparoi,Rext,Zext,t] = formeTS(numchoc)
% [Rparoi,Zparoi,Rext,Zext] = formeTS(numchoc);
%   --  input  --
% numchoc : numero de choc
%  --  outputs  --
% donne la forme de la paroi de TS (Rparoi,Zparoi)
% + la derniere surface magnetique [Rext(temps,:), Zext(temps,:)]
%
x = tsmat(numchoc,'APOLO;+Parametres;Paroi');

Rparoi = x(:,1);
Zparoi = x(:,2);

theta  = (0:15:23*15)*pi/180;
R0     = 2.42;
[y,t,void1,void2]  = tsbase(numchoc,'grho');

Rext   = R0+y.*cos(ones(size(y,1),1)*theta);
Zext   = y.*sin(ones(size(y,1),1)*theta);

Rext(:,end+1) = Rext(:,1);
Zext(:,end+1) = Zext(:,1);

