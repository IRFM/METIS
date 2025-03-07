% reference : P.G. Carolan, Plasma Physics, vol 25 number 10,pp 1065-1086 (1983)
function rz = z0noncronal(tz,lz,xli,tep,nep,taue,dsol,a)

% temps de residence
if nargin > 6
	taul = max((taue .* dsol ./ a) * ones(size(xli)), 1e14 ./ max(1e13,1e6 .* nep)); % s
else
	taul = max(taue * (1.05 - xli), 1e14 ./ max(1e13,1e6 .* nep)); % s
end
% exposant / parametre
% 
lambda = min(1,1e6 .* nep .* taul ./ 2e18);

% temperature du maximum de lz
tz_bar = mean(tz(lz == max(lz)));
lz_bar = mean(lz(lz == max(lz)));

% calcul de la temperature apparante
%t_star = tep .* (tep <= tz_bar) + tep .* (tz_bar ./ tep) .^ lambda .* (tep > tz_bar);
%t_star = tep .* (tep <= tz_bar) + tz_bar .* (tep ./ tz_bar) .^ lambda .* (tep > tz_bar);
t_star = tep .* (tep <= tz_bar) + tz_bar .* (tep > tz_bar);
% calcul du "cooling rate" modifié
%rz   = reshape(10 .^ pchip(log10(tz),log10(lz),log10(t_star(:))),size(tep));
rz   = reshape(10 .^ pchip(log10(tz),log10(lz),log10(tep(:))),size(tep));
rz   = rz .* (tep <= tz_bar) + (rz .* lambda  + (1 - lambda) .* lz_bar) .* (tep > tz_bar);



%  rref   = reshape(10 .^ pchip(log10(tz),log10(lz),log10(tep(:))),size(tep));
%  figure(21);clf;
%  subplot(2,2,1)
%  semilogy(xli,rz,'r',xli,rref,'.b');
%  subplot(2,2,2)
%  semilogy(xli,t_star,'r',xli,tep,'.b')
%  subplot(2,2,3)
%  loglog(tz,lz,'b',tep',rz','.r');
%  subplot(2,2,4);
%  plot(xli,lambda);
%  drawnow
%  keyboard;
