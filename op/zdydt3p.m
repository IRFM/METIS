% derivee temporelle sur trois points  
% ynew = y a t+dt
% yold = y a dt
% dydtold  = dydt a t
% dt = pas de temps
function dydt = zdydt3p(ynew,yold,dydtold,dt)

if any(~isfinite(dydtold))
	dydt = (ynew - yold) ./ dt;
else
	dydt = 0.5 .* (ynew - yold) ./ dt + 0.5 .* dydtold;
end	
