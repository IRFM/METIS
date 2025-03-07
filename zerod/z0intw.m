function [ww,dwdt] = z0intw(t,p,tau,w);
% this function compute an simple solution of equation dw/dt = - w(t) / tau(t) + P(t);
% warning : this method is more precise than z0ode but time must be short due to exp()
% on suppose que tau = tau0 / sqrt(P)

pn  = max(1,max(p));
tau = max(1e-6,tau);
p   = max(0,p ./ pn);
tauf = tau .* sqrt(p);
inttauinv   = sqrt(p) .* cumtrapz(t,1./tauf);
intmtauinv  =  cumtrapz(t,- 1./tau);
ww          =  (pn .*cumtrapz(t,p .* exp(inttauinv)) + w(1)) .* exp(intmtauinv);
dwdt        =  z0dxdt(ww,t);

figure(134);clf
plot(t,ww,t,dwdt)
drawnow
