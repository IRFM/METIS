function [ww,dwdt] = zdwdt0(w,p,tau,t,vp,flag)


% limite  sur tau
tau   = max(2e-6,abs(tau(:)));

% calcul de la limite en dwdt
if all(isfinite(vp(:)))
   p     = min(vp .* 10e6,max(vp ./ 10e6,abs(p)));
end


% integration pour le calcul de w
%[tw,ww] = ode23t('ztf0',t,w(1),[],t,p,tau);
%[tw,ww] = ode23('ztf0',t,w(1),[],t,p,tau);
[tw,ww] = z0ode(t,p,tau,w(1));
% calcul de dwdt
dwdt  = z0dxdt(ww,t);

if nargin > 5
	figure(16);clf;
	subplot(4,1,1)
	plot(t,dwdt,t,p,t,-p)
	subplot(4,1,2)
	plot(t,tau);
	subplot(4,1,3)
	plot(t,w,t,ww);
	set(gca,'ylim',[0,max(ww)]);
	subplot(4,1,4)
	plot(t,dwdt,'b+',t,z0dxdt(w,t),'or')
	drawnow
	%disp('in zdwdt0')
	%keyboard
end
