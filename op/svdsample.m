function [vout,nb,sd] = svdsample(t,x,vin,tt,xx,e)

addpath /usr/local/matlab5/toolboxdrfc/BSplines

% test entree
if nargin < 6
	e= 0.95;
end

% svd 
[u,s,v]=svd(vin,0);
sd=diag(s);

% seuil en energie
if e <=1
	ed = sqrt(cumsum(sd.^2));
	ed = ed ./ed(end);
	nb = min(find(ed>=e));
else
	nb =e;
end
if isempty(nb)
	nb =1;
end
if (length(e) >5)&(nb <2)
	nb =2;
end

% mise en forme
x  = x(:);
xx = xx(:);
t  = t(:);
tt = tt(:);

% boucle sur les modes
uu = zeros(length(tt),nb);
vv = zeros(length(xx),nb);
warning off
for k =1:nb
	% interpolation temporelle
	uu(:,k)= interp1(t,medfilt1(u(:,k),3),tt,'linear');
	
	% interpolation spatiale
	if x(1) <=0.01
		vp = cat(1,v(end:-1:2,k),v(:,k),0,0);
 		xp = cat(1,x(1)-x(end:-1:2),x,1.05,1.1);
 	else
 		vp = cat(1,v(end:-1:1,k),v(:,k),0,0);
 		xp = cat(1,-x(end:-1:1),x,1.05,1.1);
 	end
 	C = [0,1,0,0;0,1,0,0.05];
 	L = [0;0];
 	nbsp =fix(length(vp)/3-1);
 	[B,tauh,F0,F1,F2,IntF]=BSPLINFIT1D(xp,vp,nbsp,xx,[0 0],C,L);
 	vv(:,k)  = F0(:);
% 	vvp     = interp1(x,v(:,k),xx,'spline');
%  	dvvdxx  = pdederive(xx,vvp,0,2,1,1);
%   	dvvdxx  = dvvdxx .* (dvvdxx < 0);
%  	vvc     = cumtrapz(xx,dvvdxx);
%  	vvp     = interp1(xp,vp,xx,'spline');
%  	vvc     = vvc -vvc(end) +vvp(end);
%  	vv(:,k) = vvc ./ trapz(xx,vvc) .* trapz(xx,vvp);
end
warning on

% reconstitution
vout = uu * diag(sd(1:nb)) * vv';

% nomalisation
fact=  sqrt(sum(sd.^2)) ./ sqrt(sum(sd(1:nb).^2))
vout = vout .* fact;
