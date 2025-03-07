% calcul du piquage de la pression
function [hitb,ate,aitb,teshape,xieshape,xieshape_itb] = z0piqpe(x,fs,fsdb,nep,qp,qep,tep,kishape,qdds,kidds,xieorkie)

% decodage 
ki_expo =  imag(kishape);
kishape =  real(kishape);

ve = ones(size(x));
vt = ones(size(nep,1),1);
if kishape == 0
	%fq = fsdb;
	fq = vt * (1 + 0 .* x);
elseif kishape <= 0
	fq = max(1,qp) .^ abs(kishape);
elseif ki_expo < 0
	% grad(T) = Ct
	kishape_default  = 2/3;
	fq = vt * (x - 2/3 .* x .^ 2 + abs(ki_expo) ./ max(eps,x) ./ 30 + kishape .* x .^ 20);
	fq(:,1)  = 2 .* fq(:,2) - fq(:,3);
	%figure(21);clf;plot(x,fq);
else
	fq = vt * (1 + kishape .* x .^ ki_expo);
end

% applatissement DDS
%fq = fq + (max(fq,[],2) * ve) .* (1 + tanh(30 .* (qdds - qp))) ;
if qdds == 0
	% pas de dents de scie
elseif qdds < 0
 	%fq = fq;
	%figure(21);clf;plot(x,fq);drawnow
elseif qdds <= 1
	fq = fq .* (1 + 0.5 .* max(0,kidds - 1) .* (1 + tanh(30 .* (1 - qp)))) ;
elseif qdds <= 1.5
	fq = fq .* (1 + 0.5 .* max(0,kidds - 1) .* (1 + tanh(30 .* (1.5 - qp)))) ;
else
	fq = fq .* (1 + 0.5 .* max(0,kidds - 1) .* (1 + tanh(30 .* (ceil(qdds) - qp)))) ;
end

% prise en compte de ne (la forme donnee est kie ou xie)
if (xieorkie ~= 0) && (kishape ~= 0)
	fq =  fq .* (nep ./ (max(nep,[],2) * ve)) .^ xieorkie;
end

Qe = cumtrapz(x,qep .* (vt*x),2);
%inte = -Qe ./ fq;
inte = -Qe ./ fq./max(eps,vt*x);
inte(:,1) = 0;
te = cumtrapz(x(:,end:-1:1),inte(:,end:-1:1),2);
te = max(eps,real(te(:,end:-1:1)));
% normalisation
temoy  = trapz(x(1:end-1),te(:,1:end-1) .* (vt*x(1:end-1)),2) ./ trapz(x(1:end-1),(vt*x(1:end-1)),2);
ate    = te(:,1) ./temoy - 1;
teshape = te;

%tepl = tep - (tep(:,end-1) *ve);
%figure(61);clf;plot(x,te ./(te(:,1)*ve),'r',x,tepl ./(tepl(:,1) *ve),'b');drawnow

pe = te .* nep;
we    = trapz(x(1:end-1),pe(:,1:end-1) .* (vt*x(1:end-1)),2);

%inteh = -Qe ./ fq ./ fs .* fsdb;
inteh = -Qe ./ fq ./ fs .* fsdb./max(eps,vt*x);
inteh(:,1) = 0;
teh   = cumtrapz(x(:,end:-1:1),inteh(:,end:-1:1),2);
teh   = max(eps,real(teh(:,end:-1:1)));
%teh    = (teh - teh(:,end-1) *ve) .* (fact*ve) + tep(:,end-1) * ve;
%teh(:,end) = tep(:,end);
peh   = teh .* nep;
weh   = trapz(x(1:end-1),peh(:,1:end-1) .* (vt*x(1:end-1)),2);
hitb  = weh ./ max(eps,we);

hitb     = max(1,hitb);
hitb(~isfinite(hitb)) = 1;

aitb  = max(1,teh(:,1) ./ te(:,1));

% forme des coefficient de transport
xieshape     = fq;
xieshape_itb = fq .* fs ./ fsdb;

%ti = 1:length(vt);
%figure(61);clf;plot(ti,hitb,'b',ti,aitb,'r');drawnow
