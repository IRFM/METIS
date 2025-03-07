function [alpha,beta,teout] = z0expote(x,te)

% le fit se fait sur 80% du rayon
vt     = ones(size(te,1),1);
indok  = find(x <= 0.9);
ve     = ones(1,length(indok));
forme  = max(30,te(:,indok) - te(:,max(indok)) * ve) ./max(30,(te(:,1) - te(:,max(indok))) * ve);
xf     = x(indok);
indfit = find((xf>=0.1) & (xf <= 0.8));

err   = Inf .* vt;
alpha = NaN .* vt;
beta  = NaN .* vt;
betav = linspace(1,5,41);
for k =1:length(betav)
	betac = betav(k);
	ux    = vt * (1 - xf .^ betac);
	alphac = max(0.1,min(10,mean(log(max(eps,forme(:,indfit))) ./ log(ux(:,indfit)),2)));
	fx     = ux .^ (alphac * ve);
	% errc   = abs(sum(fx(:,indfit) - forme(:,indfit) ,2));
	errc   = sqrt(sum((fx(:,indfit) - forme(:,indfit)) .^ 2,2));
	indbest = find(errc < err);
	if ~isempty(indbest)
		beta(indbest)  = betac;
		alpha(indbest) = alphac(indbest);
		err(indbest)   = errc(indbest);

	end
	%figure(61);clf;plot(xf,forme,'b',xf,fx,'r');title(sprintf('%g %g %g %g',err,errc,alphac,betac));drawnow;pause
end

teout = ((te(:,1) - te(:,max(indok))) * ones(size(x))) .*  ...
        (1 - (vt *x) .^ (beta * ones(size(x)))) .^ (alpha*ones(size(x))) +  ...
	te(:,max(indok)) * ones(size(x));


