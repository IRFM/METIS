choc = input('# choc ? ');
%ordre  = input('nb mode svd ? ') + 1;
reflex = cgcgettrait(choc,'treflex','ref');
times = reflex.times;
pgeo   = cgcgettrait(choc,'gplasma@');
dens   = cgcgettrait(choc,'gnl@');
xest   = pgeo.plasma(:,6); 
nbar   = dens.nl(:,4) ./ pgeo.plasma(:,3)./  pgeo.plasma(:,12) ./ 2;

%[u,s,v] = svd(reflex.rx,0);
%s = diag(s);
%disp(s(1:10)');
%s(ordre:end) = 0;
%reflex.rxc = u * diag(s) * v';
%[u,s,v] = svd(reflex.nex,0);
%s = diag(s);
%disp(s(1:10)');
%s(ordre:end) = 0;
%reflex.nexc = u * diag(s) * v';

%reflex.dndr = pdederive(reflex.rx,reflex.nex,2,2,2,1);
%reflex.dndrc = pdederive(reflex.rxc,reflex.nexc,2,2,2,1);
%reflex.d2ndr2 = pdederive(reflex.rx,reflex.nex,2,2,2,2);
%reflex.d2ndr2s = cat(2,zeros(size(reflex.d2ndr2,1),1),reflex.d2ndr2(:,1:end-1)) .* ...
% 					  reflex.d2ndr2; 
rsepa = NaN .* ones(length(reflex.times),1);
indsepa = NaN .* ones(length(reflex.times),1);
nesepa = NaN .* ones(length(reflex.times),1);
dnesepadr = NaN .* ones(length(reflex.times),1);
nexest = NaN .* ones(length(reflex.times),1);
nxexp  = NaN .* ones(length(reflex.times),1);
lxexp = NaN .* ones(length(reflex.times),1);
%nbar  = NaN .* ones(length(reflex.times),1);

rx     = NaN .* ones(size(reflex.nex));
nex     = NaN .* ones(size(reflex.nex));

for k = 1: length(times)
   [rxc,inds]   =  sort(reflex.rx(k,:));
	nexc         = reflex.nex(k,inds); 
	rxc          = rxc(end:-1:1);
	nexc         = nexc(end:-1:1);
	rx(k,:)      = rxc;
	nex(k,:)     = nexc;
	
	%nbar(k) = abs(trapz(rxc,nexc,2)) ./ (max(rxc) - min(rxc));
	
	%d2ndr2(k,:)  = pdederive(rxc,nexc,2,2,2,2);
	%d2ndr2s(k,:) = cat(2,0,d2ndr2(k,1:end-1)) .* d2ndr2(k,:); 
   %d2ndr2s(k,1:9) = 0;
   %indsepa(k)  = min(find(d2ndr2s(k,:) == min(d2ndr2s(k,:))) + 1,length(d2ndr2s(k,:)));
	
	% calcul du fit exponetiel
	nbp =length(rxc) - 10;
	err = 1e300 .* ones(1,nbp);
	n0  = 1e300 .* ones(1,nbp);
	ln0  = 1e300 .* ones(1,nbp);
	for l = 1:nbp
	   indp   = 5:(l +10);
		indok  = find(nexc(indp) > 3e17);
		if rxc(l+5) < (xest(k) -0.05)
			break
		end
		if length(indok) >= 5
			pp      = polyfit(rxc(indp(indok)),log(nexc(indp(indok))),1);
			indok2  = find(nexc > 3e17);
			nfit1    = exp(polyval(pp,rxc(indp(indok))));
			nexp1    = nexc(indp(indok));
			nfit2    = exp(polyval(pp,rxc(indok2)));
			nexp2    = nexc(indok2);
			err(l)  = sum((nexp1 - nfit1 ) .^ 2) ./ length(nfit) +  ...
			          sum((nfit2 < nexp2) .* nexp2 .^ 2);
         n0(l)      = pp(2);
         ln0(l)     = pp(1);
			
		end
	end
	%figure(4);
	%plot(err(err <1e299));
	%drawnow;
	
	fprintf('.');
   iedeb       = min(min(find(rxc < xest(k))),length(err));
   vemax       = max(err(iedeb));
	indnok      = max(find(err == vemax));
	err(indnok:end) = 1e300;
	indsepa(k)  = min(max(find(err == min(err))) + 10,length(rxc));
	rsepa(k)    = rx(k,indsepa(k));
	nesepa(k)   = nex(k,indsepa(k));
	indfit      = indsepa(k)+ (0:4);
	if max(indfit) <= length(rxc)
		pp = polyfit(rxc(indfit),nexc(indfit),1);
		dnesepadr(k) = pp(1);
	end     
	lxexp(k)    = ln0(max(find(err == min(err))));
	nxexp(k)    = n0(max(find(err == min(err))));
	
	% densite sur xest
	d           = abs(rxc - xest(k));
	nexest(k)   = min(nexc(find(d == min(d))));
	
end
fprintf('\n');
%rsepa = medfilt1(rsepa,3);
%nesepa = medfilt1(nesepa,3);
%nexest = medfilt1(nexest,3);

pp          = [ -0.2867   27.7006 -622.8086];
nebord      = max(3e17,exp(polyval(pp,log(nbar))));

figure(1)
clf
zplotprof(gca,reflex.times,rx,nex,'color','k', ...
          'linestyle','none','marker','.');
nxe  = exp(nxexp * ones(1,size(rx,2)))  .* exp ((lxexp * ones(1,size(rx,2))) .*rx);
nxe(nxe > 1e20) = NaN;
zplotprof(gca,reflex.times,rx,nxe,'color','r');

xx = [xest,xest];
yy = [zeros(size(nbar)),nbar];
zplotprof(gca,times,xx,yy,'color','b');
 

figure(2);
subplot(3,1,1)
plot(times,xest,'b',times,rsepa,'r');
subplot(3,1,2)
semilogy(times,nbar,'k',times,nesepa,'r',times,nexest,'b');
subplot(3,1,3)
semilogy(times,-dnesepadr,'k.');


pp = polyfit(nbar .^ 2,nesepa,1);
figure(3)
plot(nbar,nesepa,'r.',nbar,nexest,'b.',nbar,polyval(pp,nbar .^ 2),'m', ...
     nbar,nebord,'g');
