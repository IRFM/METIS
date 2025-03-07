function zplotbench(short)


if nargin <1
	short = 0;
end


try
     load /usr/drfc/cgc/cgc_data/zineb/bench/benchdata
catch
     [out,ref,jdiff,pe,pepion] = zloadbench;
end
disp('------------------')
disp(' 0 = out')
disp(' 1 = jdiff')
disp(' 2 = pe')
disp(' 3 = pepion')
mode = input('mode ?')
switch mode
	
case 1
	out = jdiff;
	ref = jdiff{1};
case 2
	out = pe;
	ref = pe{1};
case 3
	out = pepion;
	ref = pepion{1};
end


% intervalle util
ind = 10:80;
nbt = length(ind);
%ref             = out{kref};

% information statistiques :
for k =1:length(out)
	outc              = out{k};
	if ~isempty(outc)
		cpu(k)            = (max(outc.cpu(ind)) - min(outc.cpu(ind)))/nbt;
		memory(k)         = max(outc.memory(ind));
		datation(k)       = (max(outc.date(ind)) - min(outc.date(ind)))/nbt;
		cn(k)             = outc.cn;
		amorti(k)         = outc.amorti;
		conv.mean(k)      = mean(outc.conv(ind));
		conv.std(k)       = std(outc.conv(ind));
		conv.median(k)    = median(outc.conv(ind));
		nbsplit.mean(k)   = mean(outc.nbsplit(ind));
		nbsplit.std(k)    = std(outc.nbsplit(ind));
		nbsplit.median(k) = median(outc.nbsplit(ind));
		
		ind2 = find((outc.tmin >= min(outc.t(ind))) & (outc.tmax <= max(outc.t(ind))));
		
		dt.mean(k)        = mean(outc.dt(ind2)*1e3);
		dt.std(k)         = std(outc.dt(ind2)*1e3);
		dt.median(k)      = median(outc.dt(ind2)*1e3);
		
		derr               = (2 .* (outc.li(ind) - outc.liexp(ind))./(outc.li(ind) + outc.liexp(ind))) .^ 2;
		errli.mean(k)     = round(sqrt(mean(derr))*100);
		errli.std(k)      = round(sqrt(std(derr))*100);
		errli.med(k)      = round(sqrt(median(derr))*100);
		
		derr               = (2 .* (outc.ip(ind) - outc.ipexp(ind))./(outc.ip(ind) + outc.ipexp(ind))) .^ 2;
		errip.mean(k)     = round(sqrt(mean(derr))*100);
		errip.std(k)      = round(sqrt(std(derr))*100);
		errip.med(k)      = round(sqrt(median(derr))*100);
		
		derr               = (2 .* (outc.vres(ind) - outc.vexp(ind))./(outc.vres(ind) + outc.vexp(ind))) .^ 2;
		errvs.mean(k)     = round(sqrt(mean(derr))*100);
		errvs.std(k)      = round(sqrt(std(derr))*100);
		errvs.med(k)      = round(sqrt(median(derr))*100);
		
		derr               = (2 .* (outc.qa(ind) - outc.qaexp(ind))./(outc.qa(ind) + outc.qaexp(ind))) .^ 2;
		errqa.mean(k)     = round(sqrt(mean(derr))*100);
		errqa.std(k)      = round(sqrt(std(derr))*100);
		errqa.med(k)      = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.ne(ind,1) - ref.ne(ind,1))./ (outc.ne(ind,1) + ref.ne(ind,1))) .^ 2;
		ne0.mean(k)       = round(sqrt(mean(derr))*100);
		ne0.std(k)        = round(sqrt(std(derr))*100);
		ne0.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.ne(ind,51) - ref.ne(ind,51))./ (outc.ne(ind,51) + ref.ne(ind,51))) .^ 2;
		nem.mean(k)       = round(sqrt(mean(derr))*100);
		nem.std(k)        = round(sqrt(std(derr))*100);
		nem.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.pe(ind,1) - ref.pe(ind,1))./ (outc.pe(ind,1) + ref.pe(ind,1))) .^ 2;
		pe0.mean(k)       = round(sqrt(mean(derr))*100);
		pe0.std(k)        = round(sqrt(std(derr))*100);
		pe0.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.pe(ind,51) - ref.pe(ind,51))./ (outc.pe(ind,51) + ref.pe(ind,51))) .^ 2;
		pem.mean(k)       = round(sqrt(mean(derr))*100);
		pem.std(k)        = round(sqrt(std(derr))*100);
		pem.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.pion(ind,1) - ref.pion(ind,1))./ (outc.pion(ind,1) + ref.pion(ind,1))) .^ 2;
		pi0.mean(k)       = round(sqrt(mean(derr))*100);
		pi0.std(k)        = round(sqrt(std(derr))*100);
		pi0.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.pion(ind,51) - ref.pion(ind,51))./ (outc.pion(ind,51) + ref.pion(ind,51))) .^ 2;
		pim.mean(k)       = round(sqrt(mean(derr))*100);
		pim.std(k)        = round(sqrt(std(derr))*100);
		pim.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.qdiff(ind,1) - ref.qdiff(ind,1))./ (outc.qdiff(ind,1) + ref.qdiff(ind,1))) .^ 2;
		qdiff0.mean(k)       = round(sqrt(mean(derr))*100);
		qdiff0.std(k)        = round(sqrt(std(derr))*100);
		qdiff0.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.qdiff(ind,51) - ref.qdiff(ind,51))./ (outc.qdiff(ind,51) + ref.qdiff(ind,51))) .^ 2;
		qdiffm.mean(k)       = round(sqrt(mean(derr))*100);
		qdiffm.std(k)        = round(sqrt(std(derr))*100);
		qdiffm.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.qeq(ind,1) - ref.qeq(ind,1))./ (outc.qeq(ind,1) + ref.qeq(ind,1))) .^ 2;
		qeq0.mean(k)       = round(sqrt(mean(derr))*100);
		qeq0.std(k)        = round(sqrt(std(derr))*100);
		qeq0.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.qeq(ind,51) - ref.qeq(ind,51))./ (outc.qeq(ind,51) + ref.qeq(ind,51))) .^ 2;
		qeqm.mean(k)       = round(sqrt(mean(derr))*100);
		qeqm.std(k)        = round(sqrt(std(derr))*100);
		qeqm.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.jmoy(ind,1) - ref.jmoy(ind,1))./ (outc.jmoy(ind,1) + ref.jmoy(ind,1))) .^ 2;
		jmoy0.mean(k)       = round(sqrt(mean(derr))*100);
		jmoy0.std(k)        = round(sqrt(std(derr))*100);
		jmoy0.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.jmoy(ind,51) - ref.jmoy(ind,51))./ (outc.jmoy(ind,51) + ref.jmoy(ind,51))) .^ 2;
		jmoym.mean(k)       = round(sqrt(mean(derr))*100);
		jmoym.std(k)        = round(sqrt(std(derr))*100);
		jmoym.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.jmoyeq(ind,1) - ref.jmoyeq(ind,1))./ (outc.jmoyeq(ind,1) + ref.jmoyeq(ind,1))) .^ 2;
		jmoyeq0.mean(k)       = round(sqrt(mean(derr))*100);
		jmoyeq0.std(k)        = round(sqrt(std(derr))*100);
		jmoyeq0.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.jmoyeq(ind,51) - ref.jmoyeq(ind,51))./ (outc.jmoyeq(ind,51) + ref.jmoyeq(ind,51))) .^ 2;
		jmoyeqm.mean(k)       = round(sqrt(mean(derr))*100);
		jmoyeqm.std(k)        = round(sqrt(std(derr))*100);
		jmoyeqm.med(k)        = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.li(ind,1) - ref.li(ind,1))./ (outc.li(ind,1) + ref.li(ind,1))) .^ 2;
		li.mean(k)        = round(sqrt(mean(derr))*100);
		li.std(k)         = round(sqrt(std(derr))*100);
		li.med(k)         = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.ip(ind,1) - ref.ip(ind,1))./ (outc.ip(ind,1) + ref.ip(ind,1))) .^ 2;
		ip.mean(k)        = round(sqrt(mean(derr))*100);
		ip.std(k)         = round(sqrt(std(derr))*100);
		ip.med(k)         = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.vres(ind,1) - ref.vres(ind,1))./ (outc.vres(ind,1) + ref.vres(ind,1))) .^ 2;
		vs.mean(k)        = round(sqrt(mean(derr))*100);
		vs.std(k)         = round(sqrt(std(derr))*100);
		vs.med(k)         = round(sqrt(median(derr))*100);
		
		derr              = (2 .* (outc.qa(ind,1) - ref.qa(ind,1))./ (outc.qa(ind,1) + ref.qa(ind,1))) .^ 2;
		qa.mean(k)        = round(sqrt(mean(derr))*100);
		qa.std(k)         = round(sqrt(std(derr))*100);
		qa.med(k)         = round(sqrt(median(derr))*100);
		
		amin              = min(outc.amorti,ref.amorti);
		amax              = max(outc.amorti,ref.amorti);
		
	end
end	

% affichage du tableau
switch short
case 1
	
	fprintf('--------------------------------------------------');
	fprintf('\ncn             =   ');fprintf('%8.3g  ',cn);
	fprintf('\namorti         =   ');fprintf('%8.3g  ',amorti);
	fprintf('\ncpu            =   ');fprintf('%8.3g  ',cpu);
	fprintf('\n<conv>         =   ');fprintf('%8.3g  ',conv.mean);
	fprintf('\n<nbsplit>      =   ');fprintf('%8.3g  ',nbsplit.mean);
	fprintf('\n<dt>      (ms) =   ');fprintf('%8.3g  ',dt.mean);
	fprintf('\n<D(li)>        =   ');fprintf('%8.3g  ',errli.mean);
	fprintf('\n<D(ip)>        =   ');fprintf('%8.3g  ',errip.mean);
	fprintf('\n<D(Vs)>        =   ');fprintf('%8.3g  ',errvs.mean);
	fprintf('\n<D(qa)>        =   ');fprintf('%8.3g  ',errqa.mean);
	fprintf('\n<ne0>          =   ');fprintf('%8.3g  ',ne0.mean);
	fprintf('\n<nem>          =   ');fprintf('%8.3g  ',nem.mean);
	fprintf('\n<pe0>          =   ');fprintf('%8.3g  ',pe0.mean);
	fprintf('\n<pem>          =   ');fprintf('%8.3g  ',pem.mean);
	fprintf('\n<pion0>        =   ');fprintf('%8.3g  ',pi0.mean);
	fprintf('\n<pionm>        =   ');fprintf('%8.3g  ',pim.mean);
	fprintf('\n<jmoy0>        =   ');fprintf('%8.3g  ',jmoy0.mean);
	fprintf('\n<jmoym>        =   ');fprintf('%8.3g  ',jmoym.mean);
	fprintf('\n<jmoyeq0>      =   ');fprintf('%8.3g  ',jmoyeq0.mean);
	fprintf('\n<jmoyeqm>      =   ');fprintf('%8.3g  ',jmoyeqm.mean);
	fprintf('\n<qdiff0>       =   ');fprintf('%8.3g  ',qdiff0.mean);
	fprintf('\n<qdiffm>       =   ');fprintf('%8.3g  ',qdiffm.mean);
	fprintf('\n<qeq0>         =   ');fprintf('%8.3g  ',qeq0.mean);
	fprintf('\n<qeqm>         =   ');fprintf('%8.3g  ',qeqm.mean);
	fprintf('\n<li>           =   ');fprintf('%8.3g  ',li.mean);
	fprintf('\n<ip>           =   ');fprintf('%8.3g  ',ip.mean);
	fprintf('\n<Vs>           =   ');fprintf('%8.3g  ',vs.mean);
	fprintf('\n<qa>           =   ');fprintf('%8.3g  ',qa.mean);
	fprintf('\n--------------------------------------------------\n');
		
otherwise
	
	fprintf('--------------------------------------------------');
	fprintf('\ncn             =   ');fprintf('%8.3g  ',cn);
	fprintf('\namorti         =   ');fprintf('%8.3g  ',amorti);
	fprintf('\ncpu            =   ');fprintf('%8.3g  ',cpu);
	fprintf('\nmemory         =   ');fprintf('%8.3g  ',memory);
	fprintf('\ndatation       =   ');fprintf('%8.3g  ',datation);
	fprintf('\n<conv>         =   ');fprintf('%8.3g  ',conv.mean);
	fprintf('\nmed(conv)      =   ');fprintf('%8.3g  ',conv.median);
	fprintf('\nstd(conv)      =   ');fprintf('%8.3g  ',conv.std);
	fprintf('\n<nbsplit>      =   ');fprintf('%8.3g  ',nbsplit.mean);
	fprintf('\nmed(nbsplit)   =   ');fprintf('%8.3g  ',nbsplit.median);
	fprintf('\nstd(nbsplit)   =   ');fprintf('%8.3g  ',nbsplit.std);
	fprintf('\n<dt>      (ms) =   ');fprintf('%8.3g  ',dt.mean);
	fprintf('\nmed(dt)   (ms) =   ');fprintf('%8.3g  ',dt.median);
	fprintf('\nstd(dt)   (ms) =   ');fprintf('%8.3g  ',dt.std);
	fprintf('\n<D(li)>        =   ');fprintf('%8.3g  ',errli.mean);
	fprintf('\nmed(D(li))     =   ');fprintf('%8.3g  ',errli.med);
	fprintf('\nstd(D(li))     =   ');fprintf('%8.3g  ',errli.std);
	fprintf('\n<D(ip)>        =   ');fprintf('%8.3g  ',errip.mean);
	fprintf('\nmed(D(ip))     =   ');fprintf('%8.3g  ',errip.med);
	fprintf('\nstd(D(ip))     =   ');fprintf('%8.3g  ',errip.std);
	fprintf('\n<D(Vs)>        =   ');fprintf('%8.3g  ',errvs.mean);
	fprintf('\nmed(D(Vs))     =   ');fprintf('%8.3g  ',errvs.med);
	fprintf('\nstd(D(Vs))     =   ');fprintf('%8.3g  ',errvs.std);
	fprintf('\n<D(qa)>        =   ');fprintf('%8.3g  ',errqa.mean);
	fprintf('\nmed(D(qa))     =   ');fprintf('%8.3g  ',errqa.med);
	fprintf('\nstd(D(qa))     =   ');fprintf('%8.3g  ',errqa.std);
	fprintf('\n<ne0>          =   ');fprintf('%8.3g  ',ne0.mean);
	fprintf('\nmed(ne0)       =   ');fprintf('%8.3g  ',ne0.med);
	fprintf('\nstd(ne0)       =   ');fprintf('%8.3g  ',ne0.std);
	fprintf('\n<nem>          =   ');fprintf('%8.3g  ',nem.mean);
	fprintf('\nmed(nem)       =   ');fprintf('%8.3g  ',nem.med);
	fprintf('\nstd(nem)       =   ');fprintf('%8.3g  ',nem.std);
	fprintf('\n<pe0>          =   ');fprintf('%8.3g  ',pe0.mean);
	fprintf('\nmed(pe0)       =   ');fprintf('%8.3g  ',pe0.med);
	fprintf('\nstd(pe0)       =   ');fprintf('%8.3g  ',pe0.std);
	fprintf('\n<pem>          =   ');fprintf('%8.3g  ',pem.mean);
	fprintf('\nmed(pem)       =   ');fprintf('%8.3g  ',pem.med);
	fprintf('\nstd(pem)       =   ');fprintf('%8.3g  ',pem.std);
	fprintf('\n<pion0>        =   ');fprintf('%8.3g  ',pi0.mean);
	fprintf('\nmed(pion0)     =   ');fprintf('%8.3g  ',pi0.med);
	fprintf('\nstd(pion0)     =   ');fprintf('%8.3g  ',pi0.std);
	fprintf('\n<pionm>        =   ');fprintf('%8.3g  ',pim.mean);
	fprintf('\nmed(pionm)     =   ');fprintf('%8.3g  ',pim.med);
	fprintf('\nstd(pionm)     =   ');fprintf('%8.3g  ',pim.std);
	fprintf('\n<jmoy0>        =   ');fprintf('%8.3g  ',jmoy0.mean);
	fprintf('\nmed(jmoy0)     =   ');fprintf('%8.3g  ',jmoy0.med);
	fprintf('\nstd(jmoy0)     =   ');fprintf('%8.3g  ',jmoy0.std);
	fprintf('\n<jmoym>        =   ');fprintf('%8.3g  ',jmoym.mean);
	fprintf('\nmed(jmoym)     =   ');fprintf('%8.3g  ',jmoym.med);
	fprintf('\nstd(jmoym)     =   ');fprintf('%8.3g  ',jmoym.std);
	fprintf('\n<jmoyeq0>      =   ');fprintf('%8.3g  ',jmoyeq0.mean);
	fprintf('\nmed(jmoyeq0)   =   ');fprintf('%8.3g  ',jmoyeq0.med);
	fprintf('\nstd(jmoyeq0)   =   ');fprintf('%8.3g  ',jmoyeq0.std);
	fprintf('\n<jmoyeqm>      =   ');fprintf('%8.3g  ',jmoyeqm.mean);
	fprintf('\nmed(jmoyeqm)   =   ');fprintf('%8.3g  ',jmoyeqm.med);
	fprintf('\nstd(jmoyeqm)   =   ');fprintf('%8.3g  ',jmoyeqm.std);
	fprintf('\n<qdiff0>       =   ');fprintf('%8.3g  ',qdiff0.mean);
	fprintf('\nmed(qdiff0)    =   ');fprintf('%8.3g  ',qdiff0.med);
	fprintf('\nstd(qdiff0)    =   ');fprintf('%8.3g  ',qdiff0.std);
	fprintf('\n<qdiffm>       =   ');fprintf('%8.3g  ',qdiffm.mean);
	fprintf('\nmed(qdiffm)    =   ');fprintf('%8.3g  ',qdiffm.med);
	fprintf('\nstd(qdiffm)    =   ');fprintf('%8.3g  ',qdiffm.std);
	fprintf('\n<qeq0>         =   ');fprintf('%8.3g  ',qeq0.mean);
	fprintf('\nmed(qeq0)      =   ');fprintf('%8.3g  ',qeq0.med);
	fprintf('\nstd(qeq0)      =   ');fprintf('%8.3g  ',qeq0.std);
	fprintf('\n<qeqm>         =   ');fprintf('%8.3g  ',qeqm.mean);
	fprintf('\nmed(qeqm)      =   ');fprintf('%8.3g  ',qeqm.med);
	fprintf('\nstd(qeqm)      =   ');fprintf('%8.3g  ',qeqm.std);
	fprintf('\n<li>           =   ');fprintf('%8.3g  ',li.mean);
	fprintf('\nmed(li)        =   ');fprintf('%8.3g  ',li.med);
	fprintf('\nstd(li)        =   ');fprintf('%8.3g  ',li.std);
	fprintf('\n<ip>           =   ');fprintf('%8.3g  ',ip.mean);
	fprintf('\nmed(ip)        =   ');fprintf('%8.3g  ',ip.med);
	fprintf('\nstd(ip)        =   ');fprintf('%8.3g  ',ip.std);
	fprintf('\n<Vs>           =   ');fprintf('%8.3g  ',vs.mean);
	fprintf('\nmed(Vs         =   ');fprintf('%8.3g  ',vs.med);
	fprintf('\nstd(V)         =   ');fprintf('%8.3g  ',vs.std);
	fprintf('\n<qa>           =   ');fprintf('%8.3g  ',qa.mean);
	fprintf('\nmed(qa)        =   ');fprintf('%8.3g  ',qa.med);
	fprintf('\nstd(qa)        =   ');fprintf('%8.3g  ',qa.std);
	fprintf('\n--------------------------------------------------\n');
end	



h=findobj(0,'type','figure','tag','zplotbench1');
if isempty(h)
	h=figure('tag','zplotbench1');
else
	figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	   'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',3,'units','normalized','position',[0.03 0.03 0.95 0.95]);
	   
subplot(2,1,1)
cc= get(gca,'colororder');
hold on
legstr={};
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		plot(outc.t(ind),outc.ne(ind,1),'color',cc(kc,:));
		plot(outc.t(ind),outc.ne(ind,51),'linestyle','-.','color',cc(kc,:));
		legstr{end+1} = sprintf('r = 0  , cn = %g, amorti = %g',outc.cn,outc.amorti);
		legstr{end+1} = sprintf('r = 0.5, cn = %g, amorti = %g',outc.cn,outc.amorti);
	end
	
end
xlabel('temps (s)');
ylabel('Ne (m^-^3)');
legend(legstr,-1);
title('Bench');

subplot(2,1,2)
hold on
nb =0;
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		tp = outc.t(ind) -outc.t(ind(1));
		tt =  cat(1,-tp(end:-1:2),tp);
		xne0 = xcov(outc.ne(ind,1),ref.ne(ind,1))./xcov(ref.ne(ind,1),ref.ne(ind,1));
		xnem = xcov(outc.ne(ind,51),ref.ne(ind,51))./xcov(ref.ne(ind,51),ref.ne(ind,51));
		plot(tt,xne0,'color',cc(kc,:));
		nb = nb +1;
	end
	
end
xlabel('Delta(temps) (s)');
ylabel('XCov(Ne-Neref,Neref)');
legend(legstr{1:nb},-1);
drawnow
pos =get(gca,'pos');
legend off
title('Bench');
set(gca,'pos',pos);
	   

	   
h=findobj(0,'type','figure','tag','zplotbench2');
if isempty(h)
	h=figure('tag','zplotbench2');
else
	figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	   'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',3,'units','normalized','position',[0.03 0.03 0.95 0.95]);
	   
	   
subplot(2,1,1)
cc= get(gca,'colororder');
hold on
legstr={};
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		plot(outc.t(ind),outc.pe(ind,1),'color',cc(kc,:));
		plot(outc.t(ind),outc.pe(ind,51),'linestyle','-.','color',cc(kc,:));
		legstr{end+1} = sprintf('r = 0  , cn = %g, amorti = %g',outc.cn,outc.amorti);
		legstr{end+1} = sprintf('r = 0.5, cn = %g, amorti = %g',outc.cn,outc.amorti);
	end
	
end
xlabel('temps (s)');
ylabel('Pe (Pa)');
legend(legstr,-1);
title('Bench');

subplot(2,1,2)
hold on
nb =0;
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		tp = outc.t(ind) -outc.t(ind(1));
		tt =  cat(1,-tp(end:-1:2),tp);
		xne0 = xcov(outc.pe(ind,1),ref.pe(ind,1))./xcov(ref.pe(ind,1),ref.pe(ind,1));
		xnem = xcov(outc.pe(ind,51),ref.pe(ind,51))./xcov(ref.pe(ind,51),ref.pe(ind,51));
		plot(tt,xne0,'color',cc(kc,:));
		nb = nb +1;
	end
	
end
xlabel('Delta(temps) (s)');
ylabel('XCov(Pe-Peref,Peref)');
legend(legstr{1:nb},-1);
drawnow
pos =get(gca,'pos');
legend off
title('Bench');
set(gca,'pos',pos);
	   

h=findobj(0,'type','figure','tag','zplotbench3');
if isempty(h)
	h=figure('tag','zplotbench3');
else
	figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	   'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',3,'units','normalized','position',[0.03 0.03 0.95 0.95]);
	   
	   
subplot(2,1,1)
cc= get(gca,'colororder');
hold on
legstr={};
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		plot(outc.t(ind),outc.pion(ind,1),'color',cc(kc,:));
		plot(outc.t(ind),outc.pion(ind,51),'linestyle','-.','color',cc(kc,:));
		legstr{end+1} = sprintf('r = 0  , cn = %g, amorti = %g',outc.cn,outc.amorti);
		legstr{end+1} = sprintf('r = 0.5, cn = %g, amorti = %g',outc.cn,outc.amorti);
	end
	
end
xlabel('temps (s)');
ylabel('Pion (Pa)');
legend(legstr,-1);
title('Bench');

subplot(2,1,2)
hold on
nb =0;
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		tp = outc.t(ind) -outc.t(ind(1));
		tt =  cat(1,-tp(end:-1:2),tp);
		xne0 = xcov(outc.pion(ind,1),ref.pion(ind,1))./xcov(ref.pion(ind,1),ref.pion(ind,1));
		xnem = xcov(outc.pion(ind,51),ref.pion(ind,51))./xcov(ref.pion(ind,51),ref.pion(ind,51));
		plot(tt,xne0,'color',cc(kc,:));
		nb = nb +1;
	end
	
end
xlabel('Delta(temps) (s)');
ylabel('XCov(Pion-Pionref,Pionref)');
legend(legstr{1:nb},-1);
drawnow
pos =get(gca,'pos');
legend off
title('Bench');
set(gca,'pos',pos);
	   

h=findobj(0,'type','figure','tag','zplotbench4');
if isempty(h)
	h=figure('tag','zplotbench4');
else
	figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	   'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',3,'units','normalized','position',[0.03 0.03 0.95 0.95]);
	   
	   
subplot(2,1,1)
cc= get(gca,'colororder');
hold on
legstr={};
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		plot(outc.t(ind),outc.jmoy(ind,1),'color',cc(kc,:));
		plot(outc.t(ind),outc.jmoy(ind,51),'linestyle','-.','color',cc(kc,:));
		legstr{end+1} = sprintf('r = 0  , cn = %g, amorti = %g',outc.cn,outc.amorti);
		legstr{end+1} = sprintf('r = 0.5, cn = %g, amorti = %g',outc.cn,outc.amorti);
	end
	
end
xlabel('temps (s)');
ylabel('Jmoy diff (A*m^-^2)');
legend(legstr,-1);
title('Bench');

subplot(2,1,2)
hold on
legstr={};
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		plot(outc.t(ind),outc.jmoyeq(ind,1),'color',cc(kc,:));
		plot(outc.t(ind),outc.jmoyeq(ind,51),'linestyle','-.','color',cc(kc,:));
		legstr{end+1} = sprintf('r = 0  , cn = %g, amorti = %g',outc.cn,outc.amorti);
		legstr{end+1} = sprintf('r = 0.5, cn = %g, amorti = %g',outc.cn,outc.amorti);
	end
	
end
xlabel('temps (s)');
ylabel('Jmoy eq (A*m^-^2)');
legend(legstr,-1);

drawnow
pos =get(gca,'pos');
legend off
set(gca,'pos',pos);
	   

h=findobj(0,'type','figure','tag','zplotbench5');
if isempty(h)
	h=figure('tag','zplotbench5');
else
	figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	   'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',3,'units','normalized','position',[0.03 0.03 0.95 0.95]);
	   
	   
subplot(2,1,1)
cc= get(gca,'colororder');
hold on
legstr={};
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		plot(outc.t(ind),outc.qdiff(ind,1),'color',cc(kc,:));
		plot(outc.t(ind),outc.qdiff(ind,51),'linestyle','-.','color',cc(kc,:));
		legstr{end+1} = sprintf('r = 0  , cn = %g, amorti = %g',outc.cn,outc.amorti);
		legstr{end+1} = sprintf('r = 0.5, cn = %g, amorti = %g',outc.cn,outc.amorti);
	end
	
end
xlabel('temps (s)');
ylabel('q diff');
legend(legstr,-1);
title('Bench');

subplot(2,1,2)
hold on
legstr={};
for k =1:length(out)
	outc =out{k};
	if ~isempty(outc)
		kc   = rem((k-1),size(cc,1))+1 ;
		plot(outc.t(ind),outc.qeq(ind,1),'color',cc(kc,:));
		plot(outc.t(ind),outc.qeq(ind,51),'linestyle','-.','color',cc(kc,:));
		legstr{end+1} = sprintf('r = 0  , cn = %g, amorti = %g',outc.cn,outc.amorti);
		legstr{end+1} = sprintf('r = 0.5, cn = %g, amorti = %g',outc.cn,outc.amorti);
	end
	
end
xlabel('temps (s)');
ylabel('q eq');
legend(legstr,-1);

drawnow
pos =get(gca,'pos');
legend off
set(gca,'pos',pos);
	   

	
	
	
	
