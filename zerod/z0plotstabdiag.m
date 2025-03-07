function z0plotstabdiag(post)

% script du  plot plasma stability indicators
h = findobj(0,'type','figure','tag','z0plotstabdiag');
if isempty(h)
       h=figure('tag','z0plotstabdiag');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

% vertical stability index
%  ref: P.C. de Vries et al, NF 58 (2018) 026019 and refrence within
Ka   = post.zerod.vp ./ (2*pi^2.* post.z0dinput.geo.R.*post.z0dinput.geo.a.^2);
ms = (1.47 .* (1 + exp(-(2 .* post.zerod.li +1))) ./ 2 ./ (Ka - 1.13) - 1) .* (1 + 0.6 .* (post.zerod.betaptot - 0.1));
ms(ms < 0) = 0;
ms(ms > pi) = pi;
ms(Ka < 1.13) = pi;

% disruption diagramm
% from EXS/P6-37 D.Humphreys/IAEA FEC Kyoto/October 2016
subplot(2,4,1)
q0span = linspace(1,ceil(max(post.zerod.q0)),1001);
qlim   = 1./ q0span .^ 20 - 0.5 + 2 .* q0span + 0.2 .* sin(4 .* pi .* (q0span -0.5));  
plot(post.zerod.q0,post.zerod.q95,'r',q0span,qlim,'k',[1,1],[qlim(1),ceil(max(max(qlim),max(post.zerod.q95)))],'k');
xlabel('q_0')
ylabel('q_{95}');
legend('METIS','Stability limit (disruption)','Location','best');
hold on
plot(post.zerod.q0,post.zerod.q95,'.r')
title(sprintf('METIS : %s@%d/ stability diagram', ...
		    post.z0dinput.machine,post.z0dinput.shot));
set(gca,'ButtonDownFcn','zdataplot(''extrait'');');
set(gca,'xlim',[0,ceil(min(10,max(post.zerod.q0)))],'ylim',[0,ceil(min(20,max(post.zerod.q95)))]);

subplot(2,4,2)
qedge = linspace(min(qlim),10);
lilim    = 0.12 .* qedge + 0.6;
plot(post.zerod.q95,post.zerod.li,'r',[min(qlim),min(qlim),2.7,ceil(max(post.zerod.q95))],[lilim(1),0.7,0.5,0.5],'k',qedge,lilim,'k');
ylabel('l_i')
xlabel('q_{95}');
legend('METIS','Stability limit','Location','best');
hold on
plot(post.zerod.q95,post.zerod.li,'.r')
set(gca,'ButtonDownFcn','zdataplot(''extrait'');');
set(gca,'xlim',[0,ceil(min(20,max(post.zerod.q95)))]);

subplot(2,4,3)
Kamax = max(2.5,ceil(max(Ka)*2)/2);
plot(Ka,ms,'r',[1.3,1.75,NaN,2.25,1.75],[2,0,NaN,0,2],'b-.',[1,Kamax],[0.4,0.4],'g',[1,Kamax],[0.25,0.25],'k-.',[1,Kamax],[0.15,0.15],'k')
xlabel('K')
ylabel('ms vertical stability index');
title('lower ms value = more vertical instable plasma');
legend('METIS','Experimental envelope','reasonable limit','limit with in vessel coils','ultimate limit','Location','best');
hold on
plot(Ka,ms,'.r')
set(gca,'ButtonDownFcn','zdataplot(''extrait'');');

subplot(2,4,4);
iota = linspace(min(1./(post.zerod.q95 +1)),0.5);
fi   = mean(post.zerod.negr .* post.zerod.q95);
limh = min(1,fi .* iota .* mean(post.z0dinput.geo.R) ./ mean(post.z0dinput.geo.b0) ./ 1e20);
fi   = mean(post.zerod.ip .* post.zerod.q95);
jave =  iota .* fi ./ pi ./  mean(post.z0dinput.geo.a) .^ 2 ./ mean(Ka);
pp = polyfit([1e5,2e6],[1e18,1.5e19],1);
nlimrunaway = polyval(pp,jave);
plot(post.zerod.nbar .* post.z0dinput.geo.R ./post.z0dinput.geo.b0 ./ 1e20,1./post.zerod.q95,'r',limh,iota,'k',[max(limh),0.1],[0.5,0.5],'k', ...
nlimrunaway.* mean(post.z0dinput.geo.R) ./ mean(post.z0dinput.geo.b0) ./ 1e20,iota,'k', ...
[min(nlimrunaway.* mean(post.z0dinput.geo.R) ./ mean(post.z0dinput.geo.b0) ./ 1e20),min(limh)],[min(iota),min(iota)],'k');
ylabel('\iota_{95}')
xlabel('n_{bar} R / B_t (10^{20} m^{-2} T^{-1})');
set(gca,'ylim',[0,0.6]);
legend('METIS','Stability limit','Location','best');
hold on
plot(post.zerod.nbar .* post.z0dinput.geo.R ./post.z0dinput.geo.b0 ./ 1e20,1./post.zerod.q95,'.r')
set(gca,'ButtonDownFcn','zdataplot(''extrait'');');
title('Hugill diagram');

subplot(2,2,3);
plot(post.zerod.temps,Ka,post.zerod.temps,post.zerod.li,post.zerod.temps,post.zerod.betaptot, ...
     post.zerod.temps,1 ./ post.zerod.q95,post.zerod.temps,post.zerod.q0 ./ post.zerod.q95, ...
     post.zerod.temps,post.zerod.nbar .* post.z0dinput.geo.R ./post.z0dinput.geo.b0 ./ 1e20);
xlabel('time (s)');
legend('K (elongation)','l_i','\beta_p','\iota_{95} (1/q_{95})','q_0 / q_{95}','n_{bar} R / B_t (10^{20} m^{-2} T^{-1})','Location','best')
title('limits are sketched for illustration')
set(gca,'ButtonDownFcn','zdataplot(''extrait'');');

subplot(2,2,4);
vt = ones(size(post.zerod.temps));
plot(post.zerod.temps,ms,'r',post.zerod.temps,0.4 * vt,'g',post.zerod.temps,0.25 * vt,'k-.',post.zerod.temps,0.15 * vt,'k');
xlabel('time (s)');
ylabel('ms vertical stability index');
title('lower ms value = more vertical instable plasma');
legend('METIS','reasonable limit','limit with in vessel coils','ultimate limit','Location','best');
set(gca,'ButtonDownFcn','zdataplot(''extrait'');');

