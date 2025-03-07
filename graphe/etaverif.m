% plot des IP
t = data.gene.temps;
x = param.gene.x;


lnL  = 31.474 + log((data.prof.ne .^ (-0.5)) .* data.prof.te);
etas = 0.51 * sqrt(param.phys.me * param.phys.e) / 3 / param.phys.epsi0 / param.phys.epsi0 / ...
       ((2 * pi) ^ 1.5) * 1.97 / 3.40 .* lnL .* data.prof.zeff .* ((data.prof.te/1) .^ (-1.5)) .* ...
       (2.67 + data.prof.zeff) ./ (1.13 + data.prof.zeff);


eta = (data.prof.jeff - data.source.totale.j)./(data.prof.epar);
figure('color',[1 1 1],'defaultaxesfontsize',12);
plotprof(gca,t,x,1./data.neo.eta,'linestyle','+','color',[0 1 1]);
plotprof(gca,t,x,1./data.coef.eta,'linestyle','x','color',[1 0 1]);
plotprof(gca,t,x,eta,'linestyle','o','color',[1 0 0]);
plotprof(gca,t,x,1./etas,'linestyle','-','color',[0 0 0]);
%set(gca,'yscale','log');
title(['differentes resistivite (choc # ',int2str(param.from.shot.num),')'])
xlabel('x (s)')
ylabel('1/eta (Ohm^-^1 m^-^1)');
legend('neoclass','diffusion','(J-Jni)/Epar','spitzer');

