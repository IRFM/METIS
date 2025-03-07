% plot des IP
t=data.gene.temps;
figure('color',[1 1 1],'defaultaxesfontsize',12);
plot(t,data.gene.ip,'xr',t,data.gene.ipepar,'om',t,data.gene.ipjmoy,'c+', ...
     t,data.gene.ipohm,'g',t,data.exp.ip,'b');
title(['different Ip (choc # ',int2str(param.from.shot.num),')'])
xlabel('temps (s)')
ylabel('Ip (A)');
legend('analytique','E/eta+Jni','int(Jmoy)','Ohm','mesure',-1);

