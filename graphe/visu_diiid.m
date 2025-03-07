t = data.gene.temps;
x = param.gene.x;

h = findobj(0,'type','figure','tag','fig1');
if isempty(h)
       h=figure('tag','fig1');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(2,1,1)
plot(t,data.gene.iboot/1e6,'r',t,data.gene.iidn/1e6,'b')
axis([0.5 1.5 0 0.2]);
xlabel('t (s)')
ylabel('MA')
legend('I_b_o_o_t','I_N_B_I')
title('shot 119562, non inductive current')
subplot(2,2,3)
plot(x,data.neo.jboot(1,:)/1e6,'b',x,data.neo.jboot(end,:)/1e6,'r--')
legend('t=0.5s','t=1.5s')
ylabel('MA/m^2')
xlabel('x') 
title('bootstrap current')	
subplot(2,2,4)
plot(x,data.source.idn.j(1,:)/1e6,'b',x,data.source.idn.j(end,:)/1e6,'r--')
ylabel('MA/m^2')
xlabel('x') 
title('NB current')	

h = findobj(0,'type','figure','tag','fig2');
if isempty(h)
       h=figure('tag','fig2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
subplot(1,2,1) 
plot(t,data.cons.flux,'r--',t,data.prof.psi(:,end),'r')
title('edge flux')
xlabel('t (s)')
ylabel('Weber')
legend('exp.','CRONOS')
axis([0.5 1.5 0.1 0.35])

subplot(2,2,2)
plot(t,data.prof.q(:,15),'r',t,data.prof.q(:,30),'b',t,data.prof.q(:,60),'m',[0.5 1.5],[1 1],'k')
axis([0.5 1.5 0 4])	
legend('x=0.15','x=0.3','x=0.6')
xlabel('t (s)')
title('safety factor')
subplot(2,2,4)
plot(x,data.prof.q(1,:),x,data.prof.q(26,:),x,data.prof.q(end,:),[0 1],[1 1])
legend('t=0.5s','t=1s','t=1.5s')
title('q profile')
xlabel('x')


h = findobj(0,'type','figure','tag','fig3');
if isempty(h)
       h=figure('tag','fig3');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(3,1,1)
plot(t,data.gene.paddidn/1e6,'r',t,data.gene.nbar/1e19,'b',t,data.gene.ip/1e6,'m')
axis([0.5 1.5 0 4])
title('shot 119562, CRONOS simulation')
%legend('Pnbi (MW)','nbar (10^1^9 m^-^3)','Ip (MA')
subplot(3,1,2)
plot(t,data.gene.li,'b',diiidtemp.tli,diiidtemp.li,'b--')
hold on
plot(t,data.gene.vloop,'r',diiidtemp.tvl,diiidtemp.vl,'r--')
axis([0.5 1.5 0 2])
subplot(3,1,3)
plot(t,data.gene.wdia/1e6,'r',diiidtemp.twdia,diiidtemp.wdia/1e6,'r--')
hold on
plot(t,data.gene.paddohm/1e6,'b',diiidtemp.tpoh,diiidtemp.poh/1e6,'b--')
axis([0.5 1.5 0 1])
%legend('Wdia (MJ)','','Pohm (MW)')
xlabel('t (s)')

