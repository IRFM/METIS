% script du  plot des puissance 0d
h = findobj(0,'type','figure','tag','z0plot_reference');
if isempty(h)
       h=figure('tag','z0plot_reference');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

t    = z0dinput.cons.temps;
subplot(3,1,1);	
plot(t,z0dinput.cons.ip./1e6);
ylabel('I_p (MA)')
title(sprintf('Zerod : %s@%d / Overview', ...
          z0dinput.machine,z0dinput.shot));
subplot(3,1,2);	
hl = plot(t,z0dinput.cons.pecrh./1e6,t,z0dinput.cons.picrh./1e6,t,z0dinput.cons.plh./1e6,t,real(z0dinput.cons.pnbi)./1e6 + imag(z0dinput.cons.pnbi)./1e6 );
legend('ECRH','ICRH','LH','NBI');
ylabel('MW');
subplot(3,1,3)
plot(t,real(z0dinput.cons.nbar) ./ 1e19);
ylabel('n_{bar} (10^{19} m^{-3})')
xlabel('time (s)');
edition2