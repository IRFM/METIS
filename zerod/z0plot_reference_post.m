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

t    = post.z0dinput.cons.temps;
subplot(3,1,1);	
plot(t,post.z0dinput.cons.ip./1e6);
ylabel('I_p (MA)')
title(sprintf('Zerod : %s@%d / Overview', ...
          post.z0dinput.machine,post.z0dinput.shot));
subplot(3,1,2);	
hl = plot(t,post.z0dinput.cons.pecrh./1e6,t,post.z0dinput.cons.picrh./1e6,t,post.z0dinput.cons.plh./1e6,t,real(post.z0dinput.cons.pnbi)./1e6 + imag(post.z0dinput.cons.pnbi)./1e6 );
if post.z0dinput.option.lhmode == 5
    legend('ECRH1','ICRH','ECRH2','NBI');
else
    legend('ECRH','ICRH','LH','NBI');
end
ylabel('MW');
subplot(3,1,3)
plot(t,post.z0dinput.cons.nbar ./ 1e19);
ylabel('n_{bar} (10^{19} m^{-3})')
xlabel('time (s)');
edition2