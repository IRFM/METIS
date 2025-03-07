% script du  plot de la convergence du 0d
zs   = post.zerod;
t    = zs.temps;
% precalcul

h = findobj(0,'type','figure','tag','z0plotconv');
if isempty(h)
       h=figure('tag','z0plotconv');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 4;
subplot(k,1,1);	
plot(t,100 .* zs.dw,'r');
ylabel('dW/W (%)')
title(sprintf('Zerod : %s@%d/convergence 0D in %d loops',post.z0dinput.machine,post.z0dinput.shot,zs.nb));
subplot(k,1,2);	
plot(t,100 .* zs.dpfus,'r');
ylabel('dPfus/Pfus (%)');
subplot(k,1,3);	
plot(t,100 .* zs.dini,'r');
ylabel('dIni/Ini (%) ');
subplot(k,1,4);	
plot(t,100 .* zs.diboot,'r');
ylabel('dIboot/Iboot (%)');
xlabel('time (s)')
joint_axes(h,k);
