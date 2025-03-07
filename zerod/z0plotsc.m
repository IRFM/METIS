% script du  plot des puissance 0d
h = findobj(0,'type','figure','tag','z0dsc');
if isempty(h)
       h=figure('tag','z0dsc');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

zs   = post.zerod;
ipar = max(0,zs.ipar);
    
cons = post.z0dinput.cons;
exp0d= post.z0dinput.exp0d;
t    = zs.temps;
subplot(3,1,1);	
plot(t,zs.ieccd/1e6,t,zs.ifwcd./1e6,t,zs.ilh./1e6,t,real(zs.inbicd)./1e6 + imag(zs.inbicd)./1e6,t,zs.iboot./1e6, ...
     t,zs.iohm./1e6,t,ipar./1e6,t,zs.ifus./1e6,'-.');
if post.z0dinput.option.lhmode == 5
  legend('ECCD 1','ICRH','ECCD 2','NBI','Boot','Ohm','Ip (//B)','Alpha');
else
  legend('ECCD','ICRH','LHCD','NBI','Boot','Ohm','Ip (//B)','Alpha');
end
ylabel('MA')
title(sprintf('Zerod : %s@%d / Overview', ...
          post.z0dinput.machine,post.z0dinput.shot));
subplot(3,1,2);	
hl = plot(t,zs.pecrh./1e6,t,zs.picrh./1e6,t,zs.plh./1e6,t,real(zs.pnbi)./1e6 + imag(zs.pnbi)./1e6 , ...
     t,zs.pfus./1e6,t,zs.pohm./1e6,t,zs.prad./1e6,'-.', ...
     t,zs.pbrem./1e6,'-.',t,zs.pcyclo./1e6,'-.');
if post.z0dinput.option.lhmode == 5
    legend('ECRH 1','ICRH','ECRH 2','NBI','Alpha','Ohm','Rad','Brem','cyclo');
else
    legend('ECRH','ICRH','LH','NBI','Alpha','Ohm','Rad','Brem','cyclo');
end
hold on 
plot(t,cons.picrh./1e6,'linestyle',':','color',get(hl(2),'color'));
plot(t,cons.plh./1e6,'linestyle',':','color',get(hl(3),'color'));

ylabel('MW');
subplot(3,1,3)
plot(t,zs.vloop,'b',t,zs.li,'r',t,exp0d.li,'-.r',t,zs.betap,'g',t,zs.betaptot,'g:');
legend('Vloop','li 0D','li ref','betap (th)','betap (tot)');
set(gca,'ylim',[-0.1,3]);
ylabel('V, -, -')
xlabel('time (s)');
