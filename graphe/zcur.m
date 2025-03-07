
h = findobj(0,'type','figure','tag','current');
if isempty(h)
       h=figure('tag','current');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])



ind = param.gene.kmin:(param.gene.k-1);

Ilh   = data.gene.ihyb(ind)/1e3;
jlh   = data.source.hyb.j(ind,:)/1e6;
Inbi  = data.gene.iidn(ind)/1e3;
Ip    = data.gene.ip(ind)/1e3;
Iboot = data.gene.iboot(ind)/1e3;
time  = data.gene.temps(ind);
mlh   = data.mode.hyb(ind);
Plh   = data.gene.paddhyb(ind);
Pcons = abs(data.cons.hyb(ind,:));
Ini   = Ilh + Iboot + Inbi;
if size(Pcons,2) > 1
  Pcons = sum(Pcons,2);
end
x     = param.gene.x;
subplot(2,2,1)
if mean(Inbi) == 0 & mean(Ilh) > 0
  ind2 = find(mlh == 2);
  ind3 = find(mlh == 3);
  plot(time,Ini,time,Iboot,time,Ip,time(ind2),Ilh(ind2),'ro',time(ind3),Ilh(ind3),'k*')
  legend('non ind.','bootstrap','Ip','LH (calc.)','LH (interp.)')
  ylabel('kA')
  title('current')
else
  ind2 = find(mlh == 2);
  ind3 = find(mlh == 3);
  plot(time,Ini,time,Inbi,time,Iboot,time,Ip,time(ind2),Ilh(ind2),'ro',time(ind3),Ilh(ind3),'k*')
  legend('non ind.','NBI','bootstrap','Ip','LH (calc.)','LH (interp.)')
  ylabel('kA')
  title('current')  
end

subplot(2,2,2)
if mean(Inbi) == 0 & mean(Ilh) > 0
  plot(time,Iboot./Ip*100,time,Ilh./Ip*100,time,Plh./Pcons*100)
  legend('I_b_o_o_t/I_p','I_L_H/I_p','Pabs/Pref')
  ylabel('%')
  title('')
else
   plot(time,Inbi./Ip*100,time,Iboot./Ip*100,time,Ilh./Ip*100,time,Plh./Pcons*100)
  legend('I_N_B_I/I_p','I_b_o_o_t/I_p','I_L_H/I_p','Pabs/Pref')
  ylabel('%')
  title('')
   
end

subplot(2,2,3)


contourf(x,time,jlh);
colorbar
ylabel('t (s)');
xlabel('x')
title('j_L_H current (MA/m^2)')

subplot(2,2,4)
if size(data.cons.hyb,2) > 1
plot(data.gene.temps(ind),data.gene.paddhyb(ind)/1e6,data.gene.temps(ind),sum(abs(data.cons.hyb(ind,:)),2)/1e6)
else
plot(data.gene.temps(ind),data.gene.paddhyb(ind)/1e6,data.gene.temps(ind),abs(data.cons.hyb(ind))/1e6,'--')
end

xlabel('t (s)')
ylabel('MW');
title('P_L_H')
