% script de  la figure de verification de l'execution de zineb en diffusion du courant
x=param.gene.x;
t=data.gene.temps;
if length(t)> 100
	pas =fix(length(t)/100);
	if pas ==0
		pas =1;
	end
else
	pas =1;
end

ind=param.gene.kmin:min(param.gene.kmax,length(t));
tt=t(ind);
inde = 1:pas:length(t);
t=t(inde);

figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname','times','defaultlinelinewidth',3);
subplot(4,1,1)
plot(t,data.exp.ip(inde)/1e6,'or',tt,data.gene.ip(ind)/1e6,'r',tt,data.gene.iboot(ind)/1e6,'k')
legend('I_p exp','I_p diffusion','I_b_o_o_t_s_t_r_a_p',-1)
ylabel('(MA)')
axis([0 20 0 1])

st =sprintf('choc %s #%d, t = %g:%g:%g',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin);
title(st);
grid

subplot(4,1,2)
plot(t,data.exp.vloop(inde),'ro',tt,data.gene.vres(ind),'k');
legend('Vl exp.','Vl resistive',-1)
ylabel('Vloop (V)')
axis([0 20 0 2])
grid

subplot(4,1,3)
plot(t,data.exp.li(inde),'or',tt,data.gene.li(ind),'b',tt,data.equi.li(ind),'k')
legend('li exp.','li Bpol','li equilibre',-1)
ylabel('li')
axis([0 20 0.5 2.5])
grid

subplot(4,1,4)
ind95 =  min( iround(param.gene.x,0.95));
if strcmp(param.from.machine,'TS')
   disp('q au bord')
   ind95 =length(param.gene.x);
else
   disp('q a 0.95')
end
plot(t,data.exp.qa(inde),'or',tt,data.equi.q(ind,ind95),'k')
legend('q_a exp.','q_a equilibre',-1)
if strcmp(param.from.machine,'TS')
    ylabel('q_a')
else
    ylabel('q_95')

end
grid
xlabel('temps (s)')
axis([0 20 5 15])

