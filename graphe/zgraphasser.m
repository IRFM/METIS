h = findobj(0,'type','figure','tag','asser');
if isempty(h)
       h=figure('tag','asser');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])



ind = param.gene.kmin:param.gene.k;
%ind = param.gene.kmin:108;

subplot(4,1,1)
plot(data.gene.temps(ind),data.gene.paddhyb(ind)/1e6,'r')
ylabel('P_L_H (MW)')
title(['file : [',param.gene.origine(22:end),'] dt_a_s_s_e_r=',num2str(param.cons.asser.dt,2),' s'])
text('units','normalized','position',[0.05 0.2],'string',...
     ['D_P_L_H = ',num2str(param.cons.asser.plhinc,2),' MW'])
subplot(4,1,2)
plot(data.gene.temps(ind),data.gene.nbar(ind))
ylabel('nbar (m^-^3)')
text('units','normalized','position',[0.05 0.6],'string',...
     ['D_i_n_j_e_c_t_i_o_n = ',num2str(param.cons.asser.injgazinc,2),' p/m^3'])

subplot(4,1,3)
val = data.cons.c(ind,1)-data.cons.pomp(ind);
val=val./max(val);
plot(data.gene.temps(ind),val,...
     data.gene.temps(ind),data.gene.paddfus(ind)./max(data.gene.paddfus(ind)))
legend('injection','paddfus',2)
ylabel('u.a.')
subplot(4,1,4)
%sourtot = trapz(param.gene.x,data.source.totale.ne(ind,:)'.*data.equi.vpr(ind,:)')'.*data.equi.rhomax(ind);
%fluxsor = trapz(param.gene.x,data.prof.flux.sortant.ge(ind,:)'.*data.equi.vpr(ind,:)')'.*data.equi.rhomax(ind);
%plot(data.gene.temps(ind),sourtot,'r',...
%     data.gene.temps(ind),fluxsor,'b')
%legend('source','flux sortant')
plot(data.gene.temps(ind),data.gene.ip(ind)/1e6)
ylabel('MA')
