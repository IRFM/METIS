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
indt = find(param.gene.file=='_');
param.gene.file(indt) = '-';
indt = find(param.gene.origine=='_');
param.gene.origine(indt) = '-';

%ind = param.gene.kmin:108;

subplot(4,1,1);
plot(data.gene.temps(ind),data.gene.paddhyb(ind)/1e6,'r');
ylabel('P_L_H (MW)')
set(gca,'xticklabel',[],'xgrid','on')
if isfield(param.cons.asser,'dt')
title(['file : [',param.gene.origine(22:end),'] dt_a_s_s_e_r=',num2str(param.cons.asser.dt,2),' s'])
end
text('units','normalized','position',[0.05 0.2],'string',...
     ['D_P_L_H = ',num2str(param.cons.asser.plhinc,2),' MW'])
subplot(4,1,2)
plot(data.gene.temps(ind),data.gene.nbar(ind))
set(gca,'xticklabel',[],'xgrid','on')
ylabel('nbar (m^-^3)')
text('units','normalized','position',[0.05 0.6],'string',...
     ['D_n_b_a_r = ',num2str(param.cons.asser.nbarinc,2),' m^-^2'])
title(['asser : ',param.fonction.asser])

subplot(4,1,3)
plot(data.gene.temps(ind),data.gene.paddfus(ind))
set(gca,'xticklabel',[],'xgrid','on')
legend('paddfus',2)

subplot(4,1,4)
%sourtot = trapz(param.gene.x,data.source.totale.ne(ind,:)'.*data.equi.vpr(ind,:)')'.*data.equi.rhomax(ind);
%fluxsor = trapz(param.gene.x,data.prof.flux.sortant.ge(ind,:)'.*data.equi.vpr(ind,:)')'.*data.equi.rhomax(ind);
%plot(data.gene.temps(ind),sourtot,'r',...
%     data.gene.temps(ind),fluxsor,'b')
%legend('source','flux sortant')
fgrew = (data.gene.nbar/1e20)./(data.gene.ip/1e6).*3.1415.*(data.geo.a).^2;
plot(data.gene.temps(ind),data.gene.ip(ind)/1e6,data.gene.temps(ind),fgrew(ind),'r')
set(gca,'xgrid','on')
legend('Ip (MA)','fr Greenwald')

h = findobj(0,'type','figure','tag','asser2');
if isempty(h)
       h=figure('tag','asser2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
val=data.gene.ip/1e6./(pi.*(data.geo.a.^2))*10;

plot(jeux.data.gene.ip(ind)/1e6,jeux.data.gene.nbar(ind)/1e19,'g+')
hold on

plot(data.gene.ip(ind)/1e6,data.gene.nbar(ind)/1e19,...
     data.gene.ip(ind)/1e6,val(ind),'ro')
xlabel('I_p (MA)')
ylabel('nbar (10^1^9 m^-^3)')     


