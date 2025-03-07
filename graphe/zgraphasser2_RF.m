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
tigor = data.gene.temps(ind);
Plh   = data.gene.paddhyb(ind)/1e6;
nbar  = data.gene.nbar(ind);
Pfus  = data.gene.paddfus(ind);
Ip    = data.gene.ip(ind)/1e6;
subplot(4,1,1);
plot(tigor,data.gene.paddhyb(ind)/1e6,'r');
ylabel('P_L_H (MW)')
set(gca,'xticklabel',[],'xgrid','on')
if isfield(param.cons.asser,'dt')
title(['dt_a_s_s_e_r=',num2str(param.cons.asser.dt,2),' s'])
end
text('units','normalized','position',[0.05 0.2],'string',...
     ['D_P_L_H = ',num2str(param.cons.asser.plhinc,2),' MW'])
subplot(4,1,2)
plot(data.gene.temps(ind),data.gene.nbar(ind))
set(gca,'xticklabel',[],'xgrid','on')
ylabel('nbar (m^-^3)')
text('units','normalized','position',[0.05 0.6],'string',...
     ['D_i_n_j_e_c_t_i_o_n = ',num2str(param.cons.asser.injgazinc,2),' p/m^3'])


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
fgreenwald = fgrew(ind);
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
Ipmax = val(ind);
plot(data.gene.ip(ind)/1e6,data.gene.nbar(ind)/1e19,...
     data.gene.ip(ind)/1e6,val(ind),'ro')
xlabel('I_p (MA)')
ylabel('nbar (10^1^9 m^-^3)')     

Ti0 = data.prof.ti(ind,1);
save file_igor_RF2003_simu5 tigor Plh nbar Pfus Ipmax fgreenwald Ti0 Ip 
