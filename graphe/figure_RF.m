[ilh,fgre] = zcalicible(data);

ecartip    = 100*((ilh+data.gene.iboot)-data.gene.ip)./data.gene.ip;
ng         = 1;
while ng<10
plot(data.gene.temps(ind),data.gene.paddfus(ind))
[xi,yi]    = ginput(2);
indng      = find(data.gene.temps > min(xi) & data.gene.temps < max(xi))';
eval(['ind',int2str(ng),'=indng;'])
ng=ng+1;
end
for i=1:(ng-1)
  eval(['ind=ind',int2str(i),';']);
  [a,b]      = polyfit(ecartip(ind),data.gene.paddfus(ind),1);
  nbardeb(i) = data.gene.nbar(ind(1));
  pente(i)   = a(1);
  origine(i) = a(2);
end
plot(data.gene.paddfus(ind1)-min(data.gene.paddfus(ind1)),ecartip(ind1))
hold on
plot(data.gene.paddfus(ind2)-min(data.gene.paddfus(ind2)),ecartip(ind2))
[a2,b2]      = polyfit(ecartip(ind2),data.gene.paddfus(ind2),1);
hold off


h = findobj(0,'type','figure','tag','augmentation');
if isempty(h)
       h=figure('tag','augmentation');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])


for i=1:(ng-1)
  eval(['ind=ind',int2str(i),';']);
  t    = data.gene.temps(ind);
%  fus  = data.gene.paddfus(ind)/max(data.gene.paddfus(ind));
  fus  = data.gene.paddfus(ind);
  ecip = abs(ecartip(ind))/100;
  eval(['subplot(3,3,',int2str(i),')'])
  plot(t,fus,'r',t,ecip,'b--')
  title(['n_b_a_r=',num2str(data.gene.nbar(ind(1))/1e19,3)])
  if i > 6
    xlabel('t (s)')
  end
  grid
  if i == 9
  
    legend('P_f_u_s','D_I_p')
  
  end
  axis([fix(min(t)) ceil(max(t)) 0 0.2])
end
