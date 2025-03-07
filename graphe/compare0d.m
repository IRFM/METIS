% script de plot des temperatures

chargedata
h = findobj(0,'type','figure','tag','compare');
if isempty(h)
       h=figure('tag','compare');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

ind  = param.gene.kmin:min(param.gene.k,length(t));
tt   = t(ind);



  subplot(2,1,1)
  plot(tt,data.prof.te(ind,1),'r',tt,post.zerod.te0(ind),'+');

  ylabel('Te(0) eV')
  legend('mesure','0D');
    if ~isempty(bile.ip)
    st =sprintf('shot %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d',param.from.machine, ...   
               fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
               param.gene.nbrho);
  else
    st =sprintf('shot %s simulation, t = %4.2g:%4.2g:%4.2g, nbrho = %d',param.from.machine, ...   
               param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
               param.gene.nbrho);

  end
  title(st);

  subplot(2,1,2)
  plot(tt,data.gene.taue(ind),'r',tt,post.zerod.taue(ind),'+');
  title('taue ')
  ylabel('s')
  xlabel('time (s)')


