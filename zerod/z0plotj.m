% script du  plot des courants 0d
h = findobj(0,'type','figure','tag','z0plotj');
if isempty(h)
       h=figure('tag','z0plotj');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

zs   = post.zerod;
exp0d  = post.z0dinput.exp0d;
cons = post.z0dinput.cons;
t    = zs.temps;
k    = 9;
if post.z0dinput.mode_exp == 0
   subplot(k,1,1);	
   plot(t,zs.ip ./ 1e6,'r',data.gene.temps,data.gene.ip ./ 1e6,':b');
   ylabel('Ip (MA)')
   title(sprintf('Zerod : %s@%d/ccurrent ( r -> 0D, b: -> cronos)', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.iboot./1e6,'r',data.gene.temps,data.gene.iboot./1e6,':b');
   ylabel('Iboot (MA)');
   subplot(k,1,3);	
   plot(t,zs.iohm./1e6,'r',data.gene.temps,data.gene.ipohm./1e6,':b',t,zs.irun./1e6,'g',t,zs.ipar./1e6,'c:');
   ylabel('Iohm, I runaway & I_{ // B} (MA)');
   subplot(k,1,4);	
   plot(t,zs.icd./1e6,'r',data.gene.temps,data.gene.icd./1e6,':b');
   ylabel('Icd (MA)');
   subplot(k,1,5);	
   plot(t,zs.ilh./1e6,'r',data.gene.temps,data.gene.ihyb./1e6,':b');
   if post.z0dinput.option.lhmode == 5
    ylabel('Ieccd 2 (MA)');
   else
    ylabel('Ilh (MA)');
   end
   subplot(k,1,6);	
   plot(t,zs.ieccd./1e6,'r',data.gene.temps,data.gene.ifce./1e6,':b');
   ylabel('Ieccd (MA)');
   subplot(k,1,7);	
   plot(t,zs.ifwcd./1e6,'r',data.gene.temps,data.gene.ifci./1e6,':b');
   ylabel('Ifwcd (MA)');
   subplot(k,1,8);	
   plot(t,real(zs.inbicd)./1e6 + imag(zs.inbicd)./1e6,'r',data.gene.temps,data.gene.iidn./1e6,':b');
   ylabel('Inbicd (MA)');
   subplot(k,1,9);
   if isfield(data.gene,'ifus')	
       plot(t,zs.ifus./1e6,'r',data.gene.temps,data.gene.ifus./1e6,':b');
   else
       plot(t,zs.ifus./1e6,'r');
   end
   ylabel('Ialpha (MA)');
   xlabel('time (s)');
else
   subplot(k,1,1);	
   plot(t,zs.ip ./ 1e6,'r',t,exp0d.ip ./ 1e6,':b');
   ylabel('Ip (MA)')
   title(sprintf('Zerod : %s@%d/ccurrent ( r -> 0D, b: -> experiment)', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.iboot./1e6,'r',t,exp0d.iboot./1e6,':b');
   ylabel('Iboot (MA)');
   subplot(k,1,3);	
   plot(t,zs.iohm./1e6,'r',t,exp0d.iohm./1e6,':b',t,zs.irun./1e6,'g',t,zs.ipar./1e6,'m-.');
   ylabel('Iohm (r) , I_{runaway} (g) & I_{ // B} (m) (MA)');
   subplot(k,1,4);	
   plot(t,zs.icd./1e6,'r',t,exp0d.icd./1e6,':b');
   ylabel('Icd (MA)');
   subplot(k,1,5);	
   plot(t,zs.ilh./1e6,'r',t,exp0d.ilh./1e6,':b');
   if post.z0dinput.option.lhmode == 5
      ylabel('Ieccd 2 (MA)');
   else
      ylabel('Ilh (MA)');
   end
   subplot(k,1,6);	
   plot(t,zs.ieccd./1e6,'r',t,exp0d.ieccd./1e6,':b');
   if post.z0dinput.option.lhmode == 5
      ylabel('Ieccd 1 (MA)');
   else
      ylabel('Ieccd (MA)');
   end
   subplot(k,1,7);	
   plot(t,zs.ifwcd./1e6,'r',t,exp0d.ifwcd./1e6,':b');
   ylabel('Ifwcd (MA)');
   subplot(k,1,8);	
   plot(t,real(zs.inbicd)./1e6 + imag(zs.inbicd)./1e6,'r',t,exp0d.inbicd./1e6,':b');
   ylabel('Inbicd (MA)');
   subplot(k,1,9);	
   if isfield(exp0d,'ifus')
      plot(t,zs.ifus./1e6,'r',t,exp0d.ifus./1e6,':b');
   else
      plot(t,zs.ifus./1e6,'r');
   end
   ylabel('Ialpha (MA)');
   xlabel('time (s)');

end

hh = findobj(h,'type','axes');
for l=1:k
   lim = get(hh(l),'ylim');
   set(hh(l),'ylim',[0,max(max(lim),0.1)]);
end
joint_axes(h,k);


z0plot_passive_current;