% script du  plot des temperature 0d
zs   = post.zerod;
cons = post.z0dinput.cons;
exp  = post.z0dinput.exp;
t    = zs.temps;
% precalcul
if post.z0dinput.mode_exp == 0
	vpr  = data.equi.vpr;
	nvpr = trapz(param.gene.x,vpr,2);
	%   n1m = trapz(param.gene.x,vpr .* squeeze(sum(data.impur.impur(:,:,ind1),3)),2)./ nvpr;
end 



h = findobj(0,'type','figure','tag','z0plottuk');
if isempty(h)
       h=figure('tag','z0plott');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1])

k    = 4;
if post.z0dinput.mode_exp == 0

   subplot(k,1,1);	
   plot(t,zs.tem/1e3,'r',t,zs.te0/1e3,'r--');
   ylabel('keV')
   title(sprintf('Data validation : %s@%d/ r -> 0D, b: -> cronos)', ...
             post.z0dinput.machine,post.z0dinput.shot));	
   hold on 
   plot(t,data.gene.temoy ./ 1e3,'b:',t,data.prof.te(:,1) ./ 1e3,'b+');
   legend('<Te>','Te0')

   subplot(k,1,2);	
   plot(t,zs.w/1e6,'r');
   ylabel('MJ')
   hold on
   plot(t,data.gene.wdia ./ 1e6,'b:');
   legend('Wdia')

   subplot(k,1,3);	
   plot(t,zs.vloop,'r',t,data.gene.vres,':b',t,0.*t,'g');
   ylabel('V');
   set(gca,'ylim',[max(-0.1,min(zs.vloop)),min(3,max(zs.vloop))])
   legend('Vloop')

   subplot(k,1,4);	
   plot(t,zs.nem./1e19,'r',t,data.gene.nemoy./1e19,':b');
   ylabel('10^1^9 m^-^3');
   legend('<ne>')
   xlabel('time (s)')

end
hh = findobj(h,'type','axes');
%for l=1:k
%   lim = get(hh(l),'ylim');
%   set(hh(l),'ylim',[0,max(max(lim),0.1)]);
%end
%joint_axes(h,k);
