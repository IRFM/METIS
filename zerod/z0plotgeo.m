% script du  plot de la geometrie du 0d
zs   = post.zerod;
cons = post.z0dinput.cons;
exp0d  = post.z0dinput.exp0d;
geo  = post.z0dinput.geo;
t    = zs.temps;
% precalcul

% limite de l'elongation pour la stabilite verticale 
% loi de J-L Duchateau
k95lim = 2.22 - 0.17 .* geo.R ./ geo.a;
% de iter basis et FDR :
klim   = sqrt(1.25 .* (1+ k95lim .^ 2) - 1);

h = findobj(0,'type','figure','tag','z0plotgeo');
if isempty(h)
       h=figure('tag','z0plotgeo');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 9;
if post.z0dinput.mode_exp == 0

   subplot(k,1,1);	
   plot(t,geo.R,'r',data.gene.temps,data.geo.r0 );
   ylabel('R (m)')
   title(sprintf('Zerod : %s@%d/geometry ( r -> 0D, b: -> cronos, g -> vertical stability limit)', ...
          post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,geo.a,'r',data.gene.temps,data.equi.a(:,end),':b');
   ylabel('a (m)');
   subplot(k,1,3);	
   plot(t,geo.K,'r',data.gene.temps,data.equi.e(:,end),':b',t,k95lim,'g:',t,klim,'g');
   ylabel('b/a');
   subplot(k,1,4);	
   plot(t,geo.d,'r',data.gene.temps,0.5 .*( data.equi.trh(:,end) + data.equi.trl(:,end)) ,':b');
   ylabel('triangularity');
   subplot(k,1,5);
   plot(t,geo.b0,'r',data.gene.temps,data.equi.F(:,end) ./ data.geo.r0,':b');
   ylabel('B0 (T)');
   subplot(k,1,6);	
   plot(t,zs.vp,'r',data.gene.temps,data.gene.volume,':b');
   ylabel('Vol (m^3)');
   subplot(k,1,7);	
   plot(t,zs.sp,'r',data.gene.temps,data.gene.surface,':b');
   ylabel('S (m^2)');
   subplot(k,1,8);	
   plot(t,zs.sext,'r',data.gene.temps,data.equi.vpr(:,end).* data.equi.grho(:,end),':b');
   ylabel('Sext (m^2)');
   subplot(k,1,9);	
   plot(t,zs.peri,'r',data.gene.temps,data.equi.spr(:,end) .* data.equi.grho(:,end),':b');
   ylabel('LCMS length (m)');
   xlabel('time (s)');
else
   subplot(k,1,1);	
   plot(t,geo.R,'r');
   ylabel('R (m)')
   title(sprintf('Zerod : %s@%d/geometry ( r -> 0D, b: -> experiment, g -> vertical stability limit)', ...
          post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,geo.a,'r');
   ylabel('a (m)');
   subplot(k,1,3);	
   plot(t,geo.K,'r',t,k95lim,'g:',t,klim,'g');
   ylabel('b/a');
   subplot(k,1,4);	
   plot(t,geo.d,'r');
   ylabel('triangularity');
   subplot(k,1,5);	
   plot(t,geo.b0,'r');
   ylabel('B0 (T)');
   subplot(k,1,6);	
   plot(t,zs.vp,'r',t,exp0d.vp,':b');
   ylabel('Vol (m^3)');
   subplot(k,1,7);	
   plot(t,zs.sp,'r',t,exp0d.sp,':b');
   ylabel('S (m^2)');
   subplot(k,1,8);	
   plot(t,zs.sext,'r',t,exp0d.sext,':b');
   ylabel('Sext (m^2)');
   subplot(k,1,9);	
   plot(t,zs.peri,'r');
   ylabel('LCMS length (m)');
   xlabel('time (s)');

end
   
joint_axes(h,k);


z0plot_sepa_evolution;