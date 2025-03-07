% test de gg0d
mode_expo_inte = 1;
morphing = 15;
rho_in = data.equi.rhomax * param.gene.x;
x = param.gene.x;
xa = data.equi.a ./ (data.equi.a(:,end) * ones(size(x)));
t = data.gene.temps;
tm = t * ones(size(x));
xm = ones(size(t)) * x;
raxe = data.equi.raxe;
a = data.equi.a(:,end);
kx = griddata(tm,xa,data.equi.e,tm,xm,'cubic');
dx = griddata(tm,xa, (data.equi.trh + data.equi.trl)/2,tm,xm,'cubic');
Fin = griddata(tm,xa,data.equi.F,tm,xm,'cubic');
rho_in = griddata(tm,xa,rho_in,tm,xm,'cubic');

[phi,dphidx,vpr,grho2r2,r2i,ri,C2,C3,grho,grho2,phiplasma] = gg0d(x,raxe,a,kx, ...
                                   dx,data.gene.volume,data.gene.surface,Fin,double(data.geo.R),double(data.geo.Z),rho_in,morphing);
rmx = sqrt(phi ./ pi .* ((data.equi.raxe(:,end) ./data.equi.F(:,end)) * ones(1,size(phi,2))));

grho2r2_tor  = tsplinet(rmx,grho2r2,rmx(:,end) * x);
r2i_tor      = tsplinet(rmx,r2i,rmx(:,end) * x);
ri_tor       = tsplinet(rmx,ri,rmx(:,end) * x);
C2_tor       = tsplinet(rmx,C2,rmx(:,end) * x);
C3_tor       = tsplinet(rmx,C3,rmx(:,end) * x);
vpr_tor      = tsplinet(rmx,vpr,rmx(:,end) * x);
spr_tor      = tsplinet(rmx,spr,rmx(:,end) * x);	
phi_tor      = tsplinet(rmx,phi,rmx(:,end) * x);	
%x_tor        = tsplinet(rmx,vt * x ,rmx(:,end) * x);	
grho2_tor    = tsplinet(rmx,grho2,rmx(:,end) * x);
grho_tor     = tsplinet(rmx,grho,rmx(:,end) * x);


figure(13);
clf
subplot(3,3,1)
plot(x,data.equi.vpr,'b',x,vpr_tor,'r');
subplot(3,3,2)
plot(x,(data.equi.rhomax * x),'b',x,rmx,'r');
subplot(3,3,3)
plot(x,data.equi.c2c,'b',x,C2_tor,'r');
subplot(3,3,4)
plot(x,data.equi.r2i,'b',x,r2i_tor,'r');
subplot(3,3,5)
plot(x,data.equi.grho2r2,'b',x,grho2r2_tor,'r');
subplot(3,3,6)
plot(x,data.equi.ri,'b',x,ri_tor,'r');
subplot(3,3,7)
plot(x,data.equi.grho,'b',x,grho_tor,'r');
subplot(3,3,8)
plot(x,data.equi.grho2,'b',x,grho2_tor,'r');


% test de z0curdiff
t = data.gene.temps;
%  [psiv,dpsidtv,qv,jmoyv,jresv,eparv,Fv,bpolv,liv,rmx,vpr,grho2r2,r2i,ri,C2,C3] = ...
%     z0curdiff(t,x,data.neo.eta,data.source.totale.j,data.prof.ptot, ...
%     data.equi.spr .* (data.equi.rhomax * ones(size(x))), ...
%     data.equi.e,data.equi.raxe, ...
%     data.geo.b0,data.equi.a(:,end),(data.equi.trh(:,end) + data.equi.trl(:,end))/2,...
%     data.equi.ip,data.gene.vloop,data.equi.li);
eta = griddata(tm,xa,data.neo.eta,tm,xm,'cubic');
jni = griddata(tm,xa,data.source.totale.j,tm,xm,'cubic');
ptot = griddata(tm,xa,data.prof.ptot,tm,xm,'cubic');
Sp  = data.gene.surface;
Vp  = data.gene.volume;
b0   = data.geo.b0;
ip   = data.equi.ip;
vloop = data.gene.vloop;
liin = data.equi.li;
asser = zeros(size(t));
qdds = 0;
peri = data.equi.spr(:,end);
Rsepa = double(data.geo.R);
Zsepa = double(data.geo.Z);
amorti = 0;
drmdt  = data.equi.drhomaxdt;
qeff   = data.equi.q(:,end);
psi_old = griddata(tm,xa,data.equi.psi,tm,xm,'cubic');
phi_old = griddata(tm,xa,data.equi.phi,tm,xm,'cubic');
difcurconv = morphing;
lao_change = 1;
evolution = 0;
edge_flux = data.equi.psi(:,end);
mode_expo_inte = 1;
%
[psi,dpsidt,q,jmoy,jres,epar,F,bpol,li, ...
 rmx,vpr,grho2r2,r2i,ri,C2,C3,ej,ipout,grho,grho2,jeff,spr,phi,dphidx, ...
 difcurconv,phiplasma,indice_inv,poynting] = ...
 z0curdiff(t,x,eta,jni,ptot, ...
          Sp,kx,Raxe,b0,a,dx, ....
          ip,vloop,liin,asser, ...
          Vp,Fin,qdds,peri,Rsepa,Zsepa, ...
	  amorti,drmdt,qeff,psi_old,phi_old,difcurconv,lao_change,evolution,edge_flux,mode_expo_inte);

figure(14);
clf
subplot(2,2,1)
zplotprof(gca,t,xa,data.equi.psi,'color','b');
zplotprof(gca,t,x,psi,'color','r');
subplot(2,2,2)
zplotprof(gca,t,xa,data.equi.jmoy,'color','b');
zplotprof(gca,t,x,jmoy,'color','r');
subplot(2,2,3)
zplotprof(gca,t,xa,data.equi.q,'color','b');
zplotprof(gca,t,x,q,'color','r');
subplot(2,2,4)
zplotprof(gca,t,xa,data.equi.F,'color','b');
zplotprof(gca,t,x,F,'color','r');

figure(15)
subplot(4,1,1)
plot(t,data.equi.li,'b',t,data.gene.li,'g',t,li,'r');
subplot(4,1,2)
plot(t,data.equi.psi(:,end) - data.equi.psi(1,end),'b',t,psi(:,end) -psi(1,end),'r')
subplot(4,1,3)
plot(t,data.equi.psi(:,1) - data.equi.psi(:,end),'b',t,psi(:,1) -psi(:,end),'r')
subplot(4,1,4)
plot(t,medfilt1(data.prof.dpsidt3p(:,end),3),'b',t,dpsidt(:,end),'r')

figure(16);
clf
subplot(2,3,1)
plot(xa',data.equi.vpr','b',x,vpr,'r');
subplot(2,3,2)
plot(xa',(data.equi.rhomax * x)','b',x,rmx,'r');
subplot(2,3,3)
plot(xa',data.equi.c2c','b',x,C2,'r');
subplot(2,3,4)
plot(xa',data.equi.r2i','b',x,r2i,'r');
subplot(2,3,5)
plot(xa',data.equi.grho2r2','b',x,grho2r2,'r');
subplot(2,3,6)
plot(xa',data.equi.ri','b',x,ri,'r');

