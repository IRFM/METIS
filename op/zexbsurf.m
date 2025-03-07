% script d'extraction des composantes du champs sur la DSMF
br    = squeeze(double(data.equi.BR(:,end,:)));
bz    = squeeze(double(data.equi.BZ(:,end,:)));
bphi  = squeeze(double(data.equi.BPHI(:,end,:)));
r     = squeeze(double(data.equi.R(:,end,:)));
z     = squeeze(double(data.equi.Z(:,end,:)));
r0    = data.geo.r0;
theta = unwrap(angle(r - r0*ones(1,size(r,2)) + sqrt(-1) .* z));
temps = data.gene.temps;
bpol  = sqrt(br.^2+bz.^2);
figure
subplot(2,2,1)
plotprof(gca,temps,theta,bpol,'color','r','linestyle','.');
ylabel('Bpol (T)')
subplot(2,2,2)
plotprof(gca,temps,r,z,'linestyle','.');
title('DSMF ')
xlabel('R (m)')
ylabel('Z (m)')
subplot(2,2,3)
plotprof(gca,temps,theta,br,'color','r','linestyle','.');
plotprof(gca,temps,theta,bz,'color','c','linestyle','.');
ylabel('Br et Bz (T)')
xlabel('theta')
subplot(2,2,4)
xlabel('theta')
ylabel('Bphi (T)')
plotprof(gca,temps,theta,bphi,'linestyle','.');

save(sprintf('bsurf%d',fix(param.from.shot.num)),'temps','theta','br','bz','bphi','bpol','r','z','r0');