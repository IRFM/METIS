% calcul des composantes de jpol
br = double(data.equi.BR);
bz = double(data.equi.BZ);
r = double(data.equi.R);
z = double(data.equi.Z);
rho = double(data.equi.rhoRZ);
psi = double(data.equi.psiRZ);
fdia = tsplinet(ones(size(data.equi.F,1),1) * param.gene.x ,data.equi.F,rho ./ (max(rho,[],2) * ones(1,size(rho,2))));
dfdpsi = pdederive(psi,fdia,0,2,2,1);
jr   = repmat(dfdpsi,[1,1,size(br,3)]) .* br ./ param.phys.mu0;
jz   = repmat(dfdpsi,[1,1,size(bz,3)]) .* bz ./ param.phys.mu0;

fjp_r = - repmat(fdia,[1,1,size(br,3)]) ./ r .* jz;
fjp_z =   repmat(fdia,[1,1,size(br,3)]) ./ r .* jr;


%  figure
%  nbt = length(data.gene.temps);
%  n    = ceil(sqrt(nbt));
%  m    = ceil(nbt./n);
%  
%  for k = 1:nbt
%  	subplot(n,m,k)
%  	contour(squeeze(r(k,:,:)),squeeze(z(k,:,:)),sqrt(squeeze(jr(k,:,:)) .^2 + squeeze(jz(k,:,:)) .^ 2));
%  	hold on
%  	quiver(squeeze(r(k,:,:)),squeeze(z(k,:,:)),squeeze(jr(k,:,:)),squeeze(jz(k,:,:)));
%  	set(gca,'ButtonDownFcn','zdataplot(''extrait'');');
%  end
%  
tsscenario
drawnow
[t,void] = ginput(1);
k= find(data.gene.temps >= t,1);


figure
clf
set(gcf,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
contour(squeeze(r(k,:,:)),squeeze(z(k,:,:)),sqrt(squeeze(jr(k,:,:)) .^2 + squeeze(jz(k,:,:)) .^ 2));
hold on
quiver(squeeze(r(k,:,:)),squeeze(z(k,:,:)),squeeze(jr(k,:,:)),squeeze(jz(k,:,:)));
xlabel('R (m)');
ylabel('Z (m)');
title('poloidal current iso norm and vector')

figure
clf
set(gcf,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
contour(squeeze(r(k,:,:)),squeeze(z(k,:,:)),sqrt(squeeze(fjp_r(k,:,:)) .^2 + squeeze(fjp_z(k,:,:)) .^ 2));
hold on
quiver(squeeze(r(k,:,:)),squeeze(z(k,:,:)),squeeze(fjp_r(k,:,:)),squeeze(fjp_z(k,:,:)));
xlabel('R (m)');
ylabel('Z (m)');
title('Jpol x Bphi strength iso norm and vector')




