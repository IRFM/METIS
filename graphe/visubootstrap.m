%
% coupled to CRONOS by V. Basiuk
% version v2.0, 17/12/2002
%
% last modifications
%
% * 17/12/2002 : calcul en psi pour Sauter
%

qe    = param.phys.e;
xcron = ones(size(data.gene.temps))*param.gene.x;
Rp    = data.geo.r0*ones(size(xcron(1,:)));

xrho  			= double(data.equi.rhoRZ);
boot.xrho		= xrho;
boot.corx		= pdederive(xrho,xcron,0,2,2,1);
boot.corx(:,1) = boot.corx(:,2);
ap 				= data.equi.rhomax*ones(size(xcron(1,:)));
boot.corx		= (ap.^2).*boot.corx;
%
ia 				= (1./Rp).*xrho;

ia(1) = eps;
Zi    = param.compo.z(1);
xq    = data.equi.q;
Bt    = data.geo.b0;
xTe   = data.prof.te/1e3;
xne   = data.prof.ne;
xni   = squeeze(data.impur.impur(:,:,1));
xTi   = data.prof.ti/1e3;
%
%
%
dTedr   = pdederive(xcron,xTe,0,2,2,1).*boot.corx;
dTidr   = pdederive(xcron,xTi,0,2,2,1).*boot.corx;
dnedr   = pdederive(xcron,xne,0,2,2,1).*boot.corx;
dnidr   = pdederive(xcron,xni,0,2,2,1).*boot.corx;

dTedr   = ap.*data.prof.gte/1e3.*boot.corx;
dnedr   = ap.*data.prof.gne.*boot.corx;%
dnidr   = ap.*data.prof.gni.*boot.corx;
dTidr   = ap.*data.prof.gti/1e3.*boot.corx;
dTedr   = dTedr .* (dTedr<0);
dTidr   = dTidr .* (dTidr<0);
dnedr   = dnedr .* (dnedr<0);
dnidr   = dnidr .* (dnidr<0);
%
dlnTedr = 1./(xTe+eps).*dTedr;
dlnnedr = 1./(xne+eps).*dnedr;
dlnTidr = 1./(xTi+eps).*dTidr;
dlnnidr = 1./(xni+eps).*dnidr;

%
xpe     = xne.*xTe;
xpi     = xni.*xTi;
xp      = xpe + xpi;
dlnpedr = dlnTedr + dlnnedr;
dlnpidr = dlnTidr + dlnnidr;
%
x       = (1.46*ia.^(0.5) + 2.40*ia)./(1 - ia).^(1.5);
x       = data.equi.ftrap./(1-data.equi.ftrap);
%
Dx      = 1.414*Zi + Zi^2 + x*(0.754 + 2.657*Zi + 2*Zi^2) + x.^2*(0.348 + 1.243*Zi + Zi^2);
%
L31     = x.*(0.754 + 2.21*Zi + Zi^2 + x*(0.348 + 1.243*Zi + Zi^2))./Dx;
L32     = -x*(0.884 + 2.074*Zi)./Dx;
L34     = L31;
alpha   = -1.172./(1 + 0.462*Zi);
%
% Sauter Correction for L32 in case of large aspect ratio limit (S.Schultz thesis p 92)
%
x_s     = 1 - (1 - ia).^2./sqrt(1-ia.^2)./(1 + 1.46*sqrt(ia));
L32_s   = -(0.51 + 1.31*Zi)/Zi/(1 + 0.44*Zi)*(x_s - x_s.^4)...
          +(5.95 + 3.57*Zi)/(1 + 2.70*Zi + 0.546*Zi^2)*(x_s.^2 - x_s.^4)...
          -(3.92 + 3.57*Zi)/(1 + 2.70*Zi + 0.546*Zi^2)*(x_s.^3 - x_s.^4);
%
A1                   = dlnpedr + xTi./(Zi*xTe).*dlnpidr;
A2                   = dlnTedr;
A4                   = alpha.*xTi./(Zi*xTe).*dlnTidr;
%
Jbavnorm             = -xpe.*xq./ia.*(L31.*A1 + L32.*A2 + L34.*A4); 
Jbavnorm_s           = -xpe.*xq./ia.*(L31.*A1 + L32_s.*A2 + L34.*A4); 
%Jbavnorm_s           = -xpe.*data.equi.F.*(L31.*A1 + L32_s.*A2 + L34.*A4); 
Jbavnorm_hh          = -xpe.*xq./sqrt(ia).*(0.69*dlnTedr + 2.44*dlnnedr); 
%
Btr                  = Bt * ones(size(xcron(1,:)));
boot.hir             = Jbavnorm./Btr*qe*1e3; 
boot.sau             = Jbavnorm_s./Btr*qe*1e3;
boot.hin1            = Jbavnorm_hh./Btr*qe*1e3;
boot.sau(:,1)        = 0;

boot.hir(:,1)        = 0;
[boot.nustare,boot.nustari,boot.L31,boot.L32,boot.L34,boot.alpha,boot.sau1,boot.eta,boot.spit,boot.sau2]  = sauter(data,param);
boot.sau1            = real(boot.sau1);
boot.isau            = zintsurf(boot.sau,param.gene.x,data.equi.spr,data.equi.rhomax);
boot.isau1           = zintsurf(boot.sau1,param.gene.x,data.equi.spr,data.equi.rhomax);
boot.isau2           = zintsurf(boot.sau2,param.gene.x,data.equi.spr,data.equi.rhomax);
boot.etasau          = zintsurf(1./boot.eta,param.gene.x,data.equi.spr,data.equi.rhomax);
boot.etacron         = zintsurf(1./data.neo.eta,param.gene.x,data.equi.spr,data.equi.rhomax);
boot.etasptz         = zintsurf(1./boot.spit,param.gene.x,data.equi.spr,data.equi.rhomax);
boot.surf            = zintsurf(ones(size(boot.spit)),param.gene.x,data.equi.spr,data.equi.rhomax);
boot.etasau          = boot.surf ./ boot.etasau ;
boot.etacron         = boot.surf ./ boot.etacron;
boot.etasptz         = boot.surf ./ boot.etasptz;
%
% Ion Flux to be added in the fluid limit: Hinton & Hazeltine, Hirshman and Sauter relations use the same term
%
A1_i                 = xTi./(Zi*xTe).*dlnpidr;
Jbavnorm_i           = -xpe.*xq./ia.*(L31.*A1_i + L34.*A4);
flux.ion             = Jbavnorm_i./Btr*qe*1e3; 
%
boot.hin2            = boot.hin1 + flux.ion;
%
%Taguchi's model for anisotropic electron distribution function (Ref: J. Physical Soc. Japan, 65,11 (1996) 3530)
%
h = findobj(0,'type','figure','tag','bootstrap');
if isempty(h)
       h=figure('tag','bootstrap');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1])
subplot(2,2,1)

plot(data.gene.temps,data.gene.iboot,data.gene.temps,boot.isau1,'g--',data.gene.temps,boot.isau2,'g+');
legend('CRONOS','Sauter Z=<Zeff>','Sauter Z=Z1')
ylabel('A');
xlabel('time (s)')
title(['Bootstrap current shot ',param.from.machine,' #',int2str(param.from.shot.num)])

subplot(2,2,2)

plot(data.gene.temps,boot.etacron,data.gene.temps,boot.etasau,'g--',...
     data.gene.temps,boot.etasptz,'ro');
legend('CRONOS','Sauter','SPitzer')
ylabel('Ohm m');
xlabel('time (s)')
title(['integrated Resisitivity'])

subplot(2,2,3)
plotprof('gca',data.gene.temps,param.gene.x,data.neo.jboot,'color','b');
plotprof('gca',data.gene.temps,param.gene.x,boot.sau1,'color','g','linestyle','--');
plotprof('gca',data.gene.temps,param.gene.x,boot.sau2,'color','g','linestyle','+');
title('bootstrap')
ylabel('A/m^2')
legend('CRONOS','Sauter Z=<Zeff>','Sauter Z=Z1')

subplot(2,2,4)
plotprof('gca',data.gene.temps,param.gene.x,data.neo.eta,'color','b');
plotprof('gca',data.gene.temps,param.gene.x,boot.eta,'color','g','linestyle','--');
plotprof('gca',data.gene.temps,param.gene.x,boot.spit,'color','r','linestyle','o');
title('resistivity')
ylabel('Ohm m')
legend('CRONOS','Sauter Z=Zeff','Spitzer')
set(gca,'yscale','log');

