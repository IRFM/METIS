function [yout] = struct_qlkz(flnm,pth,plt)
   
if ~exist('plt')
  plt = 0;
end

path = [pth '/' flnm];

% Read and stores inputs
debug_path = [path, '/debug/'];
yout.kthetarhos = load([debug_path, 'kthetarhos.dat']);
yout.Ane        = load([debug_path, 'Ane.dat']);
yout.Ani        = load([debug_path, 'Ani.dat']);
yout.Ati        = load([debug_path, 'Ati.dat']);
yout.Ate        = load([debug_path, 'Ate.dat']);
yout.qx         = load([debug_path, 'q.dat']);
yout.smag       = load([debug_path, 'smag.dat']);
yout.Tex        = load([debug_path, 'Te.dat']);
yout.Tix        = load([debug_path, 'Ti.dat']);
yout.gammaE     = load([debug_path, 'gammaE.dat']);
yout.Aupar      = load([debug_path, 'Aupar.dat']);
yout.Autor      = load([debug_path, 'Autor.dat']);
yout.Machpar    = load([debug_path, 'Machpar.dat']);
yout.Machtor    = load([debug_path, 'Machtor.dat']);
yout.Bo         = load([debug_path, 'Bo.dat']);
yout.Ro         = load([debug_path, 'Ro.dat']);
%yout.R0         = load([debug_path, 'R0.dat']);
yout.Nex        = load([debug_path, 'ne.dat']);
yout.ninorm     = load([debug_path, 'normni.dat']);
yout.x          = load([debug_path, 'x.dat']);
yout.rho          = load([debug_path, 'rho.dat']);
yout.Zi         = load([debug_path, 'Zi.dat']);
yout.alphax     = load([debug_path, 'alpha.dat']);
yout.Ai         = load([debug_path, 'Ai.dat']);
yout.Rmin       = load([debug_path, 'Rmin.dat']);
yout.scann      = length(yout.x);
yout.phi       = load([debug_path, 'phi.dat']);

x = 1:yout.scann;

Lambe=1-0.078.*log10(yout.Nex.*0.1)+0.15.*log10(yout.Tex);
q_ele  = 1.6022e-19;
me     = 9.1094e-31;
cthe=sqrt(2*yout.Tex*1e3*q_ele./me);
Athe=cthe./(yout.qx.*yout.Ro);
Epsilonx=yout.Rmin.*yout.x./yout.Ro;
ft=2.*(2.*Epsilonx).^(0.5)./pi; %trapped particle fraction

% Reading outputs

yout.eps = Epsilonx;
yout.gam_GB = load([path '/output/gam_GB.dat']);
yout.ome_GB = load([path '/output/ome_GB.dat']);
yout.gam_SI = load([path '/output/gam_SI.dat']);
yout.ome_SI = load([path '/output/ome_SI.dat']);

yout.ntor = load([path '/output/primitive/ntor.dat']);

try
yout.rLecirci = load([path '/output/primitive/rLecirci.dat']);
yout.iLecirci = load([path '/output/primitive/iLecirci.dat']);
yout.rLcirci = load([path '/output/primitive/rLcirci.dat']);
yout.iLcirci = load([path '/output/primitive/iLcirci.dat']);

yout.rLepiegi = load([path '/output/primitive/rLepiegi.dat']);
yout.iLepiegi = load([path '/output/primitive/iLepiegi.dat']);
yout.rLpiegi = load([path '/output/primitive/rLpiegi.dat']);
yout.iLpiegi = load([path '/output/primitive/iLpiegi.dat']);

catch
end

yout.efi_SI = load([path '/output/efi_SI.dat']);
yout.efi_cm = load([path '/output/efi_cm.dat']);
yout.efi_GB = load([path '/output/efi_GB.dat']);
yout.efe_SI = load([path '/output/efe_SI.dat']); 
yout.efe_GB = load([path '/output/efe_GB.dat']);

yout.pfi_SI = load([path '/output/pfi_SI.dat']); 
yout.pfi_GB = load([path '/output/pfi_GB.dat']); 
yout.pfe_SI = load([path '/output/pfe_SI.dat']); 
yout.pfe_GB = load([path '/output/pfe_GB.dat']); 

yout.vce_SI = load([path '/output/vce_SI.dat']); 
yout.vte_SI = load([path '/output/vte_SI.dat']); 

% Contains part of centrifugal effects i.e. exp term
yout.vci_SI = load([path '/output/vci_SI.dat']); 
yout.vti_SI = load([path '/output/vti_SI.dat']); 
yout.vri_SI = load([path '/output/vri_SI.dat']); 

yout.dfi_SI = load([path '/output/dfi_SI.dat']); 
yout.dfe_SI = load([path '/output/dfe_SI.dat']); 
yout.dfe_GB = load([path '/output/dfe_GB.dat']); 
yout.vce_GB = load([path '/output/vce_GB.dat']); 
yout.vte_GB = load([path '/output/vte_GB.dat']); 

try
yout.chie_SI = load([path '/output/chiee_SI.dat']);
yout.chieITG_SI = load([path '/output/chieeITG_SI.dat']); 
yout.chii_SI = load([path '/output/chiei_SI.dat']); 

yout.vene_SI = load([path '/output/vene_SI.dat']); 
yout.vece_SI = load([path '/output/vece_SI.dat']); 
yout.veni_SI = load([path '/output/veni_SI.dat']); 
yout.veci_SI = load([path '/output/veci_SI.dat']); 
catch
yout.chie_SI = NaN;
yout.chieITG_SI = NaN;
yout.chii_SI = NaN;

yout.vene_SI = NaN;
yout.vece_SI = NaN;
yout.veni_SI = NaN;
yout.veci_SI = NaN;
end

% Contains complete set of centrifugal terms
try
cfmat = load([path '/output/cftrans.dat']);
cfmat=reshape(cfmat,[yout.scann 7 size(cfmat,2)]);

yout.vci_tot = squeeze(cfmat(:,3,:)+cfmat(:,7,:));
yout.vti_tot = squeeze(cfmat(:,2,:)+cfmat(:,5,:));
yout.vri_tot = squeeze(cfmat(:,4,:)+cfmat(:,6,:));
yout.dfi_tot = squeeze(cfmat(:,1,:));
catch
yout.vci_tot = NaN;
yout.vti_tot = NaN;
yout.vri_tot = NaN;
yout.dfi_tot = NaN;
end

% Compute effective chi
yout.chieff_i = sum(yout.efi_SI,2)./yout.Ati(:,1)./(yout.ninorm(:,1).*yout.Nex*1e19)./(yout.Tix(:,1)*q_ele*1e3).*yout.Ro;
yout.chieff_e = yout.efe_SI./yout.Ate./(yout.Nex*1e19)./(yout.Tex*q_ele*1e3).*yout.Ro;

%save(['~/matlab/W_transport/fig/' flnm '_qlkz.mat']);

%gam=reshape(yout.gam_GB,[length(yout.Ate) 3 length(yout.kthetarhos)]);
%gam1=squeeze(gam(:,1,:));
	     
if ~(plt)
 return
end
 
figure; 
plot(yout.x,yout.dfe_SI,'-r+','LineWidth',2);
hold on
plot(yout.x,yout.vce_SI+yout.vte_SI,'--b+','LineWidth',2);
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$D_N, C_N$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

% A figure summarizing the main input parameters
figure;
subplot(221)
plot(yout.x,yout.Ati(:,1),'-rs',yout.x,yout.Tix(:,1),'-ro',yout.x,yout.Ate,'-bs',yout.x,yout.Tex,'-bo','LineWidth',2)
l1=legend('-$R\nabla{T_i}/T_i$','$T_i$','-$R\nabla{T_e}/T_e$','$T_e$')
set(l1,'Interpreter','latex')
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')

subplot(222)
plot(yout.x,yout.Ane,'-gs',yout.x,yout.Nex,'-go','LineWidth',2)
hold on
  plot(yout.x,yout.Ani(:,1),'-rs',yout.x,yout.ninorm(:,1).*yout.Nex,'-ro','LineWidth',2)
l1=legend('-$R\nabla{n_e}/n_e$','$n_e$')
set(l1,'Interpreter','latex')
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')

subplot(223)
plot(yout.x,yout.Tex./yout.Tix(:,1),'-cs',yout.x,ft,'-co','LineWidth',2)
l1=legend('$T_e/T_i$','$ft$')
set(l1,'Interpreter','latex')
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')

subplot(224)
plot(yout.x,yout.smag,'-ms',yout.x,yout.qx,'-mo',yout.x,yout.alphax,'m*','LineWidth',2)
l1=legend('s','q','$\alpha$')
set(l1,'Interpreter','latex')
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
  
figure
plot(yout.x,yout.Ani(:,1),'g',yout.x,yout.ninorm(:,1).*yout.Nex,'g--','LineWidth',2)
l1=legend('-$R\nabla{n_W}/n_W$','$n_W$')
set(l1,'Interpreter','latex')
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')

figure; 
plot(yout.x,yout.chie_SI,'-b+',yout.x,yout.chii_SI(:,1),'-r+','LineWidth',2);
hold on
plot(yout.x,yout.chieff_e,'--b+',yout.x,yout.chieff_i(:,1),'--r+','LineWidth',2);
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$\chi_{i,e} (m^{2}/s)$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

figure; 
plot(yout.x,yout.pfe_SI,'-b+','LineWidth',2);
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$\Gamma_e (m^{-2}/s)$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

figure; 
plot(yout.x,yout.dfi_SI(:,1)./yout.chii_SI(:,1),'-b+','LineWidth',2);
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$D_W / \chi_{i}$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

if 0

figure; 
subplot(1,2,1),plot(yout.x,yout.dfi_SI(:,2),'-r+','LineWidth',2);
hold on
plot(yout.x,yout.dfi_tot(:,2),'--r+','LineWidth',2);
%plot(yout.x,-(yout.vci_SI(:,2)+yout.vti_SI(:,2)+yout.vri_SI(:,2))./yout.dfi_SI(:,2),'-r+','LineWidth',2);
%plot(yout.x,-(yout.vci_tot(:,2)+yout.vti_tot(:,2)+yout.vri_tot(:,2))./yout.dfi_tot(:,2),'--r+','LineWidth',2);
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$D_W$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

subplot(1,2,2),plot(yout.x,yout.vci_SI(:,2)+yout.vti_SI(:,2)+yout.vri_SI(:,2),'-r+','LineWidth',2);
hold on
plot(yout.x,yout.vci_tot(:,2)+yout.vti_tot(:,2)+yout.vri_tot(:,2),'--r+','LineWidth',2);
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$V_W$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

figure
plot(yout.eps,yout.dfi_tot(:,2),'--r+','LineWidth',2);
l2=xlabel('$\epsilon$')
set(l2,'Interpreter','latex')
l3=ylabel('$D_W$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

end

% Spectra for checks
figure
for ii=1:size(yout.gam_GB,1)
    subplot(1,2,1),plot(yout.kthetarhos',yout.gam_GB(ii,:),'-')
    hold on
end
xlim([0 2])
box on
set(gca,'FontSize',20)
xlabel('k_\theta \rho_s','FontSize',20)
ylabel('\gamma GB','FontSize',20)

for ii=1:size(yout.gam_GB,1)
  subplot(1,2,2),  plot(yout.kthetarhos',yout.ome_GB(ii,:),'-')
  hold on
end
box on
set(gca,'FontSize',20)
xlabel('k_\theta \rho_s','FontSize',20)
ylabel('\omega GB','FontSize',20)
	 
figure; 
plot(yout.x,yout.efe_SI,'-b+','LineWidth',2)
hold on
plot(yout.x,yout.efi_SI(:,1),'-r+','LineWidth',2)
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$Q_{e,i}$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

return
	 
figure; 
plot(yout.x,(yout.efe_SI./yout.Nex./yout.Tex./1e19/1.6e-16)./yout.Ate.*yout.R0,'-b+','LineWidth',2)
hold on
plot(yout.x,(yout.efi_SI(:,1)./yout.Nex./yout.Tix(:,1)./1e19/1.6e-16)./yout.Ati(:,1).*yout.R0,'-r+','LineWidth',2)
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$\chi_{e,i}$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

figure; 
plot(yout.x,yout.efe_SI,'-b+','LineWidth',2)
hold on
plot(yout.x,yout.efi_SI(:,1),'-r+','LineWidth',2)
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$Q_{e,i}$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

figure; 
plot(yout.x,yout.pfi_SI(:,1),'-r+','LineWidth',2)
hold on
plot(yout.x,yout.pfe_SI,'-b+','LineWidth',2)
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$\gamma_{e,i}$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)

figure; 
plot(yout.x,yout.pfi_SI(:,2),'-r+','LineWidth',2)
hold on
plot(yout.x,yout.pfi_SI(:,1),'-b+','LineWidth',2)
l2=xlabel('$\rho$')
set(l2,'Interpreter','latex')
l3=ylabel('$\gamma_{W,i}$');
set(l3,'Interpreter','latex')
set(gca,'FontSize',20)


