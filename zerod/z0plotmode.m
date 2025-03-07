% constante physique (phys)
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)


% partie gradient critique
ve     = ones(size(post.profil0d.xli));
vt     = ones(size(post.profil0d.temps));
lmin   = (post.profil0d.rmx(:,end) ./ 101) *ve;
lmax   = (2 .* pi .* post.profil0d.qjli(:,end) .* post.profil0d.Raxe(:,end)) *ve;
nepd1  = pdederive(post.profil0d.xli,post.profil0d.nep,0,2,2,1);
nipd1  = pdederive(post.profil0d.xli,post.profil0d.nip,0,2,2,1);
tepd1  = pdederive(post.profil0d.xli,post.profil0d.tep,0,2,2,1);
tipd1  = pdederive(post.profil0d.xli,post.profil0d.tip,0,2,2,1);
ptotd1  = pdederive(post.profil0d.xli,post.profil0d.ptot,0,2,2,1);
rhomax = post.profil0d.rmx(:,end) * ve;
rolne  = min(1./lmin,max(1./lmax,abs(nepd1) ./ max(post.profil0d.nep,1e13) ./ rhomax)) .* post.profil0d.Raxe;
rolni  = min(1./lmin,max(1./lmax,abs(nipd1) ./ max(post.profil0d.nip,1e13) ./ rhomax)) .* post.profil0d.Raxe;
rolte  = min(1./lmin,max(1./lmax,abs(tepd1) ./ max(post.profil0d.tep,13.6) ./ rhomax)) .* post.profil0d.Raxe;
rolti  = min(1./lmin,max(1./lmax,abs(tipd1) ./ max(post.profil0d.tip,13.6) ./ rhomax)) .* post.profil0d.Raxe;
rolptot  = min(1./lmin,max(1./lmax,abs(ptotd1) ./ max(post.profil0d.ptot,13.6) ./ rhomax)) .* post.profil0d.Raxe;
shear  = pdederive(post.profil0d.xli,post.profil0d.qjli,0,2,2,1)  .* (vt * post.profil0d.xli) ./ post.profil0d.qjli;
% plot ITG TEM et ETG
% ref : 
% PRL Hoang ETG 
% EFDA-JET-PR(05)11 E.Aps
% Simulation Synthesis of Formulas Jenko et al 2001 PoP
% Garbet 
Kt             = max(eps,post.profil0d.ftrap ./ (1 - post.profil0d.ftrap));
b2             = post.profil0d.bpol .^ 2 + (post.profil0d.fdia .* post.profil0d.ri) .^ 2;
meff           = interp1(post.zerod.temps,post.zerod.meff,post.profil0d.temps,'pchip','extrap');
rhos           = 4.57e-3 .* sqrt((meff * ve) .* post.profil0d.tep ./ 1e3 ./ b2);
rolte_tem_gn   = 20 ./ 9 ./ Kt + (2/3) .* rolne + Kt ./ 2 .* (1 -  rolne) .^ 2; 
rolte_tem_gn(:,1)  = rolte_tem_gn(:,2);
%rolte_etg_cr   = (1.33 + 1.91 .* abs(shear) ./ post.profil0d.qjli) .* (1 + post.profil0d.zeff .* post.profil0d.tep ./ max(13.6,post.profil0d.tip));
rolti_itg_gn       =(2/3) .* rolni + (20/9) .* post.profil0d.tip ./ max(13.6,post.profil0d.tep) + 0.5 .* (1  - rolni + 0.25 .* rolni .^ 2) .* post.profil0d.tep ./ max(13.6,post.profil0d.tip);
%

%sans gradient de densite
% Article de Casati 
lambda_e        = 0.25 +  2/3 .* abs(shear);
rolte_tem_0gn   = lambda_e ./ post.profil0d.ftrap  .*  (4.4 + 1.1 .*  ...
                   min(1.5,post.profil0d.zeff .* post.profil0d.tep ./ max(13.6,post.profil0d.tip))); 
rolte_tem_0gn(:,1) =  rolte_tem_0gn(:,2);
% article Fourment (sans le 1/2 pour coller aux experiences avec les ITGs)
sf2             = (sign(shear) >= 0) .* (1.1 + 1.4 .* abs(shear) + 1.9 .* abs(shear) ./ post.profil0d.qjli) + ...
                  (sign(shear) < 0) .* (0.9 + 1.6 .* abs(shear)  + 9.9 .*abs(shear) ./ post.profil0d.qjli); 
rolti_itg_0gn  =  (4/3) .* (1 + post.profil0d.tip ./ max(13.6,post.profil0d.zeff .* post.profil0d.tep)) .* sf2;

% limitation
rolte_tem_gn = min(1./lmin,max(1./lmax,rolte_tem_gn));
rolti_itg_gn = min(1./lmin,max(1./lmax,rolti_itg_gn));
rolte_tem_0gn = min(1./lmin,max(1./lmax,rolte_tem_0gn));
rolti_itg_0gn = min(1./lmin,max(1./lmax,rolti_itg_0gn));

% transition (tend vers 0 quand le gradient de densite est negligeable)
trans_e = (1 + tanh(rolne - 2 .* (1 + post.profil0d.zeff .* post.profil0d.tep ./ max(13.6,post.profil0d.tip)))) ./ 2;
rolte_cr = rolte_tem_gn .* trans_e + (1 - trans_e) .* rolte_tem_0gn;
trans_i = (1 + tanh(rolni - 2 .* (1 + post.profil0d.tip ./ max(13.6,post.profil0d.zeff .* post.profil0d.tep)))) ./ 2;
rolti_cr = rolti_itg_gn .* trans_i + (1 - trans_i) .* rolti_itg_0gn;

% J. Cintrin Phd formulae
rolti_itg_citrin = max((1 + post.profil0d.tip ./ max(13.6,post.profil0d.tep)) .* (1.33 + 1.91 .* abs(shear) ./ post.profil0d.qjli) .* ...
                       (1- 3 ./ 2 .* post.profil0d.epsi) .* (1 + 0.3 .* post.profil0d.epsi .*  ...
                        pdederive(post.profil0d.epsi,post.profil0d.kx,2,2,2,1)), 0.8 .* rolni);


% ref: A.G. Peeters et al PoP 12 002505 (2005)
nua = 0.1 .* post.profil0d.nep ./ 1e19 .* post.profil0d.zeff ./ (post.profil0d.tep/1e3) .^ 2;
rolte_peeters = (0.357 .* sqrt(post.profil0d.epsi) + 0.271) ./  sqrt(post.profil0d.epsi) .* ...
                    (4.9 - 1.31 .* rolne + 2.68 .* abs(shear) + log(1 + 20 .* nua));
% P. Mantica 
rolte_GKW_fit = (0.357 .* sqrt(post.profil0d.epsi) + 0.271) ./  sqrt(post.profil0d.epsi) .* ...
                    (-1.5 - 0.1 .*  rolne + 2.5 .* abs(shear)  + 2.5 .*  log(1 + 20 .* nua) + 0.2 .* post.profil0d.tep ./ max(13.6,post.profil0d.tip));

xcronos = post.profil0d.rmx ./ rhomax;

rapte_cr = (rolte - rolte_cr) ./ max(0.1,rolte + rolte_cr);
rapti_cr = (rolti - rolti_cr) ./ max(0.1,rolti + rolti_cr);

rolte_cr = min(1./lmin,max(1./lmax,rolte_cr));
rolti_cr = min(1./lmin,max(1./lmax,rolti_cr));



% indicateur de turbulence type Alfven drift wave
% ref : E. fable et al, Nucl. Fus. 52 (2012) 063017 
switch post.z0dinput.option.gaz
    case   1
        fmass = 1;
        
    case 2
        fmass = 2;
        
    case 3
        fmass = (2 +  3 .* post.z0dinput.cons.iso) ./ (1 + post.z0dinput.cons.iso);
        
    case 4
        fmass = 4;
        
    case 5
        fmass = (1 +  3 .* real(post.z0dinput.cons.iso)) ./ (1 + real(post.z0dinput.cons.iso));
        warning('nHe3onD & nTonD not yet implemented !');
        
    case 11
        fmass = (1 +  11 .* post.z0dinput.cons.iso) ./ (1 + post.z0dinput.cons.iso);
        
end
if length(fmass) > 1
	fmass = interp1(post.zerod.temps,fmass,post.profil0d.temps,'nearest','extrap') * ve;
end
cs  = min(phys.c,sqrt(post.profil0d.tep .* 1.602176462e-19  .* post.profil0d.zeff ./1.6726485e-27 ./ fmass  +  ...
           post.profil0d.tip .* 1.602176462e-19 ./1.6726485e-27 ./ fmass ));
VA       = min(phys.c,2.18e16 .* (post.profil0d.fdia .* post.profil0d.ri) ./ sqrt(post.profil0d.nip .* fmass)); 
lnldei   = 15.2 - 0.5 .* log(post.profil0d.nep./1e20) + log(post.profil0d.tep ./1e3);
tauei    = phys.mp .* fmass ./ phys.me ./ 2 .* 1.09e16 .* (post.profil0d.tep ./1e3) .^ (3/2) ./ post.profil0d.nep ./ lnldei;
mu_hat   = (phys.me ./ phys.mp ./ fmass) .* rolptot .^ 2 .* post.profil0d.qjli .^ 2; 
beta_hat = (rolptot .* cs ./ VA) .^ 2 .* post.profil0d.qjli .^ 2;
C_hat    = 0.51 ./  tauei ./ (cs .* rolptot .* post.profil0d.ri) .* mu_hat;
dwbmi    =  cumtrapz(post.profil0d.xli,C_hat + beta_hat + mu_hat,2);
seuil_dwbm = 1;
indoff   = dwbmi < seuil_dwbm;
dwbmi(indoff) = NaN;
%  % longueur de gradient comparée à la lageur banane 
%  % Debye
%  rho_d  = 2.35e5  .* sqrt(max(13.6,post.profil0d.tep) ./ max(1e13,post.profil0d.nep)./ 1e3);
%  % 1- pour les eletrons
%  % Larmor
%  rho_e   = 1.07e-4 .*  sqrt(max(13.6,post.profil0d.tep) ./ 1e3) ./ sqrt((post.profil0d.fdia .* post.profil0d.ri) .^ 2 + post.profil0d.bpol .^2);
%  % banana 
%  rho_b   = rho_e .* post.profil0d.qjli ./ sqrt(post.profil0d.epsi);	
%  % potatoes
%  rho_pot = (rho_e .^ 2 .* post.profil0d.qjli .^ 2 .* post.profil0d.Raxe) .^ (1/3);
%  % melange au centre
%  rho_b   = rho_b .* (rho_b < post.profil0d.rmx) + rho_pot .* (rho_b >= post.profil0d.rmx);
%  % valeur central
%  rho_b(:,1) = rho_pot(:,1);
%  % securite pour la valeur minimale
%  rho_bd_e = max(rho_b,rho_d);
%  % 1- pour les ions
%  % Larmor
%  rho_i   =  4.57e-3 .* sqrt(meff * ve) .* sqrt(max(13.6,post.profil0d.tip) ./ 1e3) ./ sqrt((post.profil0d.fdia .* post.profil0d.ri) .^ 2 + post.profil0d.bpol .^2);
%  % banana 
%  rho_b   = rho_i .* post.profil0d.qjli ./ sqrt(post.profil0d.epsi);	
%  % potatoes
%  rho_pot = (rho_i .^ 2 .* post.profil0d.qjli .^ 2 .* post.profil0d.Raxe) .^ (1/3);
%  % melange au centre
%  rho_b   = rho_b .* (rho_b < post.profil0d.rmx) + rho_pot .* (rho_b >= post.profil0d.rmx);
%  % valeur central
%  rho_b(:,1) = rho_pot(:,1);
%  % securite pour la valeur minimale
%  rho_bd_i = max(rho_b,rho_d);
%  
%  % valeurs normalisees
%  rolbd_i = post.profil0d.Raxe ./ rho_bd_i;
%  rolbd_e = post.profil0d.Raxe ./ rho_bd_e ;

% zon DDS
indic_dds = zeros(size(post.profil0d.qjli));
indice_inv = interp1(post.zerod.temps,post.zerod.indice_inv,post.profil0d.temps,'nearest',0);
for k= 1:size(indic_dds,1)
	ind = round(indice_inv(k));
        if ind > 0
		indic_dds(k,1:ind) = 1;
	end
end
indic_dds(indic_dds == 0) = NaN;

if mode_plot_transport == 1

      hz =findobj(0,'type','figure','tag','itgtemetg1');
      if isempty(hz)
		hz=figure('tag','itgtemetg1','name','gradient length');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',3,'color',[1 1 1])
      subplot(2,2,1)
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.tip ./ max(13.6,post.profil0d.tep));
      ylabel('Ti/Te');
      set(gca,'xlim',[0.2 0.8]);

      subplot(2,2,3)
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rolte_cr,'color','b');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rolte_peeters,'color','m');
      %zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rolte_GKW_fit,'color','m','linestyle','-.');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rolte,'color','r');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dwbmi,'color','g');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,indic_dds,'color','c');
      set(gca,'YScale','linear');
      ylabel('R/LT_e');
      legend('TEM','TEM Peeters','R/LT_e profile','DW & BM prevailing','Sawteeth zone');
      xlabel('r/a');
      z0loglin(gca);
      set(gca,'xlim',[0.2 0.8]);

      subplot(2,2,4)
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rolti_cr,'color','b');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rolti_itg_citrin,'color','m');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rolti,'color','r');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dwbmi,'color','g');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,indic_dds,'color','c');
      set(gca,'YScale','linear');
      ylabel('R/LT_i');
      legend('ITG','ITG J. Citrin','R/LT_i profile','DW & BM prevailing','Sawteeth zone');
      xlabel('r/a');
      z0loglin(gca);
      set(gca,'xlim',[0.2 0.8]);

      subplot(2,2,2)
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rolne,'color','b');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rolni,'color','r','linestyle','-.');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,0.2 .* min(rolte , rolti),'color','c','linestyle','-');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,0.3 .* max(rolte , rolti),'color','m','linestyle','-');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,2.*ones(size(rolne)),'color','g','linestyle','-');
      legend('R/Ln_e','R/Ln_i','0.2*min(R/LT_e,R/LT_i)','0.3*max(R/LT_e,R/LT_i)','R/Ln = 2');
      set(gca,'YScale','linear');
      z0loglin(gca);
      set(gca,'xlim',[0.2 0.8]);



      hz =findobj(0,'type','figure','tag','itgtemetg2');
      if isempty(hz)
		hz=figure('tag','itgtemetg2','name','instable modes  if > 0 (images)');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',1,'color',[1 1 1])

      subplot(2,2,1)
      imagesc(post.profil0d.xli,post.profil0d.temps,rapte_cr);
      colormap('default');
      colorbar;
      set(gca,'ydir','normal');
      title('R/LT_e - R/LT_e critical (TEM, normalised)');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,2)
      imagesc(post.profil0d.xli,post.profil0d.temps,2 .* trans_e - 1);
      colormap('default');
      colorbar;
      set(gca,'ydir','normal');
      title('Flat (<0) or non Flat (>0) density (for TEM)');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,3)
      imagesc(post.profil0d.xli,post.profil0d.temps,rapti_cr);
      colormap('default');
      colorbar;
      set(gca,'ydir','normal');
      title('R/LT_i - R/LT_i critical (ITG, normalised)');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,4)
      imagesc(post.profil0d.xli,post.profil0d.temps,2 .* trans_i - 1);
      colormap('default')
      colorbar;
      set(gca,'ydir','normal');
      title('Flat (<0)or non Flat (>0) density (for ITG)');
      xlabel('x');
      ylabel('time (s)');


      hz =findobj(0,'type','figure','tag','ADW');
      if isempty(hz)
		hz=figure('tag','ADW','name','Alfven drift wave');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',3,'color',[1 1 1]);

      subplot(2,2,1)
      imagesc(post.profil0d.xli,post.profil0d.temps,log10(mu_hat));
      colormap('default');
      colorbar;
      set(gca,'ydir','normal');
      title('log_{10}(mu_{hat})');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,2)
      imagesc(post.profil0d.xli,post.profil0d.temps,log10(beta_hat));
      colormap('default');
      colorbar;
      set(gca,'ydir','normal');
      title('log_{10}(beta_{hat})');
      xlabel('x');
      ylabel('time (s)');



      subplot(2,2,3)
      imagesc(post.profil0d.xli,post.profil0d.temps,log10(C_hat));
      colormap('default')
      colorbar;
      set(gca,'ydir','normal');
      title('log_{10}(C_{hat})');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,4)
      imagesc(post.profil0d.xli,post.profil0d.temps,sign(cumtrapz(post.profil0d.xli,C_hat + beta_hat + mu_hat,2) - seuil_dwbm));
      colormap('default');
      colorbar;
      set(gca,'ydir','normal');
      title('if < 0, ITG and/or TEM ; if > 0, drift wave and ballooning mode');
      xlabel('x');
      ylabel('time (s)');
end


