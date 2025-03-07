% script pour tracer la periode des dents de scie, le seui lde declenchement des NTM et la puissance necessaire pour les declenchers.
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

% position des DDS
vt   = ones(size(post.profil0d.temps));
ve   = ones(size(post.profil0d.xli));
indice_inv = interp1(post.zerod.temps,post.zerod.indice_inv,post.profil0d.temps,'nearest','extrap');
dd         = abs(post.profil0d.qjli - 1);
mask_q1    = (dd == (min(dd,[],2) * ve));
indice_q1  = sum((vt * (1:size(post.profil0d.xli,2))) .* mask_q1,2) ./ sum(mask_q1,2);
x_q1  = sum((vt * post.profil0d.xli) .* mask_q1,2) ./ sum(mask_q1,2);

% l'effet de l'augmentation de la periode est supposee lineaire avec P_//
%ref: M.F.F. Nave  et al, Nuclear Fusion 42 (2002) 281-289 figure 6
% scaling pour la periode de dents de scies provient d'un fit experimental
%ref:  W. Park and D.A Monticello, Nuclear Fusion  vol 30  -11 (1990) p 2413-2418
% et de ref: A. Zocco arXiv:1209.6321v3
% ref Zocco et al, PPCF 55 (2013à 074005 (+ ref 2 within)
% le seuil de declenchement des NTM est donne par une loi d'echelle experimentale: 
% ref: I.T. Chapman et al, Nucl. Fusion 50 (2010) 102001 (7pp)  doi:10.1088/0029-5515/50/10/102001

% le courant necessaire pour diminuer la période est estimer à partir de la formule donnees dans:
% A. Merkolv Phd 2006 FOM - never published ? 

% resistive time for ST:
%dd   = abs(vt * (1:21) - indice_inv * ve);
%mask = dd == min(dd,[],2) * ve;
%tau_r_ST = phys.mu0 .* post.profil0d.rmx  .^ 2  ./ 1.22 ./ post.profil0d.eta;
a    = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest','extrap');
tau_r_ST = phys.mu0 .* (a .* x_q1)  .^ 2  ./ 1.22 ./ post.profil0d.eta(:,1);
%tau_r_ST = sum(tau_r_ST .* mask_q1,2) ./ max(1,sum(mask_q1,2));
% Park and Monticello scaling
tau_saw_PM = 9 .* post.profil0d.Raxe(:,1) .^ 2 .* (post.profil0d.tep(:,1) ./ 1000) .^(3/2) ./ post.profil0d.zeff(:,1) ./ 1000;
tau_saw_PM(indice_q1 <= 1) = NaN;
% Zocco 
lne            = 31.3-log(sqrt(post.profil0d.nep(:,1))./post.profil0d.tep(:,1));
fact           = post.profil0d.Raxe(:,1) ./ (post.profil0d.epsi(:,end).^(1.5));
nustare        = fact .* 6.921e-18.* post.profil0d.nep(:,1).*lne.*sum(post.profil0d.qjli .* mask_q1,2) ./ max(1,sum(mask_q1,2))./(post.profil0d.tep(:,1).^2).*post.profil0d.zeff(:,1);
%nustare        = sum(nustare .* mask,2) ./ max(1,sum(mask,2));
tau_eta        = phys.mu0 .* post.profil0d.rmx(:,end)  .^ 2  ./ 1.22 ./ post.profil0d.eta(:,1);
tau_saw_Z      = tau_eta .* (7e-4 + 0.034 .* nustare .^ (2/3) + 0.004 .* nustare - 0.08 .* nustare .^ 2);
tau_saw_Z(indice_q1 <= 1) = NaN;

% stabilisation par les suprathermique
% pour le alpha d'apres experience JET.
ptot_th = phys.e .* post.profil0d.tep .* post.profil0d.nep + ...
	  phys.e .* post.profil0d.tip .* post.profil0d.nip;
wtot_th = (3./2) .* cumtrapz(post.profil0d.xli,ptot_th.* post.profil0d.vpr,2);
wtot_th(:,1) = wtot_th(:,2);
psup_alpha  = post.profil0d.pfus;
esup_fus    = interp1(post.zerod.temps,post.zerod.esup_fus,post.profil0d.temps,'nearest','extrap');
psup_alpha  = psup_alpha ./ max(1, trapz(post.profil0d.xli,psup_alpha .* post.profil0d.vpr,2) * ve) .*  (esup_fus * ve) ./ (3/2);

psup_nbi1   = real(post.profil0d.nbinesource);
esup_nbi1   = interp1(post.zerod.temps,real(post.zerod.esup_nbi),post.profil0d.temps,'nearest','extrap');
psup_nbi1   = psup_nbi1 ./ max(1, trapz(post.profil0d.xli,psup_nbi1 .* post.profil0d.vpr,2) * ve) .*  (esup_nbi1 * ve) ./ (3/2);

psup_nbi2   = imag(post.profil0d.nbinesource);
esup_nbi2   = interp1(post.zerod.temps,imag(post.zerod.esup_nbi),post.profil0d.temps,'nearest','extrap');
psup_nbi2   = psup_nbi2 ./ max(1, trapz(post.profil0d.xli,psup_nbi2 .* post.profil0d.vpr,2) * ve) .*  (esup_nbi2 * ve) ./ (3/2);

psup_icrh  = post.profil0d.picrh;
esup_icrh  = interp1(post.zerod.temps,post.zerod.esup_icrh,post.profil0d.temps,'nearest','extrap');
psup_icrh  = psup_icrh ./ max(1, trapz(post.profil0d.xli,psup_icrh .* post.profil0d.vpr,2) * ve) .*  (esup_icrh * ve) ./ (3/2);

psup_lh  = post.profil0d.plh;
esup_lh  = interp1(post.zerod.temps,post.zerod.esup_lh,post.profil0d.temps,'nearest','extrap');
psup_lh  = psup_lh ./ max(1, trapz(post.profil0d.xli,psup_lh .* post.profil0d.vpr,2) * ve) .*  (esup_lh * ve) ./ (3/2);
	  
ptot    = ptot_th + psup_alpha + psup_nbi1 + psup_nbi2 + psup_icrh  + psup_lh;
wtot    = (3./2) .* cumtrapz(post.profil0d.xli,ptot.* post.profil0d.vpr,2);
wtot(:,1) = wtot(:,2);

fact_press = sum(wtot./ wtot_th .* mask_q1,2) ./ max(1,sum(mask_q1,2));

tau_saw_PM_fast = tau_saw_PM .* fact_press;
tau_saw_Z_fast  = tau_saw_Z .* fact_press;

% calcul de betan NTM onset
% etalonage sur les donnees du papier; definitions pas  precise !
bpol    = phys.mu0 ./ 2 ./ pi .* 15e6 ./ 2;
rho_i_r1 = 4.57e-3 .* sqrt(2.5) .* sqrt(23) ./ bpol ./ (2 .* [0.33 0.45 0.33 0.45]);
betaN_NTM_iter_50 = 2.614 .* (0.0446) .^ -0.4084 .* rho_i_r1 .^ 0.5721 .* ([1.3 1.3 1.7 1.7]) .^ 0.4204 .* (1e20 ./ 1e19) .^ 0.4948;
%betaN_NTM_iter_20 = 2.614 .* (0.0178) .^ -0.4084 .* rho_i_r1 .^ 0.5721 .* ([1.3 1.3 1.7 1.7]) .^ 0.4204 .* (1e20 ./ 1e19) .^ 0.4948;
fact_paper = mean(2.09 ./ betaN_NTM_iter_50)

meff    = interp1(post.zerod.temps,post.zerod.meff,post.profil0d.temps,'nearest','extrap');
bpol    = phys.mu0 ./ 2 ./ pi .* interp1(post.zerod.temps,post.zerod.ip,post.profil0d.temps,'nearest','extrap') ./ a;
tim     = interp1(post.zerod.temps,post.zerod.tite .*post.zerod.tem,post.profil0d.temps,'nearest','extrap');
warning off
rho_i_r1 = 4.57e-3 .* sqrt(meff) .* sqrt(post.profil0d.tip(:,1) ./ 1e3) ./ bpol ./ (a .* x_q1);
%rho_i_r1 = 4.57e-3 .* sqrt(meff) .* sqrt(tim ./ 1e3) ./ bpol ./ (a .* x_q1);

%rho_i_r1(:,1) = rho_i_r1(:,2);
%rho_i_r1 = sum(rho_i_r1 .* mask,2) ./ max(1,sum(mask,2));
nbar    = interp1(post.zerod.temps,post.zerod.nbar,post.profil0d.temps,'nearest','extrap');
sext    = interp1(post.zerod.temps,post.zerod.sext,post.profil0d.temps,'nearest','extrap');
pin     = interp1(post.zerod.temps,post.zerod.pin,post.profil0d.temps,'nearest','extrap');
%pl2h    = interp1(post.zerod.temps,post.zerod.plossl2h,post.profil0d.temps,'nearest','extrap');
pl2h    = 0.0488e6 .*  (nbar ./ 1e20) .^ 0.717 .* (post.profil0d.fdia(:,end) .* post.profil0d.ri(:,end)) .^ 0.8 .*  sext .^ 0.941; 
betaN_NTM_onset_PM_fast = fact_paper .* 2.614 .* (tau_saw_PM_fast ./ tau_r_ST) .^ -0.4084 .* rho_i_r1 .^ 0.5721 .* (pin ./ pl2h) .^ 0.4204 .* (nbar ./ 1e19) .^ 0.4948;
betaN_NTM_onset_PM      = fact_paper .* 2.614 .* (tau_saw_PM ./ tau_r_ST)      .^ -0.4084 .* rho_i_r1 .^ 0.5721 .* (pin ./ pl2h) .^ 0.4204 .* (nbar ./ 1e19) .^ 0.4948;
betaN_NTM_onset_Z_fast  = fact_paper .* 2.614 .* (tau_saw_Z_fast ./ tau_r_ST)  .^ -0.4084 .* rho_i_r1 .^ 0.5721 .* (pin ./ pl2h) .^ 0.4204 .* (nbar ./ 1e19) .^ 0.4948;
betaN_NTM_onset_Z       = fact_paper .* 2.614 .* (tau_saw_Z ./ tau_r_ST)       .^ -0.4084 .* rho_i_r1 .^ 0.5721 .* (pin ./ pl2h) .^ 0.4204 .* (nbar ./ 1e19) .^ 0.4948;
betaN_NTM_onset_PM_fast(betaN_NTM_onset_PM_fast == 0) = NaN;
betaN_NTM_onset_PM(betaN_NTM_onset_PM == 0) = NaN;
betaN_NTM_onset_Z_fast(betaN_NTM_onset_Z_fast == 0) = NaN;
betaN_NTM_onset_Z(betaN_NTM_onset_Z == 0) = NaN;
warning on

% le critere doit etre complete avec ta_saw > tau_E
% doit être en mode H
betan    = interp1(post.zerod.temps,post.zerod.betan,post.profil0d.temps,'nearest','extrap') .* 100;
modeh    = interp1(post.zerod.temps,post.zerod.modeh,post.profil0d.temps,'nearest','extrap');
taue     = interp1(post.zerod.temps,min(post.zerod.taue,post.zerod.taue_alt),post.profil0d.temps,'nearest','extrap');
NTM_onoff_PM_fast = (betan > betaN_NTM_onset_PM_fast) & (modeh > 0) & (tau_saw_PM_fast >= taue);
NTM_onoff_PM_fast(~isfinite(NTM_onoff_PM_fast)) = 0;
NTM_onoff_PM = (betan > betaN_NTM_onset_PM) & (modeh > 0) & (tau_saw_PM >= taue);
NTM_onoff_PM(~isfinite(NTM_onoff_PM)) = 0;
NTM_onoff_Z_fast = (betan > betaN_NTM_onset_Z_fast) & (modeh > 0) & (tau_saw_Z_fast >= taue);
NTM_onoff_Z_fast(~isfinite(NTM_onoff_Z_fast)) = 0;
NTM_onoff_Z = (betan > betaN_NTM_onset_Z) & (modeh > 0) & (tau_saw_Z >= taue);
NTM_onoff_Z(~isfinite(NTM_onoff_Z)) = 0;
% resume
NTM_onoff = double(NTM_onoff_PM_fast | NTM_onoff_PM |  NTM_onoff_Z_fast | NTM_onoff_Z);

% calcul du courant necessaire pour le controle:
% efficacite de generartion de courant a q=1
% utilise les reglages de METIS
% efficacite de ECCD (communication privee G. Giruzzi et
% G. Giruzzi, Nucl. Fus. 27, (1987) )
% modification de la dependance en Zeff:
% ref : Y.T. Lin-Liu, GA-A24257
aece  = post.z0dinput.option.angle_ece./180*pi;
%xece = max(0,(indice_inv - 1) ./ 20);
xece = x_q1;
%indece = fix(interp1(post.profil0d.xli,1:length(post.profil0d.xli),abs(xece),'nearest','extrap'));
indece = indice_q1;
tece   = diag(post.profil0d.tep(:,indece));
nece   = diag(post.profil0d.nep(:,indece));
zeff   = diag(post.profil0d.zeff(:,indece));
R    = interp1(post.zerod.temps,post.z0dinput.geo.R,post.profil0d.temps,'nearest','extrap');
%a    = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest','extrap');
mut  = sqrt(a .* abs(xece) .* (1 + cos(aece)) ./ (R + a .* abs(xece) .* cos(aece)));
etaece = 1e20 ./ (1 + 100 ./ (tece./1000)) .* (1 - (1 + (5+zeff) ./ 3 ./ (1+zeff)) .* ...
	(sqrt(2) .* mut) .^ ((5+zeff)./ (1+zeff))) .* 6 ./  ...
	(1 + 4 .* (1 - sqrt(2 .* a .* abs(xece) ./  (R + a .* abs(xece)))) + zeff);
ieccd_pecrh   =  (post.z0dinput.option.eccdmul .* etaece) ./ nece ./ R .* (xece >= 0);

% largeur du depot ECCD
vth2     = 2 .* 1.602176462e-19 .* tece ./9.10938188e-31;
deccd    = sqrt(vth2 .^ 2 ./ 4 ./ 2.99792458e8 .^ 4  + ...
	    min(1,abs(ieccd_pecrh)) .* vth2 ./ 2.99792458e8 .^ 2);		   

% Merkulov formula for current need for the control
iofx       = cumtrapz(post.profil0d.xli,post.profil0d.spr .*post.profil0d.jeff,2);
ir1        = diag(iofx(:,indece));
I_eccd_min = 2 .* (deccd ./ a ./ max(eps,xece)) .^ 2 .* ir1; 
P_eccd_min = I_eccd_min ./ ieccd_pecrh;
I_eccd_min_dim = max(I_eccd_min .* NTM_onoff);
P_eccd_min_dim = max(P_eccd_min .* NTM_onoff);

% NTM stabilisation 
%ref (and ref within): Laura Urso Phd, 22 April 2009, chapiter 2.8.1
% j_eccd > [0.6 1.8] * j_bs  + integration on EC deposition
% J'ai cale le résultat pour retrouve 20 MW pour ITER
% fisrt evaluation for 3/2 
dd   = abs(post.profil0d.qjli - 3/2);
mask = dd == min(dd,[],2) * ve;
xece = sum((vt * post.profil0d.xli) .* mask,2) ./ max(1,sum(mask,2));
indece = fix(interp1(post.profil0d.xli,1:length(post.profil0d.xli),abs(xece),'nearest','extrap'));
tece   = diag(post.profil0d.tep(:,indece));
nece   = diag(post.profil0d.nep(:,indece));
zeff   = diag(post.profil0d.zeff(:,indece));
%R    = interp1(post.zerod.temps,post.z0dinput.geo.R,post.profil0d.temps,'nearest','extrap');
%a    = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest','extrap');
mut  = sqrt(a .* abs(xece) .* (1 + cos(aece)) ./ (R + a .* abs(xece) .* cos(aece)));
etaece = 1e20 ./ (1 + 100 ./ (tece./1000)) .* (1 - (1 + (5+zeff) ./ 3 ./ (1+zeff)) .* ...
	(sqrt(2) .* mut) .^ ((5+zeff)./ (1+zeff))) .* 6 ./  ...
	(1 + 4 .* (1 - sqrt(2 .* a .* abs(xece) ./  (R + a .* abs(xece)))) + zeff);
ieccd_pecrh   =  (post.z0dinput.option.eccdmul .* etaece) ./ nece ./ R .* (xece >= 0);

vth2     = 2 .* 1.602176462e-19 .* tece ./9.10938188e-31;
deccd    = sqrt(vth2 .^ 2 ./ 4 ./ 2.99792458e8 .^ 4  + ...
	    min(1,abs(ieccd_pecrh)) .* vth2 ./ 2.99792458e8 .^ 2);		   
shape_eccd    = exp(-(vt*post.profil0d.xli -(abs(xece) * ve)).^ 2 ./2 ./ (deccd * ve) .^ 2);
% calcul du courant
ieccd_3o2  = 0.25 .*  trapz(post.profil0d.xli,post.profil0d.spr .* post.profil0d.jboot .* shape_eccd,2);
peccd_3o2  = ieccd_3o2 ./ ieccd_pecrh;

% second evaluation for 2 
dd   = abs(post.profil0d.qjli - 2);
mask = dd == min(dd,[],2) * ve;
xece = sum((vt * post.profil0d.xli) .* mask,2) ./ max(1,sum(mask,2));
indece = fix(interp1(post.profil0d.xli,1:length(post.profil0d.xli),abs(xece),'nearest','extrap'));
tece   = diag(post.profil0d.tep(:,indece));
nece   = diag(post.profil0d.nep(:,indece));
zeff   = diag(post.profil0d.zeff(:,indece));
%R    = interp1(post.zerod.temps,post.z0dinput.geo.R,post.profil0d.temps,'nearest','extrap');
%a    = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest','extrap');
mut  = sqrt(a .* abs(xece) .* (1 + cos(aece)) ./ (R + a .* abs(xece) .* cos(aece)));
etaece = 1e20 ./ (1 + 100 ./ (tece./1000)) .* (1 - (1 + (5+zeff) ./ 3 ./ (1+zeff)) .* ...
	(sqrt(2) .* mut) .^ ((5+zeff)./ (1+zeff))) .* 6 ./  ...
	(1 + 4 .* (1 - sqrt(2 .* a .* abs(xece) ./  (R + a .* abs(xece)))) + zeff);
ieccd_pecrh   =  (post.z0dinput.option.eccdmul .* etaece) ./ nece ./ R .* (xece >= 0);

vth2     = 2 .* 1.602176462e-19 .* tece ./9.10938188e-31;
deccd    = sqrt(vth2 .^ 2 ./ 4 ./ 2.99792458e8 .^ 4  + ...
	    min(1,abs(ieccd_pecrh)) .* vth2 ./ 2.99792458e8 .^ 2);		   
shape_eccd    = exp(-(vt*post.profil0d.xli -(abs(xece) * ve)).^ 2 ./2 ./ (deccd * ve) .^ 2);
% calcul du courant
ieccd_2o1  = 0.125 .* trapz(post.profil0d.xli,post.profil0d.spr .* post.profil0d.jboot .* shape_eccd,2);
peccd_2o1  = ieccd_2o1 ./ ieccd_pecrh;

% NTM control
I_eccd_NTM_control = max(ieccd_3o2,ieccd_2o1);
P_eccd_NTM_control = max(peccd_3o2,peccd_2o1);
%I_eccd_NTM_dim = max(I_eccd_NTM_control .* NTM_onoff);
%I_eccd_NTM_dim = max(I_eccd_NTM_control .* double((modeh > 0) & (tau_saw_PM_fast >= taue)));
I_eccd_NTM_dim = max(I_eccd_NTM_control .* double((modeh > 0)));
%P_eccd_NTM_dim = max(P_eccd_NTM_control .* NTM_onoff);
%P_eccd_NTM_dim = max(P_eccd_NTM_control .* double((modeh > 0) & (tau_saw_PM_fast >= taue)));
P_eccd_NTM_dim = max(P_eccd_NTM_control .* double((modeh > 0)));


% figure
fullscreen = get(0,'ScreenSize');
h = findobj(0,'type','figure','tag','z0plot_sawtooth_info');
if isempty(h)
       h=figure('tag','z0plot_sawtooth_info');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1],'Position',fullscreen);
colormap('hot')
k = 4;
l = 0;

l = l+1;
subplot(k,1,l)
plot(post.profil0d.temps,NTM_onoff_PM .* 0.2,'xm',post.profil0d.temps,NTM_onoff_PM_fast.* 0.4,'or', ...
     post.profil0d.temps,NTM_onoff_Z .* 0.6,'+c',post.profil0d.temps,NTM_onoff_Z_fast .* 0.8,'*b');
ylabel('NTM on/off ');
legend('Park-Monticello (thermal)','Park-Monticello (fast)','Zocco-Connor (thermal)','Zocco-Connor (fast)','location','best')
title(sprintf('METIS : %s@%d/ NTM and Sawtooth',post.z0dinput.machine,post.z0dinput.shot));
set(gca,'ylim',[0 1]);

l = l+1;
subplot(k,1,l)
plot(post.profil0d.temps,tau_saw_PM,'m',post.profil0d.temps,tau_saw_PM_fast,'r', ...
     post.profil0d.temps,tau_saw_Z,'c',post.profil0d.temps,tau_saw_Z_fast,'b', ...
     post.zerod.temps,post.zerod.taue,'k',post.zerod.temps,post.zerod.taue_alt,'k:');
ylabel('\tau_{saw}  & \tau_E (s)')
legend('Park-Monticello (thermal)','Park-Monticello (fast)','Zocco-Connor (thermal)','Zocco-Connor (fast)' , ...
       'Energy confinement time','location','best')

l = l+1;
subplot(k,1,l)
plot(post.profil0d.temps,betaN_NTM_onset_PM,'m',post.profil0d.temps,betaN_NTM_onset_PM_fast,'r', ...
     post.profil0d.temps,betaN_NTM_onset_Z,'c',post.profil0d.temps,betaN_NTM_onset_Z_fast,'b', ...
     post.zerod.temps,post.zerod.betan .* 100,'k');
ylabel('\beta_{N, NTM onset} & \beta_{N} ');
legend('Park-Monticello (thermal)','Park-Monticello (fast)','Zocco-Connor (thermal)','Zocco-Connor (fast)' , ...
       '\beta_{N} plasma','location','best')

l = l+1;
subplot(k,1,l)
semilogy(post.profil0d.temps,I_eccd_min ./ 1e6,'b',post.profil0d.temps,P_eccd_min ./ 1e6,'r', ...
     post.profil0d.temps,I_eccd_NTM_control ./ 1e6,'c',post.profil0d.temps,P_eccd_NTM_control ./ 1e6,'m');
ylabel('I_{ECCD} (MA) & P_{ECCD} (MW)');
legend(sprintf('I_{ECCD,min}; ST control = %6.3g (MA)',I_eccd_min_dim ./ 1e6),sprintf('P_{ECCD,min}; ST control = %6.3g (MW)',P_eccd_min_dim ./ 1e6), ...
       sprintf('I_{ECCD}; NTM control = %6.3g (MA)',I_eccd_NTM_dim ./ 1e6),sprintf('P_{ECCD}; NTM control = %6.3g (MW)',P_eccd_NTM_dim ./ 1e6), ...
       'location','best');

xlabel('time (s)');
z0loglin(gca);

joint_axes(h,k);
edition2
set(h,'Position',fullscreen);
drawnow