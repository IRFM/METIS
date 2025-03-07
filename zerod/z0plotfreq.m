% script pedagogique pour les frequences
function z0plotfreq(post)
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
% en Hz; pas de pulsation
% profil de densites
ve = ones(size(post.profil0d.xli));
n1m  = interp1(post.zerod.temps,post.zerod.n1m,post.profil0d.temps,'linear','extrap');
nDm  = interp1(post.zerod.temps,post.zerod.nDm,post.profil0d.temps,'linear','extrap');
nTm  = interp1(post.zerod.temps,post.zerod.nTm,post.profil0d.temps,'linear','extrap');
nhep = post.profil0d.nhep;
nep = post.profil0d.nep;
nHp  = post.profil0d.n1p .* ((max(0,n1m - nDm - nTm) ./ n1m) * ve);
nDp  = post.profil0d.n1p .* ((nDm ./ n1m) * ve);
nTp  = post.profil0d.n1p .* ((nTm ./ n1m) * ve);
nzp  = post.profil0d.nzp;
if isfield(post.profil0d,'nwp')
	nwp  = post.profil0d.nwp;
else
	nwp  = 0 .* nzp;		
end
if ~isfield(post.z0dinput.option,'Sn_fraction')
   Sn_fraction = 0;
else
   Sn_fraction =  post.z0dinput.option.Sn_fraction;
end
switch post.z0dinput.option.gaz
    case 5
        nBp   = zeros(size(nTp));
        nhe3p = nhep; 
        nhep  = post.z0dinput.option.frhe0 .* nep;
   case 11
        nBp   = nTp;
        nHp   = nHp + nTp;
        nTp   = zeros(size(nTp));
    otherwise
        nBp   = zeros(size(nTp));
        nhe3p = zeros(size(nhep));        
end

% precalcul
x    = cat(2,-post.profil0d.xli(end:-1:2),post.profil0d.xli);
t    = post.profil0d.temps;
ux   = (1  - x .^ 2);
ve   = ones(size(x));
vt   = ones(size(t));
qp   = post.profil0d.qjli;
R    = post.profil0d.Raxe;
a    = post.profil0d.Raxe .* post.profil0d.epsi;
rh   = R - a;
rl   = R + a;
rr   = cat(2,rh(:,end:-1:2),rl);
Raxe = cat(2,post.profil0d.Raxe(:,end:-1:2),post.profil0d.Raxe);
xa   = cat(2,-a(:,end:-1:2),a);
xa   = cat(2,-a(:,end:-1:2),a);
qpp  = max(0.5,cat(2,qp(:,end:-1:2),qp));
% 
fdia = cat(2,post.profil0d.fdia(:,end:-1:2),post.profil0d.fdia);
btor = fdia ./ rr;
psi  = cat(2,post.profil0d.psi(:,end:-1:2),post.profil0d.psi);
phi  = cat(2,post.profil0d.phi(:,end:-1:2),post.profil0d.phi);
% computation of Bpol 
% in equatorial plane
bpol = - pdederive(x,psi,2,2,2,1) ./ (pdederive(x,Raxe,2,2,2,1) + a(:,end) * ve) ./ rr;
btot = sqrt(btor .^ 2 + bpol .^ 2);
% other profiles
nep  = cat(2,post.profil0d.nep(:,end:-1:2),post.profil0d.nep);
nip  = cat(2,post.profil0d.nip(:,end:-1:2),post.profil0d.nip);
tep  = cat(2,post.profil0d.tep(:,end:-1:2),post.profil0d.tep);
tip  = cat(2,post.profil0d.tip(:,end:-1:2),post.profil0d.tip);
zeff = cat(2,post.profil0d.zeff(:,end:-1:2),post.profil0d.zeff);
% extend density profiles 
nHp    = cat(2,nHp(:,end:-1:2),nHp);
nDp    = cat(2,nDp(:,end:-1:2),nDp);
nTp    = cat(2,nTp(:,end:-1:2),nTp);
nBp    = cat(2,nBp(:,end:-1:2),nBp);
nhep   = cat(2,nhep(:,end:-1:2),nhep);
nhe3p   = cat(2,nhe3p(:,end:-1:2),nhe3p);
nimpp  = cat(2,nzp(:,end:-1:2),nzp);
nmaxp  = post.z0dinput.option.rimp .* nimpp;
nwp    = cat(2,nwp(:,end:-1:2),nwp);

% base: Wesson
fpe = 8.98.* sqrt(nep);
%
meff = interp1(post.zerod.temps,post.zerod.meff,post.profil0d.temps,'linear','extrap') * ones(size(x));
Zlog = ((double(post.z0dinput.option.gaz == 4) + 1) .^ 2 .* (nep ./ nip) .* zeff) .^ (1/4);
fpi = sqrt(Zlog .^ 2 .* nip .* phys.e .^ 2 ./ meff ./ phys.ua ./ phys.epsi0) ./ 2 ./ pi;
fpi_alt = 0.21 .* sqrt(nip ./ meff);
% fce
fce = phys.e ./ phys.me .* btot ./ (2 .* pi);
% fci for each ion
fci_H   = phys.e ./ phys.mp .* btot ./ (2 .* pi);
fci_D   = phys.e ./ (2 .* phys.ua) .* btot ./ (2 .* pi);
fci_T   = phys.e ./ (3 .* phys.ua) .* btot ./ (2 .* pi);
fci_B   = phys.e ./ (11 .* phys.ua) .* btot ./ (2 .* pi);
fci_He3 = (2 .* phys.e) ./ (3 .* phys.ua) .* btot ./ (2 .* pi);
fci_He4 = (2 .* phys.e) ./ (4 .* phys.ua) .* btot ./ (2 .* pi);
[name_imp,a_imp] = z0ion_name(post.z0dinput.option.zimp);
fci_imp = (post.z0dinput.option.zimp .* phys.e) ./ (a_imp .* phys.ua) .* btot ./ (2 .* pi);
[name_max,a_max] = z0ion_name(post.z0dinput.option.zmax);
fci_max = (post.z0dinput.option.zmax .* phys.e) ./ (a_max .* phys.ua) .* btot ./ (2 .* pi);
zave = z0wavez(tep);
fci_W_ave = (zave .* phys.e) ./ (183.84 .* phys.ua) .* btot ./ (2 .* pi);
zave = z0snavez(tep);
fci_Sn_ave = (zave .* phys.e) ./ (118.71 .* phys.ua) .* btot ./ (2 .* pi);
% lower hybride for main ion
switch post.z0dinput.option.gaz
case 1
  flh_main = sqrt(1 ./ (1 ./ (fci_H .^ 2 + fpi .^ 2) + 1 ./ (fce .* fci_H)));
  fci_main = fci_H;
  n_main   = nHp;
  a_main   = 1;
  z_main   = 1;
case 2
  flh_main = sqrt(1 ./ (1 ./ (fci_D .^ 2 + fpi .^ 2) + 1 ./ (fce .* fci_D))); 
  fci_main = fci_D;
  n_main   = nDp;
  a_main   = 2;
  z_main   = 1;
case 3
  iso = post.z0dinput.cons.iso * ve;
  A_DT     =  mean(mean((2 + 3 .* iso) ./ (1 + iso)));
  fci_DT   =  phys.e ./ (A_DT .* phys.ua) .* btot;
  fci_main = fci_DT;
  flh_main = sqrt(1 ./ (1 ./ (fci_DT .^ 2 + fpi .^ 2) + 1 ./ (fce .* fci_DT))); 
  n_main   = nDp + nTp;
  a_main   = A_DT;
  z_main   = 1;
case 4
  flh_main = sqrt(1 ./ (1 ./ (fci_He4 .^ 2 + fpi .^ 2) + 1 ./ (fce .* fci_He4))); 
  fci_main = fci_He4;
  n_main   = nhep;
  a_main   = 4;
  z_main   = 2;
case 5
  flh_main = sqrt(1 ./ (1 ./ (fci_He3 .^ 2 + fpi .^ 2) + 1 ./ (fce .* fci_He3))); 
  fci_main = fci_He3;
  n_main   = nhep;
  a_main   = 3;
  z_main   = 2;
case 11
  flh_main = sqrt(1 ./ (1 ./ (fci_B .^ 2 + fpi .^ 2) + 1 ./ (fce .* fci_B))); 
  fci_main = fci_B;
  n_main   = nBp;
  a_main   = 11;
  z_main   = 5;
end
% upper hybrid
f_uh = sqrt(fce .^ 2 + fpe .^ 2);
% FCI minoritaire (hyrbid ion-ion dit Bushsbaum)
switch post.z0dinput.option.mino
case 'H'
  z_mino = 1;
  a_mino = 1;
  n_mino = nHp;
  fci_mino = fci_H;
case 'D'
  z_mino = 1;
  a_mino = 2;
  n_mino = nDp;
  fci_mino = fci_D;
case 'T'
  z_mino = 1;
  a_mino = 3;
  n_mino = nTp;
  fci_mino = fci_T;
case {'He','He4' }
  z_mino = 2;
  a_mino = 4;
  n_mino = nhep;
  fci_mino = fci_He4;
case 'He3'
  z_mino = 2;
  a_mino = 4;
  n_mino = nhep .* post.z0dinput.option.cmin;
  fci_mino = fci_He3;
case {'B','B11'}
  z_mino = 5;
  a_mino = 11;
  n_mino = nBp;
  fci_mino = fci_B;
otherwise
  error('unknown minority ion');
end
f_ion_ion_h = sqrt((fci_main .* fci_mino .* (1 + (n_mino .* a_mino) ./ (n_main .* a_main))) ./ ...
              ((a_mino .* z_main) ./ (a_main .* z_mino) + (n_mino .* z_mino) ./ (n_main .* z_main)));


% le plot
fullscreen = get(0,'ScreenSize');
hz =findobj(0,'type','figure','tag','frequences');
if isempty(hz)
  	  hz=figure('tag','frequences','name','Frequencies in plasma');
else
  	  figure(hz);
end
clf
set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
collist = get(gca,'ColorOrder');
zplotprof(gca,post.profil0d.temps,x,fpe,'color','b');
zplotprof(gca,post.profil0d.temps,x,fpi,'color','r');
zplotprof(gca,post.profil0d.temps,x,fce,'color','c');
zplotprof(gca,post.profil0d.temps,x,fci_mino,'color','m');
zplotprof(gca,post.profil0d.temps,x,fci_main,'color','b','linestyle','none','marker','*');
zplotprof(gca,post.profil0d.temps,x,flh_main,'color','k');
zplotprof(gca,post.profil0d.temps,x,f_uh,'color','k','linestyle','-.');
zplotprof(gca,post.profil0d.temps,x,f_ion_ion_h,'color','b','linestyle','--');
zplotprof(gca,post.profil0d.temps,x,fci_imp,'color','m','linestyle','none','marker','x');
zplotprof(gca,post.profil0d.temps,x,fci_max,'color','c','linestyle','none','marker','o');
zplotprof(gca,post.profil0d.temps,x,fci_W_ave,'color','b','linestyle','none','marker','.');
zplotprof(gca,post.profil0d.temps,x,fci_Sn_ave,'color','r','linestyle','none','marker','.');
hold on
plot(x,ones(size(x)).* post.z0dinput.option.freq .* 1e6,'g--');
plot(x,ones(size(x)).* post.z0dinput.option.freq .* 2e6,'g-.');
plot(x,ones(size(x)).* post.z0dinput.option.freq .* 3e6,'g:');
plot(x,ones(size(x)).* post.z0dinput.option.freqlh .* 1e9,'g');
set(gca,'yscale','log');
legend('f_{PE}','f_{PI}','f_{CE}','f_{CI,minority}','f_{CI,main}','f_{LH,main}','f_{UH}','f_{ion-ion}', ...
       sprintf('f_{CI,%s}',name_imp),sprintf('f_{CI,%s}',name_max),'f_{CI,W@<Z>}','f_{CI,Sn@<Z>}', ...
       'f_{ICRH,harm=1}','f_{ICRH,harm=2}','f_{ICRH,harm=3}','f_{LHCD}');
xlabel('r/a signed (<0 @ HFS and >0 @ LFS, equatorial plane)');
ylabel('Hz');
title(sprintf('METIS : %s@%d/Frequencies ', ...
          post.z0dinput.machine,post.z0dinput.shot));
edition2
z0loglin(gca);
drawnow
set(hz,'Position',fullscreen)

