function [times,rhofit,te,teptot,tesc,tefit,zeffm,zeff,epar,eta0,jboot,jmoy,data]=coher_te(choc,mode_zeff,plotfl)

if nargin <2
	mode_zeff = 1;
end

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

% acces aux donnees
data=cgcgettrait(choc,'tprof');
% profil de j estime
if fix(choc) < 28500
	[rhophi,rhog,d,vpr,grho2,r2i,ri,grho2r2,psi,phi,fdia, ...
    ptot,jmoy,q,bpol,ftrap,b2,b2i,sp,dpsidrho,tau,sqrtg, ...
    rmoy,r2moy,r2tau2,r3invtau3,r3invtau] = equilTS(data.rmaj,data.amin,data.ip,data.beli,data.wdia,-data.btor);
   rho = rhog ./ (max(rhog')'*ones(1,size(rhog,2))); 
   jmoy = tsplinet(rho,jmoy,ones(size(data.times,1),1) * data.rhofit); 
   q = tsplinet(rho,q,ones(size(data.times,1),1) * data.rhofit); 
   ptot = tsplinet(rho,ptot,ones(size(data.times,1),1) * data.rhofit); 
   ftrap = ftrap .* (ftrap <1) + (ftrap > 1);
   ftrap = tsplinet(rho,ftrap,ones(size(data.times,1),1) * data.rhofit); 
   ri = tsplinet(rho,ri,ones(size(data.times,1),1) * data.rhofit); 
   r2i = tsplinet(rho,r2i,ones(size(data.times,1),1) * data.rhofit); 
   rmoy = tsplinet(rho,rmoy,ones(size(data.times,1),1) * data.rhofit); 
else
	
end

% estimation du zeff
if mode_zeff == 0
	zeffm = data.zeff;
else
	zeffm = mode_zeff .* zeffscaling(data.nbar,data.ploss,data.ip,data.amin,data.rmaj,data.itor,data.zmain(1));
        zeffm = real(zeffm);
        zeffm = max(min(zeffm,6),data.zmain(1));
end
% profil de zeff
expo       = -0.4

if expo ~= 0
   zeff       = data.nefit .^ expo;
   zeff(end)  = 2 .* zeff(end-1) - zeff(end-2);
   g          = trapz(data.rhofit,zeff,2);
   alpha      = (zeffm - data.zmain(1)) ./ (g - 1e21.^expo);
   gamma      = (g .* data.zmain(1)  - zeffm .* 1e21.^expo) ./(g - 1e21.^expo); 
   comp       = ones(1,size(data.rhofit,2));
   zeff       = (alpha * comp) .* zeff + (gamma * comp);
else
   zeff = zeffm *ones(1,size(data.rhofit,2));
end

% champ electrique
vloop = medfilt1(data.vs,5);
comp  = ones(1,size(data.rhofit,2));
epar  = (vloop * comp) ./ (2*pi) .* rmoy .* r2i; 

epar  = epar .* (epar >0);

% rapport d'aspect 
ep    = data.amin ./ data.rmaj;

% estimation de te 
teptot     = ptot ./ (data.nefit + data.nifit) .* 1e3; 

% estimation loi d'echelle
te0   =4 .* 1.23 .* sqrt(data.amin .* data.btor ./ data.rmaj ./ data.nl) .* data.ip  .^ -0.29;
if data.zmain(1)  == 2
   te0   = 1.8 .* zeffm .* sqrt(data.amin .* data.btor ./ data.rmaj ./ data.nl) .* data.ip  .^ -0.29;
else
   te0   = 4.92 .* sqrt(data.amin .* data.btor ./ data.rmaj ./ data.nl) .* data.ip  .^ -0.29;
end
alpha = 0.71 .* (data.rmaj ./ data.amin) .^ 1.02 .* (data.ip ./ data.btor) .^ 0.77;
beta  = 5.82 .* data.amin .^ 2.41 .* data.rmaj .^ -1.1 .*  ...
                data.btor .^ 0.47 .* data.ip .^ -0.79 .* alpha .^ 0.84;
comp  = ones(1,size(data.rhofit,2));
vt    = ones(size(te0,1),1);
tesc  =  (te0 * comp) .* (1 - (vt * data.rhofit) .^ (2 .* (alpha * comp))) .^ (beta *comp) .* 1e3;

% calcul de taue
te   = tesc;  
lnL  = 31.474 + log((data.nefit .^ (-0.5)) .* te);
ind  = find(~isfinite(lnL) | (lnL <10));
if ~isempty(ind)
	lnL(ind) = 10 .* ones(1,length(ind));	
end

taue      = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.me))  .* ...
                ((phys.e .* te) .^ (3/2) ./ data.nefit ./ lnL);
                
comp       = ones(1,size(data.rhofit,2));
nus  = (ep * comp) .^ (-3/2) .* (data.rmaj * comp) .* q ./ taue ./ sqrt(te .* phys.e ./ phys.me);                  
                
% resistivite sans dependance en te
c    = 0.56 .* (3 - zeff) ./ (3 + zeff) ./ zeff;  
fi   = ftrap ./ (1 + (0.58 + 0.2 .* zeff).*nus);
ff   = (1 - fi) .* (1 - c .* fi);               
eta0 = 1.65e-9 .* lnL .* zeff .* (1 + 0.27 .* (zeff-1)) ./ ( 1 + 0.47 .* (zeff - 1)) ./ ff .* (1e-3 ^(-3/2)) ;

% test 
ind  = find(eta0 <1e-9); 
if ~isempty(ind)
	eta0(ind) = 1e-9 * ones(1,length(ind));
end

% courant de bootstrap (estimation ti = te et pe = 2/3 ptot)
c1    = (4 + 2.6 .* ftrap)./ (1 + 1.02 .* sqrt(nus) + 1.07 * nus) ./ (1 + 1.07 .* (ep*comp) .^ (3/2) .* nus); 
jboot = - (data.rmaj * comp) .* ftrap .* (2/3 .* ptot .* phys.e .* 1e3) ./ (2.4 + 5.4 .* ftrap + 2.6 .* ftrap .^2) .* ...
        c1 .* pdederive(data.rhofit,ptot,0,2,2,1) ./ ptot;
jboot(:,end) = zeros(size(jboot,1),1);
        
% estimation de te
te = (eta0 .* (jmoy - jboot) ./ epar) .^(2/3);
te(:,end) = 2 * te(:,end-1) - te(:,end-2); 
% validation du profil
dtedr = pdederive(data.rhofit,real(te),0,2,2,1);
ind = find(dtedr >0);
%te = cumtrapz(data.rhofit,dtedr .* (dtedr <=0),2);
%te = te - te(:,end) * comp;
if ~isempty(ind)
	te(ind) = NaN *ones(1,length(ind));
end
% suppresion des points complexe
ind = find(imag(te));
if ~isempty(ind)
	te(ind) = NaN *ones(1,length(ind));
end

% validation temporelle 
ind = find(data.ploss > max(data.poh));
if ~isempty(ind)
	te(ind,:) = NaN * ones(length(ind),size(te,2));
end

% autre variables
tefit = data.tefit.*1e3;
times = data.times;
rhofit = data.rhofit;

% plot 
if nargin > 2
	figure;
        subplot(2,1,1)
	plotprof(gca,times,rhofit,te,'color','r');
	plotprof(gca,times,rhofit,tefit,'color','b','linestyle','none','marker','o');
	plotprof(gca,times,rhofit,teptot,'color','g');
	plotprof(gca,times,rhofit,tesc,'color','k');
        plotprof(gca,times,data.rhote(:,1:12),data.te(:,1:12).*1e3,'linestyle','none','marker','+','color','c');
        plotprof(gca,times,data.rhote(:,13:end),data.te(:,13:end).*1e3,'linestyle','none','marker','x','color','m');
	xlabel('r/a')
	ylabel('Te (eV)')
	title(sprintf('Verification de Te (choc # %.1f)',choc));
	legend('resistif','fit','Wdia','scalling','Thomson','ECE');
        subplot(2,1,2);
        plot(times,data.ip.*10,times,data.nl,times,data.ploss);
        xlabel('temp (s)')
        legend('ip','nl','ploss');
end	

function zeffsc = zeffscaling(nbar,ptot,ip,a,R0,itor,gaz)
% zeffsc = zeffscaling(nbar,ptot,ip,a,R0,itor,gaz);
%
%
if gaz == 2
  alpha = [-0.6907    0.1467    0.3994   -2.0644    0.9234   -0.1183];

  c0 = 3.0160;

  a0 = c0;
  a1 = alpha(1);
  a2 = alpha(2);
  a3 = alpha(3);
  a4 = alpha(4);
  a5 = alpha(5);
  a6 = alpha(6);

  zeffsc = a0 * nbar.^a1 .* ptot.^a2 .* ip.^a3 .* a.^a4 .* R0.^a5 .* itor.^a6;
end
if gaz == 1

  alpha = [-0.3931    0.1116   3.6174];
  c0 = 185.2499;

  a0 = c0;
  a1 = alpha(1);
  a2 = alpha(2);
  a3 = alpha(3);

  zeffsc = a0 * nbar.^a1 .* ptot.^a2 .* (a./R0).^a3;


end


