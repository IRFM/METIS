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

% plot taue (rho)
wx = (3/2) .* cumtrapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.ptot,2);
px_th = phys.e .* (post.profil0d.nep .* post.profil0d.tep +post.profil0d.nip .* post.profil0d.tip);
wx_th = (3/2) .* cumtrapz(post.profil0d.xli,post.profil0d.vpr .* px_th,2);

qx = cumtrapz(post.profil0d.xli,post.profil0d.vpr .* (post.profil0d.source_el + post.profil0d.source_ion),2);
qx(:,1) = qx(:,2);

taue_fit_1 = zeros(size(qx,1),1);
taue_fit_0 = zeros(size(qx,1),1);
for k=1:size(qx,1)
	pp = polyfit(post.profil0d.xli(2:end),wx_th(k,2:end)./qx(k,2:end),2);
        taue_fit_1(k) = polyval(pp,1);
        taue_fit_0(k) = polyval(pp,0);
end
figure
subplot(3,1,1)
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,wx_th,'color','b');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,wx,'color','r');
subplot(3,1,2)
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,wx_th./qx,'color','b');
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,wx./qx,'color','r');
subplot(3,1,3)
plot(post.zerod.temps,post.zerod.taue,'b',post.zerod.temps,post.zerod.taue_alt,'c', ...
post.profil0d.temps,taue_fit_1,'r',post.profil0d.temps,taue_fit_0,'m');
set(gca,'ylim',[0,10])

