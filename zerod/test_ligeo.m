% script pour le test du calcul de li avec varoition de geometrie
zs     = post.zerod;
cons   = post.z0dinput.cons;
geo    = post.z0dinput.geo;
exp0d  = post.z0dinput.exp0d;
R  = geo.R;
a  = geo.a;
K =  geo.K;
Bt = geo.b0;
mu0 = 4e-7 .* pi; 
RRref = zs.RR;
RR    = zs.RR;
taufilt  =  mu0 .* R .* 0.1 ./ RRref ./ 4;
temps    = cons.temps;
qa       = zs.qa;
ip       = zs.ip;
li       = cons.li;
ioh      = zs.iohm;
vloop    = zs.vloop;
tauli = min(1e6,max(1e-6,zs.tauip));
% il faut tenir compte de la proportiond des sources dans la constante de temps
fract = max(0.01,min(1,abs(medfilt1(zs.vloop ./ RR,3)./ip)));
%rm       = a .* sqrt(K);
rm        = a .* (1+1./4.*(-1+K)-3./64.*(-1+K).^2+5./256.*(-1+K).^3-175./16384.*(-1+K).^4);
% derivee temporellle de rm
[rm_void,drmdt]  = zdwdt0(rm, rm ./ taufilt,taufilt,temps,NaN);
% derivee temporelle de dphi/dt
dphidt    = 2 .* pi .* rm .* Bt .* drmdt;
vsc       = drmdt .* (R .* mu0 .* ip) ./ (4 .* pi .^ 2 .* rm .^ 2);
vrc       = dphidt ./ qa;
lext     = 2 .* (log(8 .* R ./ a) - 7/4);
ltot     = lext + li;
[R_void,dRdt]  = zdwdt0(R, R ./ taufilt,taufilt,temps,NaN);
[lext_void,dlextdt]  = zdwdt0(lext, lext ./ taufilt,taufilt,temps,NaN);
[ltot_void,dltotdt]  = zdwdt0(ltot, ltot ./ taufilt,taufilt,temps,NaN);
[wbp_void,dwbpdt]  = zdwdt0(zs.wbp, zs.wbp ./ taufilt,taufilt,temps,NaN);
dwbpdt   = dwbpdt ./ zs.wbp;
   % la correction de courant est deja prise en compte dans jli
   %licor    = - ((4 ./ mu0) .* (vloop - RR .* ioh - 0.* mu0 ./ 2 .* R .* li .* diohdt) ./ ipv - li .* dRdt)  ./ R .* tauli;
licor    = - ((4 ./ mu0) .* (vloop - RR .* ioh) ./ ip  - ltot .* dRdt -  R .* dlextdt )  ./ R;
lict    = li .* dRdt ./ R ;
lice    = R .* dlextdt ./ R;
[l_void,dlidt]  = zdwdt0(li, li ./ taufilt,taufilt,temps,NaN);
[l_void,dli0ddt]  = zdwdt0(zs.li,zs.li ./ taufilt,taufilt,temps,NaN);

figure(51);
clf
x = linspace(-1,1);
subplot(2,2,1)
pp = polyfit(dlidt-dli0ddt,lict,1)
plot(dlidt-dli0ddt,lict,'o',x,polyval(pp,x))
subplot(2,2,2)
pp = polyfit(dlidt-dli0ddt,lice,1)
plot(dlidt-dli0ddt,lice,'o',x,polyval(pp,x))
subplot(2,2,3)
pp = polyfit(dlidt-dli0ddt,drmdt,1)
plot(dlidt-dli0ddt,drmdt,'o',x,polyval(pp,x))
subplot(2,2,4)
pp = polyfit(dlidt-dli0ddt,dwbpdt,1)
plot(dlidt-dli0ddt,dwbpdt,'o',x,polyval(pp,x))

