% compute betap and li as in MECS
function [time,betap_th_mecs,betap_tot_mecs,li_mecs,ip] = betali_mecs(post)

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

time = post.profil0d.temps;
ip   = interp1(post.zerod.temps,post.zerod.ip,time,'linear','extrap');
peri = interp1(post.zerod.temps,post.zerod.peri,time,'linear','extrap');
wth  =  interp1(post.zerod.temps,post.zerod.wth,time,'linear','extrap');
w    =  interp1(post.zerod.temps,post.zerod.w,time,'linear','extrap');
Rp   =  interp1(post.zerod.temps,post.z0dinput.geo.R,time,'linear','extrap');

Bpa      = phys.mu0 .* ip ./ peri;
Pth      = phys.e .* (post.profil0d.nep  .* post.profil0d.tep + post.profil0d.nip .* post.profil0d.tip);
betap_th_mecs = trapz(post.profil0d.xli,Pth .* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2) ./ ...
                      (Bpa .^ 2 ./ 2 ./ phys.mu0);

betap_tot_mecs = w ./ wth .* betap_th_mecs;

li_efit = trapz(post.profil0d.xli,post.profil0d.bpol .^ 2 .* post.profil0d.vpr,2) ./ ...
          trapz(post.profil0d.xli,post.profil0d.vpr,2) ./ Bpa .^ 2;
 
li_mecs = trapz(post.profil0d.xli,post.profil0d.bpol .^ 2 .* post.profil0d.vpr,2) .* ...
          2 ./ (phys.mu0 .^ 2 .* Rp .* ip .^ 2);
     
if nargout > 0
  return
end
     
h = findobj(0,'type','figure','tag','z0dsc_libeta');
if isempty(h)
       h=figure('tag','z0dsc_libeta');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
 
subplot(3,1,1)
plot(post.zerod.temps,post.zerod.betap,time,betap_th_mecs);
legend('METIS','MECS');
ylabel('\beta_{p,th}');
title(sprintf('METIS : %s@%d / \\beta_p and l_i various definitions', ...
          post.z0dinput.machine,post.z0dinput.shot));

subplot(3,1,2)
plot(post.zerod.temps,post.zerod.betaptot,time,betap_tot_mecs);
legend('METIS','MECS');
ylabel('\beta_{p,total}');

subplot(3,1,3)
plot(post.zerod.temps,post.zerod.li,time,li_mecs);
legend('METIS li(3)','MECS');
ylabel('l_i');
xlabel('time (s)')

joint_axes(gcf,3);
edition2
     