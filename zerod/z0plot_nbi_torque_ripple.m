% calcul du ripple pour NBI
zs=post.zerod;option=post.z0dinput.option;cons=post.z0dinput.cons;geo=post.z0dinput.geo;profli=post.profil0d;
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

bpol     = (geo.a ./ geo.R) ./ zs.q95 .* geo.b0;

% case of Hydrogen NBI in DT plasma
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
   gas_nbi = -1; 
else
   gas_nbi = option.gaz; 
end


switch gas_nbi
    case -1
        Ab = 1 * ones(size(ftnbi));
    case 11
        Ab = mean(1 .* (1-ftnbi) + 11 .* ftnbi);
    case {3,5}
        Ab = mean(2 .* (1-ftnbi) + 3 .* ftnbi);
    otherwise
        Ab = mean(2 .* (1-ftnbi) + 1 .* ftnbi);
end
fracrip = [0,0.1,0.2,0.4,0.8];
sripnbi = NaN .* (zs.pnbi * fracrip);
for k = 1:length(fracrip)
  pnbi_rip = real(zs.pnbi) .* real(zs.frnbi) .* fracrip(k);   
  einj = real(zs.esup_nbi) .* option.einj ./ max(1,real(zs.pnbi)) ./ max(eps,real(zs.taus_nbi));
  %einj = option.einj;
  inbirip = - phys.e .* (real(cons.pnbi) - real(zs.pnbi) .* real(zs.frnbi)  + pnbi_rip) ./ max(1,einj) .* sqrt(2 .* einj ./ phys.ua ./ real(Ab) ./ (1 - real(zs.mu0_nbi) .^ 2));
  if option.nb_nbi >1
	      pnbi_rip = imag(zs.pnbi) .* imag(zs.frnbi) .* fracrip(k);   
	      einj2 = imag(zs.esup_nbi) .* option.einj2  ./  max(1,imag(zs.pnbi)) ./ max(eps,imag(zs.taus_nbi));
	      %einj2 = option.einj2;
              inbirip = inbirip - phys.e .* (imag(cons.pnbi) - imag(zs.pnbi).* imag(zs.frnbi) + pnbi_rip) ./ max(1,einj2) .* sqrt(2 .* einj2 ./ phys.ua ./ imag(Ab)./ (1 - imag(zs.mu0_nbi) .^ 2));             
  end
  sripnbi(:,k)   = 2 .* pi .* inbirip .* bpol .* geo.R .^ 2;
end

h = findobj(0,'type','figure','tag','z0plot_nbi_torque_ripple');
if isempty(h)
       h=figure('tag','z0plot_nbi_torque_ripple');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
colormap('hot')

subplot(3,1,1);
plot(zs.temps,zs.snbi,'k',zs.temps,sripnbi);
%xlabel('time (s)');
ylabel('torque (N.m)')
legend('NBI col.','no ripple loss (only first orbit)','10% power lost','20% power lost','40% power lost','80% power lost');
title(sprintf('METIS : %s@%d / NBI torque contribution estimation', ...
          post.z0dinput.machine,post.z0dinput.shot));


% recalcul des champs non sauvegardes
[wrad,snbi,sicrh,sfus,sripth,sriplh,sripicrh,sturb,fact,wrot,slh,profli] = z0rot3(zs,option,cons,geo,profli);

% source du au champ electrique // induit. cette source est negligeable mais brise la symetrie co et contre courant (ref: Kim 1993)
% sacling pour taurotmul
if option.taurotmul  == 0
	% il semble que le temps de confinement de la rotation ne depasse jamais le temps de confinement de l'energie
	%option.taurotmul = min(1,max(0.1,zs.nbar ./ zs.nsat));
	%figure(21);clf;plot(zs.temps,option.taurotmul);drawnow
	% reference : P.C. de Vries PPCF 48 (2006) p 1693
	% reference : P.C. de Vries NF 48 (2008) 
	tau_rot = min(zs.taue,zs.tauii);
else
	tau_rot = option.taurotmul .* zs.taue;
end
vion_epar = zs.meff .* (phys.me ./ phys.mp) .*  (zs.iohm + zs.irun) ./ (zs.nem .* zs.sp) ./ (phys.e .* zs.zeff);
s_epar = fact  ./ max(1e-6,tau_rot) .* vion_epar;

%  h = findobj(0,'type','figure','tag','z0plot_nbi_torque_ripple2');
%  if isempty(h)
%         h=figure('tag','z0plot_nbi_torque_ripple2');
%  else
%         figure(h);
%  end
%  clf
%  set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
%  	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
%  colormap('hot')
subplot(3,1,2)
plot(zs.temps,snbi,'r',zs.temps,slh,'b',zs.temps,sturb,'c',zs.temps,s_epar,'m',zs.temps,sfus,'g');
hold on
plot(zs.temps,sripnbi(:,end - 1),'r-.',zs.temps,sripicrh,'k-.',zs.temps,sriplh,'b-.',zs.temps,sripth,'c-.',zs.temps,-zs.sn0fr,'-.g');
%xlabel('time (s)');
ylabel('torque (N.m)')
legend('NBI col.','LH wave','Intrinsic (turbulent + neoclassic)','E_{par}','Fast alpha','NBI ripple (40%)','ICRH ripple','LH ripple','Thermal ripple (JxB)','Neutral friction');
%title(sprintf('METIS : %s@%d / torque contribution estimation', ...
%          post.z0dinput.machine,post.z0dinput.shot));

subplot(3,1,3)
plot(zs.temps,zs.wrad ./ 1e3,profli.temps,profli.omega(:,1) ./ 1e3,profli.temps,profli.omega(:,end) ./ 1e3,profli.temps,max(profli.omega,[],2) ./ 1e3,'-.',profli.temps,min(profli.omega,[],2) ./ 1e3,'-.');
ylabel('Plasma rotation (10^{3} rad/s)')
legend('Volume averaged','central value','edge value','maximum value','minimum value');
xlabel('time (s)');
joint_axes(h,3);


