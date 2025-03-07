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


% script pour tracer la fonction  <s/q> * qa
vt    = ones(size(post.profil0d.qjli,1),1);
ve    = ones(size(post.profil0d.xli));
s     = pdederive(post.profil0d.xli,post.profil0d.qjli,2,2,2,1) ./ post.profil0d.qjli .* (vt * post.profil0d.xli);
Fs    = trapz(post.profil0d.xli,s .* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2); 
soq   = s ./ post.profil0d.qjli;
Fsoq  = post.profil0d.qjli(:,end) .* trapz(post.profil0d.xli,soq .* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2); 
iotaq = (post.profil0d.qjli(:,end) * ve) ./ post.profil0d.qjli;
Fiota = trapz(post.profil0d.xli,iotaq .* post.profil0d.vpr,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2); 


rm      = interp1(post.zerod.temps,post.zerod.rm,post.profil0d.temps,'linear','extrap');
meff    = interp1(post.zerod.temps,post.zerod.meff,post.profil0d.temps,'linear','extrap') * ones(size(post.profil0d.xli));
zeff    = interp1(post.zerod.temps,post.zerod.zeff,post.profil0d.temps,'linear','extrap');
nem     = interp1(post.zerod.temps,post.zerod.nem,post.profil0d.temps,'linear','extrap');
nim     = interp1(post.zerod.temps,post.zerod.nim,post.profil0d.temps,'linear','extrap');
rho_e   = 1.07e-4 .*  sqrt(max(13.6,post.profil0d.tep) ./ 1e3) ./ sqrt((post.profil0d.fdia .* post.profil0d.ri) .^ 2 + post.profil0d.bpol .^2);
rho_b   = rho_e .* post.profil0d.qjli ./ sqrt(post.profil0d.epsi);	
rho_pot = (rho_e .^ 2 .* post.profil0d.qjli .^ 2 .* post.profil0d.Raxe) .^ (1/3);
rho_b   = rho_b .* (rho_b < post.profil0d.rmx) + rho_pot .* (rho_b >= post.profil0d.rmx);
rho_b(:,1) = rho_pot(:,1);
% debyes
rho_d  = 2.35e5  .* sqrt(max(13.6,post.profil0d.tep) ./ max(1e13,post.profil0d.nep)./ 1e3);
rho_bd = max(rho_b,rho_d);
rho_bd = post.profil0d.ftrap .* rho_bd + (1 - post.profil0d.ftrap) .* rho_e;
%nblarmor_l =  rm ./ max(eps,trapz(post.profil0d.xli,rho_bd,2));  
nblarmor_e =  0.95 .* rm ./ max(eps,trapz(post.profil0d.xli(1:end-1),rho_bd(:,1:end-1),2)); 


rho_i   = 4.57e-3  .* sqrt(meff) .*  sqrt(max(13.6,post.profil0d.tip) ./ 1e3) ./ sqrt((post.profil0d.fdia .* post.profil0d.ri) .^ 2 + post.profil0d.bpol .^2);
rho_b   = rho_i .* post.profil0d.qjli ./ sqrt(post.profil0d.epsi);	
rho_pot = (rho_i .^ 2 .* post.profil0d.qjli .^ 2 .* post.profil0d.Raxe) .^ (1/3);
rho_b   = rho_b .* (rho_b < post.profil0d.rmx) + rho_pot .* (rho_b >= post.profil0d.rmx);
rho_b(:,1) = rho_pot(:,1);
% debyes
rho_d  = 2.35e5  .* sqrt(max(13.6,post.profil0d.tep) ./ max(1e13,post.profil0d.nep)./ 1e3);
rho_bd = max(rho_b,rho_d);
rho_bd = post.profil0d.ftrap .* rho_bd + (1 - post.profil0d.ftrap) .* rho_i;
%nblarmor_i =  rm ./ max(eps,trapz(post.profil0d.xli,rho_bd,2));  
nblarmor_i =  0.95 .* rm ./ max(eps,trapz(post.profil0d.xli(1:end-1),rho_bd(:,1:end-1),2)); 

nblarmor_tot =  nem ./ (nem + nim) .* nblarmor_e + sqrt(phys.mp ./ phys.me) .* nim ./ (nem + nim) .* nblarmor_i;

if mode_plot_transport == 1

	hz =findobj(0,'type','figure','tag','soq');
	if isempty(hz)
		  hz=figure('tag','soq','name','<s/q>  confinement factor');
	else
		  figure(hz);
	end
	clf
	set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
		'defaultlinelinewidth',1,'color',[1 1 1])

	subplot(3,1,1)
	plot(post.profil0d.temps,Fsoq);
	legend('<s/q> * q_{edge} confinement factor','location','best');
	ylabel('<s/q> * q_{edge}');
	xlabel('time (s)');
	title(sprintf('METIS : %s@%d/confinement factor', ...
		    post.z0dinput.machine,post.z0dinput.shot));

	subplot(3,1,2)
	plot(post.profil0d.temps,Fiota);
	legend('Rotational transformation factor','location','best');
	ylabel('<q_{edge}/q>');
	%xlabel('time (s)');

	subplot(3,1,3)
	plot(post.profil0d.temps,Fs);
	legend('Volume averaged magnetic shear','location','best');
	ylabel('<s>');
	xlabel('time (s)');

	joint_axes(hz,3);
	edition2

	hz =findobj(0,'type','figure','tag','insulation');
	if isempty(hz)
		  hz=figure('tag','insulation','name','insulation factor');
	else
		  figure(hz);
	end
	clf
	set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
		'defaultlinelinewidth',1,'color',[1 1 1])


	subplot(3,1,1)
	%semilogy(post.profil0d.temps,nblarmor_tot,post.profil0d.temps,nblarmor_e,post.profil0d.temps,nblarmor_i);
	%legend('total weighted','electrons','ions','location','best');
	semilogy(post.profil0d.temps,nblarmor_e,post.profil0d.temps,nblarmor_i);
	legend('electrons','ions','location','best');
	hy = ylabel('number Larmor radius or orbit width in the plasma');
	title(sprintf('METIS : %s@%d/confinement factor', ...
		    post.z0dinput.machine,post.z0dinput.shot));
	z0loglin(gca);

	subplot(3,1,2)
	plot(post.profil0d.temps,(post.profil0d.tep(:,1) - post.profil0d.tep(:,end -1)) ./ nblarmor_e);
	legend('core electron insulation factor','location','best');
	ylabel('eV / layer');

	subplot(3,1,3)
	plot(post.profil0d.temps,(post.profil0d.tip(:,1) - post.profil0d.tip(:,end - 1)) ./ nblarmor_i);
	legend('core ion insulation factor','location','best');
	ylabel('eV / layer');
	xlabel('time (s)');

	%  subplot(4,1,4)
	%  plot(post.profil0d.temps,nblarmor_e ./ nblarmor_i);
	%  legend('ratio electron / ion','location','best');
	%  ylabel('-');

	joint_axes(hz,3);
	edition2
	set(hy,'String',sprintf('number Larmor radius or \norbit width in the plasma'));
end