% script de plot de nustar, rhostar, betan
if exist('mode_question_nustar','var')
    mode_plot_nustar = 0;
    mode_plot_transport = 0;
    mode_plot_edge = 0;
    mode_plot_stab = 0;
    mode_plot_edu = 0;
    K = menu('What plots do you want ?','nu* / rho*','transport related quantities (ITG, TEM, <s/q>, ...)', ...
             'SOl related quantities (SOL regime and lambda_q)','Stability diagrams','EDU plot (Frequencies an Coulomb logarithm)','all graphs','cancel');
    switch K
    case 1
	  mode_plot_nustar = 1;
    case 2
	  mode_plot_transport = 1;
    case 3
	  mode_plot_edge = 1;
    case 4
	  mode_plot_stab = 1;
    case 5
	  mode_plot_edu = 1;
    case 6
	  mode_plot_nustar = 1;
	  mode_plot_transport = 1;
	  mode_plot_edge = 1;
	  mode_plot_stab = 1;
	  mode_plot_edu = 1;
    otherwise 
      return
    end
else
    if ~exist('mode_plot_nustar','var')
	  mode_plot_nustar = 1;
    end
    if ~exist('mode_plot_transport','var')
	  mode_plot_transport = 1;
    end
    if ~exist('mode_plot_edge','var')
	  mode_plot_edge = 1;
    end
    if ~exist('mode_plot_stab','var')
	  mode_plot_stab = 1;
    end
    if ~exist('mode_plot_edu','var')
	  mode_plot_edu = 1;
    end
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

% palsma de fond
vt  =  ones(size(post.profil0d.temps));
switch post.z0dinput.option.gaz
    case 1
        zj = 1 * vt;
        aj = 1 * vt;
    case 2
        zj = 1 * vt;
        aj = 2 * vt;
    case 3
        zj = 1 * vt;
        aj = pchip(post.z0dinput.cons.temps,(2 + 3 .* post.z0dinput.cons.iso) ./ (1 + post.z0dinput.cons.iso),post.profil0d.temps);
    case 4
        zj = 2 * vt;
        aj = 4 * vt;
    case 5
         zj = pchip(post.z0dinput.cons.temps,(1 + 2 .* real(post.z0dinput.cons.iso)) ./ (1 + real(post.z0dinput.cons.iso)),post.profil0d.temps);
         aj = pchip(post.z0dinput.cons.temps,(1 + 3 .* real(post.z0dinput.cons.iso)) ./ (1 + real(post.z0dinput.cons.iso)),post.profil0d.temps);    
         warning('nHe3onD & nTonD not yet implemented !');
    case 11
         zj = pchip(post.z0dinput.cons.temps,(1 + 5 .* post.z0dinput.cons.iso) ./ (1 + post.z0dinput.cons.iso),post.profil0d.temps);
         aj = pchip(post.z0dinput.cons.temps,(1 + 25 .* post.z0dinput.cons.iso) ./ (1 + post.z0dinput.cons.iso),post.profil0d.temps);       
end

% profil de B0
B0 = sqrt(post.profil0d.fdia .^2 .* post.profil0d.r2i + post.profil0d.bpol .^ 2);

ve = ones(size(post.profil0d.xli));
aj   = aj * ve;
zj   = zj * ve;
amin = (post.profil0d.Raxe(:,end) .* post.profil0d.epsi(:,end)) * ve;

% calcul de rhostar
vte = sqrt(2 .* phys.e .* post.profil0d.tep ./ phys.me);
vti = sqrt(2 .* phys.e .* post.profil0d.tip ./ phys.mp ./ aj);
wce = phys.e * B0 ./ phys.me; 
wci = zj .* phys.e .*  B0./ phys.mp ./ aj;
%
rhoe = vte ./ wce;
rhostare = rhoe ./ amin;
rhoi = vti ./ wci;
rhostari=rhoi./ amin;

rhostare(:,1) =NaN;
rhostari(:,1) =NaN;

% debye length
debye_length = sqrt(2 .* phys.epsi0 ./ phys.e ./(post.profil0d.nep ./ post.profil0d.tep + ...
                    post.profil0d.nep .* post.profil0d.zeff ./ post.profil0d.tip));
%ld = 2.35e5 .* sqrt(post.profil0d.tep./ 1e3 ./ post.profil0d.nep);
%keyboard
ldstar   = debye_length ./ amin;
ldstar(:,1) =NaN;


% calcul de nustare
[nustare,nustari,tei]=z0fnustar(post.profil0d.tep,post.profil0d.tip,post.profil0d.nep,post.profil0d.qjli, ...
                                post.profil0d.Raxe,post.profil0d.Raxe .* post.profil0d.epsi,post.profil0d.zeff,zj);
nustare(:,1) =NaN;
nustari(:,1) =NaN;


% les champ en Rmax
rloc         = post.profil0d.Raxe .* ( 1 + post.profil0d.epsi);
ve           = ones(size(post.profil0d.xli));
btor         = post.z0dinput.option.signe .* (post.profil0d.fdia ./ rloc);
grho         = abs((post.profil0d.rmx(:,end) * ve) ./ max(eps,pdederive(post.profil0d.xli,rloc,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(post.profil0d.xli,post.profil0d.psi,0,2,2,1)./ rloc ./ (post.profil0d.rmx(:,end) * ve);
btot         = sqrt(btor .^ 2 + bpol .^ 2);
aloc         = post.profil0d.Raxe .* post.profil0d.epsi;
% largeur d'orbite injection de neutre
ral    = 7.2e-3.* sqrt(post.z0dinput.option.einj ./ 1e3) ./ btot; % en m
% patato
dp1     = post.profil0d.Raxe .* (2 .* post.profil0d.qjli .* ral ./ post.profil0d.Raxe) .^ (2/3);
% banana
dp2     = sqrt(post.profil0d.epsi) .* ral .* post.profil0d.qjli;
dp      = dp2 .* (dp2 < aloc) + dp1 .* (dp2 >= aloc);
% seul l'injection contre courant envoie les particules vers l'exterieur
dp     = (post.z0dinput.option.angle_nbi <= 0) .* dp + ral;
width_idn1 = max(eps,dp ./ amin);
width_idn1(:,1) = NaN;

% largeur d'orbite injection de neutre
ral    = 7.2e-3.* sqrt(post.z0dinput.option.einj2 ./ 1e3) ./ btot; % en m
% patato
dp1     = post.profil0d.Raxe .* (2 .* post.profil0d.qjli .* ral ./ post.profil0d.Raxe) .^ (2/3);
% banana
dp2     = sqrt(post.profil0d.epsi) .* ral .* post.profil0d.qjli;
dp      = dp2 .* (dp2 < aloc) + dp1 .* (dp2 >= aloc);
% seul l'injection contre courant envoie les particules vers l'exterieur
dp     = (post.z0dinput.option.angle_nbi2 <= 0) .* dp + ral;
width_idn2 = max(eps,dp ./ amin);
width_idn2(:,1) = NaN;


% largeur d'orbite alpha
ral    = 7.2e-3.* sqrt(3.55e6 ./ 1e3) ./ btot; % en m
% seul l'injection contre courant envoie les particules vers l'exterieur
% patato
dp1     = post.profil0d.Raxe .* (2 .* post.profil0d.qjli .* ral ./ post.profil0d.Raxe) .^ (2/3);
% banana
dp2     = sqrt(post.profil0d.epsi) .* ral .* post.profil0d.qjli;
dp      = dp2 .* (dp2 < aloc) + dp1 .* (dp2 >= aloc) + ral;
width_alpha = max(eps,dp ./ amin);
width_alpha(:,1) = NaN;


% ion thermique
switch post.z0dinput.option.gaz
    case 1
        mion = 1;
        zion = 1;
    case 2
        mion = 2;
        zion = 1;
    case 3
        mion = 3; % take T for losses
        zion = 1;
    case 4
        mion = 4;
        zion = 2;
    case 5
        mion = 3; % take He3 for losses
        zion = 2;
    case 11
        mion = 11; % take B11 for losses
        zion = 5;
    otherwise
        error('unknown main ion');
end

rho_i   = 4.57e-3 ./  zion .* sqrt(mion) .*  sqrt(max(13.6,post.profil0d.tip) ./ 1e3) ./ btot;
%          sqrt((post.profil0d.fdia .* post.profil0d.ri) .^ 2 + post.profil0d.bpol .^2);

rho_b   = rho_i .* post.profil0d.qjli ./ sqrt(post.profil0d.epsi);	
rho_pot = (rho_i .^ 2 .* post.profil0d.qjli .^ 2 .* post.profil0d.Raxe) .^ (1/3);
rho_b   = rho_b .* (rho_b < post.profil0d.rmx) + rho_pot .* (rho_b >= post.profil0d.rmx);
rho_b(:,1) = rho_pot(:,1);
rho_d  = 2.35e5  .* sqrt(max(13.6,post.profil0d.tep) ./ max(1e13,post.profil0d.nep)./ 1e3);
rho_bd = max(rho_b,rho_d);
width_ion = max(eps,rho_bd ./ amin);
width_ion(:,1) = NaN;


% Coulomb logarithm 
% ref Hager et Chang,  PoP 23  042503 (2016)
% + Wesson
Zlog = ((double(post.z0dinput.option.gaz == 4) + 1) .^ 2 .* (post.profil0d.nep ./ post.profil0d.nip) .* post.profil0d.zeff) .^ (1/4);
lambda_ee = 14.9  - 0.5 .* log(post.profil0d.nep ./ 1e20) + log(post.profil0d.tep ./ 1e3);
lambda_ei = 15.2  - 0.5 .* log(post.profil0d.nep ./ 1e20) + log(post.profil0d.tep ./ 1e3);
lambda_ii = 30 - log(Zlog .^ 3 .* sqrt(post.profil0d.nip) ./ (post.profil0d.tip ./ 1e3) .^ (3/2));

% electron mean free path 
% Wesson 
lpm_e    = 1.44e23 .* (post.profil0d.tep ./ 1e3) .^ 2 ./ post.profil0d.nep ./ lambda_ei;
lpm_e_star = lpm_e ./ amin;

% selon le plot demande
if mode_plot_nustar == 1
      % le plot
      hz =findobj(0,'type','figure','tag','nustar1');
      if isempty(hz)
		hz=figure('tag','nustar1','name','dimensionless parameters (profiles)');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',1,'color',[1 1 1])

      subplot(2,2,1)
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,nustare,'color','b');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,nustari,'color','r');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,ones(size(nustari)),'color','g');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli(2:end),post.profil0d.epsi(:,2:end) .^ (-3/2),'color','g');
      set(gca,'yscale','log');
      ylabel('\nu ^{*}_e (b) & \nu ^{*}_i (r)');
      xlabel('x');
      text(0.1,0.1,'banana','color','g');
      text(0.1,3,'plateau','color','g');
      text(0.5,3 .* max(max(post.profil0d.epsi(:,2:end) .^ (-3/2))),'Pfirsch-Schluter','color','g');
      z0loglin(gca);

      subplot(2,2,2)
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rhostare,'color','b');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,rhostari,'color','r');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,ldstar,'color','g');
      %zplotprof(gca,post.profil0d.temps,post.profil0d.xli,ones(size(rhostari)),'color','g');
      set(gca,'yscale','log');
      ylabel('\rho ^{*}_e (b) ; \rho ^{*}_i (r) & \lambda_{D}^{*} (g)');
      xlabel('x');
      z0loglin(gca);

      subplot(2,2,3)
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,width_idn1,'color','b');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,width_idn2,'color','c','linestyle','-.');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,width_alpha,'color','r');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,width_ion,'color','m');
      x1001 = linspace(0,1,1001);
      lim1001 = interp1(post.profil0d.xli',(eps + ones(size(post.profil0d.temps)) * (1 - post.profil0d.xli))',x1001,'linear','extrap')';
      zplotprof(gca,post.profil0d.temps,x1001,lim1001,'color','g');
      %zplotprof(gca,post.profil0d.temps,post.profil0d.xli,eps + ones(size(post.profil0d.temps)) * (1 - post.profil0d.xli),'color','g');
      set(gca,'yscale','log');
      ylabel('normalised orbit width');
      legend('NBI 1','NBI 2','\alpha','main ion');
      set(gca,'ylim',[0.01,1]);
      xlabel('x');
      text(0.1,0.1,'confined','color','g');
      text(0.5,1,'lost','color','g');
      z0loglin(gca);

      subplot(2,2,4)
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.qjli,'color','b');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,ones(size(post.profil0d.qjli)),'color','g');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,2 .* ones(size(post.profil0d.qjli)),'color','g');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,1.5 .* ones(size(post.profil0d.qjli)),'color','g');
      ylabel('q');
      xlabel('x');


      % le plot
      hz =findobj(0,'type','figure','tag','nustar2');
      if isempty(hz)
		hz=figure('tag','nustar2','name','dimensionless parameters (averaged)');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',1,'color',[1 1 1])
      %  qmoy = pchip(post.profil0d.temps,trapz(post.profil0d.xli,post.profil0d.spr .* post.profil0d.qjli,2) ./ ...
      %               trapz(post.profil0d.xli,post.profil0d.spr,2),post.zerod.temps);
      %  [nustarem,nustarim,teim]=z0fnustar(post.zerod.tem,post.zerod.tem .* post.zerod.tite,post.zerod.tem,qmoy, ...
      %                                  post.z0dinput.geo.R,post.z0dinput.geo.a,post.zerod.zeff, ...
      %                              (post.z0dinput.option.gaz == 4) + 1);
      nustarem = trapz(post.profil0d.xli(2:end),post.profil0d.vpr(:,2:end) .* nustare(:,2:end),2) ./ ...
		    trapz(post.profil0d.xli(2:end),post.profil0d.vpr(:,2:end),2);
      nustarim = trapz(post.profil0d.xli(2:end),post.profil0d.vpr(:,2:end) .* nustari(:,2:end),2) ./ ...
		    trapz(post.profil0d.xli(2:end),post.profil0d.vpr(:,2:end),2);
      rhostarem = trapz(post.profil0d.xli(2:end),post.profil0d.vpr(:,2:end) .* rhostare(:,2:end),2) ./ ...
		    trapz(post.profil0d.xli(2:end),post.profil0d.vpr(:,2:end),2);
      rhostarim = trapz(post.profil0d.xli(2:end),post.profil0d.vpr(:,2:end) .* rhostari(:,2:end),2) ./ ...
		    trapz(post.profil0d.xli(2:end),post.profil0d.vpr(:,2:end),2);
      subplot(3,1,1)
      semilogy(post.profil0d.temps,nustarem,'b',post.profil0d.temps,nustarim,'r');
      legend('<\nu ^{*}_e>','<\nu ^{*}_i>');
      %xlabel('time (s)');
      title(sprintf('Zerod : %s@%d/energy ',post.z0dinput.machine,post.z0dinput.shot));
      z0loglin(gca);
      subplot(3,1,2)
      semilogy(post.profil0d.temps,rhostarem,'b',post.profil0d.temps,rhostarim,'r');
      legend('<\rho ^{*}_e>','<\rho ^{*}_i>');
      %xlabel('time (s)');
      z0loglin(gca);
      subplot(3,1,3)
      plot(post.zerod.temps,post.zerod.betap,'b',post.zerod.temps,post.zerod.betaptot,'c',post.zerod.temps,post.zerod.betan .* 100,'r');
      legend('\beta _{P(th)}','\beta _{P}','\beta _{N} (%)');
      xlabel('time (s)');

      joint_axes(hz,3);

      % le plot
      hz =findobj(0,'type','figure','tag','nustar3');
      if isempty(hz)
		hz=figure('tag','nustar3','name','dimensionless parameters (images)');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',1,'color',[1 1 1])

      subplot(2,2,1)
      imagesc(post.profil0d.xli(2:end),post.profil0d.temps,log(nustare(:,2:end))./log(10));			
      colormap('default')
      colorbar
      set(gca,'ydir','normal')
      title('log( \nu ^{*}_e)');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,2)
      imagesc(post.profil0d.xli(2:end),post.profil0d.temps,log(rhostare(:,2:end))./log(10));
      colormap('default')
      colorbar
      set(gca,'ydir','normal')
      title('log(rho ^{*}_e) ');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,3)
      imagesc(post.profil0d.xli(2:end),post.profil0d.temps,log(nustari(:,2:end))./log(10));
      colormap('default')
      colorbar
      set(gca,'ydir','normal')
      title('log( \nu ^{*}_i)');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,4)
      imagesc(post.profil0d.xli(2:end),post.profil0d.temps,log(rhostari(:,2:end))./log(10));
      colormap('default')
      colorbar
      set(gca,'ydir','normal')
      title('log(rho ^{*}_i) ');
      xlabel('x');
      ylabel('time (s)');
end
if mode_plot_transport == 1
      hz =findobj(0,'type','figure','tag','nustar4');
      if isempty(hz)
		hz=figure('tag','nustar4','name','dimensionless parameters (images)');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',1,'color',[1 1 1])

      subplot(2,2,1)
      imagesc(post.profil0d.xli(2:end),post.profil0d.temps,1./post.profil0d.qjli(:,2:end));
      colormap('default')
      colorbar
      set(gca,'ydir','normal')
      title('\iota = 1/q');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,2)
      imagesc(post.profil0d.xli(2:end),post.profil0d.temps,post.profil0d.omega(:,2:end));
      colormap('default')
      colorbar
      set(gca,'ydir','normal')
      title('\omega (rad/s)');
      xlabel('x');
      ylabel('time (s)');


      marge = ones(size(post.profil0d.temps)) * (1 - post.profil0d.xli);

      subplot(2,2,3)
      imagesc(post.profil0d.xli(2:end),post.profil0d.temps,tanh(max(width_idn1(:,2:end),width_idn2(:,2:end)) - marge(:,2:end)));
      colormap('default')
      colorbar
      set(gca,'ydir','normal')
      title('width_{nbi} (Fast ions lost if >0)');
      xlabel('x');
      ylabel('time (s)');

      subplot(2,2,4)
      imagesc(post.profil0d.xli(2:end),post.profil0d.temps,tanh(width_alpha(:,2:end) - marge(:,2:end)));
      colormap('default')
      colorbar
      set(gca,'ydir','normal')
      title('width_{alpha} (Fast alpha lost if >0)');
      xlabel('x');
      ylabel('time (s)');
end

% le plot
if mode_plot_edu == 1
      hz =findobj(0,'type','figure','tag','nustar5');
      if isempty(hz)
		hz=figure('tag','nustar5','name','Coulomb logarithms');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',1,'color',[1 1 1]);

      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,lambda_ee,'color','b');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,lambda_ei,'color','r');
      zplotprof(gca,post.profil0d.temps,post.profil0d.xli,lambda_ii,'color','k');
      legend('\lambda_{ee}','\lambda_{ei}','\lambda_{ii}');
      ylabel('Coulomb logarithms');
      xlabel('x');
end

% plot suplementaires
if mode_plot_stab == 1
    z0plotstabdiag(post);
end
if mode_plot_edu == 1
    z0plotfreq(post);
end
z0plotmode;
z0plotsoq;
if mode_plot_edge
    z0plot_SOL_regime;
end