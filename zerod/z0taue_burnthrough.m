% reference : H.T. Kim et al, NF 52 (2012) 103016
%             28 GHz START-UP SYSTEM ON MAST,V. Shevchenko et al, 34th EPS Conference on Plasma Phys. Warsaw, 2 - 6 July 2007 ECA Vol.31F, P-4.159 (2007)
%             A.D. Mac Donald, microwave breakdown in gases New York 1966
% post.z0dinput.option.berror  = 1e-3;
% post.z0dinput.option.L_eddy = 0;
% post.z0dinput.option.R_eddy = 0;
% taue = z0taue_burnthrough(post.z0dinput.option,post.z0dinput.geo,post.z0dinput.cons,post.zerod,post.zerod.tauthl);
function [taue,f,ieddy,nec,tref,prf,flux_bord_cor,dt_breakdown,ip_c,x_ioniz,gas_puff_equi,f_wth] = z0taue_burnthrough(option,geo,cons,zerod,tauref,ieddy_mem,flux_bord_cor_mem)

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

%securite
tauref = max(1e-6,tauref);

% isotopic composition for option.gaz == 5
if option.gaz == 5
    nHe3onD = real(cons.iso);
    nTonD   = imag(cons.iso);
    warning('nHe3onD & nTonD not yet used !');
else
    nHe3onD = zeros(size(cons.iso));
    nTonD   = real(cons.iso);
end
cons.iso = real(cons.iso);



if option.L_eddy == 0
 	Led = phys.mu0 * geo.R .* (log(8 *  2) - 3/2);
else
	Led = ones(size(cons.temps)) .* option.L_eddy;
end 
if option.R_eddy == 0
	Red =  2 .* pi .* geo.R .* 1e-5;
else
	Red =  ones(size(cons.temps)) .* option.R_eddy;
end
% cas du couplage avec FREEBIE
if isfield(cons,'b_error')
	berror = cons.b_error;
elseif option.berror == 0
	berror = eps;
else
	berror = ones(size(cons.temps)) .* option.berror;
end

tau_eddy    = Led ./ Red;
% la structure passive voie la meme variation de flux que le bord du plasma
vloop       = zerod.vloop;


if option.breakdown < 0
	v0 = abs(option.breakdown);
	% loop voltage limitation to the value @ breakdown that is the maximum allows for a given tokamak
	if option.vloop == 0
	    vloop = max(min(vloop,abs(option.breakdown)),-abs(option.breakdown));
	end
else
	v0 = zerod.vloop(1);
end
%  if (nargin >= 6)
%      % rien
%  elseif ~isfield(option,'evolution')
if ~isfield(option,'evolution')
    vloop(1) = v0;
elseif option.evolution == 0
    vloop(1) = v0;
end
if isfield(cons,'ieddy')
    % case METIS coupled to FREEBIE
    ted   = cons.temps;
    ieddy = cons.ieddy;
    %disp('g0')
elseif (nargin >= 6) && (length(ieddy_mem) == length(cons.temps))
    % case evolution computation of confinement time
    ted   = cons.temps;
    ieddy = ieddy_mem;
    %disp('g1')
elseif (nargin >= 6)
    % case evolution
    ted   = cons.temps;
    ieddy = ones(size(ted)) * ieddy_mem;
    [ted(2:end),ieddy(2:end)] = z0ode(cons.temps(2:end),vloop(2:end) ./ Led(2:end),tau_eddy(2:end),ieddy_mem);
    %disp('g2')
elseif option.I_eddy == 0 
    [ted,ieddy] = z0ode(cons.temps,vloop ./ Led,tau_eddy,0);
    %disp('g3')
elseif option.I_eddy == -1
   [ted,ieddy] = z0ode(cons.temps,vloop ./ Led,tau_eddy,-max(v0 ./ Red(1),zerod.ip(1)));
   %disp('g3')
else
    [ted,ieddy] = z0ode(cons.temps,vloop ./ Led,tau_eddy,max(v0 ./ Red(1),zerod.ip(1)));
    %disp('g4')
end
% champ du au courant induit dans les structures
beddy       = option.B_eddy .* phys.mu0 .* abs(ieddy) ./ (2 .* pi .* geo.R / 2);
% champ vertical (en supposant le plasma en mode limiter en non le cannal initial sans contact)
bv = phys.mu0 .* zerod.ip ./ (4 .* pi .* geo.R) .* (log(8 .* geo.R ./((zerod.sp/pi).^0.5)) + zerod.betaptot + zerod.li ./ 2 - 3/2);
% poloidal field at plasma surface
bpol = phys.mu0 .* zerod.ip ./ zerod.sp;
% passage d'un regime  � l'autre
f  = max(0,tanh(bpol ./ sqrt(beddy .^ 2  + berror .^ 2) - 3));
f_plot = f;
%figure(21);plot(cons.temps,f);drawnow

if isfield(option,'initiation_only') && ~isempty(option.initiation_only)
  if option.initiation_only == 1
      f = zeros(size(f));
  end
end
if nargout == 0
      f = zeros(size(f));
end
%figure(21);clf;plot(cons.temps,bv ./ sqrt(beddy .^ 2  + berror .^ 2));set(gca,'xlim',[0 3]);drawnow
%calcul du petit rayon effectif au demarrage du plasma
ap = min(geo.a,max(geo.a ./ 10,phys.mu0 .* zerod.ip ./ sqrt(beddy .^ 2  + berror .^ 2 + bv .^ 2) ./ 2 ./ pi ./ (1 -  (zerod.betaptot + zerod.li ./ 2 -1) .* geo.a ./ geo.R)));
if ap(1) == geo.a(1)
  ap(1) = min(geo.a(1),ap(2));
end
x_channel = ap ./ geo.a;
%figure(21);clf;plot(cons.temps,x_channel);set(gca,'xlim',[0 3]);drawnow
% utilisation de la definition originale
% for confinement
lf        = 0.25 .* ap .* geo.b0 ./ sqrt( beddy .^ 2  + berror .^ 2 + bv .^ 2);
% for breaking down before any plasma current
% take minor radius as proxy for effective radius
lf0        = 0.25 .* geo.a .* geo.b0 ./ sqrt( beddy .^ 2  + berror .^ 2);
%figure(21);clf;plot(cons.temps,lf,'b',cons.temps,lf0,'r');set(gca,'xlim',[0 3]);drawnow


% correction du champ electrique au bord du plasma
mu0   = 4 .* pi .* 1e-7;
L_p   = mu0 * geo.R .* (log(8 *  geo.R ./ ((zerod.sp/pi).^0.5)) - 2)  + mu0 * geo.R .* zerod.li / 2;
if isfield(option,'PSI_eddy')
     M_pv  = option.PSI_eddy .* sqrt(L_p .* Led);
else
     M_pv  = 1 .* sqrt(L_p .* Led);
end
dpsidt_edge_cor  = - z0dxdt(M_pv .* ieddy,cons.temps) ./ 2 ./ pi;
dpsidt_edge_cor(1) = dpsidt_edge_cor(2);
%figure(22);clf;plot(cons.temps,dpsidt_edge_cor);set(gca,'xlim',[0 3]);drawnow
if isfield(cons,'ieddy')
    % case METIS coupled to FREEBIE -> FREEBIE provides directly the rigth value of cons.flux
    flux_bord_cor = zeros(size(cons.temps));
    %disp('h0')
elseif (nargin >= 7) && (length(flux_bord_cor_mem) == length(cons.temps))
    % case evolution computation of confinement time
    flux_bord_cor = flux_bord_cor_mem;
    %disp('h1')
elseif (nargin >= 7)
    % case evolution
    flux_bord_cor    = cumtrapz(cons.temps,dpsidt_edge_cor,1);
    flux_bord_cor    = flux_bord_cor_mem - flux_bord_cor(2) + flux_bord_cor;
    %disp('h2')
else
    flux_bord_cor    = cumtrapz(cons.temps,dpsidt_edge_cor,1);
    %disp('h3')
end

%  % these E. Tsitrone
%  lclim = pi .* geo.R .* zerod.qa;
%  lcpol = pi .* geo.R;
%  lcx = sqrt(zerod.peri .^ 2  + (pi .* geo.R .* option.lcx .* zerod.qa) .^ 2);  
%  switch option.configuration
%  case 0
%  	lc = lcpol;
%  case 1
%  	lc = lclim;
%  case 2
%  	lc  = zerod.xpoint .* lcx + (~zerod.xpoint) .* lcpol;
%  case 3
%  	lc  = zerod.xpoint .* lcx + (~zerod.xpoint) .* lclim;
%  otherwise
%  	lc  = lcx;
%  end
%  lf          = max(lc,lf0);
% temps de propagation le long des lignes de champs
cs          = sqrt(zerod.tem .* phys.e .* (zerod.tite + zerod.zeff)./ phys.mp ./ zerod.meff);
tau_par     = lf ./ cs;
%tau_par     = min(tau_par,tauref);

% modification pour utiliser le modele de la reference
% ref NF H.T. Kim et al, NF 52 (2012) 103016 (15p)
DBohm   = zerod.tem ./ 16 ./ geo.b0; 
tau_per = ap .^ 2 ./ 2 ./ DBohm;
%figure(21);plot(cons.temps,tau_par,'b',cons.temps,tau_per,'r',cons.temps,tauref,'g');drawnow
tau_par = 1 ./ (1 ./ tau_par + 1 ./ tau_per); 
tau_par     = min(tau_par,tauref);

% convergence pour le point initial
prf = zeros(size(tau_par));
for k=1:31
    % densit� critique pour l'amorcage
    % densite initiale
    na      = option.p_prefill ./ phys.k ./ option.temp_vac;
    tref    = option.temp_vac .* phys.k ./ phys.e + tau_par .* prf ./ 3 ./ na ./ zerod.vp ./ phys.e;
    gamma   = 2.e-3 .* tref .^ (3/2);
    nec     = gamma ./ (1 + gamma) .* na;


    % absorbtion of ECRH 
    % ref: Theory of electron cyclotron absorption in magnetized plasma, M. Bornatici, PoP 1982 vol 24 number 6 p 629 to 638.
    % ref :Young Soon Bae and A.C. England, Journal of Korean Physical Society, vol 51 number 4 2007 p 1313-1319
    % n_par = 0.5 + first harmonique & mode O or mode X optimized for the machine
    wce = phys.e .* geo.b0 ./ phys.me;
    lce = 2 .* pi .* phys.c ./ wce;
    nc    = phys.epsi0 .* geo.b0 .^ 2 ./ phys.me;
    alpha = max(nec,zerod.nem) ./ nc;
    n_par = 0.5;
    % kept only main ndependances
    opdepth = pi .^ 2 .* geo.R ./ lce .* phys.e .* zerod.tem ./ (phys.me .* phys.c .^ 2) .* alpha;
    % ajout de la fraction absorb�e avnt que le plasma confine ne soit etabli
    if option.lhmode == 5
      prf_new = (1 - f) .* (cons.pecrh + cons.plh) .* (1 - exp(- opdepth));
    else
      prf_new = (1 - f) .* cons.pecrh .* (1 - exp(- opdepth));
    end
    %  figure(21);
    %  clf
    %  subplot(4,1,1)
    %  plot(cons.temps,alpha);
    %  subplot(4,1,2)
    %  plot(cons.temps,opdepth);
    %  subplot(4,1,3)
    %  plot(cons.temps,(1 - exp(- opdepth)));
    %  subplot(4,1,4)
    %  plot(cons.temps,prf);
    %  drawnow
    %  %keyboard
    err = sqrt(sum((prf_new-prf) .^ 2) ./ sum(prf .^ 2 + prf_new .^ 2));
    if k == 1
      prf = prf_new;
    else
      prf = 0.3 .* prf_new + 0.7 .* prf;
    end
    if err < eps
      break;
    end
end

% nouveau model pour la modification du champ electrique � l'amorcage
% ref : M. Hasegawa et al, Plasma and fusion research vol 2 (2007) p 007-
% pour LH = mode de cavite tout le volume de la chambre est pris en compte
% pour EC, volume de la couche est donn� dans :
% S. K . Borowski et al, Fusion and Technology 1984 vol 6 p7-
% formule 7
delta_R = 1.03e-19 .* zerod.nem ./ geo.b0 .* geo.R ./ 2;
vol_EC  = 4 .* pi .* geo.a .* geo.K .*  delta_R .* geo.R;
if isfield(option,'VV_volume') && (option.VV_volume > 0)
     vol_EC = min(option.VV_volume,vol_EC);
else
     vol_EC = min(zerod.vp,vol_EC);
end

if option.lhmode == 5
  erf_ec  = sqrt(2 .* (cons.plh + cons.pecrh) ./ phys.epsi0 ./ vol_EC);
  erf_lh  = zeros(size(erf_ec));
else
  if isfield(option,'VV_volume') && (option.VV_volume > 0)
      erf_lh  = sqrt(2 .* cons.plh ./ phys.epsi0 ./ option.VV_volume);
  else
      erf_lh  = sqrt(2 .* cons.plh ./ phys.epsi0 ./ max(zerod.vp));
  end
  erf_ec  = sqrt(2 .* cons.pecrh ./ phys.epsi0 ./ vol_EC);
end
if isfield(option,'VV_volume') && (option.VV_volume > 0)
    erf_ic  = sqrt(2 .* cons.picrh ./ phys.epsi0 ./ option.VV_volume);
else
    erf_ic  = sqrt(2 .* cons.picrh ./ phys.epsi0 ./ max(zerod.vp));
end
%
lambda_lh = phys.c ./ option.freqlh ./ 1e9;
lambda_ec = phys.c ./ (phys.e .* geo.b0 ./ phys.me);
lambda_ic = phys.c ./ option.freq ./ 1e6;
%% formules sont en Torr
%%p_prefill_torr = option.p_prefill ./ 133.322368;
lpath      = 133.322368 .* 0.556 .* sqrt(zerod.meff ./ 2) ./ option.p_prefill ./ zerod.zeff;
eeff_lh    = erf_lh ./ sqrt(1 + (lpath ./ lambda_lh) .^ 2);
eeff_ec    = erf_ec ./ sqrt(1 + (lpath ./ lambda_ec) .^ 2);
eeff_ic    = erf_ic ./ sqrt(1 + (lpath ./ lambda_ic) .^ 2);
% negligeable
%prf     = phys.epsi0 .* eeff .^ 2 ./ 2 .* zerod.vp;
% tension par tour minimum	
switch option.gaz
    case 4
        % For He, ref: K.T.A.L. Burm, Contrib. Plasma Phys. 47, no 3, 177-182 (2007)
        vloop_min = max(eps,2 .* pi .* geo.R .* 25.4 .* option.p_prefill ./ max(eps,log(2.10 .* option.p_prefill .* lf0)));
    case 5
        % assume ratio given by iso
        vloop_min_he = max(eps,2 .* pi .* geo.R .* 25.4 .* option.p_prefill ./ max(eps,log(2.10 .* option.p_prefill .* lf0)));
        vloop_min_h  = max(eps,2 .* pi .* geo.R .* 93.8  .* option.p_prefill ./ max(eps,log(3.83 .* option.p_prefill .* lf0)));
        vloop_min    = (vloop_min_h + cons.iso .* vloop_min_he) ./ (1 + cons.iso);
    otherwise
        % for option.gaz == 11 assume to start the plasma in hydrogen
        vloop_min = max(eps,2 .* pi .* geo.R .* 93.8  .* option.p_prefill ./ max(eps,log(3.83 .* option.p_prefill .* lf0)));
end
% effet RF sur le champ min
% ref: M. Hasegawa e al, Plasma and fusion reasearch, vol 2 (2007) p 7-
e0         = vloop_min ./ (2 .* pi .* geo.R);
fcoll      = 1.78e7 .* option.p_prefill .* zerod.zeff;
factor_rf  = min(1,0.5 .* (eeff_ec ./ e0) .^ 2 + (0.5 .* (eeff_lh ./ e0) .^ 2 + 0.5 .* (eeff_ic ./ e0) .^ 2) .* ...
             fcoll .^ 2 ./ (wce .^ 2 + fcoll .^ 2));

% tension par tour minimum
vloop_min_plot = vloop_min;
vloop_min = vloop_min .* (1 + sqrt(1 - factor_rf)) ./ 2; 
%figure(21);clf;subplot(2,1,1);semilogy(cons.temps,vloop_min,cons.temps,zerod.vloop);subplot(2,1,2);plot(cons.temps,lf);drawnow
% condition pour le breakdown
% if ip is provided, v0 is sustained up to breakdown
if option.vloop == 0
    ind_go    = find((max(abs(v0),abs(zerod.vloop)) >= vloop_min),1);
else
    ind_go    = find((abs(zerod.vloop) >= vloop_min),1);
end
if ~isempty(ind_go)
	f(1:ind_go) = 0;
	f(ind_go:end) = max(f(ind_go:end),sqrt(eps));
else
	f(:) = 0;
end
% calcul de taue en sortie
taue        = tau_par  .* (1 - f) + tauref .* f;

% duree pour le break down
ind_break = find(f>= sqrt(eps),1);
if isempty(ind_break)
    dt_breakdown = 0;
else
    dt_breakdown = ted(ind_break);
end


% densite de courant critique au breakdown
% ref: H.-T. Kim et al, NF 52 '(2012) 103016
jc    = gamma ./ 0.002 .* 382.5 .* zerod.vloop  ./ (2 .* pi .* geo.R);
% temperature must be kept over e0
e0 = z0eioniz_div(option.temp_vac.* phys.k ./ phys.e,zerod.ne0);
Vioniz  = zerod.wth ./ ((3/2) .* phys.e .* e0 .* zerod.ne0);
%x_ioniz = max(0.1,sqrt(min(1,Vioniz ./ zerod.vp)));
%x_ioniz = max(0.1,min(x_channel,sqrt(min(1,Vioniz ./ zerod.vp))));
x_ioniz = max(0.1,min(x_channel,sqrt(min(1,Vioniz ./ zerod.vp))));
ip_c  = zerod.sp .* jc .* x_ioniz .^2;
%fprintf('ip_c = %g\n',ip_c(1));
%figure(21);clf;plot(cons.temps,zerod.wth);drawnow

% equivalent gas puffing during breakdown due to pumping from plasma
n0a_static   = zerod.nem .* zerod.vp ./ max(1e-6,zerod.taup); % flux out from plasma if steady state
if isfield(option,'VV_volume') && (option.VV_volume > 0)
    n0_vac_available = max(0,na .* option.VV_volume  - zerod.nem .* zerod.vp.* x_ioniz .^ 2) ./ option.VV_volume; %  pressure evolution 
else
    n0_vac_available = max(0,na .* max(zerod.vp)  - zerod.nem .* zerod.vp.* x_ioniz .^ 2) ./ max(zerod.vp); %  pressure evolution 
end
seffective  = zerod.sext .* x_ioniz;
gas_puff_equi = max(0,n0_vac_available .* sqrt(option.temp_vac .* phys.k ./ phys.ua ./ zerod.meff) .* seffective - n0a_static);
%figure(21);clf;subplot(2,1,1);plot(cons.temps,gas_puff_equi,cons.temps,n0a_static);subplot(2,1,2);plot(cons.temps,n0_vac_available,cons.temps,na);

%Effect of neutral charge exchange on plasma energy contents and temperature
%W = 3/2* (ne*te + ni*ti + n0_cx * T_cx);
%we assume T_cx = T_i
tesec = max(zerod.tem,option.temp_vac .* phys.k ./ phys.e);
[svi1s,svi2s,svcx,svrec,sii,sss,Ass] = z0sectionh(tesec,zerod.tite .* tesec);
n0_in = zerod.n0a  ./ zerod.vp ./ max(nec,zerod.nem) ./ svcx;
f_wth = 1 + (n0_in .* max(0.1,zerod.tite) .* tesec ./ ...
      (max(nec,zerod.nim) .* max(0.1,zerod.tite) .* tesec + ...
      max(nec,zerod.nem) .* tesec));
%figure(21);clf;plot(cons.temps,f_wth);set(gca,'ylim',[0,10]);drawnow;


if nargout > 0
   return
end
fullscreen = get(0,'ScreenSize');
hz =findobj(0,'type','figure','tag','breakdown_burnthrough');
if isempty(hz)
  	  hz=figure('tag','breakdown_burnthrough','name','Breakdown & Burnthrough');
else
  	  figure(hz);
end
clf
set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

%
k  = 6;
subplot(k,1,1);
semilogy(cons.temps,tauref,'r',cons.temps,tau_par,'b');
z0loglin(gca);
title(sprintf('METIS : %s@%d/Breakdown @ t = %g s', ...
          option.machine,option.shot,dt_breakdown));
legend('tau_L (s)','tau_{E,breakdown & burn through} (s)');
subplot(k,1,2)
plot(cons.temps,f_plot,cons.temps,x_ioniz);
legend('Confimement indicator: 0= open field line -> 1= close field line', ...
       'Width of ionised plasma normalized to limited plasma minor radius');
subplot(k,1,3)
plot(cons.temps,vloop_min,'r',cons.temps,vloop_min_plot,'.b',cons.temps,abs(vloop),'k',dt_breakdown,0,'og');
z0loglin(gca);
legend('V_{Townsend,RF} (V)','V_{Townsend} (V)','Vloop (V)','Breakdown time');
set(gca,'ylim',[eps,max(abs(vloop))])
subplot(k,1,4)
plot(cons.temps,prf ./ 1e6,cons.temps,ip_c ./ 1e3 );
legend('Absorbed EC (MW)','Initial current @ breakdown (kA)');
subplot(k,1,5)
semilogy(cons.temps,nec ./ 1e19,cons.temps,gas_puff_equi ./ 1e19);
z0loglin(gca);
legend('Initial electron density @ breakdown (10^{19} m^{-3})','Plasma pump in (10^{19} e^-/s)');
subplot(k,1,6)
plot(cons.temps,beddy,cons.temps,bv,cons.temps,ones(size(cons.temps)) .* berror);
z0loglin(gca);
legend('B_{eddy}','B_v','B_{error}');
xlabel('time (s)');
joint_axes(hz,k);
edition2
set(hz,'Position',fullscreen)
drawnow


%figure(21);clf;plot(ted,taue);drawnow
%  figure(21);clf;
%  subplot(2,1,1)
%  plot(ted,nec,ted,zerod.nem);
%  subplot(2,1,2)
%  plot(ted,vloop_min,ted,vloop);
%  drawnow
