% cette fonction donne la solution pour 1 appel
% option -> option du zerod
%    option.gaz  = 1 -> H, 2 -> D , 3 -> DT, 4 -> He , 5 -> D-He3 et 11 -> p-B11'
%    option.modeh =  0 -> pas de mode H possible, 1-> mode h possible (autre valeur pour amelioration intermediaire)
%    option.zmax = charge de l'impurete assurant le rayonnement
%    option.rw   = coefficient de reflexion des parois pour le rayonnement cyclotronique
%    option.fwcd = 1 si chauffage FCI en mode FWCD
%    option.transitoire   =  1 -> traite les transitoire; 0 -> stationnaire uniquement (d/dt = 0)
%    option.vloop  = 0 -> ip donne ; 1 -> vloop = 0
%    option.xece     = 0..1 : position radiale d'injection
%    option.angle_ece  = 0..180 degres: angle equivalent d'injection
%    option.synergie = effet de la synergie sur le courant genere par eccd (0 = pas d'effet)
%    option.sens     = 1 -> co courant, -1 -> contre courant
%    option.angle_nbi  = angle d'injection des neutres -90..90 degres (<0 -> contre courant, >0 co courant, 0 = pas de courant)
%    option.einj       = energie des neutres en eV
%    option.lhmode     = choix de la formule pour le calcul de eta LH
%    option.etalh      = efficacite de generation de courant de l'hybride dans le mode constant
%    option.itb      =  seuil de declenchement des ITB en unite de beta +li/2 , si >= 1000 pas d'itb
%     cons.ip    = courant plasma (A)
%     cons.nbar  = densite lineique (m^ -3)
%     cons.picrh = puissance FCI (W)
%     cons.plh   = puissance LH (W)
%     cons.pnbi  = puissance NBI (W)
%     cons.pecrh = puissance FCE (W)
%     cons.zeff  = zeff du plasma sans les alpha
%     cons.iso   = consigne pour le rapport isotopique
%     cons.temps = vecteur temps du zerod
% geo    -> geometrie
%     geo.R     = grand rayon (m)
%     geo.a     = petit rayon (m)
%     geo.K     = elongation du plasma sur la DSMF (b/a)
%     geo.b0    = champ toroidal (T)
%  donnees optionnelles :
%     geo.vp    = volume du plasma (m^3)
%     geo.sp    = surface d'une section poloidale du plasma (m^ 2)
%     geo.sext  = surface externe du plasma (m^ 2)
% zs     -> strcuture de donnees du zerod

function [zs,profli] = zero1t(option,cons,geo,zs,amorti,profli)

if nargin == 0
   zs = zerod_0dinfo;
   return 
end

%save('sat1','option','cons','geo','zs','amorti','profli');
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

% energy of fusion produced particle
fus.dhe3_he4   = 3.71e6; % (eV) 
fus.dhe3_p     = 14.64e6; % (eV)
fus.ddn_he3    = 0.82e6; % (eV)
fus.ddp_t      = 1.01e6; % (eV)
fus.ddp_p      = 3.02e6; % (ev)
fus.dt_he4     = 3.56e6; % (eV)
                
% flag pour debug
fwr = option.plotconv;


if nargin < 4
   mode  = 0;  % mode init
   zs    = [];
   amorti = 0.5;
elseif isempty(zs)
   % mode init
   mode  = 0;
else
   mode  = 1;
end

if nargin < 6
	profli.xli = linspace(0,1,21);
end

% dans zero1t c'est donnees ne sont pas utilisees
zs.dw      = 0;
zs.dpfus   = 0;
zs.dini    = 0;
zs.diboot  = 0;

% separation  nbar et gaspuff
gas_puff_mem = imag(cons.nbar);  % en e-/s
cons.nbar    = max(eps,real(cons.nbar));

% activation ou desactivation de IDN2
if option.nb_nbi == 1
  cons.pnbi = real(cons.pnbi);
end

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

% valeurs utiles pour les impuretees
zu1 = (option.zimp + option.rimp .* option.zmax);
zu2 = (option.zimp .^ 2 + option.rimp .* option.zmax .^ 2);

if (option.W_effect == 1) && isfield(profli,'nep') && isfield(profli,'nzp')
    [profli.nwp,zs.nwm,zu1w,zu2w] = z0acctungsten(option,geo,cons,zs,profli);
    zu1 = zu1 + zu1w; 
    zu2 = zu2 + zu2w;
else
    profli.nwp = zeros(length(cons.temps),21);
    zs.nwm     = zeros(length(cons.temps),1);
end

% donnees de volume et de surface
if (mode == 0) || (option.init_geo == 1)
    if fwr == 1; disp('zgeo0');end
    if isappdata(0,'EQUILIBRIUM_EXP')  
        equi_ext  = getappdata(0,'METIS_EXTERNAL_CURDIF_EQUI');
        vp        = equi_ext.vp;
        sp        = equi_ext.sp;
        sext      = equi_ext.sext;
        peri      = equi_ext.peri;
        geo       = equi_ext.geo;
        zs.xpoint = double(equi_ext.xpoint);
        Rsepa     = equi_ext.Rsepa;
        Zsepa     = equi_ext.Zsepa;      
    elseif isfield(profli,'Rsepa') && isfield(profli,'Zsepa') && ~isempty(profli.Rsepa) && ~isempty(profli.Zsepa)
        [vp,sp,sext,peri,geo,zs.xpoint] = zgeo0(geo,profli.Rsepa,profli.Zsepa,option.evolution);
        Rsepa = profli.Rsepa;
        Zsepa = profli.Zsepa;
        zs.xpoint = double(zs.xpoint);
    else
        [vp,sp,sext,peri] = zgeo0(geo);
        Rsepa = [];
        Zsepa = [];
        zs.xpoint = zeros(size(cons.temps));
    end
    zs.sp = sp;
    zs.vp = vp;
    zs.sext = sext;
    zs.peri = peri;
end

if isfield(profli,'Rsepa') & isfield(profli,'Zsepa')
    Rsepa = profli.Rsepa;
    Zsepa = profli.Zsepa;
    
else
    Rsepa = [];
    Zsepa = [];
    if isfield(zs,'modeh')
        switch option.configuration
            case {0,1}
                zs.xpoint = zeros(size(zs.xpoint));
            case {2,3}
                zs.xpoint = zs.modeh;
                if (option.Kappa_xpoint > 0) || (option.delta_xpoint ~= 0) || ...
                   (option.R_HFS_xpoint > 0) || (option.R_LFS_xpoint > 0)
                     if option.R_HFS_xpoint > 0
                        xflag_hfs = ((geo.R - geo.a) >= option.R_HFS_xpoint);
                     else
                        xflag_hfs = (ones(size(geo.R)) > 0);
                     end
                     if option.R_LFS_xpoint > 0
                        xflag_lfs = ((geo.R + geo.a) <= option.R_LFS_xpoint);
                     else
                        xflag_lfs = (ones(size(geo.R)) > 0);
                     end                                       
                     zs.xpoint = double(xflag_hfs & xflag_lfs & (geo.K >= option.Kappa_xpoint) & ...
                                       ((geo.d .* sign(option.delta_xpoint)) >= abs(option.delta_xpoint)));
                end
            case 4
                zs.xpoint = ones(size(zs.xpoint));
            otherwise
                error('option not implemented');
        end
    end
end


% masse effective
if mode == 0
    switch option.gaz
        case 1
            zs.meff    = ones(size(cons.temps));
        case 2
            zs.meff    = 2 .* ones(size(cons.temps));
        case 11
            zs.meff    = (11 .* cons.iso + 1) ./ (cons.iso + 1);            
        case {3,5}
            zs.meff    = (3 .* cons.iso + 2) ./ (cons.iso + 1);
        otherwise
            zs.meff    = 4 .* ones(size(cons.temps));
    end
else
    switch option.gaz
        case {1,2,3}
            zs.meff = ((zs.n1m - zs.nDm - zs.nTm) + 2 .* zs.nDm + 3 .* zs.nTm) ./ zs.n1m;
        case 5
            % in this He code for He3
            zs.meff = ((zs.n1m - zs.nDm - zs.nTm) + 2 .* zs.nDm + 3 .* (zs.nTm + zs.nhem)) ./ (zs.n1m + zs.nhem);
        case 11
            % in this case Tritium field code for boron 11 !
            zs.meff = ((zs.n1m - zs.nDm) + 2 .* zs.nDm + 11 .* zs.nTm) ./ (zs.n1m + zs.nTm);
        otherwise
            zs.meff    = 4 .* ones(size(cons.temps));
    end
end
% valeur initiale pour debuter la convergence
if mode == 0
   zs.pfus    = zeros(size(cons.temps));
   zs.ifus    = zeros(size(cons.temps));
   zs.xfus    = zeros(size(cons.temps));
   zs.jxfus   = zeros(size(cons.temps));
   zs.j0fus   = zeros(size(cons.temps));
   zs.salpha  = zeros(size(cons.temps));
   zs.salpha_he4 = zeros(size(cons.temps));
   zs.salpha_p   = zeros(size(cons.temps));
   zs.salpha_t   = zeros(size(cons.temps));
   zs.salpha_he3 = zeros(size(cons.temps));
   zs.salpha_n   = zeros(size(cons.temps));
   zs.RR      = 0.1 .* 2 .* pi .* geo.R ./ cons.ip;
   zs.pohm    = zs.RR .* (cons.ip./2) .^ 2;
   zs.iohm    = cons.ip./2;
   zs.vloop   = zs.pohm ./ zs.iohm;
   zs.pbrem   = zeros(size(cons.temps));
   zs.prad    = zeros(size(cons.temps));
   zs.pradsol = zeros(size(cons.temps));
   zs.pcyclo  = zeros(size(cons.temps));
   zs.zeff    = real(cons.zeff);
   switch option.gaz
       case 4
            zs.zeffsc  = 3 .* ones(size(cons.zeff));
       case 5
            zs.zeffsc  = 2.5 .* ones(size(cons.zeff)); 
       case 11
            zs.zeffsc  = 5 .* ones(size(cons.zeff));
       otherwise
            zs.zeffsc  = 2 .* ones(size(cons.zeff));
   end           
   zs.ip      = cons.ip;
   zs.ifwcd   = zeros(size(cons.temps));
   zs.ieccd   = zeros(size(cons.temps));
   zs.inbicd  = zeros(size(cons.temps));
   zs.iboot   = zeros(size(cons.temps));
   zs.asser   = zeros(size(cons.temps));
   zs.ane     = ones(size(cons.temps));
   zs.nem     = cons.nbar;
   zs.nim     = 0.8 .* cons.nbar;
   zs.nebord  = 1e17 * ones(size(cons.temps));
   zs.nhe4m   = cons.nbar .* option.frhe0;
   zs.nDm     = cons.nbar ./ sqrt(2);
   zs.n1m     = cons.nbar ./ sqrt(2);
   zs.nhem    = cons.nbar ./ 3;
   zs.nbar    = cons.nbar;
   zs.qa      = 3.2 .* ones(size(cons.temps));
   zs.q95     = 3 .* ones(size(cons.temps));
   zs.qmin    = ones(size(cons.temps));
   zs.q0      = ones(size(cons.temps));
   zs.ate     = ones(size(cons.temps));
   zs.tem     = 1e3 * ones(size(cons.temps));
   zs.te0     = 2e3 * ones(size(cons.temps));
   zs.tebord  = option.min_te_LCFS * ones(size(cons.temps));
   zs.tite    = 0.5 .*ones(size(cons.temps));
   zs.betap   = 0.3 .* ones(size(cons.temps));
   zs.betaptot= 0.3 .* ones(size(cons.temps));
   zs.piqj    = 2 .* ones(size(cons.temps));
   zs.w       = ones(size(cons.temps));
   zs.wth     = ones(size(cons.temps));
   zs.hitb    = ones(size(cons.temps));
   zs.aitb    = ones(size(cons.temps));
   zs.hmhd    = ones(size(cons.temps));
   zs.ecrit_nbi      = option.einj .* ones(size(cons.temps));
   zs.ecrit_he       = 1e6 .* ones(size(cons.temps)); 
   zs.ecrit_he4_DHe3 = 1e6 .* ones(size(cons.temps));
   zs.ecrit_p_DHe3   = 1e6 .* ones(size(cons.temps));
   zs.ecrit_he3_DDn  = 1e6 .* ones(size(cons.temps));
   zs.ecrit_t_DDp    = 1e6 .* ones(size(cons.temps));
   zs.ecrit_p_DDp    = 1e6 .* ones(size(cons.temps));
   zs.ecrit_he4_DT   = 1e6 .* ones(size(cons.temps)); 
   zs.li      = cons.li;
   zs.tauip    = geo.R ./ 4;
   zs.xnbi    = zeros(size(cons.temps));
   zs.xeccd    = zeros(size(cons.temps));
   zs.wrad    = zeros(size(cons.temps));
   zs.d0       = zeros(size(cons.temps));
   zs.frnbi    = ones(size(cons.temps));
   if option.nb_nbi == 2
	  zs.frnbi    =  zs.frnbi + sqrt(-1) .* ones(size(cons.temps));
   end
   zs.firstorb_nbi    = zeros(size(cons.temps));
   zs.wrad     = zeros(size(cons.temps));
   zs.frloss_icrh = zeros(size(cons.temps));
   zs.fwcorr    = ones(size(cons.temps));
   zs.modeh    = zeros(size(cons.temps));
   zs.zmszl    = ones(size(cons.temps));
   zs.plhrip    = zeros(size(cons.temps));
   zs.n0a    = zeros(size(cons.temps));
   zs.rres    = geo.R;
   zs.xres    = zeros(size(cons.temps));
   zs.rm      = geo.a .* sqrt(geo.K);
   zs.drmdt   = zeros(size(cons.temps));
   zs.efficiency    = 5e18 .* ones(size(cons.temps));
   zs.fracmino   = 0.25 .* ones(size(cons.temps));
   zs.ppedmax   = Inf .* ones(size(cons.temps));
   zs.pioniz    = zeros(size(cons.temps));
   zs.pioniz_i  = zeros(size(cons.temps));
   zs.irun    = zeros(size(cons.temps));
   zs.difcurconv = option.morphing .* ones(size(cons.temps));
   zs.disrup    = zeros(size(cons.temps));
   zs.harm       = ones(size(cons.temps));
   zs.nmino      = 1e13 .* ones(size(cons.temps));
   zs.einj_nbi_icrh = option.einj * ones(size(cons.temps));
   zs.pnbi_icrh     = real(cons.pnbi);
   zs.ialign     = ones(size(cons.temps));
   zs.frac_pellet    = zeros(size(cons.temps));
   zs.indice_inv    = zeros(size(cons.temps));
   zs.sn0fr      = zeros(size(cons.temps));
   zs.wrot      = zeros(size(cons.temps));
   zs.nimpm     = ones(size(cons.temps));
   zs.ane_actual = zs.ane;
   zs.dsol = geo.a ./ 100;
   % breakdown variable 
   zs.eddy_current  = zeros(size(cons.temps));
   zs.flux_edge_cor = zeros(size(cons.temps));
end

%save('sat2','option','cons','geo','zs','amorti','profli');

% effet des pertes d'injection
zs.picrh = cons.picrh .* (1 - zs.frloss_icrh);
zs.plh   = cons.plh;
zs.pnbi  = real(cons.pnbi) .* real(zs.frnbi) + sqrt(-1) .* imag(cons.pnbi) .* imag(zs.frnbi);
zs.pecrh = cons.pecrh;


% DENSITE NATURELLE
% calcul de la densite naturelle en utilisant
% ref : K. Borarass et al, NF 42 (2002) p 1251-1256
% ref : F. Miskane  and X. Garbet PoF (2000) vol 7,10 p 4197-
% ref : J. Bucalossi et al, 30th EPS 
% ref : L. Horton et al PPCF 38 (1196) A269-A280 
%
if (mode == 1)

     nbar_cons_mem = cons.nbar;
     switch option.natural
     case 1
	 cons.nbar = max(zs.nbar_nat,cons.nbar);
     case 2
	 cons.nbar = zs.nbar_nat;
     end


     if (option.berror > 0)
     
		% new variable for ITPA scaling ITPA20
		if isfield(profli,'psi') && isfield(profli,'vpr')
		  psin  = (profli.psi - profli.psi(:,1) * ones(size(profli.xli))) ./  ...
			  ((profli.psi(:,end) - profli.psi(:,1)) * ones(size(profli.xli)));
		  d95    = tsplinet(psin,profli.dx,0.95*ones(size(psin,1),1));
		  Ka95prof = cumtrapz(profli.xli,profli.vpr,2) ./ max(eps,2*pi^2.*profli.Raxe .* (geo.a * profli.xli) .^ 2); 
		  Ka95   = tsplinet(psin,Ka95prof,0.95*ones(size(psin,1),1));
		else
		  Ka95 = zs.vp ./ (2*pi^2.*geo.R.*geo.a.^2);
		  d95  = geo.d;
		end
     
                tau_bohm = max(1e-6,(zs.sp ./ pi) .* geo.b0 ./ max(option.min_te_LCFS,zs.tem) .* 2 ./ (1 + zs.zeff .* zs.tite));
                [void1,void2,tauthl_ref] = ...
                zscale0(zs.nbar./ 1e20,zs.ip ./ 1e6 ,geo.b0,zs.ploss./ 1e6, ...
                        max(zs.ploss,zs.pin - zs.dwdt) ./ 1e6,geo.R,geo.K,geo.a, zs.meff, zs.zeff,zs.vp,zs.sp, ...
	                zs.q95,zs.sext,option.scaling,option.l2hscaling,zs.pohm./1e6,zs.nem,zs.tite,zs.nim ./ zs.nem, ...
                    tau_bohm,(zs.prad + zs.pbrem + zs.pcyclo) ./ 1e6,geo.d,Ka95,d95,cons.temps);
                %
		[void1,fact_confinement,void_ieddy,ne_ini,void4,prf_ini,void_flux_bord_cor,void_dt,ip_c_break,x_ioniz,gas_puff_equi,f_wth] = ...
		            z0taue_burnthrough(option,geo,cons,zs,tauthl_ref, ...
		                               zs.eddy_current,zs.flux_edge_cor);
		if (option.neasser >= 1)
		    cons.nbar = max(1e13,min(zs.nbar_nat,ne_ini)) .* double(fact_confinement == 0) + cons.nbar .* double(fact_confinement ~= 0);
		else
		    cons.nbar = max(1e13,min(zs.nbar_nat,ne_ini)) .* (1 - fact_confinement) + cons.nbar .* fact_confinement;
		end
		if any(gas_puff_mem > 0)
		    % modify gas puff reference only if this reference is used
		    gas_puff_mem = gas_puff_mem + gas_puff_equi;
		end 
		% current_eddy is  false in evolution mode, missign memory of inductive effect
		%ip_mem = cons.ip;
                %%%cons.ip = cons.ip - min(0.99 .* cons.ip,option.C_eddy .* zs.eddy_current);
                cons.ip = cons.ip .* (1 - min(0.99, option.C_eddy .* tanh(zs.eddy_current./ max(cons.ip,zs.ip + option.C_eddy .* zs.eddy_current))));
                % correction for Ip = 0 at breakdown
		if cons.ip(1) == 0
		        %cons.temps(fact_confinement == 0)
		        %cons.temps(cons.ip == 0)
			cons.ip(fact_confinement == 0) = ip_c_break(fact_confinement == 0);
			cons.ip(cons.ip == 0)          = ip_c_break(cons.ip == 0);
		        %disp('here')
		end
		% signe OK pour le controle en flux -> cons.flux = cons.flux - flux_bord_cor;
                %flux_mem  = cons.flux;
		cons.flux = cons.flux - zs.flux_edge_cor;
		
		% profile shape control
		zs.hitb(x_ioniz < (1 + 2 .* sqrt(eps))) =  1;
		zs.hitb(x_ioniz < 1) =  max(1 + sqrt(eps),zs.hitb(x_ioniz < 1));

		
%  		  figure(21);clf;subplot(2,1,1);
%  		  plot(cons.temps,cons.flux,'b',cons.temps,flux_mem,'g',cons.temps,zs.flux_edge_cor,'r');
%   		  if length(cons.temps) > 4
%  		    set(gca,'xlim',[0 3]);
%  		  end
%  		  subplot(2,1,2);
%   		  plot(cons.temps,cons.ip,'b',cons.temps,ip_mem,'g',cons.temps,zs.ip,'k:',cons.temps,zs.eddy_current,'r');
%   		  if length(cons.temps) > 4
%  		    set(gca,'xlim',[0 3]);
%  		  end
%  		  drawnow
                %keyboard
                %figure(21);hold on;plot(cons.temps,current_eddy,'b',cons.temps,cons.ip,'r');drawnow
                %% gaz puff due to plasma volume increase
                %gas_puff_mem = gas_puff_mem + (2/3) .* option.p_prefill .* max(0,z0dxdt(zs.vp,cons.temps)) ./ (300 .* 1.3806503e-23);
                %figure(21);clf;plot(cons.temps,gas_puff_mem);set(gca,'xlim',[0 3]);drawnow
     else
		fact_confinement = ones(size(cons.nbar));
     end


     % scaling
     ulh = 0.25;
     fh  = zs.modeh;
     zs.nbar_nat   = min(zs.negr,max(1e13, 1e20 .* (geo.b0 ./ zs.q95 ./ geo.R) .^ 0.6  .*  (ulh + (1 - ulh) .* fh) .* option.fnbar_nat));

     %figure(21);clf;
     %subplot(2,1,1) 
     %plot(cons.temps,nbar_cons_mem,'b',cons.temps,zs.negr,'g',cons.temps,zs.nbar_nat,'r',cons.temps,ne_ini,'m',cons.temps,cons.nbar,'k');
     %subplot(2,1,2) 
     %plot(cons.temps,fact_forme)
     %set(gca,'xlim',[0 1])
     %drawnow


else
    zs.nbar_nat = cons.nbar;
    fact_confinement = ones(size(cons.nbar));
    x_ioniz = ones(size(cons.nbar));
    f_wth   = ones(size(cons.nbar));
end     

% external data fraction of line of field confined:
if isappdata(0,'TAUE_EXP')
    tau_ext  = getappdata(0,'TAUE_EXP');
    if isfield(tau_ext,'temps') && ~isempty(tau_ext.temps)
        if isfield(tau_ext,'fraction_closed_lines') && ~isempty(tau_ext.fraction_closed_lines)
            fact_confinement = max(0,min(1,import_external_tau(tau_ext.temps,tau_ext.fraction_closed_lines,fact_confinement,cons.temps)));
        end
    end
end

% effet du ripple
% attention au bug -> ne pas lire le nom de la machine dans le workspace

if option.rip == 1 
    switch option.gaz
        case {2,3,5}
            nmaj = zs.nDm;
        case 4
            nmaj = zs.nhem;
        case 11
            nmaj = zs.n1m;
            fprintf('q');
        otherwise
            nmaj = zs.n1m;
    end
    if fwr == 1; disp('zripts0');end
    [zs.picrh,zs.plh,zs.priptherm,zs.einj_icrh,zs.einj_lh] = zripts0(cons.picrh .* (1 - zs.frloss_icrh),cons.plh,zs.nbar,geo.a, ...
        geo.R,geo.b0,zs.betaptot,zs.li, ...
        zs.ip,zs.tem,zs.ate,zs.nem,zs.ane,nmaj,option.cmin,option.mino,option.freq, ...
        zs.tite,zs.qa,zs.qmin,zs.meff,zs.tebord,zs.nebord,zs.vp,zs.rres);
    % securite si valeur hors limites
    zs.picrh = min(cons.picrh .* (1 - zs.frloss_icrh),max(0,zs.picrh));
    zs.plh = min(cons.plh,max(0,zs.plh));
    zs.plhrip  = zs.plh - cons.plh;
    if option.fwcd ~= 0
        zs.picrh = cons.picrh .* (1 - zs.frloss_icrh);
    else
        switch option.mino
            case {'He3','He4'}
                zs.picrh = cons.picrh .* (1 - zs.frloss_icrh);
        end
    end
    if isempty(zs.einj_icrh)
        zs.einj_icrh =2 .* zs.te0;
    end
elseif (mode == 0)  || (option.evolution == 0)
   zs.einj_icrh = 2 .* zs.te0;
   zs.priptherm = zeros(size(cons.temps));
   zs.einj_lh = zeros(size(cons.temps));
   zs.plhrip    = zeros(size(cons.temps));
end

% calcul des premieres informations 
if (option.berror > 0) && (mode ~= 0)
	zs.pfus  = fact_confinement .* zs.pfus;
	zs.picrh = fact_confinement .* zs.picrh;
	zs.plh   = fact_confinement .* zs.plh;
	zs.pnbi  = fact_confinement .* zs.pnbi;
        % toutes les ondes se comporte comme ecrh a ce stade
	zs.pecrh   = fact_confinement .* zs.pecrh + (1 - fact_confinement) .* prf_ini;
end
zs.pin   = zs.pfus + zs.pohm + zs.picrh + zs.plh + real(zs.pnbi) + imag(zs.pnbi) + zs.pecrh;
switch option.ploss_exp
    case 'no_prad'
        zs.ploss =  zs.pin;
    case 'max_power'
        if isfield(profli,'qe') &&isfield(profli,'qi')
            zs.ploss  =  max(option.pth_min,max(profli.qe + profli.qi,[],2));
        else
            zs.ploss =  zs.pin;
        end
    case 'max(pel)+max(pion)'
        if isfield(profli,'qe') &&isfield(profli,'qi')
            zs.ploss  =  max(option.pth_min,max(profli.qe,[],2) + max(profli.qi,[],2));
        else
            zs.ploss =  zs.pin;
        end
    otherwise
        zs.ploss =  zs.pin - (zs.pbrem + zs.pcyclo + option.fprad .* zs.prad + zs.pioniz);
end

% partie thermique et suprthermique
if (mode == 1)
    if fwr == 1; disp('zsupra0');end
    if option.evolution == 1
        mem_pfus_th  = zs.pfus_th;
        mem_esup_fus = zs.esup_fus;
        switch option.gaz
            case 5
                [zs.esup_fus_he4_DHe3, zs.pfus_th_he4_DHe3]  = ...
                    zsupra0(cons.temps,zs.pfus_he4_DHe3,zs.taus_he,zs.taue,zs.ecrit_he4_DHe3, fus.dhe3_he4, 4, 0, zs.esup_fus_he4_DHe3);
                [zs.esup_fus_p_DHe3,   zs.pfus_th_p_DHe3]    = ...
                    zsupra0(cons.temps,zs.pfus_p_DHe3,  zs.taus_he,zs.taue,zs.ecrit_p_DHe3,   fus.dhe3_p,   1, 0, zs.esup_fus_p_DHe3);
                [zs.esup_fus_he3_DDn,  zs.pfus_th_he3_DDn]    = ...
                    zsupra0(cons.temps,zs.pfus_he3_DDn, zs.taus_he,zs.taue,zs.ecrit_he3_DDn,  fus.ddn_he3,  3, 0, zs.esup_fus_he3_DDn);
                [zs.esup_fus_t_DDp,    zs.pfus_th_t_DDp]      = ...
                    zsupra0(cons.temps,zs.pfus_t_DDp,   zs.taus_he,zs.taue,zs.ecrit_t_DDp,    fus.ddp_t,    3, 0, zs.esup_fus_t_DDp);
                [zs.esup_fus_p_DDp,    zs.pfus_th_p_DDp]      = ...
                    zsupra0(cons.temps,zs.pfus_p_DDp,   zs.taus_he,zs.taue,zs.ecrit_p_DDp,    fus.ddp_p,    1, 0, zs.esup_fus_p_DDp);
                [zs.esup_fus_he4_DT,   zs.pfus_th_he4_DT]      = ...
                    zsupra0(cons.temps,zs.pfus_he4_DT,  zs.taus_he,zs.taue,zs.ecrit_he4_DT,   fus.dt_he4,   4, 0, zs.esup_fus_he4_DT);
                zs.esup_fus = zs.esup_fus_he4_DHe3 + zs.esup_fus_p_DHe3 + zs.esup_fus_he3_DDn + zs.esup_fus_t_DDp + zs.esup_fus_p_DDp + zs.esup_fus_he4_DT;
                zs.pfus_th  = zs.pfus_th_he4_DHe3 + zs.pfus_th_p_DHe3 + zs.pfus_th_he3_DDn + zs.pfus_th_t_DDp + zs.pfus_th_p_DDp + zs.pfus_th_he4_DT;
            case 11
                [esup_fus_2p3,pfus_th_2p3] = zsupra0(cons.temps,zs.pfus,zs.taus_he,zs.taue,zs.ecrit_he,2.3e6,4,0,zs.esup_fus);
                [esup_fus_4,pfus_th_4] = zsupra0(cons.temps,zs.pfus,zs.taus_he,zs.taue,zs.ecrit_he,4e6,4,0,zs.esup_fus);
                zs.esup_fus  = 2/3 .* esup_fus_2p3 + 1/3 .* esup_fus_4;
                zs.pfus_th   = 2/3 .* pfus_th_2p3 + 1/3 .* pfus_th_4;
            otherwise
                [zs.esup_fus,zs.pfus_th]   = zsupra0(cons.temps,zs.pfus,zs.taus_he,zs.taue,zs.ecrit_he,3.56e6,4,0,zs.esup_fus);
        end
        zs.pfus_th(1:2) = mem_pfus_th(1:2);
        zs.esup_fus(1:2) = mem_esup_fus(1:2);
        
        mem_pnbi_th  = zs.pnbi_th;
        mem_esup_nbi = zs.esup_nbi;
        [zs.esup_nbi,zs.pnbi_th]   = zsupra0(cons.temps,real(zs.pnbi),real(zs.taus_nbi),zs.taue,real(zs.ecrit_nbi), ...
            option.einj,zs.meff,0,real(zs.esup_nbi));
        if option.nb_nbi == 2
            [esup_nbi2,pnbi_th2]   = zsupra0(cons.temps,imag(zs.pnbi),imag(zs.taus_nbi),zs.taue,imag(zs.ecrit_nbi), ...
                option.einj2,zs.meff,0,imag(mem_esup_nbi));
            zs.esup_nbi = zs.esup_nbi + sqrt(-1) .* esup_nbi2;
            zs.pnbi_th  = zs.pnbi_th  + sqrt(-1) .* pnbi_th2;
        end
        zs.pnbi_th(1:2) = mem_pnbi_th(1:2);
        zs.esup_nbi(1:2) = mem_esup_nbi(1:2);
    else
        switch option.gaz
            case 5
                [zs.esup_fus_he4_DHe3, zs.pfus_th_he4_DHe3]  = ...
                    zsupra0(cons.temps,zs.pfus_he4_DHe3,zs.taus_he,zs.taue,zs.ecrit_he4_DHe3, fus.dhe3_he4, 4);
                [zs.esup_fus_p_DHe3,   zs.pfus_th_p_DHe3]    = ...
                    zsupra0(cons.temps,zs.pfus_p_DHe3,  zs.taus_he,zs.taue,zs.ecrit_p_DHe3,   fus.dhe3_p,   1);
                [zs.esup_fus_he3_DDn,  zs.pfus_th_he3_DDn]    = ...
                    zsupra0(cons.temps,zs.pfus_he3_DDn, zs.taus_he,zs.taue,zs.ecrit_he3_DDn,  fus.ddn_he3,  3);
                [zs.esup_fus_t_DDp,    zs.pfus_th_t_DDp]      = ...
                    zsupra0(cons.temps,zs.pfus_t_DDp,   zs.taus_he,zs.taue,zs.ecrit_t_DDp,    fus.ddp_t,    3);
                [zs.esup_fus_p_DDp,    zs.pfus_th_p_DDp]      = ...
                    zsupra0(cons.temps,zs.pfus_p_DDp,   zs.taus_he,zs.taue,zs.ecrit_p_DDp,    fus.ddp_p,    1);
                [zs.esup_fus_he4_DT,   zs.pfus_th_he4_DT]      = ...
                    zsupra0(cons.temps,zs.pfus_he4_DT,  zs.taus_he,zs.taue,zs.ecrit_he4_DT,   fus.dt_he4,   4);
                zs.esup_fus = zs.esup_fus_he4_DHe3 + zs.esup_fus_p_DHe3 + zs.esup_fus_he3_DDn + zs.esup_fus_t_DDp + zs.esup_fus_p_DDp + zs.esup_fus_he4_DT;
                zs.pfus_th  = zs.pfus_th_he4_DHe3 + zs.pfus_th_p_DHe3 + zs.pfus_th_he3_DDn + zs.pfus_th_t_DDp + zs.pfus_th_p_DDp + zs.pfus_th_he4_DT;
            case 11
                 [esup_fus_2p3,pfus_th_2p3] = zsupra0(cons.temps,zs.pfus,zs.taus_he,zs.taue,zs.ecrit_he,2.3e6,4);
                 [esup_fus_4,pfus_th_4]     = zsupra0(cons.temps,zs.pfus,zs.taus_he,zs.taue,zs.ecrit_he,4e6,4);
                 zs.esup_fus  = 2/3 .* esup_fus_2p3 + 1/3 .* esup_fus_4;
                 zs.pfus_th   = 2/3 .* pfus_th_2p3 + 1/3 .* pfus_th_4;
                                            
           otherwise
                [zs.esup_fus,zs.pfus_th]   = zsupra0(cons.temps,zs.pfus,zs.taus_he,zs.taue,zs.ecrit_he,3.56e6,4);
        end
        [zs.esup_nbi,zs.pnbi_th]   = zsupra0(cons.temps,real(zs.pnbi),real(zs.taus_nbi),zs.taue,real(zs.ecrit_nbi), ...
            option.einj,zs.meff);
        %[zs.esup_nbi,zs.pnbi_th]   = zsupra0(cons.temps,zs.pnbi,zs.taus_nbi,zs.taue,zs.ecrit_nbi,option.einj,zs.meff);
        if option.nb_nbi == 2
            [esup_nbi2,pnbi_th2]   = zsupra0(cons.temps,imag(zs.pnbi),imag(zs.taus_nbi),zs.taue,imag(zs.ecrit_nbi), ...
                option.einj2,zs.meff);
            
            zs.esup_nbi = zs.esup_nbi + sqrt(-1) .* esup_nbi2;
            zs.pnbi_th  = zs.pnbi_th  + sqrt(-1) .* pnbi_th2;
        end
    end
    if option.cmin == 0 || option.fwcd ~= 0
        einjfci = zs.te0;
        zs.esup_icrh = 0 .* zs.te0;
        zs.pion_icrh = 0 .* zs.te0;
        zs.taus_icrh = 0 .* zs.te0;
        zs.picrh_th  = zs.picrh;
        zs.pel_icrh  = zs.picrh;
        zs.einj_icrh = zs.te0;
        zs.ecrit_icrh = 1e3 .* zs.te0;
        frpar_icrh    = zeros(size(cons.temps));
        j_fast_mino   = zeros(length(cons.temps),length(profli.xli));
        i_fast_mino   = zeros(size(cons.temps));
        puissance_fw  = zeros(size(cons.temps));
    else
%         switch option.mino
%             case 'He3'
%                 ag = 3;
%             case 'T'
%                 ag = 3;
%             case 'He4'
%                 ag = 4;
%             case 'D'
%                 ag = 2;
%             case 'B'
%                 ag = 11;
%            otherwise
%                 ag = 1;
%         end
        if fwr == 1; disp('z0icrh');end
        switch option.icrh_model
            case 'Dumont-Vu'
                if option.evolution == 1
                    
                    mem_picrh_th  = zs.picrh_th;
                    mem_esup_icrh = zs.esup_icrh;
                    
                    [zs.picrh_th,zs.pel_icrh,zs.pion_icrh,zs.esup_icrh,zs.einj_icrh, ...
                        zs.ecrit_icrh,teff0,zs.taus_icrh,zs.frloss_icrh,zs.rres,zs.xres,zs.fracmino,zs.harm,zs.nmino, ...
                        Pe,PM,esupra_icrh_par,j_fast_mino,i_fast_mino,puissance_fw] = ...
                        z0icrh_trang(zs,geo,cons,option,profli);
                    frpar_icrh          = max(0,min(1,zs.esup_icrh ./ max(eps,esupra_icrh_par)));
                    zs.picrh_th(1:2)    = mem_picrh_th(1:2);
                    zs.esup_icrh(1:2)   = mem_esup_icrh(1:2);
                    profli.picrh        = Pe + PM;
                    profli.picrh_ion    = PM;
                    
                else
                    [zs.picrh_th,zs.pel_icrh,zs.pion_icrh,zs.esup_icrh,zs.einj_icrh, ...
                        zs.ecrit_icrh,teff0,zs.taus_icrh,zs.frloss_icrh,zs.rres,zs.xres,zs.fracmino,zs.harm,zs.nmino, ...
                        Pe,PM,esupra_icrh_par,j_fast_mino,i_fast_mino,puissance_fw] = ...
                        z0icrh_trang(zs,geo,cons,option,profli);
                    frpar_icrh          = max(0,min(1,zs.esup_icrh ./  max(eps,esupra_icrh_par)));
                    profli.picrh        = Pe + PM;
                    profli.picrh_ion    = PM;
                end
            otherwise
                if option.evolution == 1
                    
                    mem_picrh_th  = zs.picrh_th;
                    mem_esup_icrh = zs.esup_icrh;
                    
                    [zs.picrh_th,zs.pel_icrh,zs.pion_icrh,zs.esup_icrh,zs.einj_icrh, ...
                        zs.ecrit_icrh,teff0,zs.taus_icrh,zs.frloss_icrh,zs.rres,zs.xres,zs.fracmino,zs.harm,zs.nmino] = ...
                        z0icrh(zs,geo,cons,option,profli);
                    
                    zs.picrh_th(1:2) = mem_picrh_th(1:2);
                    zs.esup_icrh(1:2) = mem_esup_icrh(1:2);
                    
                else
                    [zs.picrh_th,zs.pel_icrh,zs.pion_icrh,zs.esup_icrh,zs.einj_icrh, ...
                        zs.ecrit_icrh,teff0,zs.taus_icrh,zs.frloss_icrh,zs.rres,zs.xres,zs.fracmino,zs.harm,zs.nmino] = ...
                        z0icrh(zs,geo,cons,option,profli);
                end
        end
    end
    if option.lhmode ~= 5
        % formule de F. Imbeaux
        we           = (4/3 * 1.602176462e-19) .* zs.nem .* (1 + zs.ane) .* zs.te0 ./ (zs.ate + zs.ane +1) .* zs.vp;
        nemi         = zs.nem .* (1+ zs.ane) .* (3/4) .^ zs.ane;
        temi         = zs.te0 .*  (3/4) .^ zs.ate;
        b0mi         = geo.b0 .* geo.R ./ (geo.R + geo.a/2);
        zs.esup_lh   = min(0.2,max(0,1.8e-6 .* ((b0mi .^ 2 ./ (4*pi*1e-7) ./ nemi./ (temi  .* 1.602176462e-19)) .^ (3/2) -43)))  ...
            .* we .* min(1,zs.plh ./(1+zs.ip)); % critere de saturation heursitique
    else
        zs.esup_lh   = zeros(size(zs.plh));
    end
    zs.plh_th    = zs.plh; % on neglige le temps de thermalisation
    
else
    zs.dwdt      = zeros(size(cons.temps));
    zs.dwthdt      = zeros(size(cons.temps));
    %
    zs.pfus_th   = zs.pfus;
    zs.picrh_th  = zs.picrh;
    zs.pel_icrh  = 0.5 .* zs.picrh;
    zs.pion_icrh = 0.5 .* zs.picrh;
    zs.pel_icrh  = zs.picrh;
    zs.pnbi_th   = zs.pnbi;
    zs.plh_th    = zs.plh;
    zs.esup_fus  = zeros(size(cons.temps));
    zs.esup_nbi  = zeros(size(cons.temps));
    zs.esup_icrh = zeros(size(cons.temps));
    zs.esup_lh   = zeros(size(cons.temps));
    zs.ecrit_he  = 1e5 .* ones(size(cons.temps));
    zs.ecrit_nbi = option.einj ./ 2 .* ones(size(cons.temps));
    zs.ecrit_icrh = 1e3 .* zs.te0;
    zs.nhem       = zeros(size(cons.temps));
    zs.taus_icrh  = eps .* ones(size(cons.temps));
    frpar_icrh    = zeros(size(cons.temps));
    j_fast_mino   = zeros(length(cons.temps),length(profli.xli));
    i_fast_mino   = zeros(size(cons.temps));
    puissance_fw  = zeros(size(cons.temps));
    if option.gaz == 5
        zs.pfus_th_he4_DHe3 = zs.pfus;
        zs.pfus_th_p_DHe3   = zs.pfus;
        zs.pfus_th_he3_DDn  = zs.pfus;
        zs.pfus_th_t_DDp    = zs.pfus;
        zs.pfus_th_p_DDp    = zs.pfus;
        zs.pfus_th_he4_DT   = zs.pfus;
        zs.esup_fus_he4_DHe3 = zeros(size(cons.temps));
        zs.esup_fus_p_DHe3   = zeros(size(cons.temps));
        zs.esup_fus_he3_DDn  = zeros(size(cons.temps));
        zs.esup_fus_t_DDp    = zeros(size(cons.temps));
        zs.esup_fus_p_DDp    = zeros(size(cons.temps));
        zs.esup_fus_he4_DT   = zeros(size(cons.temps));
    end
end

switch option.ploss_exp
    case 'no_prad'
        zs.pth   = (zs.pfus_th + zs.pohm + zs.picrh_th + zs.plh_th + real(zs.pnbi_th) + imag(zs.pnbi_th) + zs.pecrh);
    case 'max_power'
        if isfield(profli,'qe') &&isfield(profli,'qi')
            zs.pth  =  max(option.pth_min,max(profli.qe + profli.qi,[],2));
        else
            zs.pth   = (zs.pfus_th + zs.pohm + zs.picrh_th + zs.plh_th + real(zs.pnbi_th) + imag(zs.pnbi_th) + zs.pecrh);
        end
    case 'max(pel)+max(pion)'
        if isfield(profli,'qe') &&isfield(profli,'qi')
            zs.pth  =  max(eps,max(profli.qe,[],2) + max(profli.qi,[],2));
        else
            zs.pth   = (zs.pfus_th + zs.pohm + zs.picrh_th + zs.plh_th + real(zs.pnbi_th) + imag(zs.pnbi_th) + zs.pecrh);
        end
    otherwise
        zs.pth   = (zs.pfus_th + zs.pohm + zs.picrh_th + zs.plh_th + real(zs.pnbi_th) + imag(zs.pnbi_th) + zs.pecrh) - ...
            (zs.pbrem + zs.pcyclo + option.fprad .* zs.prad + zs.pioniz);
end
zs.pth  =  max(option.pth_min,zs.pth);

pion_icrh = zeros(size(zs.picrh_th));

% limite pour dwthdt
plim          = max(option.pth_min,min(zs.vp .* 1e6,max(0.8 .* zs.pth,0.2 .* abs(zs.pin))));

% puissance ion/electron
if isfield(profli,'nbishape_ion')
	rapnbi1      = trapz(profli.xli,profli.vpr .* real(profli.nbishape_ion),2) ./ ...
		       max(1,trapz(profli.xli,profli.vpr .* (real(profli.nbishape_el) + real(profli.nbishape_ion)),2));
	rapnbi2      = trapz(profli.xli,profli.vpr .* imag(profli.nbishape_ion),2) ./ ...
		       max(1,trapz(profli.xli,profli.vpr .* (imag(profli.nbishape_el) + imag(profli.nbishape_ion)),2));
	zs.pion_nbi = rapnbi1 .* real(zs.pnbi_th) + sqrt(-1) .* rapnbi2 .* imag(zs.pnbi_th);
else
	zs.pion_nbi = zfract0(real(zs.ecrit_nbi),option.einj) .* real(zs.pnbi_th) + ...
                      sqrt(-1) .* zfract0(imag(zs.ecrit_nbi),option.einj2) .* imag(zs.pnbi_th);
end
%figure(21);clf;plot(cons.temps,real(zs.pion_nbi),cons.temps,imag(zs.pion_nbi));drawnow

zs.pel_nbi  = max(0,real(zs.pnbi_th) - real(zs.pion_nbi)) + sqrt(-1) .* max(0,imag(zs.pnbi_th) - imag(zs.pion_nbi));
if isfield(profli,'tep') && isfield(profli,'salf')
    switch option.gaz
        case 5
            [ecrit_alpha,zs.pel_fus,zs.pion_fus,pfus_shape_el,profli.pfus_ion] = z0fusDHe3_el_ion(zs,profli,option,zu1,zu2);
            profli.pfus  = pfus_shape_el + profli.pfus_ion;
            zs.ecrit_he4_DHe3 = ecrit_alpha.he4_DHe3;
            zs.ecrit_p_DHe3   = ecrit_alpha.p_DHe3;
            zs.ecrit_he3_DDn  = ecrit_alpha.he3_DDn;
            zs.ecrit_t_DDp    = ecrit_alpha.t_DDp;
            zs.ecrit_p_DDp    = ecrit_alpha.p_DDp;
            zs.ecrit_he4_DT   = ecrit_alpha.he4_DT; 
        case 11
            % the function z0fus_el_ion has been updated for pB11
            [zs.ecrit_he,zs.pel_fus,zs.pion_fus,pfus_shape_el,profli.pfus_ion] = z0fus_el_ion(zs,profli,option,zu1,zu2);
            profli.pfus  = pfus_shape_el + profli.pfus_ion;
       otherwise
            [zs.ecrit_he,zs.pel_fus,zs.pion_fus,pfus_shape_el,profli.pfus_ion] = z0fus_el_ion(zs,profli,option,zu1,zu2);
            profli.pfus  = pfus_shape_el + profli.pfus_ion;
    end
else
    switch option.gaz
        case 5
            zs.pion_fus  =   min(1 - eps,(option.alpha_channeling + zfract0(zs.ecrit_he4_DHe3, fus.dhe3_he4))) .* zs.pfus_th_he4_DHe3 ...
                           + min(1 - eps,(option.alpha_channeling + zfract0(zs.ecrit_p_DHe3  , fus.dhe3_p  ))) .* zs.pfus_th_p_DHe3   ...
                           + min(1 - eps,(option.alpha_channeling + zfract0(zs.ecrit_he3_DDn , fus.ddn_he3 ))) .* zs.pfus_th_he3_DDn  ...
                           + min(1 - eps,(option.alpha_channeling + zfract0(zs.ecrit_t_DDp   , fus.ddp_t   ))) .* zs.pfus_th_t_DDp    ...
                           + min(1 - eps,(option.alpha_channeling + zfract0(zs.ecrit_p_DDp   , fus.ddp_p   ))) .* zs.pfus_th_p_DDp    ... 
                           + min(1 - eps,(option.alpha_channeling + zfract0(zs.ecrit_he4_DT  , fus.dt_he4  ))) .* zs.pfus_th_he4_DT  ;
            zs.pel_fus  = max(0,zs.pfus_th - zs.pion_fus);
        case 11
            zs.pion_fus = min(1 - eps,(option.alpha_channeling + (2./3) .* zfract0(zs.ecrit_he,2.3e6) + (1/3) .* zfract0(zs.ecrit_he,4e6))) .* zs.pfus_th;
            zs.pel_fus  = max(0,zs.pfus_th - zs.pion_fus);
       otherwise
            zs.pion_fus = min(1 - eps,(option.alpha_channeling + zfract0(zs.ecrit_he,3.56e6))) .* zs.pfus_th;
            zs.pel_fus  = max(0,zs.pfus_th - zs.pion_fus);
    end
end
if fwr == 1; disp('zfract0');end
if mode == 1
   if isappdata(0,'ICRH_SHAPE_EXP') & isfield(profli,'qjli');
	% il faut aussi modifier dans zero1t
	icrhexp = getappdata(0,'ICRH_SHAPE_EXP');

	fpicrh_el  = max(1,interp1_ex(icrhexp.temps,icrhexp.pel,zs.temps,'nearest','extrap'));
	fpicrh_el  = pchip(icrhexp.x,fpicrh_el,profli.xli);
   	indnok = find(any(~isfinite(fpicrh_el),2));
   	indok  = find(all(isfinite(fpicrh_el),2));
   	fpicrh_el(indnok,:) = ones(length(indnok),1) * mean(fpicrh_el(indok,:),1);

	fpicrh_ion  = max(1,interp1_ex(icrhexp.temps,icrhexp.pion,zs.temps,'nearest','extrap'));
	fpicrh_ion  = pchip(icrhexp.x,fpicrh_ion,profli.xli);
   	indnok = find(any(~isfinite(fpicrh_ion),2));
   	indok  = find(all(isfinite(fpicrh_ion),2));
   	fpicrh_ion(indnok,:) = ones(length(indnok),1) * mean(fpicrh_ion(indok,:),1);

	fpicrh_fw  = max(1,interp1_ex(icrhexp.temps,icrhexp.pfw,zs.temps,'nearest','extrap'));
	fpicrh_fw  = pchip(icrhexp.x,fpicrh_fw,profli.xli);
   	indnok = find(any(~isfinite(fpicrh_fw),2));
   	indok  = find(all(isfinite(fpicrh_fw),2));
   	fpicrh_fw(indnok,:) = ones(length(indnok),1) * mean(fpicrh_fw(indok,:),1);


	rapicrh    = trapz(profli.xli,profli.vpr .* fpicrh_ion,2) ./ ...
		        max(1,trapz(profli.xli,profli.vpr .* (fpicrh_el + fpicrh_ion + fpicrh_fw),2));
	zs.pion_icrh = rapicrh .* zs.picrh_th;
        zs.pion  = zs.pion_fus + zs.pion_icrh  + real(zs.pion_nbi) + imag(zs.pion_nbi);
        zs.pel_icrh  =   zs.picrh_th - zs.pion_icrh;

   elseif option.cmin == 0
      zs.pion  = zs.pion_fus + 0.5 .* zs.picrh_th .* (option.fwcd == 0) + real(zs.pion_nbi) + imag(zs.pion_nbi);
      % 0.5 car pas d'information
      zs.pion_icrh =   0.5 .* zs.picrh_th .* (option.fwcd == 0);
      zs.pel_icrh  =   zs.picrh_th - zs.pion_icrh;

   else
      zs.pion  = zs.pion_fus + zs.pion_icrh .* (option.fwcd == 0) + real(zs.pion_nbi) +imag(zs.pion_nbi);
      zs.pion_icrh =   zs.pion_icrh .* (option.fwcd == 0);
      zs.pel_icrh  =   zs.picrh_th - zs.pion_icrh;

   end
else
   zs.pion  = zs.pion_fus  + (2/3) .* zs.picrh_th .* (option.fwcd == 0) + real(zs.pion_nbi) +imag(zs.pion_nbi) ;
   zs.pion_icrh =   (2/3) .* zs.picrh_th .* (option.fwcd == 0);
   zs.pel_icrh  =   zs.picrh_th - zs.pion_icrh;
end

% complement au electrons
if option.cx_ion ~= 0
    zs.pion = zs.pion - zs.pioniz_i;
    zs.pel   = max(0,(zs.pfus_th + zs.pohm + zs.picrh_th + zs.plh_th + real(zs.pnbi_th) + imag(zs.pnbi_th) + zs.pecrh) ...
            		- zs.pion - (zs.pbrem + zs.pcyclo + option.fprad .* zs.prad) - zs.pioniz);
else
    zs.pel   = max(0,(zs.pfus_th + zs.pohm + zs.picrh_th + zs.plh_th + real(zs.pnbi_th) + imag(zs.pnbi_th) + zs.pecrh) ...
            		- zs.pion - (zs.pbrem + zs.pcyclo + option.fprad .* zs.prad) - zs.pioniz);
end
% calcul du terme correctif important pendant les rampup
% effet sur l'aimant toroidal
% (P. Helander, Collisional transport in magnetized plasma, Cambridge University Press ,2002,  p 164)
if isfield(profli,'fdia')

    grr              = z0dxdt(profli.Raxe,cons.temps) ./ profli.Raxe - ...
        ((zs.drmdt .* zs.rm) * ones(size(profli.xli))) ./ profli.Raxe .^ 2 ./ 2;
    inte             = (profli.fdia - (profli.fdia(:,end) * ones(size(profli.xli)))) .* profli.r2i .* ...
        (z0dxdt(profli.fdia,cons.temps) - profli.fdia .* grr);

    if option.evolution == 1
        switch option.dwdt_method
            case 'freebie'
                inte = ones(size(cons.temps)) * mean(inte,1);
        end
    end
    switch option.dwdt_method
        case 'none'
            zs.dwmagtordt = zeros(size(cons.temps));
            zs.wmagtor    = zeros(size(cons.temps));
        otherwise
            zs.dwmagtordt    = trapz(profli.xli,profli.vpr .* inte,2) ./ (4 * pi * 1e-7);
            zs.wmagtor       = cumtrapz(cons.temps,zs.dwmagtordt);
    end
else
    zs.dwmagtordt = zeros(size(cons.temps));
    zs.wmagtor    = zeros(size(cons.temps));
end

% geometrical effect of pth (plasma compression)
if isfield(profli,'vpr') && isfield(profli,'rmx') && isfield(profli,'Raxe')  && (option.adiabatic == 1)
  pth_mem = zs.pth;
  % hange in plasma minor radius and shapping
  dshapedt = (z0dxdt(profli.vpr(:,end),cons.temps) ./ profli.vpr(:,end) -  ...
              z0dxdt(profli.rmx(:,end),cons.temps) ./ profli.rmx(:,end));
  % correction problem derivative at first time
  dshapedt(1) = dshapedt(2);
  % 
  %zs.pth   = zs.pth - (z0dxdt(profli.vpr(:,end),cons.temps) ./ profli.vpr(:,end) -  ...
  %                     z0dxdt(profli.rmx(:,end),cons.temps) ./ profli.rmx(:,end))  .* zs.wth;
  % tanh -> this is just a perturbation, must be less than 1
  %zs.pth   = zs.pth - tanh(dshapedt)  .* zs.wth;                    
  % change in major radius (invariant concervation during fast motion)
  dRaxedt = z0dxdt(profli.Raxe,cons.temps);
  Radiat   = trapz(profli.xli,profli.vpr .* dRaxedt ./ profli.Raxe,2) ./ trapz(profli.xli,profli.vpr,2);
  % correction problem derivative at first time
  Radiat(1) = Radiat(2);
  %figure(21);clf;plot(cons.temps,Radiat,'b',cons.temps,dshapedt,'r',cons.temps,-dshapedt -  4/3 .* Radiat,'k');drawnow;                 
  % tanh -> this is just a perturbation, must be less than 1
  %zs.pth  = zs.pth - 4/3 .* tanh(Radiat) .* fact_confinement  .* zs.wth;
  padiabatic  = tanh(dshapedt + 4/3 .* Radiat .* fact_confinement)  .* zs.wth;
  zs.pth  = zs.pth - padiabatic;
  %figure(21);clf;plot(cons.temps,pth_mem,'b',cons.temps,zs.pth,'r',cons.temps,- padiabatic,'k');drawnow;                 
end

% correction de dwdt
if (option.transitoire == 1) && (mode ~=0)
	switch option.dwdt_method
	case 'none'
	    zs.ploss = max(option.pth_min,real(zs.ploss));
	    zs.pth   = max(option.pth_min,real(zs.pth));
	otherwise 
	    zs.pth   = max(option.pth_min,zs.pth - zs.dwthdt);
	    zs.ploss = max(option.pth_min,zs.ploss - zs.dwdt);
	end
else
	zs.ploss = max(option.pth_min,real(zs.ploss));
	zs.pth   = max(option.pth_min,real(zs.pth));
end

% securite
zs.pion  = max(option.pth_min,zs.pion);
zs.pel   = max(option.pth_min,zs.pel);
zs.pin   = max(option.pth_min,real(zs.pin));

% appel du module de loi d'echelle
if fwr == 1; disp('zscale0');end
%
tau_bohm = max(1e-6,(zs.sp ./ pi) .* geo.b0 ./ max(option.min_te_LCFS,zs.tem) .* 2 ./ (1 + zs.zeff .* zs.tite));
% nombre de rayon de larmor dans le plasma (ou en nombre de longueur de Debye, on prend le plus grand, c'est une securite a faible densite)	
%rho_i  = max(4.57e-3 .* sqrt(zs.meff) .* sqrt(max(option.min_te_LCFS,zs.tem .* zs.tite) ./ 1e3) ./ geo.b0, ...
%		     2.35e5  .* sqrt(max(option.min_te_LCFS,zs.tem) ./ max(1e13,zs.nem)));
% effet des pieges
%  if isfield(profli,'ftrap') 
%  	%% INVERSE ASPECT RATIO = r/R (PROFILE)
%  	epsi = profli.epsi;
%  
%  	%% PARTICLE MASS (KG) AND CHARGE (C) AND VELOCITY (M/S)
%  	%mass = ((((zs.n1m - zs.nDm - zs.nTm) + 2 .* zs.nDm  + 3 .* zs.nTm) ./ zs.n1m) * ones(size(profli.xli)))  .* profli.n1p + ...
%          %		4 .* profli.nhep  + ((option.zimp+1) * 2 + option.rimp .* (option.zmax+1) .* 2) .* profli.nzp;
%  	%mass   = mass ./ profli.nip;
%  	
%  	%% ORBIT WIDTH (M)
%  	%btot       = sqrt(profli.fdia .^2 .* profli.r2i + profli.bpol .^ 2);
%  	%rho_im = rho_i
%  	%rho_ip       = max(4.57e-3 .* sqrt(mass) .* sqrt(max(option.min_te_LCFS,profli.tip) ./ 1e3) ./ btot, ...
%  	%	     2.35e5  .* sqrt(max(option.min_te_LCFS,profli.tep) ./ max(1e13,profli.nep)));
%  	dr_banana  = (epsi).^(-0.5) .* profli.qjli;
%  	dr_banana(:,1) = 0;
%  	% calcul de l'effet moyen des pieges
%   	%figure(52);clf;plot(profli.xli,(1 - profli.ftrap) .*  rho_i + profli.ftrap .* max((dr_banana - rho_i) ./ 2,rho_i));drawnow;
%  	rho_i =  rho_i .* trapz(profli.xli,((1 - profli.ftrap) + profli.ftrap .* max(dr_banana ./ 2,1)),2);
%  	%rho_i = trapz(profli.xli,profli.vpr .* rho_i,2) ./ zs.vp;
%  	%figre(51);clf;plot(cons.temps,rho_i,cons.temps,rho_im);drawnow;
%  end

%nblarmor = max(1,zs.rm ./ rho_i);
% figure(151);clf;plot(cons.temps, tau_bohm) ;drawnow 
%

% ae = ni/ne;
ae = zs.nim ./ max(1e13,zs.nem); 

% new variable for ITPA scaling ITPA20
if isfield(profli,'psi') && isfield(profli,'vpr')
  psin  = (profli.psi - profli.psi(:,1) * ones(size(profli.xli))) ./  ...
          ((profli.psi(:,end) - profli.psi(:,1)) * ones(size(profli.xli)));
  d95    = tsplinet(psin,profli.dx,0.95*ones(size(psin,1),1));
  Ka95prof = cumtrapz(profli.xli,profli.vpr,2) ./ max(eps,2*pi^2.*profli.Raxe .* (geo.a * profli.xli) .^ 2); 
  Ka95   = tsplinet(psin,Ka95prof,0.95*ones(size(psin,1),1));
else
  Ka95 = zs.vp ./ (2*pi^2.*geo.R.*geo.a.^2);
  d95  = geo.d;
end

[zs.wrlw,zs.plossl2h,zs.tauthl,zs.tauh,zs.tauhe_l,zs.tauhe_h,tau_na,tau_nas,tau_ped_sc] = ...
         zscale0(zs.nbar./ 1e20,zs.ip ./ 1e6 ,geo.b0,zs.ploss./ 1e6, ...
         max(zs.ploss,zs.pin - zs.dwdt) ./ 1e6,geo.R,geo.K,geo.a, zs.meff, zs.zeff,zs.vp,zs.sp, ...
	 zs.q95,zs.sext,option.scaling,option.l2hscaling,zs.pohm./1e6,zs.nem,zs.tite,ae,tau_bohm, ...
     (zs.prad + zs.pbrem + zs.pcyclo) ./ 1e6,geo.d,Ka95,d95,cons.temps);

% using scaling for pped to get difference between H-mode and L-mode 
switch  option.scaling
    case 19
        switch option.ploss_exp
            case {'max_power','max(pel)+max(pion)'}
                p_scaling_pped =  zs.ploss;
            otherwise
                p_scaling_pped =  zs.pin;
        end
        %
        iob = (2/5) .* geo.R .* zs.ip ./ 1e6  .* geo.K .^ 2 ./  geo.a .^ 2 ./ (1 + geo.K .^ 2);
        pped   = max(0,4.53138e3 .* (min(0.5,abs(geo.d)) + 0.0034) .^ 0.435509 .* (p_scaling_pped  ./ 1e6).^ 0.121836 .* iob .^ 1.51649);
        wped      = abs(option.fpped)  .*  3/2 .*  pped .* zs.vp;
        zs.tauh   = zs.tauthl + wped ./ max(option.pth_min,zs.pth);
end

% difference de seuil selon la configuration
zs.plossl2h = (zs.xpoint  + (1 - zs.xpoint) .* option.fpl2h_lim) .* zs.plossl2h;

% effet du changement de gaz sur le seuil
if option.pl2h_mass_charge == 1
    % reference: R. Behn et all, PPCF 57 (2015) 025007 and references within.
    %plossl2h_mem = zs.plossl2h;
    switch option.gaz
        case {5,11}
             zs.plossl2h = zs.plossl2h .* (zs.meff ./ 2) .^(-1.1) .* zs.zeff .^ 1.6;
       otherwise
            zs.plossl2h = zs.plossl2h .* (zs.meff ./ 2) .^(-1.1) .* (1 + double(option.gaz == 4)) .^ 1.6;
            %figure(21);clf;plot(cons.temps,zs.plossl2h,'b',cons.temps,plossl2h_mem,'r');drawnow
    end
end

if (mode ~= 0) &&  isfield(profli,'nip') && (option.scaling == 9)
    % mode J_pol = 0 , cf. J. Garcia and G. Giruzzi, PRL 104 p 2055003 (2010)
    % calcul du profil de pression
    % securite si inversion du courant (limite de pression)
    % indok    = find(all(profli.jli(:,1:end-1) > 0,2));
    pprim_jg = profli.jli .* profli.ri .* (profli.qjli > abs(option.qdds));
    warning off
    dpdx_jg  = pprim_jg .* pdederive(profli.xli,profli.psi,0,2,2,1); % attention au facteur +/- 2*pi !
    warning on
    dpdx_jg(:,1) = 0;
    pth_jg   = cumtrapz(profli.xli,dpdx_jg,2);
    pth_jg   = pth_jg - pth_jg(:,end) * ones(size(profli.xli));
    % pression de bord
    w0  =1.602176462e-19 .* (3/2) .* (profli.nip(:,end) .* profli.tip(:,end) + profli.nep(:,end) .* profli.tep(:,end)) .* zs.vp;
    % limitation de l'energie maximum
    % l'equilibre tient compte de la pression totale + rotation (diminue la pression cinetique percue)
    wcore    = min(3/2  .* trapz(profli.xli,profli.vpr .* pth_jg,2) - max(0,zs.w - zs.wth) + 2 .* zs.wrot,10 .* zs.wbp);
    %
    %whmode    = wcore + wcore;
    whmode    = w0 + max(0,wcore);
    %wlmode    = wcore;
    wlmode    = w0 + max(0,wcore - zs.wbp);
    % temps de confinement
    %pcor = z0dxdt(zs.li,cons.temps) .* (4 .* pi .* 1e-7) .* geo.R .* zs.ip .^ 2 ./ 4;
    pcor = abs(zs.dwmagtordt + zs.dwbpdt);
    %pcor = zs.dwmagtordt + zs.dwbpdt;
    %pcor = cat(1,pcor(2:end),pcor(end));
    %pcor(pcor < 0) = 0;
    pcor(1) = 0;
    tauh   = max(tau_bohm,whmode ./ max(option.pth_min,zs.pth + pcor));
    tauthl = max(tau_bohm,wlmode ./ max(option.pth_min,zs.pth + pcor));
    %zs.tauh(indok)  = tauh(indok);
    %zs.tauthl(indok)  = tauthl(indok);
    zs.tauh  = tauh;
    zs.tauthl  = tauthl;
    %figure(21);clf;plot(cons.temps,zs.wth,'b',cons.temps,whmode,'r',cons.temps,wlmode,'m',cons.temps,w0,'g',cons.temps,tau_bohm .* zs.pth,'c');drawnow

elseif (option.fpped < 0) && (mode ~= 0) && isfield(profli,'nip')
    
    % based on alpha limit from s-alpha balooning stability limit
    % magnetic shear
    vt   = ones(size(profli.psi,1),1);
    ve   = ones(1,size(profli.psi,2));
    % psin = (profli.psi - profli.psi(:,1) * ve) ./ ((profli.psi(:,end) - profli.psi(:,1)) * ve);
    % shear = pdederive(profli.xli,profli.qjli,1,2,2,1) ./ pdederive(profli.xli,psin,1,2,2,1)./ profli.qjli .* (vt * profli.xli) ./ psin;
    shear = pdederive(profli.xli,profli.qjli,1,2,2,1) ./ profli.qjli .* (vt * profli.xli) ;
    shear(:,1) = 0;
    shear_limited = max(-0.6,shear);
    % alpha critique from Freidberg book
    %lookup table :
    alpha_crit = salpha_lookup(shear_limited);
    switch option.coef_shape
        case {'stiff_limited','stiff_limited+neo'}
            %alpha_crit = alpha_crit;
        otherwise
            alpha_crit = abs(option.fstiff) .* alpha_crit;
    end
    % compute pressure gradient (from dP/dr to dP/dx -> drop one r)
    dptotdx_crit = -alpha_crit .* profli.bpol .^ 2 ./ (2 .* (4*pi*1e-7) .*  profli.epsi);
    ptot_crit    = cumtrapz(profli.xli,dptotdx_crit,2);
    % flip integral
    ptot_crit    = ptot_crit - ptot_crit(:,end) * ve;
    % add integration constant -> the contribution is added later
    %x ptot_crit    = ptot_crit + profli.ptot(:,end) * ve;
    %  compute energy content
    wtot_crit    = (3/2) .* cumtrapz(profli.xli,ptot_crit .* profli.vpr,2);
    wtot_core    = (3/2) .* cumtrapz(profli.xli,(ptot_crit - ptot_crit(:,end-1) * ve) .* profli.vpr,2);
    % correction from non thermal content assuming all fast particles are inside the core plasma
    wcore_alpha  =  max(1,wtot_core(:,end-1) - (zs.w - zs.wth));
    wlmode_alpha =  max(1,wtot_crit(:,end)   - (zs.w - zs.wth));
    switch option.coef_shape
        case {'alpha','alpha+neo'}
           wcore  = wcore_alpha;
           wlmode = wlmode_alpha; 
         otherwise
            % nouvelle implementation du mode stiff a piedestal fixe
            % nombre de rayon de larmor dans le plasma (ou en nombre de longueur de Debye, on prend le plus grand, c'est une securite a faible densite)
            % on enleve la dependance en ti pour simuler le mode ion chaud
            % la dependance en Ti/Te n'apparait pas dans les etudes addimentionelles
            % du fait de la valeur de la constante et l'absence de dependance en Ti/Te, la dependnance en rho_e semble plus pertinente
            % De plus la dependance en masse des scaling est faible, par contre la dependance en nombre de particules est importante (effet de dilution)
            % la dependance en volume integre la dependance en shaping
            % la largeur banane semble plus pertinante puisque il y a une dependance en Ip
            % de plus cela integre la dependance en rapport d'aspect
            % La reponse est meilleur en utilisant rhomax, ce qui est logique puisque l'invariant est relie au flux toroidal
            % pour eviter le problemes de decouplage ? basse densite :
            %teff    = (profli.tip  + profli.zeff .* profli.tep) ./ (1 + profli.zeff); % selon la vitesse du son
            %rho_e   = 1.07e-4 .*  sqrt(max(option.min_te_LCFS,teff) ./ 1e3) ./ sqrt((profli.fdia .* profli.ri) .^ 2 + profli.bpol .^2);
            rho_e   = 1.07e-4 .*  sqrt(max(option.min_te_LCFS,profli.tep) ./ 1e3) ./ sqrt((profli.fdia .* profli.ri) .^ 2 + profli.bpol .^2);
            rho_b   = rho_e .* profli.qjli ./ sqrt(profli.epsi);
            rho_pot = (rho_e .^ 2 .* profli.qjli .^ 2 .* profli.Raxe) .^ (1/3);
            rho_b   = rho_b .* (rho_b < profli.rmx) + rho_pot .* (rho_b >= profli.rmx);
            rho_b(:,1) = rho_pot(:,1);
            % debyes
            rho_d  = 2.35e5  .* sqrt(max(option.min_te_LCFS,profli.tep) ./ max(1e13,profli.nep)./ 1e3);
            %figure(51);clf;plot(cons.temps,rho_b,'r',cons.temps,rho_d,'b');drawnow
            rho_bd = max(rho_b,rho_d);
            % cas du modele stiff2 avec particule passantes et piegees
            if option.fstiff < 0
                if option.isotope_stiff ~= 0
                    % reference is deuterium
                    rho_bd = profli.ftrap .* rho_bd .*  ((zs.meff / 2) * ones(size(profli.xli))) .^ (-option.isotope_stiff) + (1 - profli.ftrap) .* rho_e;
                else
                    rho_bd = profli.ftrap .* rho_bd + (1 - profli.ftrap) .* rho_e;
                end
            end
            %figure(21);clf;plot(profli.xli,rho_pot,'b',profli.xli,rho_bd,'r');drawnow;
            %les 2 derniers points compte pour la collisionalite en mode l mais pas en mode h
            nblarmor_h =  0.95 .* zs.rm ./ max(eps,trapz(profli.xli(1:end-1),rho_bd(:,1:end-1),2));
            nblarmor_l =  zs.rm ./ max(eps,trapz(profli.xli,rho_bd,2));
            %nblarmor =  (profli.Raxe(:,end) + geo.a -  profli.Raxe(:,1))./ max(eps,trapz(profli.xli,rho_bd,2));
            %nblarmor =  profli.Raxe(:,1) ./ max(eps,trapz(profli.xli,rho_bd,2));
            %nblarmor =  geo.a ./ max(eps,trapz(profli.xli,rho_bd,2));
            wcore = abs(option.fstiff)  .* 1.602176462e-19 .* nblarmor_h .* ...
                (3 / 2) .* trapz(profli.xli(1:end-1),profli.vpr(:,1:end-1) .* ...
                (profli.nip(:,1:end-1) + profli.nep(:,1:end-1)) .* ...
                (ones(size(profli.vpr,1),1) * ((profli.xli(end-1) - profli.xli(1:end-1)) ./profli.xli(end-1))),2);
            wlmode = abs(option.fstiff)  .* 1.602176462e-19 .* nblarmor_l .* ...
                (3 / 2) .* trapz(profli.xli,profli.vpr .* (profli.nip + profli.nep) .* ...
                (ones(size(profli.vpr,1),1) * (profli.xli(end) - profli.xli(1:end))),2);
            % account for scaling selection
            switch option.coef_shape
                case {'stiff_limited','stiff_limited+neo'}
                  wcore  = min(wcore,wcore_alpha);
                  wlmode = min(wlmode,wlmode_alpha);
            end
    end
    % selon la puissance utilis?e dans le scaling du temps de confinement
    switch option.ploss_exp
        case {'max_power','max(pel)+max(pion)'}
		p_scaling_pped =  zs.ploss;
	otherwise
		p_scaling_pped =  zs.pin;
	end
	if option.usepped_scl == 1
	    switch option.scaling
	    case 13
		  % use fit of pedestal database augmented with prediction from MHD stability studies for ITER and DEMO1
		  % Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
		  % Pped_{sc} = 4.53138 * (d + 0.0034) ^ 0.435509 * ptot ^ 0.121836 * iob ^ 1.51649
		  iob = (2/5) .* geo.R .* zs.ip ./ 1e6  .* geo.K .^ 2 ./  geo.a .^ 2 ./ (1 + geo.K .^ 2);
		  pped   = max(0,4.53138e3 .* (min(0.5,abs(geo.d)) + 0.0034) .^ 0.435509 .* (p_scaling_pped  ./ 1e6).^ 0.121836 .* iob .^ 1.51649); 
                  wped      = abs(option.fpped)  .*  3/2 .*  pped .* zs.vp;

	    otherwise
		  pped   = max(0,(2/3) .* zs.modeh .* tau_ped_sc .* zs.pth ./ zs.vp);
                  wped      = abs(option.fpped)  .*  3/2 .*  pped .* zs.vp;
	    end
	elseif option.usepped_scl == 2     
        switch option.scaling
            case 13
                % use fit of pedestal database augmented with prediction from MHD stability studies for ITER and DEMO1
                % Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
                % Pped_{sc} = 4.53138 * (d + 0.0034) ^ 0.435509 * ptot ^ 0.121836 * iob ^ 1.51649
                iob = (2/5) .* geo.R .* zs.ip ./ 1e6 .* geo.K .^ 2 ./  geo.a .^ 2 ./ (1 + geo.K .^ 2);
                pped   = max(0,4.53138e3 .* (min(0.5,abs(geo.d)) + 0.0034) .^ 0.435509 .* (p_scaling_pped ./ 1e6) .^ 0.121836 .* iob .^ 1.51649);
                wped      = abs(option.fpped)  .*  3/2 .*  pped .* zs.vp;
                
            otherwise
                pped   = min(zs.pped,max(0,(2/3) .* zs.modeh .* tau_ped_sc .* zs.pth ./ zs.vp));
                wped      = abs(option.fpped)  .*  3/2 .*  pped .* zs.vp;
        end
    elseif option.usepped_scl == 3
        
        % use fit of pedestal database augmented with prediction from MHD stability studies for ITER and DEMO1
        % Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
        % Pped_{sc} = 4.53138 * (d + 0.0034) ^ 0.435509 * ptot ^ 0.121836 * iob ^ 1.51649
        iob = (2/5) .* geo.R .* zs.ip ./ 1e6 .* geo.K .^ 2 ./  geo.a .^ 2 ./ (1 + geo.K .^ 2);
        pped   = max(0,4.53138e3 .* (min(0.5,abs(geo.d)) + 0.0034) .^ 0.435509 .* (p_scaling_pped ./ 1e6) .^ 0.121836 .* iob .^ 1.51649);
        wped      = abs(option.fpped)  .*  3/2 .*  pped .* zs.vp;
        
    elseif option.usepped_scl == 4
        
        % use fit of pedestal database augmented with prediction from MHD stability studies for ITER and DEMO1
        % Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
        % Pped_{sc} = 4.53138 * (d + 0.0034) ^ 0.435509 * ptot ^ 0.121836 * iob ^ 1.51649
        iob = (2/5) .* geo.R .* zs.ip ./ 1e6 .* geo.K .^ 2 ./  geo.a .^ 2 ./ (1 + geo.K .^ 2);
        pped   =  min(zs.pped,max(0,4.53138e3 .* (min(0.5,abs(geo.d)) + 0.0034) .^ 0.435509 .* (p_scaling_pped ./ 1e6) .^ 0.121836 .* iob .^ 1.51649));
        wped      = abs(option.fpped)  .*  3/2 .*  pped .* zs.vp;
    else
        % model reference :F.D. Halpern et al, PoP 15 (2008) p 062505
        [wped,pped] = z0pped_onjun(zs,profli,geo,cons,option);
        wped        = abs(option.fpped)  .*  wped;
        pped        = abs(option.fpped)  .*  pped;
    end
                    
	%
	w0  =1.602176462e-19 .* (3/2) .* (profli.nip(:,end) .* profli.tip(:,end) + profli.nep(:,end) .* profli.tep(:,end)) .* zs.vp;
	%	
	whmode     = wcore + wped;
	wlmode     = (wlmode + w0) .* double(~zs.modeh) + double(zs.modeh) .* wcore;
    
    switch option.coef_shape
        case {'alpha','alpha+neo'}
            zs.tauh   = max(tau_bohm,whmode ./ max(option.pth_min,zs.pin));
            zs.tauthl = max(tau_bohm,wlmode ./ max(option.pth_min,zs.pin));
        otherwise
            zs.tauh   = max(tau_bohm,whmode ./ max(option.pth_min,zs.pth));
            zs.tauthl = max(tau_bohm,wlmode ./ max(option.pth_min,zs.pth));
    end

%figure(21);plot(cons.temps,wped,cons.temps,wcore,cons.temps,wlmode,cons.temps,w0);drawnow

elseif option.dilution == 1

	% effet en helium	
	% ref : F. Ryter et al, NF 49 (2009) pp 062003
	% ref : D.C. McDonald et al PPCF 46 (2004) pp 419-534
	% effet de la diminution du nombre d'ions
	zscale = (zs.nem + zs.n1m + zs.nhem + zs.nimpm .* (1 + option.rimp) + zs.nwm) ./ (2 .* zs.nem);
	zs.tauh     = zscale .* zs.tauh;
	zs.tauthl   = zscale .* zs.tauthl;
end


% cas breakdown
if (option.berror > 0) && (mode ~= 0)
    fconf_taup = zs.tauthl;
    if option.evolution == 1
        [zs.tauthl,void_f,eddy_current,void_nec,void_tref,void_prf,flux_edge_cor] = z0taue_burnthrough(option,geo,cons,zs,zs.tauthl, ...
            zs.eddy_current(2),zs.flux_edge_cor(2));
        
        zs.eddy_current(3:end)   = eddy_current(3:end);
        zs.flux_edge_cor(3:end)  = flux_edge_cor(3:end);
    else
        [zs.tauthl,void_f,zs.eddy_current,void_nec,void_tref,void_prf,zs.flux_edge_cor] = z0taue_burnthrough(option,geo,cons,zs,zs.tauthl);
    end
    fconf_taup = min(1,zs.tauthl ./ fconf_taup);
end

% limitation avec la loi neo alcator
% heating, confinement and extrapolation to reactors 
% Energy for Future Centuries: Prospects for Fusion Power as a Future Energy Source
% J. Ongena, G. Van Oost
% Fusion Science and Technology / Volume 57 / Number 2T / February 2010 / Pages 3-15
% Proceedings of the Ninth Carolus Magnus Summer School on Plasma and Fusion Energy Physics
if mode ~= 0
    %figure(21);clf;plot(cons.temps,zs.tauh,'r',cons.temps,zs.tauthl,'b',cons.temps,tau_na,'g',cons.temps,tau_nas,'-.g');
    %set(gca,'xlim',[0 2]);
    %drawnow
    switch option.tau_limitation
    case 'On'
        zs.tauh     = min(tau_na,zs.tauh);
        zs.tauthl   = min(tau_na,zs.tauthl);
    case 'Saturate'
        zs.tauh     = min(tau_nas,zs.tauh);
        zs.tauthl   = min(tau_nas,zs.tauthl);
    end
end

% calcul du mode H/L
zs.plossl2h = max(1e4,real(zs.plossl2h));
switch option.plhthr
case '2*pion'
	zs.plhthr = max(option.pth_min,max(min(2 .* zs.pion,zs.pin - zs.dwdt),zs.ploss./3));
case 'P_LCFS'
	zs.plhthr = max(option.pth_min,min(zs.pin - zs.prad - zs.pbrem -zs.pcyclo - zs.pioniz,zs.pin - zs.dwdt));
case 'P_LCFS_dwdt'
	zs.plhthr = max(option.pth_min,zs.pin - zs.prad - zs.pbrem -zs.pcyclo - zs.pioniz - zs.dwdt);
otherwise
	zs.plhthr = max(option.pth_min,max(min(zs.ploss,zs.pin - zs.dwdt),zs.ploss./3));
end
if mode == 0
   zs.modeh  = (zs.plhthr >= (zs.plossl2h + 1e6 .*  option.l2hmul));
elseif (option.l2hscaling < 0) && isfield(profli,'web')
        % D. Testa Nucl. Fusion 46 (2006) 562???579
        gitg     = sqrt(1.602176462e-19 .* (profli.tip + profli.zeff .* profli.tep)./1.6726485e-27 ./  ...
                   (zs.meff * ones(size(profli.xli)))) ./ ...
                   (profli.rmx(:,end) * ones(size(profli.xli)));
	srotg    = max(exp(eps),profli.web ./ gitg);
	zs.modeh = (srotg(:,end-1) > abs(option.l2hscaling));
elseif (option.l2hscaling == 6) && isfield(profli,'fdia')
       % transition si jtheta = 0
       %ref :E. R Solano, NF52 (2012) p 114017
       fdiap    = diff(profli.fdia(:,end-2:end-1),1,2) .* diff(profli.fdia(:,end-1:end),1,2) ;
       % le changement de signe doit etre suffisant, mais il faut rejeter le bruit.
       fdiap_th = mean(abs(fdiap)/pi) .* ones(size(fdiap)); 
       %figure(21);clf;plot(fdiap,'b');hold on;plot(-fdiap_th,'r');drawnow
       zs.modeh = (fdiap < -fdiap_th);
else
   if length(zs.modeh) == 1
      hdecal     = 0;
   else
      hdecal     = cat(1,0,zs.modeh(1:end-1));
   end
   % hysteresis
   physteresis   = option.hysteresis .* zs.pin + (1 - option.hysteresis) .* zs.plhthr;
   %zs.modeh      = (zs.plhthr >= (zs.plossl2h + 1e6 .*option.l2hmul)) | (hdecal & (zs.pin >= (zs.plossl2h + 1e6 .* option.l2hmul)));
   zs.modeh      = (zs.plhthr >= (zs.plossl2h + 1e6 .*option.l2hmul)) | (hdecal & (physteresis >= (zs.plossl2h + 1e6 .* option.l2hmul)));
end
if option.modeh == 0
   zs.modeh      = zeros(size(zs.modeh));
elseif  option.modeh == 2
   zs.modeh      = ones(size(zs.modeh));
end
if isfinite(option.ton_modeh) && isfinite(option.toff_modeh)
   	zs.modeh(cons.temps < option.ton_modeh) = 0;
   	zs.modeh(cons.temps > option.toff_modeh) = 0;
elseif isfinite(option.ton_modeh)
   	zs.modeh(cons.temps < option.ton_modeh) = 0;
elseif isfinite(option.toff_modeh)
   	zs.modeh(cons.temps > option.toff_modeh) = 0;
end

if (option.berror > 0) && (mode ~= 0)
	zs.modeh(fact_confinement < 1) = 0;
end


% not useful but add numerical convergence problem
%  if (option.disrup > 0) && (mode ~= 0) && 0
%      if option.disrup == 1
%        zs.modeh(zs.disrup ~= 0) = 0;
%      elseif option.disrup == 2
%        ind_disrup = find(zs.disrup ~= 0,1);
%        zs.modeh(ind_disrup:end) = 0;    
%      end
%  end

zs.modeh = double(zs.modeh);

% transition entre taul et tauh
if option.l2hslope > 0
	f = (zs.plhthr - (zs.plossl2h + 1e6 .*  option.l2hmul)) ./ max(option.pth_min,zs.plossl2h + 1e6 .*  option.l2hmul);
	f = min(1,max(0,f ./ option.l2hslope));
	%figure(151);clf;plot(f);drawnow
	zs.tauh = zs.tauthl + f .* option.l2hslope .* (zs.tauh - zs.tauthl) +  (1 - option.l2hslope) .* (zs.tauh - zs.tauthl);
elseif (option.l2hslope < 0) && (mode ~= 0) &&  isfield(profli,'nip')
	% effet de la collisionalite sur le piedestal
	% Ce model ne rend pas compte des experiences sur JET
	% reference : M.Z. Tokar , On Greenwald density limit in H-mode, PoP 16 (2009) , p 020704
	%lbc      =  15.2 - 0.5 .* log(profli.nep(:,end - 1) ./ 1e20) + log(profli.tep(:,end - 1)./1e3);
	%tei      =  1.09e16 .* (profli.tep(:,end - 1)./1e3) .^ (3/2) ./ profli.nip(:,end - 1) ./ lbc ./ profli.zeff(:,end - 1) .^ 2;
	%lambda_c =  sqrt(2 .*  1.602176462e-19  .* profli.tep(:,end - 1) ./ 9.10938188e-31) .* tei;
	%lambda_c =  3 .* profli.tep(:,end - 1) .^ 2 ./ (2 .* sqrt(2 .* pi) .* lbc .* 1.602176462e-19 .^ 4 .* profli.zeff(:,end - 1) .* profli.nip(:,end - 1)); 
	%psr      =  profli.qjli(:,end - 1) .* profli.Raxe(:,end - 1) ./ max(eps,lambda_c);
	%factlh   =  (psr - 0.222) ./ 0.222;
	%f        =  1 - max(0,tanh(abs(option.l2hslope) .* factlh));
	%zs.tauh  = (1 - f) .* zs.tauthl + f .* zs.tauh;
	%figure(151);clf;plot(cons.temps,f,cons.temps,zs.modeh);drawnow;
	
	% utilisation du scaling lineaire en fraction de Greewald de Cordey
	% reference : J.G Cordey 28th EPS 2001 P 3.11
	% la dependance lineaire est observee dans P. Dumortier PPCF 44 2002 p 1845-1861, ainsi que le dependnace en traingularite.
	% la dependance en piquage est observee dans A.W. Leonard 27th EPS 2000
	fhh = 0.84 + 0.18 .* max(0,geo.d)  - 0.13 .*  zs.nbar ./ max(1e13,zs.negr)  +  0.51 .* max(0,(zs.nbar ./ max(1e13,profli.nep(:,end-1)) - 1));
	%tauhmem = zs.tauh;
	zs.tauh  = max(zs.tauthl , ((zs.modeh .* fhh) + (~zs.modeh)) .* zs.tauh);
	%figure(151);clf;plot(cons.temps,((zs.modeh .* fhh) + (~zs.modeh)) );drawnow;
	%figure(151);clf;plot(cons.temps,zs.tauh,cons.temps,tauhmem);drawnow;
end

% si tauhemul est donne
if option.tauhemul > 0
	zs.tauhe_l  = option.tauhemul .* zs.tauthl;
	zs.tauhe_h  = option.tauhemul .* zs.tauh;
end
% securite 
if (option.berror > 0) && (mode ~= 0)
	zs.tauthl   = min(1000,max(max(eps,1e-6 .* fact_confinement),real(zs.tauthl)));
else
	zs.tauthl   = min(1000,max(1e-6,real(zs.tauthl)));
end
zs.tauh     = min(1000,max(2e-6,real(zs.tauh)));
zs.tauhe_l  = min(1000,max(2e-6,real(zs.tauhe_l)));
zs.tauhe_h  = min(1000,max(2e-6,real(zs.tauhe_h)));


% effet de li su le confinement
if option.HH_li > 0
	HH_li = (zs.li ./ option.HH_li) .^ (2/3);
else
	HH_li = ones(size(cons.hmore));
end 
% traingularity effect
if option.HH_delta ~= 0
    HH_delta = ((1 + geo.d) ./ (1 + 0.35)) .^ option.HH_delta;
else
  	HH_delta = ones(size(cons.hmore));
end
%figure(21);plot(cons.temps,HH_delta);drawnow

% effet de l'injection de gaz
if option.HH_gas_puff > 0
    HH_gas_puff = max(eps,1 - tanh(option.HH_gas_puff .* zs.pioniz ./ max(option.pth_min,zs.pin)));
    %figure(21);clf;plot(cons.temps,HH_gas_puff);drawnow
else
    HH_gas_puff = ones(size(cons.hmore));
    %disp('here')
end
% confinement
taue_itb    = (zs.hitb-1) .* zs.tauthl;
tauhe_itb   = (zs.hitb-1) .* zs.tauhe_l;
zs.taue     = zs.hmhd .* HH_delta .* HH_li .* HH_gas_puff .* cons.hmore .* (zs.modeh .* zs.tauh  + (~zs.modeh) .* zs.tauthl)  + taue_itb;
if (option.berror > 0) && (mode ~= 0)
    zs.taue     = 	max(max(2e-6 .* fact_confinement,eps),zs.taue);
else
    zs.taue     = 	max(2e-6,zs.taue);
end


if (option.tauhemul < 0) && (mode ~= 0) && (option.Recycling < 1)
        % modele qui prend en compte le confinement reel et le recyclage dans le divertor
	% ref Stangeby section 6.7 
        % ref originale : D. Reiter et al, PPCF vol 33 (1991) p 1579-1600
	zs.tauhe    = abs(option.tauhemul) .* zs.taue + option.Recycling ./ (1 - option.Recycling) .* zs.taup;

%          figure(21);
%          clf
%          plot(cons.temps,zs.tauhe ./ zs.taue,'r',cons.temps,zs.taup ./ zs.taue,'b');
%          set(gca,'ylim',[0 30]);
%          title(sprintf('r/(1-r) = %g',option.Recycling ./ (1 - option.Recycling)));
%          drawnow
else
	zs.tauhe    = max(1e-3,zs.hmhd .* HH_delta .* HH_li .* cons.hmore .* (zs.modeh .* zs.tauhe_h + (~zs.modeh) .* zs.tauhe_l) + tauhe_itb);
end
% vieux scaling
zs.wrlw   = (zs.taue ./ zs.tauthl) .* zs.wrlw;


% betan pour la mhd (anti oscillation estime avec le jeux de donneees pre
if (option.transitoire == 1) && (mode ~=0)
	betan0 = zs.w .* (1.6.*pi./3) .* geo.a ./ zs.vp ./ geo.b0 ./ zs.ip;
else
	tau0     = max(1e-3,zs.hitb .* cons.hmore .* (zs.modeh .* zs.tauh    + (~zs.modeh) .* zs.tauthl));
	w0       = tau0 .*  zs.pth;
	betan0   = w0.* (1.6.*pi./3) .* geo .a ./ zs.vp ./ geo.b0 ./ zs.ip;
end


if isappdata(0,'TE_EXP') && isappdata(0,'TI_EXP') && isfield(profli,'nip') && isfield(profli,'tip') &&  ...
        isfield(option,'exp_shape') && (option.exp_shape == 0)
    zs.wth  = max(eps,1.602176462e-19 .* (3/2) .* trapz(profli.xli,(profli.tep .* profli.nep + profli.tip .* profli.nip) .* profli.vpr,2));
    taue_before = zs.taue;
    zs.taue = max(1e-6,zs.wth ./ max(option.pth_min,zs.pth));
    if any((taue_before / 10) > zs.taue) || any((taue_before * 10) < zs.taue)
        fprintf('t');
        indbad = find(((taue_before / 10) > zs.taue) | ((taue_before * 10) < zs.taue));
        zs.taue(indbad) = taue_before(indbad);
    end
elseif option.evolution == 1
    % premiere estimation sans transitoire
    zs.wth(end-1:end)    = zs.taue(end-1:end)  .* zs.pth(end-1:end);
    if (option.disrup ~= 0) && (mode ~= 0)
        wth_loc = zs.wth;
        if option.disrup == 1
            wth_loc(zs.disrup ~= 0) = zs.vp(zs.disrup ~= 0);
        elseif option.disrup == 2
            wth_loc(end-1:end) = zs.vp(end-1:end);
        elseif option.disrup == 3
            if isfield(profli,'qe') &&isfield(profli,'qi')
                mask = (profli.qe + profli.qi) >=0;
                fconf_disrup = (trapz(profli.xli, profli.vpr .* mask,2) ./ trapz(profli.xli, profli.vpr ,2));
                fconf_disrup(zs.disrup == 0) = 1;
                fconf_disrup = max(0,min(1, 1 + cumtrapz(cons.temps,fconf_disrup - 1)));
                fconf_disrup(1:2) = 1;
                zs.wth       = fconf_disrup .* zs.wth;
                zs.taue      = fconf_disrup .* zs.taue;
                %figure(21);clf;plot(fconf_disrup);drawnow
             end
        end
        zs.wth(end-1:end)       = wth_loc(end-1:end);
        zs.taue(end-1:end)      = zs.wth(end-1:end) ./ max(eps,zs.pth(end-1:end));
    end
else
    % premiere estimation sans transitoire
    zs.wth    = zs.taue  .* zs.pth;
    %figure(21);clf;plot(cons.temps,zs.wth);
    if (option.disrup ~= 0) && (mode ~= 0)
        wth_mem_ = zs.wth;
        taue_mem_ = zs.taue;
        if option.disrup == 1
            zs.wth(zs.disrup ~= 0) = zs.vp(zs.disrup ~= 0);
            zs.taue(zs.disrup ~= 0) = zs.wth(zs.disrup ~= 0) ./ max(eps,zs.pth(zs.disrup ~= 0));
        elseif option.disrup == 2
            ind_disrup = find(zs.disrup ~= 0,1);
            zs.wth(ind_disrup:end) = zs.vp(ind_disrup:end);
            zs.taue(ind_disrup:end) = zs.wth(ind_disrup:end) ./ max(eps,zs.pth(ind_disrup:end));
        elseif option.disrup == 3
            if isfield(profli,'qe') &&isfield(profli,'qi')
                mask = (profli.qe + profli.qi) >=0;
                fconf_disrup = (trapz(profli.xli, profli.vpr .* mask,2) ./ trapz(profli.xli, profli.vpr ,2));
                fconf_disrup(zs.disrup == 0) = 1;
                fconf_disrup = max(0,min(1, 1 + cumtrapz(cons.temps,fconf_disrup - 1)));
                zs.wth       = fconf_disrup .* zs.wth;
                zs.taue      = fconf_disrup .* zs.taue;
                %%figure(21);clf;plot(fconf_disrup);drawnow
             end
        end
        %  	    figure(21);clf
        %  	    subplot(2,1,1);
        %  	    plot(cons.temps,zs.wth,'.r',cons.temps,wth_mem_,'b');
        %  	    subplot(2,1,2);
        %  	    plot(cons.temps,zs.taue,'.r',cons.temps,taue_mem_,'b');
        %  	    drawnow
        %
        
    end
end




% fast electron due to runaway electrons addedt to LH terms
if option.runaway ~= 0
    % assuming v// ~ c
    esup_run = zs.irun .* phys.me .* phys.c ./ phys.e;
else
    esup_run = zeros(size(cons.temps));
end
%figure(21);plot(cons.temps,esup_run);drawnow

zs.w      = zs.wth + (zs.esup_fus + zs.esup_icrh + real(zs.esup_nbi) + imag(zs.esup_nbi) + zs.esup_lh + esup_run);
zs.betan  = zs.w .* (1.6.*pi./3) .* geo.a ./ zs.vp ./ geo.b0 ./ zs.ip;

% nan et imag
z0dnanimag
%save('sat3','option','cons','geo','zs','amorti','profli');

% pression en haut du piedestal
%if (option.fpped >= 0) | (mode == 0) | ~isfield(profli,'nip')
zs.pped   = max(0,(2/3) .* zs.modeh .* (zs.tauh - zs.tauthl) .* zs.pth ./ zs.vp);
% pour le cas iter std, il faut que la moitier de la difference dans le piedestal
switch option.scaling
case {0,1,6,7}
	if (option.fpped >= 0)
		zs.pped = zs.pped ./ 2;
	end
end

% selon la puissance utilis?e dans le scaling du temps de confinement
switch option.ploss_exp
case {'max_power','max(pel)+max(pion)'}
	p_scaling_pped =  zs.ploss;
otherwise
	p_scaling_pped =  zs.pin;
end

if option.usepped_scl == 1
    switch option.scaling
        case 13
            % use fit of pedestal database augmented with prediction from MHD stability studies for ITER and DEMO1
            % Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
            % Pped_{sc} = 4.53138 * (d + 0.0034) ^ 0.435509 * ptot ^ 0.121836 * iob ^ 1.51649
            iob = (2/5) .* geo.R .* zs.ip ./ 1e6  .* geo.K .^ 2 ./  geo.a .^ 2 ./ (1 + geo.K .^ 2);
            zs.pped   = max(0,4.53138e3 .* (min(0.5,abs(geo.d)) + 0.0034) .^ 0.435509 .* (p_scaling_pped  ./ 1e6).^ 0.121836 .* iob .^ 1.51649);
            
        otherwise
            zs.pped   = max(0,(2/3) .* zs.modeh .* tau_ped_sc .* zs.pth ./ zs.vp);
    end
elseif option.usepped_scl == 2
    switch option.scaling
        case 13
            % use fit of pedestal database augmented with prediction from MHD stability studies for ITER and DEMO1
            % Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
            % Pped_{sc} = 4.53138 * (d + 0.0034) ^ 0.435509 * ptot ^ 0.121836 * iob ^ 1.51649
            iob = (2/5) .* geo.R .* zs.ip ./ 1e6 .* geo.K .^ 2 ./  geo.a .^ 2 ./ (1 + geo.K .^ 2);
            zs.pped   = max(0,4.53138e3 .* (min(0.5,abs(geo.d)) + 0.0034) .^ 0.435509 .* (p_scaling_pped ./ 1e6) .^ 0.121836 .* iob .^ 1.51649);
            
        otherwise
            zs.pped   = min(zs.pped,max(0,(2/3) .* zs.modeh .* tau_ped_sc .* zs.pth ./ zs.vp));
    end
elseif option.usepped_scl == 3
    
    % use fit of pedestal database augmented with prediction from MHD stability studies for ITER and DEMO1
    % Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
    % Pped_{sc} = 4.53138 * (d + 0.0034) ^ 0.435509 * ptot ^ 0.121836 * iob ^ 1.51649
    iob = (2/5) .* geo.R .* zs.ip ./ 1e6 .* geo.K .^ 2 ./  geo.a .^ 2 ./ (1 + geo.K .^ 2);
    zs.pped   = max(0,4.53138e3 .* (min(0.5,abs(geo.d)) + 0.0034) .^ 0.435509 .* (p_scaling_pped ./ 1e6) .^ 0.121836 .* iob .^ 1.51649);
    
elseif option.usepped_scl == 4
    
    % use fit of pedestal database augmented with prediction from MHD stability studies for ITER and DEMO1
    % Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
    % Pped_{sc} = 4.53138 * (d + 0.0034) ^ 0.435509 * ptot ^ 0.121836 * iob ^ 1.51649
    iob = (2/5) .* geo.R .* zs.ip ./ 1e6 .* geo.K .^ 2 ./  geo.a .^ 2 ./ (1 + geo.K .^ 2);
    zs.pped   =  min(zs.pped,max(0,4.53138e3 .* (min(0.5,abs(geo.d)) + 0.0034) .^ 0.435509 .* (p_scaling_pped ./ 1e6) .^ 0.121836 .* iob .^ 1.51649));
    
elseif option.usepped_scl == 5
        % model reference :F.D. Halpern et al, PoP 15 (2008) p 062505
        [~,zs.pped] = z0pped_onjun(zs,profli,geo,cons,option);
end
if option.fpped > 0
	zs.pped  = option.fpped .* zs.pped;
end
switch option.hmore_pped
    case {1,4}
        zs.pped  = cons.hmore .* zs.pped;
    case {2,5}
        zs.pped  = min(1,cons.hmore) .* zs.pped;
end

% limitation de la pression en haut du piedestal (limite balooning + experimentale)
if option.fpped > 0
    switch option.hmore_pped
        case 1
            zs.pped    = min(zs.ppedmax .* abs(option.fpped) .* cons.hmore,zs.pped);
        case 2
            zs.pped    = min(zs.ppedmax .* abs(option.fpped) .* min(cons.hmore,1),zs.pped);
        case 3
            zs.pped    = min(zs.ppedmax .* abs(option.fpped) .* max(cons.hmore,1),zs.pped);
        case {4,5,6}
            zs.pped    = min(zs.ppedmax,zs.pped);
        otherwise
            zs.pped    = min(zs.ppedmax .* abs(option.fpped),zs.pped);
    end
end
% valeur minimale de pped
if mode ~= 0
    pped_min = 1.1 .* 1.602176462e-19 .* (zs.tebord .* zs.nebord + zs.tibord .* zs.nibord);
    zs.pped = max(pped_min,zs.pped);
end

% anti oscillation
if (option.ode_pped == 1) && (mode ~= 0)
  % estimation of tau_W_ped
  wped      = 3/2 .* zs.pped .* zs.vp;
  if isfield(profli,'qe') &&isfield(profli,'qi')
	pow_ped  = profli.qe(:,end-1) + profli.qi(:,end-1);
        pow_ped  = pow_ped .* (pow_ped > 0) + max(option.pth_min,min(zs.pin - zs.prad - zs.pbrem -zs.pcyclo - zs.pioniz,zs.pin - zs.dwdt)) .* (pow_ped <= 0);
  else
	pow_ped  = max(option.pth_min,min(zs.pin - zs.prad - zs.pbrem -zs.pcyclo - zs.pioniz,zs.pin - zs.dwdt));
  end
  if option.evolution == 1
      tau_W_ped = min(zs.taue,wped ./ pow_ped) .* (0.5 + 0.5 .* zs.modeh(end));
      pow_ped   = wped ./ tau_W_ped;
      [wped,dwdt]  = zdwdt0(wped(2:end),pow_ped(2:end) .* zs.modeh(end),tau_W_ped(2:end),cons.temps(2:end),zs.vp(2:end));
      %fprintf('Pped ode = %g  & Pped = %g\n',wped(3) ./ zs.vp(3) .* (2/3),zs.pped(3));
      zs.pped(3:end) = min(zs.pped(3:end),wped(3:end) ./ zs.vp(3:end) .* (2/3));
      %fprintf('mode_h ode = %g  & mode_h = %g\n',double(zs.pped(3) > pped_min(3)),zs.modeh(3));
      %fprintf('Tau_W_ped ode = %g  & Tau_E = %g\n',tau_W_ped(3),zs.taue(3));
      zs.modeh(3:end) = double(zs.pped(3:end) > pped_min(3:end));    
  else
      tau_W_ped = min(zs.taue,wped ./ pow_ped) .* (0.5 + 0.5 .* zs.modeh);
      pow_ped   = wped ./ tau_W_ped;
      [wped,dwdt]  = zdwdt0(wped,pow_ped .* zs.modeh,tau_W_ped,cons.temps,zs.vp);
      %figure(21);clf;subplot(3,1,1);plot(cons.temps,wped ./ zs.vp .* (2/3),'b',cons.temps,zs.pped,'r');
      zs.pped = min(zs.pped,wped ./ zs.vp .* (2/3));
      %subplot(3,1,2);plot(cons.temps,double(zs.pped > pped_min),'b.',cons.temps,zs.modeh,'r');
      %subplot(3,1,3);plot(cons.temps,tau_W_ped,'b',cons.temps,zs.taue,'r',cons.temps,zs.taup,'g');
      %set(gca,'ylim',[0,max(zs.taue)]);drawnow
      zs.modeh = double(zs.pped > pped_min);
  end
elseif option.evolution ~= 1
  try
    zs.pped  = sgolayfilt(zs.pped,1,3);
  catch   
    fprintf('T');
    zs.pped  = (cat(1,zs.pped(1),zs.pped(1:end-1)) + zs.pped + cat(1,zs.pped(2:end),zs.pped(end))) ./ 3;
  end
  % pas de pression de piedestal en mode L
  zs.pped  = zs.pped .* zs.modeh;
else
  % pas de pression de piedestal en mode L
  zs.pped  = zs.pped .* zs.modeh;
end
 
%figure(21);clf;plot(cons.temps,zs.pped,cons.temps,zs.ppedmax);drawnow

% insertion model imitant les elms
% la variation de densite est negligee (~qq % max, typique 3% de Nped)
dwelms = 0 .* zs.w; 
edge_density_elm_modulation = 0 .* zs.nebord;

if (option.dwow_elm ~= 0)  && (mode ~= 0)
    % evolution temporelle de la pression de reference 
    pped_mem  = zs.pped;
    if option.dwow_elm == 1
        pped_mhd_lim     = zs.ppedmax .* abs(option.fpped);
        pped_after_crash = zs.pped;
    elseif option.dwow_elm < 0
        switch option.gaz
            case 4
                zi = 2;
            otherwise
                % main gaz still have Z=1 in configuration gaz = 5 and 11
                zi = 1;
        end
        if isfield(profli,'zeff')
            zeff95 = profli.zeff(:,end-1);
        else
            zeff95 = zs.zeff;
        end
        dwow_elm         =  abs(option.dwow_elm) .* z0scalingelmsize(zs.teped,zs.tiped,zs.neped, ...
            zs.q95,geo.R,geo.a,zeff95,zi);
        pped_mhd_lim     = zs.pped ./ (1 - dwow_elm);
        pped_after_crash = zs.pped;
    else
        pped_mhd_lim     = zs.pped ./ (1 - option.dwow_elm);
        pped_after_crash = zs.pped;
    end
    pped_mhd_lim = max(pped_mhd_lim,pped_after_crash .* 1.01);
    
    [pped,void]  = zdwdt0((3/2) .* zs.pped .* zs.vp,zs.pth,option.tau_elm_factor .* zs.taue,cons.temps,zs.vp);
    pped = (2/3) .* pped ./ zs.vp;

     if (option.peeling == 2) && (isfield(profli,'jli'))
        jped_lim = zs.ip ./ zs.sp;
        jmax     = max(profli.jli(:,end-2:end),[],2);
        fpeeling = max(1,jmax ./ jped_lim);
	pped_after_crash = pped_after_crash ./ fpeeling;
	ind_elm  = find((jmax > jped_lim),1);
        
    elseif (option.peeling == 1) && (isfield(profli,'jli'))
        jped_lim = zs.ip ./ zs.sp;
        jmax     = max(profli.jli(:,end-2:end),[],2);
        fpeeling = max(1,jmax ./ jped_lim);
	pped_after_crash = pped_after_crash ./ fpeeling;
	ind_elm  = find((pped > pped_mhd_lim) | (jmax > jped_lim),1);
        
    else
	ind_elm   = find(pped > pped_mhd_lim,1);
    end
    nbmax_elm = length(cons.temps);
    while (~isempty(ind_elm)) && (nbmax_elm > 0)
        if ind_elm == length(cons.temps)
	    break;
        end
	nbmax_elm = nbmax_elm - 1;
        [pped_new,void]  = zdwdt0((3/2) .* pped_after_crash(ind_elm) .* zs.vp(ind_elm), ...
                                  zs.pth(ind_elm:end),option.tau_elm_factor .* zs.taue(ind_elm:end),cons.temps(ind_elm:end),zs.vp(ind_elm:end));
        pped(ind_elm:end) = (2/3) .* pped_new ./ zs.vp(ind_elm:end);
	ind_elm   = find(pped > pped_mhd_lim ,1);
    end

    % recopie du resultat
    dwelms  = pped - zs.pped;
    zs.pped = pped .* zs.modeh;
    edge_density_elm_modulation = min(0,max(-0.9,(pped_after_crash - zs.pped) ./ max(1,pped_mhd_lim))) .* zs.modeh;

%      figure(21);
%      clf;
%      subplot(2,1,1)
%      plot(cons.temps,pped,'b',cons.temps,pped_mhd_lim,'g',cons.temps,pped_mem,'r',cons.temps,pped_after_crash,'c');
%      legend('pped','crash limit','without elm','value after crash');
%      subplot(2,1,2) 
%      plot(cons.temps,dwelms,cons.temps,z0dxdt(dwelms,cons.temps))
%      drawnow
    
    % correction wth et w
    zs.wth = zs.wth + dwelms .* zs.vp .* zs.modeh;
    zs.w   = zs.w + dwelms .* zs.vp .* zs.modeh;

end
%d_dwelms_dt = z0dxdt(dwelms,cons.temps);

% limite MHD
zs.ppedmax   = z0maxpped(cons.temps,profli);
ppedla      = (2/3) .* zs.wth ./ zs.vp .* (zs.taue - taue_itb) ./ max(eps,zs.taue);
% suite publication de C.E .MAggi NF  47 (2007) p 535-551
zs.ppedmax  = min((2/3) .* ppedla, max(zs.ppedmax, ppedla ./ 4));

%  figure(136);clf
%  plot(cons.temps,zs.ppedmax,'ob',cons.temps,zs.pped,'r');
%  drawnow

wf = zs.w;
wf(~isfinite(wf))= 0;
wthf = zs.wth;
wthf(~isfinite(wthf))= 0;
% cas de donnees experimentales
if isappdata(0,'TE_EXP') & isappdata(0,'TI_EXP') & isfield(zs,'temps');
    if option.transitoire == 1
        pth = max(option.pth_min,zs.pth + zs.dwthdt);
    else
        pth = zs.pth;
    end
    if fwr == 1; disp('zdwdt0');end
    [void,dwthdt]  = zdwdt0(wthf,pth,zs.taue,cons.temps,zs.vp);
    dwthdt        = max(-plim,min(plim,dwthdt));
    if option.transitoire == 0
        zs.dwthdt(:)  = 0;
    elseif option.evolution == 1
        
        zs.wth(1:2) =  wthf(1:2);
        
        switch option.dwdt_method
            case {'implicit','mixed','freebie'}
                zs.dwthdt(end-1:end)     = dwthdt(end-1:end);
            case 'none'
                zs.dwthdt(:)  = 0;
            otherwise
                % rien
        end
    else
        switch option.dwdt_method
            case 'none'
                zs.dwthdt(:)  = 0;
            otherwise
                zs.dwthdt     = dwthdt;
        end
    end
    % fast electron due to runaway electrons addedt to LH terms
    if option.runaway ~= 0
        % assuming v// ~ c
        esup_run = zs.irun .* phys.me .* phys.c ./ phys.e;
    else
        esup_run = zeros(size(cons.temps));
    end
    zs.w          = zs.wth + (zs.esup_fus + zs.esup_icrh + real(zs.esup_nbi) + imag(zs.esup_nbi) + zs.esup_lh + esup_run);
    zs.dwdt       = zs.dwthdt + z0dxdt(zs.esup_fus + zs.esup_icrh + real(zs.esup_nbi) + imag(zs.esup_nbi) + zs.esup_lh + esup_run,cons.temps);
    if option.evolution == 1
        switch option.dwdt_method
            case 'freebie'
                zs.dwdt       = zs.dwthdt + z0dxdt_freebie(zs.esup_fus + zs.esup_icrh + real(zs.esup_nbi) + imag(zs.esup_nbi) + zs.esup_lh,cons.temps);
        end
    end
    
elseif mode ~= 0
    if option.transitoire == 1
        pth = max(option.pth_min,zs.pth + zs.dwthdt);
    else
        pth = zs.pth;
    end
    if fwr == 1; disp('zdwdt0');end
    [zs.wth,dwthdt]  = zdwdt0(wthf,pth,zs.taue,cons.temps,zs.vp);
    
    %  	figure(21);
    %  	clf
    %  	subplot(3,1,1);
    %  	plot(cons.temps,zs.wth,'-+',cons.temps,wthf,'-o');
    %  	ylabel('wth')
    %  	subplot(3,1,2)
    %  	plot(cons.temps,dwthdt,'-+',cons.temps,zs.dwthdt,'-o');
    %  	ylabel('dwthdt')
    %  	subplot(3,1,3)
    %  	plot(cons.temps,zs.taue,'-o');
    %  	ylabel('taue')
    %  	drawnow
    
    dwthdt        = max(-plim,min(plim,dwthdt));
    if option.transitoire == 0
        zs.dwthdt(:)  = 0;
    elseif option.evolution == 1
        
        zs.wth(1:2) =  wthf(1:2);
        
        switch option.dwdt_method
            case {'implicit','mixed'}
                zs.dwthdt(end-1:end)     = dwthdt(end-1:end);
            case 'none'
                zs.dwthdt(:)  = 0;
            otherwise
                % rien
        end
    else
        switch option.dwdt_method
            case 'none'
                zs.dwthdt(:)  = 0;
            otherwise
                zs.dwthdt     = dwthdt;
        end
    end
    % fast electron due to runaway electrons addedt to LH terms
    if option.runaway ~= 0
        % assuming v// ~ c
        esup_run = zs.irun .* phys.me .* phys.c ./ phys.e;
    else
        esup_run = zeros(size(cons.temps));
    end
    zs.w          = zs.wth + (zs.esup_fus + zs.esup_icrh + real(zs.esup_nbi) + imag(zs.esup_nbi)  + zs.esup_lh + esup_run);
    zs.dwdt       = zs.dwthdt + z0dxdt(zs.esup_fus + zs.esup_icrh + real(zs.esup_nbi) + imag(zs.esup_nbi) + zs.esup_lh + esup_run,cons.temps);
    if option.evolution == 1
        switch option.dwdt_method
            case 'freebie'
                %zs.dwdt = zs.dwthdt + polyder(polyfit(cons.temps,zs.esup_fus + zs.esup_icrh + real(zs.esup_nbi) + imag(zs.esup_nbi) + zs.esup_lh,1));
                zs.dwdt       = zs.dwthdt + z0dxdt_freebie(zs.esup_fus + zs.esup_icrh + real(zs.esup_nbi) + imag(zs.esup_nbi) + zs.esup_lh,cons.temps);
        end
    end
    
else
    zs.dwdt    = zeros(size(cons.temps));
    zs.dwthdt  = zeros(size(cons.temps));
end

% puissance equivalente
zs.pw     = zs.w ./ zs.taue;

% piguages% palsma de fond
switch option.gaz
    case 1
        A = ones(size(cons.temps));
    case 2
        A = 2 .* ones(size(cons.temps));
    case 11
        A = (1  + 11 .* cons.iso) ./ (1 + cons.iso);       
    case {3,5}
        A = (2  + 3 .* cons.iso) ./ (1 + cons.iso);
    case 4
        A = 4 .* ones(size(cons.temps));
end
% A       = 1 + (option.gaz == 4) ;
ngr     = 1e20 .* (zs.ip /1e6) ./ (pi.* geo.a .^ 2);
zs.nsat = min(ngr,0.06e20 .* (zs.ip ./ 1e6) .* geo.R .* sqrt(A) ./ geo.K ./ geo.a .^ (5/2));
switch option.ane
    case 1
        % cas plat
        zs.ane = 0.01 .* ones(size(zs.nsat));
    case 2
        % cas n  = n0 ./sqrt(q)
        zs.ane =  max(0.1,min(10,((4/3) .* zs.li + 0.25)  - 1));
    case 3
        % scaling nu_eff (H. Wiesen)
        nu_eff = 1e-14 .* geo.R .* zs.zeff .* zs.nem ./ max(30,zs.tem) .^ 2;
        ane_h = max(1,min(2,(1.28 - 0.17 .* log(nu_eff)))) - 1;
        ane_l = max(0.1,min(10,((4/3) .* zs.li + 0.25)  - 1));
        switch option.modeh
            case 0
                zs.ane = ane_h;
            otherwise
                zs.ane = zs.modeh .* ane_h + (~zs.modeh) .* ane_l;
        end
    case 4
        % cas donne en entree
        zs.ane = option.vane .* ones(size( zs.ane)) - 1;
        if any(zs.ane < 0)
            ind_neg = find(zs.ane < 0);
            if mode ~= 0
                ane_min = (zs.nebord + 1e13) ./ zs.nem - 1;
                zs.ane(ind_neg) = max(ane_min(ind_neg),zs.ane(ind_neg));
            else
                zs.ane(ind_neg) = 0.1;
            end
        end
    case {5,6,7}
        % calcule a partir du profil
        if isfield(profli,'nep')
            zs.ane = profli.nep(:,1) ./ trapz(profli.xli,profli.nep .* profli.vpr,2) .* ...
                trapz(profli.xli,profli.vpr,2) - 1;
        else
            zs.ane = ones(size( zs.ane));
        end
    case 10
        
        % en mode H scaling formule 5 de C. Angioni ref : NF 47 (2007) p 1326-1335
        % sinon en mode L comme le default
        % le piquage en mode L
        ane_l  = max(0.1,min(10, 0.3536 .* zs.nsat ./ zs.nbar));
        fnbis  = real(zs.pnbi) ./ (1.602176462e-19 .* option.einj) + imag(zs.pnbi) ./ (1.602176462e-19 .* option.einj2);
        % le facteur 2 est remplace par une meilleur approximation car disponible : (1 + zs.tite .* zs.nim ./ zs.nem)
        % utilisation de Te0 au lieu de Te2
        fnbis  = (1 + zs.tite .* zs.nim ./ zs.nem) .* fnbis .* (1.602176462e-19 .* zs.te0) ./ max(option.pth_min,zs.pin) .* (zs.ate + 1 - 0.37);
        ngr     = 1e20 .* (zs.ip /1e6) ./ (pi.* geo.a .^ 2);
        ane_h  = 0.253 - 0.499 .* (zs.nbar ./ ngr)  + 2.094 .* fnbis + 0.117 .* geo.R;
        ane_h   = max(0.01,min(10,ane_h));
        zs.ane = option.ane_factor .* (zs.modeh .* ane_h + (~zs.modeh) .* ane_l);
        
    case 11
        
        % changement de confinement Tokamak , Wesson p 176
        % le piquage en mode L :
        % observation : Ohmic energy confinement saturation and core toroidal rotation reversal in Alcator C- Mod plasmas
        % J. E. Rice et al, Physics of Plasmas (1994-present) 19, 056106 (2012); doi: 10.1063/1.3695213
        % Th?se Maxim Irishkin
        % finamelament fit de la base L-mode (ne0 <  10^21)
        fnbis  = real(zs.pnbi) ./ (1.602176462e-19 .* option.einj) + imag(zs.pnbi) ./ (1.602176462e-19 .* option.einj2);
        snbis  = 1 + fnbis ./ (zs.nem .* zs.vp);
        ane_l   = max(0.1,min(10,1.34 .* (zs.nem ./ 1e20) .^ (-0.1) .* zs.li .^ 0.2  .* snbis .^ 0.2 - 1));
        
        
        % en mode H scaling formule 5 de C. Angioni ref : NF 47 (2007) p 1326-1335
        % le facteur 2 est remplace par une meilleur approximation car disponible : (1 + zs.tite .* zs.nim ./ zs.nem)
        % utilisation de Te0 au lieu de Te2
        fnbis  = (1 + zs.tite .* zs.nim ./ zs.nem) .* fnbis .* (1.602176462e-19 .* zs.te0) ./ max(option.pth_min,zs.pin) .* (zs.ate + 1 - 0.37);
        ngr     = 1e20 .* (zs.ip /1e6) ./ (pi.* geo.a .^ 2);
        ane_h  = 0.253 - 0.499 .* (zs.nbar ./ ngr)  + 2.094 .* fnbis + 0.117 .* geo.R;
        ane_h   = max(0.01,min(10,ane_h));                 
        zs.ane = option.ane_factor .* (zs.modeh .* ane_h + (~zs.modeh) .* ane_l);
        
    case 12
        % changement de confinement Tokamak , Wesson p 176
        % le piquage en mode L :
        % observation : Ohmic energy confinement saturation and core toroidal rotation reversal in Alcator C- Mod plasmas
        % J. E. Rice et al, Physics of Plasmas (1994-present) 19, 056106 (2012); doi: 10.1063/1.3695213
        % Th?se Maxim Irishkin
        % finamelament fit de la base L-mode (ne0 <  10^21)
        fnbis  = real(zs.pnbi) ./ (1.602176462e-19 .* option.einj) + imag(zs.pnbi) ./ (1.602176462e-19 .* option.einj2);
        snbis  = 1 + fnbis ./ (zs.nem .* zs.vp);
        ane_l   = max(0.1,min(10,1.34 .* (zs.nem ./ 1e20) .^ (-0.1) .* zs.li .^ 0.2  .* snbis .^ 0.2 - 1));
        
        
        % SPARC scaling J. Plasma Phys. (2020), vol. 86, 865860502 
        % similar toC. Angioni ref : NF 47 (2007) p 1326-1335
        % le facteur 2 est remplace par une meilleur approximation car disponible : (1 + zs.tite .* zs.nim ./ zs.nem)
        % utilisation de Te0 au lieu de Te2
        fnbis  = (1 + zs.tite .* zs.nim ./ zs.nem) .* fnbis .* (1.602176462e-19 .* zs.te0) ./ max(option.pth_min,zs.pin) .* (zs.ate + 1 - 0.37);
        ngr     = 1e20 .* (zs.ip /1e6) ./ (pi.* geo.a .^ 2);
        %    ane_h  = 0.253 - 0.499 .* (zs.nbar ./ ngr)  + 2.094 .* fnbis + 0.117 .* geo.R;
        
        nu_eff  = 0.2.*zs.zeff.*(zs.nem*1e-19).*geo.R./(zs.tem*1e-3).^2;
        beta_angioni=4.02e-3.*((zs.tem.*(1 + zs.tite ))*1e-3).*(zs.nem*1e-19)./geo.b0.^2;
        ane_h   = 0.347 - 0.117.*log(nu_eff)-4.03.*beta_angioni + 0.5*1.331.*fnbis;
        ane_h   = ane_h - 0.15;
        
        ane_h   = max(0.01,min(10,ane_h));
        zs.ane = option.ane_factor .* (zs.modeh .* ane_h + (~zs.modeh) .* ane_l);
        
    otherwise
        % changement de confinement Tokamak , Wesson p 176
        % le piquage en mode L
        ane_l   = max(0.1,min(10, 0.3536 .* zs.nsat ./ zs.nbar));
        ngr     = 1e20 .* (zs.ip /1e6) ./ (pi.* geo.a .^ 2);
        ane_h   = max(0.1,min(10, 0.23 .* ngr ./ zs.nbar));
        zs.ane = zs.modeh .* ane_h + (~zs.modeh) .* ane_l;
end

% cas des dent de scie resolu dans le temps
indice_inv_min = 1.5;
if (option.qdds < 0) && isfield(zs,'indice_inv') && any(zs.indice_inv > indice_inv_min) && isfield(profli,'nep')
    
    ane_crash = profli.nep(:,1) ./ trapz(profli.xli,profli.nep .* profli.vpr,2) .* ...
        trapz(profli.xli,profli.vpr,2) - 1;
    % temps de confimenent global de la matiere = 3 * taue
    taup_loc = zs.taue ./ 7 .* min(1,zs.indice_inv);
    % seul l'evolution relative est utilisees
    if option.evolution == 1
        sane = cat(1,0,diff(ane_crash) - diff(zs.ane)) ./ zs.taue;
        [tntot,ane_evol] = z0ode(cons.temps,sane,max(taup_loc)*ones(size(cons.temps)),0);
    else
        
        ind_comp = find(zs.indice_inv > indice_inv_min);
        if max(ind_comp) == length(zs.indice_inv)
            ind_comp = ind_comp(1:end-1);
        end
        if min(ind_comp) == 1
            ind_comp = ind_comp(2:end);
        end
        dane = ane_crash(ind_comp - 1) - ane_crash(ind_comp) + zs.ane(ind_comp - 1) - zs.ane(ind_comp);
        %dane_c = ane_crash(ind_comp) - ane_crash(ind_comp + 1) + zs.ane(ind_comp) - zs.ane(ind_comp + 1);
        %ind_wrong = find((dane >=0) | (dane_c < 0) | (dane >= dane_c));
        ind_wrong = find(dane >=0);
        dane(ind_wrong) =[];
        ind_comp(ind_wrong) =[];
        ane_evol = zeros(size(ane_crash));
        sane     = zeros(size(ane_crash));
        %sane(ind_comp) = dane ./ zs.taue(ind_comp);
        for lm = 1:length(ind_comp)
            %  	        if lm == length(ind_comp)
            %  		  [tntot,ane_one_dds] = z0ode(cons.temps(ind_comp(lm):end),sane(ind_comp(lm):end),zs.taue(ind_comp(lm):end),dane(lm));
            %  		  ane_evol(ind_comp(lm):end) = ane_one_dds;
            %  		else
            %  		  ane_evol(ind_comp(lm):(ind_comp(lm + 1) - 1)) = linspace(dane(lm),0,length(ind_comp(lm):(ind_comp(lm + 1) - 1)));
            %  		end
            taup_crash = taup_loc(ind_comp(lm)) * ones(size(cons.temps(ind_comp(lm):end)));
            [tntot,ane_one_dds] = z0ode(cons.temps(ind_comp(lm):end),sane(ind_comp(lm):end),taup_crash,dane(lm) + ane_evol(ind_comp(lm)));
            ane_evol(ind_comp(lm):end) = ane_one_dds;
        end
    end
    ane_mem_ =zs.ane;
    zs.ane = max(0.5,zs.ane + ane_evol);
    %figure(21);plot(cons.temps,zs.ane,'r',cons.temps,ane_mem_,'b',cons.temps,ane_crash,'g');drawnow
    %figure(21);plot(cons.temps,ane_evol);drawnow
end

% loi empirique par rapport a une analyse analytique + effet mode h + effet Te > Ti
if mode == 0
   zs.ape  = 1./3.*(4.*zs.qa.^3+(-1-8.*log(zs.qa)).*zs.qa.^2+(-4+2.*log(zs.qa)).*zs.qa+1)./ ...
               ((2.*log(zs.qa)-3).*zs.qa.^2+4.*zs.qa-1) + tanh(1./ real(zs.tite));
   zs.ape  = zs.ape .* (~zs.modeh) +  zs.ape .* zs.modeh .* max(0.1,min(1,zs.tauthl./max(zs.tauh,eps)));
   zs.ate  = max(0.1,min(10,(zs.ape - zs.ane)));
end

% nan et imag
z0dnanimag
%save('sat4','option','cons','geo','zs','amorti','profli');

% cas de l'asservissement sur nbar
if (option.neasser >= 1) & (mode ~= 0) & isfield(profli,'nep')
   % donnees derivees des consignes
   nemi    = trapz(profli.xli,profli.nep,2);
   nem     = trapz(profli.xli,profli.nep.* profli.vpr,2) ./ trapz(profli.xli,profli.vpr,2);
   rapbar  = nem ./ nemi;
   nem_1   = cons.nbar(1) .* rapbar(1);
   nemref  = cons.nbar .* rapbar;
   % prise en compte du recyclage au niveau de la LCFS
   switch option.configuration
      case {0,1}
	fn0a    = option.fn0a ;
      case {2,3}
	fn0a    = (zs.xpoint .* option.fn0a_div + (~zs.xpoint) .* option.fn0a);
      otherwise
	fn0a    = option.fn0a_div;
   end
   tauref = min(1e3,max(min(zs.tauhe,zs.taup) ./ max(eps,1 - fact_confinement .* option.Recycling .* fn0a),1e-6));
   snem   = nemref ./ tauref .* zs.vp + (real(zs.pnbi_th) ./ option.einj + imag(zs.pnbi_th) ./ option.einj2) + option.eta_gas_puff .* gas_puff_mem;
   if option.evolution == 1
        ntot         = nem .* zs.vp;
	[tntot,ntot(2:end)] = z0ode(cons.temps(2:end),snem(2:end),tauref(2:end),ntot(2));
   else
	[tntot,ntot] = z0ode(cons.temps,snem,tauref,nem_1 .* zs.vp(1));
   end
   zs.nem  = ntot ./ zs.vp;
   zs.nbar = zs.nem ./ rapbar;
   zs.ne0  = profli.nep(:,1);

%     figure(22);
%     clf;
%     subplot(2,1,1);
%     plot(cons.temps,zs.nem,'r',cons.temps,nemref,'b');
%     set(gca,'xlim',[0 3]);
%     subplot(2,1,2);
%     semilogy(cons.temps,tauref,'r',cons.temps,zs.taup,'b',cons.temps,zs.tauhe,'g');
%     %plot(cons.temps,factor_fulling)
%     set(gca,'xlim',[0 3]);
%     pause(0.1);



elseif isfield(profli,'nep')
   zs.nbar = cons.nbar;
   nemi    = trapz(profli.xli,profli.nep,2);
   nem     = trapz(profli.xli,profli.nep.* profli.vpr,2) ./ trapz(profli.xli,profli.vpr,2);
   rapbar  = nem ./ nemi;
   zs.nem  = cons.nbar .*  rapbar;
   zs.ne0  = profli.nep(:,1);

elseif (option.neasser >= 1) & (mode ~= 0)
   % donnees derivees des consignes
   ne0ref = cons.nbar .* 2 .* gamma(zs.ane + 1.5) ./ gamma(zs.ane+1) ./ sqrt(pi); 
   nemref = ne0ref ./ (1 + zs.ane);
   tauref = min(1e3,max(min(zs.tauhe,2 .* zs.taue),min(gradient(zs.temps))));
   snem   = nemref ./ tauref .* zs.vp + (real(zs.pnbi_th) ./ option.einj + imag(zs.pnbi_th) ./ option.einj2);
   if option.evolution == 1
      ntot         = nemref .* zs.vp;
      [tntot,ntot(2:end)] = z0ode(cons.temps(2:end),snem(2:end),tauref(2:end),ntot(2));
   else
      [tntot,ntot] = z0ode(cons.temps,snem,tauref,nemref(1) .* zs.vp(1));
   end
   zs.nem = ntot ./ zs.vp;
   zs.ne0 = (1 + zs.ane) .* zs.nem;
   zs.nbar = zs.ne0 ./ (2 .* gamma(zs.ane + 1.5) ./ gamma(zs.ane+1) ./ sqrt(pi));
else
   % donnees derivees des consignes
   zs.nbar = cons.nbar;
   zs.ne0 = cons.nbar .* 2 .* gamma(zs.ane + 1.5) ./ gamma(zs.ane+1) ./ sqrt(pi);
   zs.nem = zs.ne0 ./ (1 + zs.ane);
end

% densite de greenwald
zs.negr   = 1e20 .* (zs.ip /1e6) ./ (pi.* geo.a .^ 2);

% selon la configuration
% model de nebord de resulat pour iter de Mahdavi et al, Physics of plasmas 10 (2003) p 3984-... en mode H
% et G. Porter er al J. of Nucl. Mater (1999) vol 266-269 p 917- pour le mode L 
% il y a un probleme avec le fit de Porter : ce comportement ne rend pas compte de JET et ASDEX
% on reprend le fit ITER CDA (meme publication)
% en point x mode H et mode L
% attention en mode H nebord depend aussi de la frequence des elms
%  nebord_x  = (zs.modeh) .*  cnb .* min(zs.negr,zs.nbar) .^ 2 + ...
%  	    (~zs.modeh) .* 0.00236 .* min(zs.negr,zs.nbar) .^ 1.08  .* geo.K .^ 1.11 .* geo.b0 .^ 0.78;
switch option.nea_model
    case 'Eich'
        % ref : T. Eich et al, NF 58 (2018) 034001
        nebord_xh = zs.negr .* min(0.95,5.9 .* 2 .* (geo.R ./ geo.a) .^ (-2/7) .* ((1 + geo.K .^ 2) ./ 2) .^ (-6/7) .*  ...
            zs.plhthr .^ (-11/70));
        %figure(21);plot(cons.temps,nebord_xh,'b',cons.temps,zs.negr,'r');drawnow
        nebord_x  = (zs.modeh) .*  nebord_xh .* (1 + edge_density_elm_modulation) + ...
            (~zs.modeh) .* 3.4485e-13.* min(zs.negr,zs.nbar) .^ 1.6;
    case 'fixed ratio'
        nebord_x  = min(zs.negr,zs.nbar) / 3;
        
    otherwise
        cnb       = max(1e-21,5e-21 - 6.7e-24 .* zs.tebord);
        nebord_x  = (zs.modeh) .*  cnb .* min(zs.negr,zs.nbar) .^ 2 .* (1 + edge_density_elm_modulation) + ...
            (~zs.modeh) .* 3.4485e-13.* min(zs.negr,zs.nbar) .^ 1.6;
end
% en limiteur toroidal (F. Clairet)
nebord_lim = 1e-21 .* min(zs.negr,zs.nbar) .^ 2 .* zs.qa .* geo.R;
% en limiteur poloidal (Wesson Tokamak, loi multi machine)
nebord_pol = 5e-21 .* min(zs.negr,zs.nbar) .^ 2 ; 

%  figure(151);clf
%  semilogy(cons.temps,zs.nbar,cons.temps,zs.negr,cons.temps,nebord_x,cons.temps,nebord_lim,cons.temps,nebord_pol)
%  legend('nbar','negr','x point','toroidal','poloidal');
%  drawnow
switch option.configuration
case 0
	zs.nebord = nebord_pol;
case 1
	zs.nebord = nebord_lim;
case 2
	zs.nebord = zs.xpoint .* nebord_x + (~zs.xpoint) .* nebord_pol;
case 3
	zs.nebord = zs.xpoint .* nebord_x + (~zs.xpoint) .* nebord_lim;
otherwise
	zs.nebord = nebord_x;
end

%zs.nebord = min((option.modeh) .* 5e-21 .* cons.nbar .^ 2 + ...
%            (~option.modeh) .* 1e-21 .* cons.nbar .^ 2 .* zs.qa .* geo.R, ...
%	    0.99 .* cons.nbar);
% model de nebord de resulat pour iter de Mahdavi et al, Physics of plasmas 10 (2003) p 3984-...
% et du sacling en mode limiteur toroidal pour TS (longueur de connexion)
% sinon il faut utiliser 5e-21 .* cons.nbar .^2
% normalement depend de ti_bord (mais pas d'information sur tibord)
%  cnb       = max(1e-21,5e-21 - 6.7e-24 .* zs.tebord);
%  zs.nebord = min((option.modeh) .* ( zs.modeh .* cnb  + (~zs.modeh) .* 5e-21) .* ...
%              min(zs.negr,zs.nbar) .^ 2 + ...
%              (~option.modeh) .* 1e-21 .* min(zs.negr,zs.nbar) .^ 2 .* zs.qa .* geo.R, ...
%  	    0.99 .* zs.nbar);

% fraction de fuelling de coeur
edge_gnbi = real(zs.pnbi) ./ (option.einj .* 1.602176462e-19) ./ (4.41e-4 .* 6.02214199e23) + ...
            imag(zs.pnbi) ./ (option.einj2 .* 1.602176462e-19) ./ (4.41e-4 .* 6.02214199e23);
edge_gpellet  = zs.n0a./ max(eps,1 - zs.frac_pellet) .* zs.frac_pellet ./ (4.41e-4 .* 6.02214199e23);
edge_grecycle = zs.n0a./ max(eps,1 - zs.frac_pellet) ./ (4.41e-4 .* 6.02214199e23);
edge_gtot     = max(1,edge_grecycle + edge_gpellet + edge_gnbi);
edge_gcore    = edge_gpellet + edge_gnbi;
edge_eta_c    = edge_gcore ./ edge_gtot; 
edge_eta_c(~isfinite(edge_eta_c)) = 0;
edge_eta_c    = max(0,min(1,edge_eta_c));


% cas glacon (observation TS : la densite de bord diminue avec le fuelling de coeur)
%zs.nebord = zs.nebord .* max(0.1,(1 - zs.frac_pellet)); 
if option.nea_factor < 0
    zs.nebord = abs(option.nea_factor) .* zs.nbar; 
else
    zs.nebord = option.nea_factor .* zs.nebord ./ (1 + 0.18 .* edge_eta_c); 
end

% si profils lu en entree
if isappdata(0,'NE_EXP') & isfield(profli,'nep');	 
	zs.nebord = profli.nep(:,end);
end   


% sercurite
zs.nebord = min(zs.nebord,0.99 .* min(zs.nbar,cons.nbar));
zs.nebord = max(1e13,real(zs.nebord));


% composition du plasma et  calcul du zeff interne et meff
%  pour le calcul de la densite de H (minoritaire) dans un plasam de He
nhem_mem = zs.nhem;
% accumulation des cendres de He
if option.gaz == 5
    % D-He3 reaction case
    % nhem encode here for He3 density, He4 density is encode by
%     % option.frhe0 .* zs.nem and it is not added to zs.nhem.
    % at this stage nimp is assume known as nT
    if option.transitoire == 1
        nhe4m0 = option.frhe0 .* zs.nem + ...
                  max(0,zs.esup_fus_he4_DHe3) ./ (fus.dhe3_he4 .* 1.6022e-19) .* 2 ./ zs.vp + ...
                  max(0,zs.esup_fus_he4_DT)   ./ (fus.dt_he4   .* 1.6022e-19) .* 2 ./ zs.vp;
        nHm0   = 0 .* zs.nem + ...
                max(0,zs.esup_fus_p_DHe3) ./ (fus.dhe3_p .* 1.6022e-19) .* 2 ./ zs.vp + ...
                max(0,zs.esup_fus_p_DDp)  ./ (fus.ddp_p  .* 1.6022e-19) .* 2 ./ zs.vp;
        nTm0   = 0 .* zs.nem + max(0,zs.esup_fus_t_DDp) ./ (fus.ddp_t .* 1.6022e-19) .* 2 ./ zs.vp;
        if option.evolution == 1   
            % assuming fusion produced fast ions tritium and proton have same particle confinement time (zs.tauhe) of that of alpha particle 
            zs.nhe4m(2:end) = znalpha0(cons.temps(2:end),zs.salpha_he4(2:end),zs.tauhe(2:end), ...
                                       max(1,zs.nhe4m(2:end) - nhe4m0(2:end) ),zs.vp(2:end)) ./ zs.vp(2:end) + ...
                              option.frhe0 .* zs.nem(2:end) + ...
                              max(0,zs.esup_fus_he4_DHe3(2:end)) ./ (fus.dhe3_he4 .* 1.6022e-19) .* 2 ./ zs.vp(2:end) + ...
                              max(0,zs.esup_fus_he4_DT(2:end))   ./ (fus.dt_he4   .* 1.6022e-19) .* 2 ./ zs.vp(2:end);
            zs.nHm(2:end)   = znalpha0(cons.temps(2:end),zs.salpha_p(2:end),zs.tauhe(2:end), ...
                                       max(1,zs.nHm(2:end) - nHm0(2:end) ),zs.vp(2:end)) ./ zs.vp(2:end) +...
                              max(0,zs.esup_fus_p_DHe3(2:end))   ./ (fus.dhe3_p .* 1.6022e-19) .* 2 ./ zs.vp(2:end) + ...
                              max(0,zs.esup_fus_p_DDp(2:end))    ./ (fus.ddp_p  .* 1.6022e-19) .* 2 ./ zs.vp(2:end);            
            zs.nTm(2:end)   = znalpha0(cons.temps(2:end),zs.salpha_t(2:end),zs.tauhe(2:end), ...
                                       max(1,zs.nTm(2:end) - nTm0(2:end) ),zs.vp(2:end)) ./ zs.vp(2:end) + ...
                              max(0,zs.esup_fus_t_DDp(2:end))    ./ (fus.ddp_t  .* 1.6022e-19) .* 2 ./ zs.vp(2:end);
            %          nhem1  = zs.salpha .* zs.tauhe ./ (zs.vp+eps) + option.frhe0 .* zs.nem + ...
            %                     max(0,zs.esup_fus) ./ (3.56e6 .* 1.6022e-19).* 2 ./ zs.vp;
            
            %        disp([zs.nhem(end-1:end),mem_nhem(end-1:end),nhem1(end-1:end)]./1e19)
            
            %          figure(21);clf
            %          subplot(4,1,1)
            %          plot(cons.temps,zs.nhem,'r+-',cons.temps,mem_nhem,'b-o',cons.temps,nhem1,'*-g');
            %          subplot(4,1,2)
            %          plot(cons.temps,zs.salpha,'-o');
            %          subplot(4,1,3)
            %          plot(cons.temps,zs.tauhe,'-o');
            
        else
            zs.nhe4m= znalpha0(cons.temps,zs.salpha_he4,zs.tauhe, ...
                               max(1,zs.nhe4m - nhe4m0),zs.vp) ./ zs.vp + ...
                      option.frhe0 .* zs.nem + ...
                      max(0,zs.esup_fus_he4_DHe3) ./ (fus.dhe3_he4 .* 1.6022e-19) .* 2 ./ zs.vp + ...
                      max(0,zs.esup_fus_he4_DT)   ./ (fus.dt_he4   .* 1.6022e-19) .* 2 ./ zs.vp;
            zs.nHm  = znalpha0(cons.temps,zs.salpha_p,zs.tauhe, ...
                               max(1,zs.nHm - nHm0 ),zs.vp) ./ zs.vp +...
                      max(0,zs.esup_fus_p_DHe3)   ./ (fus.dhe3_p .* 1.6022e-19) .* 2 ./ zs.vp + ...
                      max(0,zs.esup_fus_p_DDp)    ./ (fus.ddp_p  .* 1.6022e-19) .* 2 ./ zs.vp;            
            zs.nTm  = znalpha0(cons.temps,zs.salpha_t,zs.tauhe, ...
                               max(1,zs.nTm - nTm0 ),zs.vp) ./ zs.vp + ...
                      max(0,zs.esup_fus_t_DDp)    ./ (fus.ddp_t  .* 1.6022e-19) .* 2 ./ zs.vp;
        end
    else
        zs.zeff  = max(zs.zeff,1.1);
        switch option.mino
            case 'H'
                zs.nDm   = ((zs.nem .* (1 - 2 .* option.frhe0) - zs.nTm) .* zu2 - ...
                    (zs.nem .* zs.zeff .* (1 - 4 .* option.frhe0) - zs.nTm) .* zu1) ./ ...
                    ((1+ option.cmin + 2 .* cons.iso) .* zu2 - (1+ option.cmin + 4 .* cons.iso) .* zu1);
                nHm   = option.cmin .* zs.nDm; % not used in METIS
            otherwise
                if ~isfield(zs,'nTm')
                    zs.nTm  = zeros(size(zs.nem));
                    %disp('why 2')
                end
                zs.nDm   = ((zs.nem .* (1 - 2 .* option.frhe0) - zs.nTm) .* zu2 - ...
                    (zs.nem .* zs.zeff .* (1 - 4 .* option.frhe0) - zs.nTm) .* zu1) ./ ...
                    ((1 + 2 .* cons.iso) .* zu2 - (1 + 4 .* cons.iso) .* zu1);
                nHm   = zeros(size(zs.nDm));
        end
        zs.nhem =   cons.iso .* zs.nDm;
        %         zs.n1m   = nHm + zs.nDm + zs.nTm;
        %         zs.nimpm = (zs.nem .* (1 - 2 .* option.frhe0) - zs.nDm - nHm - nTm - 2 .* nhem) ./ zu1; % W is encode in zu1
        %         zs.nim   = zs.n1m + zs.nem .* option.frhe0 + zs.nhem + zs.nimpm .* (1 + option.rimp) + zs.nwm;
        %         zs.ni0   = zs.ne0 .* zs.nim ./ zs.nem;
        %         ae       = zs.nim ./ zs.nem;
        
        zs.nhe4m  = zs.salpha_he4 .* zs.tauhe  ./ (zs.vp+eps) + option.frhe0 .* zs.nem + ...
                    max(0,zs.esup_fus_he4_DHe3) ./ (fus.dhe3_he4.* 1.6022e-19).* 2 ./ zs.vp + ...
                    max(0,zs.esup_fus_he4_DT)   ./ (fus.dt_he4.* 1.6022e-19).* 2 ./ zs.vp;
        zs.nHm    = zs.salpha_p .* zs.tauhe  ./ (zs.vp+eps) + ...
                    max(0,zs.esup_fus_p_DHe3) ./ (fus.dhe3_p.* 1.6022e-19).* 2 ./ zs.vp + ...
                    max(0,zs.esup_fus_p_DDp)  ./ (fus.ddp_p.* 1.6022e-19).* 2 ./ zs.vp;
        zs.nTm    = zs.salpha_t .* zs.tauhe  ./ (zs.vp+eps) + ...
                    max(0,zs.esup_fus_t_DDp)  ./ (fus.ddp_t.* 1.6022e-19).* 2 ./ zs.vp;
    end
    

elseif option.gaz == 11
    % helium accumulation should be added has for DT (*3 taken into account in fusion module) ?
    emean_alpha = (4e6 + 2 .* 2.3e6) ./ 3;
    if option.transitoire == 1
        if option.evolution == 1
            mem_nhem = zs.nhem;
            nhem0  = option.frhe0 .* zs.nem(2:end) + max(0,zs.esup_fus(2:end)) ./ (emean_alpha .* 1.6022e-19) .* 2 ./ zs.vp(2:end);
            zs.nhem(2:end) =  znalpha0(cons.temps(2:end),zs.salpha(2:end),zs.tauhe(2:end), ...
                max(1,zs.nhem(2:end) - nhem0),zs.vp(2:end)) ./ zs.vp(2:end) + ...
                option.frhe0 .* zs.nem(2:end) +  ...
                max(0,zs.esup_fus(2:end)) ./ (emean_alpha.* 1.6022e-19).* 2 ./ zs.vp(2:end);
            
            %          nhem1  = zs.salpha .* zs.tauhe ./ (zs.vp+eps) + option.frhe0 .* zs.nem + ...
            %                     max(0,zs.esup_fus) ./ (3.56e6 .* 1.6022e-19).* 2 ./ zs.vp;
            
            %        disp([zs.nhem(end-1:end),mem_nhem(end-1:end),nhem1(end-1:end)]./1e19)
            
            %          figure(21);clf
            %          subplot(4,1,1)
            %          plot(cons.temps,zs.nhem,'r+-',cons.temps,mem_nhem,'b-o',cons.temps,nhem1,'*-g');
            %          subplot(4,1,2)
            %          plot(cons.temps,zs.salpha,'-o');
            %          subplot(4,1,3)
            %          plot(cons.temps,zs.tauhe,'-o');
            
        else
            nhem0  = option.frhe0 .* zs.nem + max(0,zs.esup_fus) ./ (3.56e6 .* 1.6022e-19) .* 2 ./ zs.vp;
            zs.nhem  = znalpha0(cons.temps,zs.salpha,zs.tauhe,max(1,zs.nhem - nhem0),zs.vp) ./ zs.vp + option.frhe0 .* zs.nem +  ...
                max(0,zs.esup_fus) ./ (emean_alpha .* 1.6022e-19).* 2 ./ zs.vp;
        end
    else

        zs.nhem  = zs.salpha .* zs.tauhe ./ (zs.vp+eps) + option.frhe0 .* zs.nem + ...
            max(0,zs.esup_fus) ./ (emean_alpha.* 1.6022e-19).* 2 ./ zs.vp;
    end
    %%  zs.nhem  = zeros(size(zs.vp)) + option.frhe0 .* zs.nem;
elseif all(cons.iso <= 1e-16) || (option.gaz~= 3)
    % works also for gaz = 11
    zs.nhem  = zeros(size(zs.vp)) + option.frhe0 .* zs.nem;
elseif option.transitoire == 1
    if option.evolution == 1
        mem_nhem = zs.nhem;
        nhem0  = option.frhe0 .* zs.nem(2:end) + max(0,zs.esup_fus(2:end)) ./ (3.56e6 .* 1.6022e-19) .* 2 ./ zs.vp(2:end);
        zs.nhem(2:end) =  znalpha0(cons.temps(2:end),zs.salpha(2:end),zs.tauhe(2:end), ...
            max(1,zs.nhem(2:end) - nhem0),zs.vp(2:end)) ./ zs.vp(2:end) + ...
            option.frhe0 .* zs.nem(2:end) +  ...
            max(0,zs.esup_fus(2:end)) ./ (3.56e6 .* 1.6022e-19).* 2 ./ zs.vp(2:end);
        
        %          nhem1  = zs.salpha .* zs.tauhe ./ (zs.vp+eps) + option.frhe0 .* zs.nem + ...
        %                     max(0,zs.esup_fus) ./ (3.56e6 .* 1.6022e-19).* 2 ./ zs.vp;
        
        %        disp([zs.nhem(end-1:end),mem_nhem(end-1:end),nhem1(end-1:end)]./1e19)
        
        %          figure(21);clf
        %          subplot(4,1,1)
        %          plot(cons.temps,zs.nhem,'r+-',cons.temps,mem_nhem,'b-o',cons.temps,nhem1,'*-g');
        %          subplot(4,1,2)
        %          plot(cons.temps,zs.salpha,'-o');
        %          subplot(4,1,3)
        %          plot(cons.temps,zs.tauhe,'-o');
        
    else
        nhem0  = option.frhe0 .* zs.nem + max(0,zs.esup_fus) ./ (3.56e6 .* 1.6022e-19) .* 2 ./ zs.vp;
        zs.nhem  = znalpha0(cons.temps,zs.salpha,zs.tauhe,max(1,zs.nhem - nhem0),zs.vp) ./ zs.vp + option.frhe0 .* zs.nem +  ...
            max(0,zs.esup_fus) ./ (3.56e6 .* 1.6022e-19).* 2 ./ zs.vp;
    end
else
    
    zs.nhem  = zs.salpha .* zs.tauhe ./ (zs.vp+eps) + option.frhe0 .* zs.nem + ...
        max(0,zs.esup_fus) ./ (3.56e6 .* 1.6022e-19).* 2 ./ zs.vp;
end

% securite nhem
% should W be included ? Does W encolde in zu1?
%nhem_max = min((zs.nem  - zu1 .* zs.nimpm) ./ 2.1,zs.nem ./ max(2.1,zs.zeff));
if option.gaz == 5
    nhem_max = min((zs.nem  - zu1 .* zs.nimpm) ./ 2.001, zs.nem .* (zs.zeff - 1) ./ 2.001);
elseif option.gaz == 11
    if ~isfield(zs,'nTm')
        zs.nTm   = cons.iso .* (zs.nem .* ((1 - 2 .* option.frhe0) .* zu2 - (zs.zeff - 4 .* option.frhe0) .* zu1 ) ./ ...
               ((1 + 5.* cons.iso + option.natural_nD_o_nH) .* zu2 - (1 + 25 .*  cons.iso + option.natural_nD_o_nH) .* zu1));
           %disp('why');
    end
    nhem_max = min((zs.nem  - zu1 .* zs.nimpm - 5 .* zs.nTm) ./ 2.1, zs.nem .* (zs.zeff - 1) ./ 2.1);
    
elseif option.gaz == 4
    nhem_max = min((zs.nem  - zu1 .* zs.nimpm) ./ 2.001, zs.nem .* (zs.zeff - 1) ./ 2.001);
else
    nhem_max = min((zs.nem  - zu1 .* zs.nimpm) ./ 2.1, zs.nem .* (zs.zeff - 1) ./ 2.1);
end
nhem_max = max(0,nhem_max);
if any(zs.nhem > nhem_max)
  fprintf('He');
end
zs.nhem = min(zs.nhem,nhem_max);

if isfield(profli,'nip') && (option.gaz ~= 4)
    
    % le contenu en impuret provient du calcul de zeff
    vpp       = trapz(profli.xli,profli.vpr,2);
    zs.nimpm  = trapz(profli.xli,profli.nzp .* profli.vpr,2) ./ vpp;
    % recalcul consitant de HDT
    zs.n1m    =  max( 1e13,zs.nem - zu1 .* zs.nimpm -  2 .* zs.nhem);
    switch option.gaz
        case 1
            zs.nDm   = 1e-2 .* zs.n1m;
            zs.nTm   = zeros(size(zs.n1m));
        case 11
            zs.n1m    = trapz(profli.xli,profli.n1p .* profli.vpr,2) ./ vpp;
            nHm      = zs.n1m ./ (1 + option.natural_nD_o_nH);
            zs.nDm   = option.natural_nD_o_nH .* nHm; %assume non deuterium depleted hydrogen
            zs.nTm   = cons.iso .* nHm;  % this is boron density in this case
        case 5

            % at this stage nhem is known as nimpm, nhem is the density of helium-3
            % nTm and nHm is computed elsewhere 
            zs.n1m = trapz(profli.xli,profli.n1p .* profli.vpr,2) ./ vpp;
            switch option.mino
                case 'H'
                    zs.nDm = zs.n1m ./ (1 + option.cmin) - zs.nTm;   % nTm is computed else where                  
                otherwise
                    zs.nDm = zs.n1m  - zs.nTm;   % nTm is computed elsewhere                                     
            end           
            zs.nhem = zs.nDm .* cons.iso;         
        case 2
            zs.nDm   = zs.n1m;
            zs.nTm   = zeros(size(zs.n1m));           
            switch option.mino
                case 'H'
                    zs.nDm = zs.nDm ./ (1 + option.cmin);
            end
        case 3
            
            switch option.mino
                case 'H'
                    zs.nDm = zs.n1m ./ (1 + option.cmin + cons.iso);
                    zs.nTm = zs.nDm .* cons.iso;
                otherwise
                    zs.nDm   = zs.n1m ./ (1 + cons.iso);
                    zs.nTm   = zs.n1m - zs.nDm;
            end
            
        case 4
            %  		switch option.mino
            %  		case 'H'
            %  			zs.nDm = zs.n1m ./ (1 + option.cmin + cons.iso);
            %  			zs.nTm = zs.nDm .* cons.iso;
            %  		otherwise
            zs.nDm   = zs.n1m ./ (1 + cons.iso);
            zs.nTm   = zs.n1m - zs.nDm;
            %		end
    end
    % densite ionique
    switch option.gaz
        case 5
            zs.nim = zs.n1m + zs.nhe4m  + zs.nhem + zs.nimpm .* (1 + option.rimp) + zs.nwm;             
        case 11
            % in this case nTm stored boron density
            zs.nim = zs.n1m + zs.nhem + zs.nimpm .* (1 + option.rimp) + zs.nwm + zs.nTm; 
        otherwise
            zs.nim = zs.n1m + zs.nhem + zs.nimpm .* (1 + option.rimp) + zs.nwm;
    end
    zs.ni0 = profli.nip(:,1);
    ae     = zs.nim ./ zs.nem;

elseif option.gaz == 4
   	zgaz = 2;
 	zs.zeff  = max(zs.zeff,1.1 .* zgaz);
	% concentration en HDT
	switch option.mino
	case 'H'
		ahdt     = max(0.001,option.cmin + sqrt(eps)/2);
	otherwise
		ahdt     = 0.001;
	end
	zs.n1m   = max(sqrt(eps) .* zs.nem,min(zs.nem ./ 2,nhem_mem) .* ahdt);
	switch option.mino
	case 'H'
		nDTm   = sqrt(eps) .* zs.n1m;
	otherwise
		nDTm   = zs.n1m;
	end
	zs.nDm   = nDTm ./ (1 + cons.iso);
	zs.nTm   = nDTm - zs.nDm;
        % densite de Helium espece principale
	%          (zu2 - zu1) n1~   (-zu2 + zu1 zeff) ne~
        % nhe~ = - --------------- - ---------------------
        %         2 (zu2 - 2 zu1)      2 (zu2 - 2 zu1)
	zs.nhem   = max(1e13,(zs.nem .* (zu1 .* zs.zeff - zu2) + zs.n1m .* (zu2 - zu1)) ./ 2 ./ (2 .* zu1 - zu2));
	% securite nhem
	zs.nhem = min(zs.nhem,zs.nem ./ 2);

	% impurete
	zs.nimpm =  (zs.nem - 2 .* zs.nhem - zs.n1m) ./ zu1;
	% densite ionique
	zs.nim = zs.n1m + zs.nhem + zs.nimpm .* (1 + option.rimp) + zs.nwm;
	zs.ni0 = zs.ne0 .* zs.nim ./ zs.nem;
	ae     = zs.nim ./ zs.nem;

elseif option.gaz == 5
    zgaz = (1 + 2.* cons.iso)./(1 + cons.iso);
    zs.zeff  = max(zs.zeff,1.1 .* zgaz);
    % zs.nHm and zs.nTm are calculated by 
    if option.mino == 'H'
        zs.nHm = max(zs.nHm, zs.nem*option.cmin);
    end
    zs.nDm   = ((zs.nem - 2 .* zs.nhe4m - zs.nTm - zs.nHm ) .* zu2 - ...
                (zs.nem .* zs.zeff - 4 .* zs.nhe4m - zs.nTm- zs.nHm ) .* zu1) ./ ...
               ((1 + 2 .* cons.iso) .* zu2 - (1 + 4 .* cons.iso) .* zu1);
    zs.nhem  = cons.iso .* zs.nDm;
    zs.n1m   = zs.nHm + zs.nDm + zs.nTm;
    zs.nimpm = (zs.nem - zs.n1m - 2 .* zs.nhem - 2 .* zs.nhe4m) ./ zu1 ; % W is encode in zu1
	zs.nim   = zs.n1m + zs.nhem + zs.nhe4m + zs.nimpm .* (1 + option.rimp) + zs.nwm;
	zs.ni0   = zs.ne0 .* zs.nim ./ zs.nem;
	ae       = zs.nim ./ zs.nem;
elseif option.gaz == 11
   	zgaz = 1;
	zs.zeff  = max(zs.zeff,1.1 .* zgaz);
    zs.nimpm = (zs.nem .* cons.zeff - zs.nem  .* (1 +  2 .* option.frhe0)) ./ (zu2-zu1); % first estimation
    ne_free  = zs.nem - 2 .* zs.nhem - zu1 .* zs.nimpm;
    nHm      = ne_free ./ (1 + option.natural_nD_o_nH + 5 * cons.iso);
    zs.nDm   = option.natural_nD_o_nH .* nHm; %assume non deuterium depleted hydrogen
    zs.nTm   = cons.iso .* nHm; % in this case it is the boron density
    zs.n1m   = nHm + zs.nDm;
	% densite ionique
	zs.nim = zs.n1m + zs.nhem + zs.nimpm .* (1 + option.rimp) + zs.nwm + zs.nTm;
	zs.ni0 = zs.ne0 .* zs.nim ./ zs.nem;
	ae     = zs.nim ./ zs.nem;

else
   	zgaz = 1;
	zs.zeff  = max(zs.zeff,1.1 .* zgaz);
	zs.nimpm = (zs.nem .* zs.zeff - zs.nem - 2 .* zs.nhem) ./ (zu2-zu1);
	zs.n1m   = max(1e13,zs.nem - 2  .* zs.nhem - zu1 .* zs.nimpm);

	switch option.gaz
	case 1
		zs.nDm   = 1e-2 .* zs.n1m;
		zs.nTm   = zeros(size(zs.n1m));

	case 2
		zs.nDm   = zs.n1m;
		zs.nTm   = zeros(size(zs.n1m));

		switch option.mino
		case 'H'
			zs.nDm = zs.nDm ./ (1 + option.cmin);
		end
	case 3

		switch option.mino
		case 'H'
			zs.nDm = zs.n1m ./ (1 + option.cmin + cons.iso);
			zs.nTm = zs.nDm .* cons.iso;
		otherwise
			zs.nDm   = zs.n1m ./ (1 + cons.iso);
			zs.nTm   = zs.n1m - zs.nDm;
		end

	end
	% densite ionique
	zs.nim = zs.n1m + zs.nhem + zs.nimpm .* (1 + option.rimp) + zs.nwm ;
	zs.ni0 = zs.ne0 .* zs.nim ./ zs.nem;
	ae     = zs.nim ./ zs.nem;

end

% calcul de la rotation toroidale
if (mode ~= 0) & (option.transitoire ~= 0) 
	[zs.wrad,zs.snbi,void_sicrh,void_sfus,void_sripth, ...
         void_sriplh,void_sripicrh,void_sturb,void_fact,zs.wrot,zs.slh,profli] = ...
                                        z0rot3(zs,option,cons,geo,profli);	

end

%test_electro
%save('sat5','option','cons','geo','zs','amorti','profli');

% donnees derivees du confinement
% solution de l'equipartition
if mode  == 0
   zs.te0     = max(100,min(option.te_max,(2/3) .* zs.wth ./ 1.602176462e-19 ./ zs.vp ./ 2 ./ max(1e13,zs.ne0) .* (1 + zs.ane +zs.ate)));
   zs.tite    = ones(size(zs.te0));
   zs.pei     = zeros(size(zs.te0));
   zs.tauee   = zs.taue;
   zs.tauii   = zs.taue;
end
if fwr == 1; disp('zeqie0');end
if (option.berror > 0) && (mode ~= 0)
      wth_eqei = zs.wth ./ f_wth;
else
      wth_eqei = zs.wth;
end

switch option.gaz
    case 5
      frhe0_loc = option.frhe0;
      nboronm_loc = zeros(size(zs.nTm));
    case 11
      frhe0_loc = 0;
      nboronm_loc = zs.nTm;        
    otherwise
      frhe0_loc = 0;
      nboronm_loc = zeros(size(zs.nTm));      
end
[zs.tite,zs.tauee,zs.tauii,zs.pei,zs.tauei,profli]= zeqie0prof(zs.vp,geo.a,geo.K,zs.te0,zs.nem,zs.pel,zs.pion,wth_eqei,ae,zs.ane,zs.ate,zs.meff, ...
					     cons.temps,zs.tite,zs.pei,zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,...
					     option.zimp,option.rimp,option.zmax, ...
                         zs.pped,zs.tebord,zs.nebord,option.xiioxie + sqrt(-1) .* option.xiioxie_ped, ...
                         option.grad_ped,zs.hitb,option.qdds,zs.indice_inv,profli,fact_confinement, ...
                         option.coef_shape,option.hollow,option.exp_shape,option.kishape, ...
                         option.Sn_fraction,option.te_max,option.extended_qei,frhe0_loc,nboronm_loc,option.min_te_LCFS);

if isfield(profli,'tep')
	zs.te0 = profli.tep(:,1);
	zs.tem = trapz(profli.xli,profli.tep .* profli.vpr,2) ./ trapz(profli.xli,profli.vpr,2);
	zs.ate = (zs.te0 ./ zs.tem) - 1;
	zs.ape  = zs.ane + zs.ate;

else
	if mode == 0
   		te0    = (2/3) .* zs.wth ./ 1.602176462e-19 ./ zs.vp ./ (zs.ne0 + zs.tite .* zs.ni0) .* (1 + zs.ane + zs.ate);
	else
   		te0    = (2/3) .* zs.wth ./ zs.fwcorr;
	end
	% ajout d'une aide a la convergence te0 doit rester < 100 keV et  > 100 eV
	te0    = max(100,min(option.te_max,te0));
	zs.tem    = te0 ./ (1 + zs.ate);
end

% puissance conduite a la separatrice
pl        = max(option.pth_min,zs.pin - zs.prad - zs.pbrem - zs.pcyclo - zs.pioniz);
%figure(21) ;plot(cons.temps,pl);drawnow;
% ref : S.K. Erents, NF vol 28 , 1988,  p 1209
% fraction perdue en volume dans la sol par rayonnement:
switch option.sol_rad
case 'coupled'
	fesol   = max(0,min(1, (zs.pradsol + max(0,1 - option.fprad) .* zs.prad) ./ max(option.pth_min,pl)));
otherwise
	fesol   = max(0,min(1, zs.pradsol ./ max(option.pth_min,pl)));
end
%
% these E. Tsitrone
lclim = pi .* geo.R .* zs.qa;
lcpol = pi .* geo.R;
%lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* 5 .* zs.qa) .^ 2);   % le 5 pour tenir compte du point X
%lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* 7 .* zs.qa) .^ 2);   % le 7 pour tenir compte du point X (valeur etalonnee sur A.S Kukushkin et al NF 43 p 716-723)
lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* option.lcx .* zs.qa) .^ 2);  
switch option.configuration
case 0
	lc = lcpol;
case 1
	lc = lclim;
case 2
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lcpol;
case 3
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lclim;
otherwise
	lc  = lcx;
end
%lc  = zs.modeh .* lcx + (~zs.modeh) .* lclim;


switch option.lambda_scale
case 0
    if option.sol_lscale  == 0
        zs.dsol        = geo.a ./ 100;
    elseif option.sol_lscale  > 0
        zs.dsol        = geo.a .* option.sol_lscale;
    else
        zs.dsol        = - geo.R .* option.sol_lscale;
    end
case 1
    % Goldston model for lambda_q Nucl. Fusion 52 (2012) 013009 (7pp)
    % revue par T. Eich et al. EPS2012
    switch option.gaz
        case {5,11}
            zs.dsol = option.factor_scale .* 5671 .* real(pl) .^ (1/8) .* (1 + (sqrt(zs.sp ./ pi) ./ geo.a) .^ 2) .^ (5/8) .*  ...
                geo.a .^ (17/8) .* geo.b0 .^ (1/4) ./ zs.ip .^(9/8) ./ geo.R .* ...
                (2.* zs.meff ./ zs.zeff);
        otherwise
            zs.dsol = option.factor_scale .* 5671 .* real(pl) .^ (1/8) .* (1 + (sqrt(zs.sp ./ pi) ./ geo.a) .^ 2) .^ (5/8) .*  ...
                geo.a .^ (17/8) .* geo.b0 .^ (1/4) ./ zs.ip .^(9/8) ./ geo.R .* ...
                (2.* zs.meff ./ (1 + double(option.gaz == 4)));
    end
case 2
    if option.sol_lscale  == 0
        dsol_L        = geo.a ./ 100;
    elseif option.sol_lscale  > 0
        dsol_L        = geo.a .* option.sol_lscale;
    else
        dsol_L        = - geo.R .* option.sol_lscale;
    end
	% Goldston model for lambda_q Nucl. Fusion 52 (2012) 013009 (7pp)
        % revue par T. Eich et al. EPS2012
    switch option.gaz
        case {5,11}
             dsol_H = option.factor_scale .* 5671 .* real(pl) .^ (1/8) .* (1 + (sqrt(zs.sp ./ pi) ./ geo.a) .^ 2) .^ (5/8) .*  ...
                geo.a .^ (17/8) .* geo.b0 .^ (1/4) ./ zs.ip .^(9/8) ./ geo.R .* ...
                (2.* zs.meff ./ zs.zeff);
       otherwise
            dsol_H = option.factor_scale .* 5671 .* real(pl) .^ (1/8) .* (1 + (sqrt(zs.sp ./ pi) ./ geo.a) .^ 2) .^ (5/8) .*  ...
                geo.a .^ (17/8) .* geo.b0 .^ (1/4) ./ zs.ip .^(9/8) ./ geo.R .* ...
                (2.* zs.meff ./ (1 + double(option.gaz == 4)));
    end
	% combinaison : Hmode & Xpoint versus  Lmode ou limiter
	flhx    = (zs.modeh .* zs.xpoint);
	zs.dsol = dsol_H .* flhx + dsol_L .* double(~flhx);
case 3 

	% F;D. Halpern et al , Nuclear Fusion 53 (2013) 122001
        % le zeff vient de la formule complete pour cs
        % the width is generally sufficient large and not need any factor.
	if isfield(zs,'tibord')
		dsol_L = 7.22e-8 .* zs.qa .^ (8/7) .* geo.R .^(5/7) .* geo.b0 .^ (-4/7) .*  ...
                 	 zs.tebord .^ (-2/7) .* zs.nebord .^ (2/7) .* (zs.zeff + zs.tibord ./ max(eps,zs.tebord)) .^ (1/7);
	else
		dsol_L = 7.22e-8 .* zs.qa .^ (8/7) .* geo.R .^(5/7) .* geo.b0 .^ (-4/7) .*  ...
                 	 zs.tebord .^ (-2/7) .* zs.nebord .^ (2/7) .* (zs.zeff + 1) .^ (1/7);
    end
    % the scaling is validated up to 0.12 m. There no point with greater width.
    % this limit prevent to run in unphysical regime.
    dsol_L = min(min(0.12,geo.a ./ 4),dsol_L);
    % Goldston model for lambda_q Nucl. Fusion 52 (2012) 013009 (7pp)
    switch option.gaz
        case {5,11}
            % revue par T. Eich et al. EPS2012
            dsol_H = option.factor_scale .* 5671 .* real(pl) .^ (1/8) .* (1 + (sqrt(zs.sp ./ pi) ./ geo.a) .^ 2) .^ (5/8) .*  ...
                geo.a .^ (17/8) .* geo.b0 .^ (1/4) ./ zs.ip .^(9/8) ./ geo.R .* ...
                (2.* zs.meff ./ zs.zeff);
        otherwise
            % revue par T. Eich et al. EPS2012
            dsol_H = option.factor_scale .* 5671 .* real(pl) .^ (1/8) .* (1 + (sqrt(zs.sp ./ pi) ./ geo.a) .^ 2) .^ (5/8) .*  ...
                geo.a .^ (17/8) .* geo.b0 .^ (1/4) ./ zs.ip .^(9/8) ./ geo.R .* ...
                (2.* zs.meff ./ (1 + double(option.gaz == 4)));
    end
    % combinaison : Hmode & Xpoint versus  Lmode ou limiter
    flhx    = (zs.modeh .* zs.xpoint);
    zs.dsol = dsol_H .* flhx + dsol_L .* double(~flhx);
    
    %figure(21);clf;plot(cons.temps,dsol_L,'b',cons.temps,dsol_H,'r',cons.temps,zs.dsol,'.g');
case 4
    % D. Brunner et al 2018 Nucl. Fusion 58 094002, equation 4
    % lambda (mm) = 0.91 <Pressure_atm> ^ -0.48
    % without ITB
    pressure_atm = (2/3) .* zs.w  ./ zs.vp ./ 101.325e3 ./ max(1,zs.hitb);
    zs.dsol = 0.91e-3 .* pressure_atm .^ -0.48;

end

% calcul de la largeur minimale pour la SOL compte tenu des orbites
% project du rayon de larmor perpendiculaire a la direction toroidal
rloc_      = geo.R + geo.a;
rb0_       = geo.R .* geo.b0;
if isfield(zs,'tibord')
	ral_       = 4.57e-3 .* sqrt(zs.meff) .* sqrt(zs.tibord ./ 1e3) ./ (rb0_ ./ rloc_); % en m
else
	ral_       = 4.57e-3 .* sqrt(zs.meff) .* sqrt(zs.tebord ./ 1e3) ./ (rb0_ ./ rloc_); % en m
end
% largeur d'orbite 
if isfield(zs,'qeff')
	dp_     = sqrt(geo.a./geo.R) .* ral_ .* zs.qeff;
	%dp_     = (2 .* zs.qeff .* ral_ ./ rloc_) .^ 2/3;
else
	dp_     = sqrt(geo.a./geo.R) .* ral_ .* 5;
	%dp_     = (2 .* 5 .* ral_ ./ rloc_) .^ 2/3;
end
% Debye
ld_  = 2.35e5 .* sqrt(zs.tebord ./ 1e3 ./ zs.nebord);
% largeur minimale de la SOL
% securite pour limiter la valeur basse 
lambda_min  = max(ral_ ,max(dp_ ,ld_));
%figure(21);clf;plot(cons.temps,lambda_min,'r',cons.temps,zs.dsol,'b');drawnow
zs.dsol = min(geo.a/2,max(zs.dsol,lambda_min));



% calcul avec divertor
% cas du divertor (pour ITER prend Tlim = 15 eV)
% A.S Kukushkin et al NF 43 p 716-723
telim_x = (1 + 0.18 .* edge_eta_c) .* 15 .* ones(size(zs.tebord));
% ref S.K. Erents et al Nuc. Fus. 40 vol 30 p .06
tebord_x = max(telim_x,1.18e-7 .* telim_x .^ 0.2 .* (zs.nebord .* lc .* (1 - 0.5 .* fesol)).^ 0.4);
tebord_x = min(zs.te0 ./ 2,max(option.min_te_LCFS,real(tebord_x)));
% temperature sur le mur ou le divertor
%nusol     = max(10,min(100,1e-16 .* zs.nebord .* lc ./ zs.tebord  .^ 2));
%zs.telim  = zs.tebord ./ 2.3e-3 .* (1 - 0.5 .* fesol) .^ 2 ./ nusol .^ 2;
%zs.telim  = max(1,real(zs.telim));
nelim_x  = 0.5 .* zs.nebord .* tebord_x ./ max(1,telim_x);
nelim_x  = max(1,real(nelim_x));

% model a 2 points
switch option.sol_model
case '2_points'
	[tebord_x,nelim_x,telim_x] = z0convergence_2points_dic(option,cons,geo,zs,profli);
end

if isfield(profli,'tep')
	%figure(21);clf;plot(cons.temps,zs.nebord,'r',cons.temps, (4 .* telim_x .* nelim_x) ./ (option.fmom .* (1 + zs.tibord ./ zs.tebord)) ./ tebord_x,'b.');hold on;drawnow 
%	figure(21);clf;plot(cons.temps,tebord_x,'r',cons.temps,telim_x,'b.');hold on;drawnow 
end

% calcul pout le limiteur
% cf. These de E. Tsitrone p 46-53
% S.K. Erents et al  NF  vol 28 p 1216 (1988)
%Agz         = 0.5;
%Agz         = 2 .* pi .*  geo.R  .* zs.dsol;
%zs.tebord  = (zs.meff .* 1.6726485e-27 ./ 2) .^ (1/3) ./ 1.602176462e-19 .* (pl ./ Agz ./ zs.nebord) .^ (2/3);
% d'apres mesure des sondes de langmuir
% rappport de longueur de decroissance
% compris entre 0.5 et 1
lqoln = 0.62;
% valeur standart de gamma (ajustable entre 5 et 8)
gtr   = 7;
% adapte de Cohen et al ; PPCF 1987
%taup  = 1.3e14 .* geo.R .* geo.a .^ 2 .* geo.K ./ (zs.nbar .^ 0.8);
% a partir de nbord
% formule p 164 de "the plasma boundary of magnetic fusion devices" de Stangeby
% avec taup = taue /3;
%zs.tebord  = tebar .* 3 ./ 7 .* max(0.1,(1 - (zs.prad + zs.pbrem  + zs.pcyclo) ./ zs.pin)) .* taupstaue ./ lqoln ;
%zs.tebord   = min(zs.tem,taup ./ zs.nem ./ zs.vp ./ gtr .* pl ./ lqoln ./ 1.602176462e-19);
switch option.gaz
    case 1
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 1);
    case 11
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
                            (1 + 11 .* cons.iso) .* (1 + cons.iso));
    case 2
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 2);
    case {3,5}
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
                             (2 + 3.* cons.iso) .* (1 + cons.iso));
    case 4
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 4);
end
% limiteur poloidal
%wlim            = zs.peri;
% limteur toroidal (TS)
switch option.configuration
case {0,2}
	wlim            = zs.peri;
otherwise
	wlim            = 2 .* pi .* geo.R ./ zs.qa;
end
zs.tebord  = (pl ./ zs.nebord ./ wlim ./2 ./ gtr ./ zs.dsol ./ 1.602176462e-19 ./cs0) .^ (2/3);	
zs.telim   = zs.tebord; % la temperature est constante le long de la ligne de champ
zs.tebord  = max(option.min_te_LCFS,zs.tebord);
zs.nelim   = 0.5 .* zs.nebord;


% selon zs.modeh
switch option.configuration
case 4
	zs.tebord  = tebord_x;
	zs.telim   = telim_x;
	zs.nelim   = nelim_x;
case {2,3}
	zs.tebord  = zs.xpoint .* tebord_x + (~zs.xpoint)  .* zs.tebord;
	zs.telim   = zs.xpoint .* telim_x  + (~zs.xpoint)  .* zs.telim;
	zs.nelim   = zs.xpoint .* nelim_x  + (~zs.xpoint)  .* zs.nelim;	
otherwise
	% rien a remplir , c'est deja fait
end

% si profils lu en entree
if isappdata(0,'TE_EXP') && isfield(profli,'tep')	 
	zs.tebord = profli.tep(:,end);
end   

% changing  minimal tebord value during breakdown 
if ((option.berror > 0) || any(fact_confinement < 1)) && (mode ~= 0) 
   % tebord must be low during breakdown   
   zs.tebord = option.temp_vac .* 1.3806503e-23./ 1.602176462e-19 .* (1 - fact_confinement) + fact_confinement .* zs.tebord;
   if option.berror == 0
       % now lower limit if the internal breakdown model is switch on
       zs.tebord =  max(option.min_te_LCFS,zs.tebord);
   end
end

% securite
tebordmax  = max(option.min_te_LCFS,zs.wth ./ 1.602176462e-19 ./ zs.nem ./ zs.vp ./ 3);
%figure(21);plot(cons.temps,zs.tebord,'b',cons.temps, tebordmax,'r');drawnow
zs.tebord  = min(tebordmax,option.fte_edge .* real(zs.tebord));
if option.te_edge_fixed > 0
    zs.tebord  = min(tebordmax,option.te_edge_fixed);
end
if isfield(profli,'tep')
	zs.tebord  = min(real(profli.tep(:,end-1)),real(zs.tebord));	
end
zs.telim  = min(zs.tebord,real(zs.telim));

% fit :Pacher et al 1992 
% caclul de la densite de puissance sur le divertor
zs.plim    = max(option.pth_min,zs.pin - zs.prad - zs.pbrem -zs.pcyclo - zs.pradsol - zs.pioniz); 
% nouvelle etalonage pour tenir compte des dernier resultat (Kukushkin NF 43 2003 p 716-723)
zs.peakdiv = 5e29 ./ zs.nebord .^ 1.82 .* (zs.plim ./ zs.sext) .^ 2.37 .* zs.q95 .^ 0.52 .* geo.R .^ 0.33; 
%zs.peakdiv = 1.02e30 ./ zs.nebord .^ 1.82 .* (zs.plim ./ zs.sext) .^ 2.37 .* zs.q95 .^ 0.52 .* geo.R .^ 0.33; 
%zs.peakdiv = min(1.02e30 ./ zs.nebord .^ 1.82 .* (zs.plim ./ zs.sext) .^ 2.37 .* zs.q95 .^ 0.52 .* geo.R .^ 0.33, ...
%                 4.44e13 ./ zs.nebord .^ 0.78 .* (zs.plim ./ zs.sext) .^ 1.56 .* zs.q95 .^ 0.56 .* geo.R .^ 0.56); 
% cette loi se mord la queue (meme htypothese que tebord)
%zs.telim  = max(1,min(zs.tebord,zs.tebord .* 1.8e30 .* (pl./zs.sext) .^ (8/7) .* lc .^ (-6/7) ./ zs.nebord .^ 2));
%McCracken in Plasma Physics de Dendy p 418-420
%zs.telim  = max(1,min(1e5,5.05e16 .^ 2 .* zs.tebord .^ 3 ./ (zs.nebord .* lc) .^ 2));
if ~isfield(zs,'taup')
  zs.taup = zs.taue;
end

% calcul d'une source de zeff
% temperature de "sputtering"
% divertor
switch option.gaz
case   1
	switch option.zeff
	case 4 
		% case W
		ya0 = z0yield('H','W',zs.telim);
	case 3 
		% case Be
		ya0 = z0yield('H','Be',zs.telim);	
	otherwise
		% case C
		ya0 = z0yield('H','C',zs.telim);
	end
	cs          = sqrt(zs.telim .* 1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 1 );
	
case 2
	switch option.zeff
	case 4 
		% case W
		ya0 = z0yield('D','W',zs.telim);
	case 3 
		% case Be
		ya0 = z0yield('D','Be',zs.telim);	
	otherwise
		% case C
		ya0 = z0yield('D','C',zs.telim);
	end
	cs          = sqrt(zs.tebord .* 1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 2 );

case 3
	switch option.zeff
	case 4 
		% case W
		ya0 = z0yield('D','W',zs.telim) ./ (1 + cons.iso) + z0yield('T','W',zs.telim) .* cons.iso ./ (1 + cons.iso);
	case 3 
		% case Be
		ya0 = z0yield('D','Be',zs.telim) ./ (1 + cons.iso) + z0yield('T','Be',zs.telim) .* cons.iso ./ (1 + cons.iso);
	otherwise
		% case C
		ya0 = z0yield('D','C',zs.telim) ./ (1 + cons.iso) + z0yield('T','C',zs.telim) .* cons.iso ./ (1 + cons.iso);
	end
	cs          = sqrt(zs.tebord .* 1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ((2 + 3 .* cons.iso ) ./ (1 + cons.iso)));
    
case 5
	switch option.zeff
        % because there is no data for He3, here we could use He4 for approximation.
	case 4 
		% case W
		ya0 = z0yield2(zs.telim,zs.telim,'D','W') ./ (1 + cons.iso) + z0yield2(zs.telim,zs.telim,'He3','W') .* cons.iso ./ (1 + cons.iso);
	case 3 
		% case Be
		ya0 = z0yield2(zs.telim,zs.telim,'D','Be') ./ (1 + cons.iso) + z0yield2(zs.telim,zs.telim,'He3','Be') .* cons.iso ./ (1 + cons.iso);
	otherwise
		% case C
		ya0 = z0yield2(zs.telim,zs.telim,'D','C') ./ (1 + cons.iso) + z0yield2(zs.telim,zs.telim,'He3','C') .* cons.iso ./ (1 + cons.iso);
	end
	cs          = sqrt(zs.tebord .* 1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ((2 + 3 .* cons.iso ) ./ (1 + cons.iso)));

case 11
	switch option.zeff
	case 4 
		% case W
		ya0 = z0yield2(zs.telim,zs.telim,'H','W') ./ (1 + cons.iso) + z0yield2(zs.telim,zs.telim,'B','W') .* cons.iso ./ (1 + cons.iso);
	case 3 
		% case Be
		ya0 = z0yield2(zs.telim,zs.telim,'H','Be') ./ (1 + cons.iso) + z0yield2(zs.telim,zs.telim,'B','Be') .* cons.iso ./ (1 + cons.iso);
	otherwise
		% case C
		ya0 = z0yield2(zs.telim,zs.telim,'H','C') ./ (1 + cons.iso) + z0yield2(zs.telim,zs.telim,'B','C') .* cons.iso ./ (1 + cons.iso);
	end
	cs          = sqrt(zs.tebord .* 1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ((1 + 11 .* cons.iso ) ./ (1 + cons.iso)));


case 4
	switch option.zeff
	case 4 
		% case W
		ya0 = z0yield('He','W',zs.telim);
	case 3 
		% case Be
		ya0 = z0yield('He','Be',zs.telim);	
	otherwise
		% case C
		ya0 = z0yield('He','C',zs.telim);
	end
	cs          = sqrt(zs.tebord .* 1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 4 );

end


% flux de particule vers le divertor
%  if option.sol_lscale  == 0
%      dsol        = geo.a ./ 100;
%  elseif option.sol_lscale  > 0
%      dsol        = geo.a .* option.sol_lscale;
%  else
%      dsol        = - geo.R .* option.sol_lscale;
%  end
wlim        = 2 .* pi .* geo.R ./ zs.qa;
phipar      = wlim .* zs.nebord .* cs .*  zs.dsol;
phiimp      = ya0 .* phipar;
% selon attenuation
switch option.configuration
case {2,3}
	phiimp = (zs.xpoint .* 0.1 + (~zs.xpoint)) .* phiimp;
case {4}
	phiimp = 0.1 .* phiimp;
otherwise
	phiimp = 1 .* phiimp;
end
% le contenu du plasma (densite moyenne)
nimpm =  zs.tauhe .* phiimp ./ zs.vp;
zgazl  = max(real(cons.zeff),1 + (option.gaz == 4));

% keep charge of main ion even for pB11 and DHe3
switch option.zeff

case 9
    % from I. Erofeev et al 2017 Nucl. Fusion 57 126067, expressed in
    % Greenwlad fraction
    zs.zeffsc  = 1 + (option.gaz == 4) + 0.0016 .* (ngr ./ zs.nbar) .^ 3;
    zs.zeffsc  = min((option.zimp + option.rimp .* option.zmax),zs.zeffsc);
    
case {7,8}
	% scaling law d'apres J. G. Cordey rapport JET-P(85) 28 + modification
	% czeff      =  17.5;
	% zs.zeffsc  =  1 + (option.gaz == 4) + (option.zimp + option.rimp .* option.zmax) .* tanh(zs.negr ./ zs.nbar ./ czeff) .* geo.K .^ 0.3;   
	%czeff_ori      =  zs.sext ./ zs.sp ./ 2 ./ (1 + (~zs.xpoint) ./ 2) ;
	fcz            =  max(0.1,min(10,zs.taup ./zs.taue));
	czeff      =  zs.sext ./ zs.sp  ./ 2 ./ (1 + (~zs.xpoint) ./ 2) ./ fcz;

	%figure(21);clf;plot(cons.temps,czeff_ori,'b',cons.temps,czeff,'r',cons.temps,fcz);drawnow

	zs.zeffsc  =  1 + (option.gaz == 4) + (option.zimp + option.rimp .* option.zmax) .* tanh(zs.negr ./ zs.nbar ./ czeff);   

case 6 
	% sorte de regulation
	prad          = max(0.1 .* zs.pin,zs.prad + zs.pradsol + zs.pbrem) ./ 1e6;
	zs.zeffsc     = (1 + (option.gaz == 4)) + prad ./ 1.6 .* 4.5  ./ (zs.nbar./1e20) .^ 1.89 ./ ...
	                zs.sext .^ 0.94 .* option.zmax .^ 0.12;
	zs.zeffsc     = 0.7 .* zs.zeffsc  + 0.3 .* zs.zeffsc;
case 5  	
	% zeff pour TS 
	% calcul du teref
	if isfield(profli,'tep')
		teref  = trapz(profli.xli,profli.tep .^ 2,2) ./ max(option.min_te_LCFS,trapz(profli.xli,profli.tep,2)) ./ 1e3;
	else
		teref  = zs.tem .* (zs.ate + 0.5) ./ 1e3;
	end
        zs.zeffsc     =  1 + (option.gaz == 4) + (1079 .* teref .^ 1.17 ./ max(1e17,zs.nem) .^ 0.17 +  ...
	                 7.68e9  .* teref  .^ 0.92 ./ max(1e17,zs.nem) .^ 0.51) ./ 2;

case 4 
	zs.zeffsc     =   zgazl + 74 .* (74 - zgazl) .* nimpm ./ zs.nem;
case 3 
	zs.zeffsc     =   zgazl + 4 .* (4 - zgazl) .* nimpm ./ zs.nem;
otherwise
	zs.zeffsc     =   zgazl + 6 .* (6 - zgazl) .* nimpm ./ zs.nem;
end
% ajout effet tungstene
if isfield(profli,'tep')
    if  option.Sn_fraction > 0
        zs.zeffsc     =   zs.zeffsc + trapz(profli.xli, (1 - option.Sn_fraction) .* z0wavez(profli.tep) .^ 2 .* profli.nwp .* profli.vpr,2) ./ ...
			              trapz(profli.xli,profli.nep .* profli.vpr,2) + ...
                          trapz(profli.xli, option.Sn_fraction .* z0snavez(profli.tep) .^ 2 .* profli.nwp .* profli.vpr,2) ./ ...
			              trapz(profli.xli,profli.nep .* profli.vpr,2);
    else
        zs.zeffsc     =   zs.zeffsc + trapz(profli.xli,z0wavez(profli.tep) .^ 2 .* profli.nwp .* profli.vpr,2) ./ ...
            trapz(profli.xli,profli.nep .* profli.vpr,2);
    end
else
        if  option.Sn_fraction > 0
            zs.zeffsc     =   zs.zeffsc + (1 - option.Sn_fraction) .* zs.nwm ./ zs.nem  .* z0wavez(zs.tem) .^ 2 + ...
                              option.Sn_fraction .* zs.nwm ./ zs.nem  .* z0snavez(zs.tem) .^ 2;
        else
            zs.zeffsc     =   zs.zeffsc + zs.nwm ./ zs.nem  .* z0wavez(zs.tem) .^ 2;
        end
end
% securite
zs.zeffsc     =   min(option.zimp  - 0.1 ,max((1 + (option.gaz == 4)) + 0.1,real(zs.zeffsc)));

%figure(21);clf;plot(cons.temps,zs.zeffsc);drawnow

% nan et imag
z0dnanimag
%save('sat6','option','cons','geo','zs','amorti','profli');

% calcul de la generation de courant
ir           = zs.ip - zs.ifwcd - zs.ieccd - real(zs.inbicd) - imag(zs.inbicd) - zs.iboot;
ir(zs.asser > 0) = 0;

%  if isfield(zs,'pnbi_th') && all(isfinite(zs.pnbi_th))
%    pnbi_ref = (1 - option.cur_nbi_time) .*  zs.pnbi +  option.cur_nbi_time .* zs.pnbi_th;
%  else
%    pnbi_ref = zs.pnbi;
%  end
%figure(21);clf;plot(cons.temps,real(zs.pnbi),'b',cons.temps,real(pnbi_ref),'.r',cons.temps,imag(zs.pnbi),'c',cons.temps,imag(pnbi_ref),'.m');drawnow

if fwr == 1; disp('zicd0');end
switch option.icrh_model
  case 'Dumont-Vu'
    option_loc = option;
    option_loc.fwcd = 1 + (option.fabs_fw <= 0);
    [zs.ilh,zs.ifwcd,zs.ieccd,zs.inbicd,zs.ecrit_nbi,zs.taus_nbi,zs.etalh0,zs.etalh1,zs.etalh, ...
             zs.xnbi,zs.piqnbi,zs.frnbi,zs.mu0_nbi,zs.xeccd,xlhout,dlhout,profli,zs.firstorb_nbi] =  ...
             zicd0(cons.temps,zs.plh,puissance_fw,zs.pecrh,zs.pnbi,zs.te0 ./ (1 + zs.ate),zs.nem,zs.zeff, ...
                   geo.R,geo.a,geo.K,zs.RR,zs.vloop,zs.ip,zs.ate,zs.qa,zs.betaptot + zs.li/2, ...
                   option,zs.hmhd,zs.ane,zs.nebord,zs.tebord,cons.ftnbi,zs.meff,zs.taue,geo.b0,zs.nbar,zs.qmin,zs.d0,zs.modeh, ...
		     zs.nDm,zs.nTm,zs.n1m,zs.nhem,cons.xece,zs.efficiency,zs.nimpm,zu1,zu2,zs.inbicd(1),profli,option.Sn_fraction);
  otherwise
    [zs.ilh,zs.ifwcd,zs.ieccd,zs.inbicd,zs.ecrit_nbi,zs.taus_nbi,zs.etalh0,zs.etalh1,zs.etalh, ...
             zs.xnbi,zs.piqnbi,zs.frnbi,zs.mu0_nbi,zs.xeccd,xlhout,dlhout,profli,zs.firstorb_nbi] =  ...
             zicd0(cons.temps,zs.plh,zs.picrh,zs.pecrh,zs.pnbi,zs.te0 ./ (1 + zs.ate),zs.nem,zs.zeff, ...
                   geo.R,geo.a,geo.K,zs.RR,zs.vloop,zs.ip,zs.ate,zs.qa,zs.betaptot + zs.li/2, ...
                   option,zs.hmhd,zs.ane,zs.nebord,zs.tebord,cons.ftnbi,zs.meff,zs.taue,geo.b0,zs.nbar,zs.qmin,zs.d0,zs.modeh, ...
		     zs.nDm,zs.nTm,zs.n1m,zs.nhem,cons.xece,zs.efficiency,zs.nimpm,zu1,zu2,zs.inbicd(1),profli,option.Sn_fraction);
end
% effet du  condinement au breakdown
if (option.berror > 0) && (mode ~= 0)
    zs.ilh    = fact_confinement .* zs.ilh;
    zs.ifwcd  = fact_confinement .* zs.ifwcd;
    zs.ieccd  = fact_confinement .* zs.ieccd;
    zs.inbicd = fact_confinement .* zs.inbicd;
    zs.ifus   = fact_confinement .* zs.ifus;
end

zs.icd   =  zs.ilh + zs.ifwcd + zs.ieccd + real(zs.inbicd) + imag(zs.inbicd) + zs.ifus;
zs.xeccd = ones(size(cons.temps)) .* zs.xeccd;
%  if (mode == 0) & (option.wlh > 0)
%     zs.dlh = ones(size(xlhout));
%     zs.xlh = zeros(size(xlhout));
%  elseif mode == 0
%     zs.dlh = real(dlhout);
%     zs.xlh = real(xlhout);
%  end
if mode == 0
	zs.xlh = real(xlhout);
	zs.dlh = real(dlhout);
end

% calcul du profil de ne
[profli.nep,profli.n0,profli.s0,profli.n0m,profli.s0m, ...
 profli.dn,profli.vn,profli.ware,zs.n0a,profli.spellet,profli.ge,zs.frac_pellet,zs.sn0fr,zs.taup] = ...
 z0profne(cons,zs,profli,option,geo);
 
% effet du  condinement au breakdown
if (option.berror > 0) && (mode ~= 0)
	zs.taup = zs.taup .* fconf_taup;
end
% le piquage reel du profil de densite n'est pas toujours celui demande du fait des autres contraintes
if isfield(profli,'vpr')
    zs.ane_actual = profli.nep(:,1) ./ (trapz(profli.xli,profli.nep .* profli.vpr,2) ./ max(eps,trapz(profli.xli,profli.vpr,2))) - 1;
end
 
% calcul de la puissance perdue du fait de la ionization
if isfield(profli,'nip') && isfield(profli,'tip') && (option.cx_ion ~= 0)
      fcxion = fraction_cx_ion(profli.tep(:,end),profli.tip(:,end),profli.nep(:,end),profli.nip(:,end));
      %figure(21);clf;plot(fcxion);drawnow
      ticx   = profli.tip(:,end);
else
      fcxion  = zeros(size(zs.taup));
      ticx    = option.min_te_LCFS .* ones(size(zs.taup));
end
% cette expression prend aussi en compte la puissance necessaire pour ioniser les glacons
if isfield(profli,'spellet') && isfield(profli,'vpr')
    spe = trapz(profli.xli,profli.vpr .* profli.spellet,2);
    if option.eioniz > 0
        pioniz_e  = (1 - fcxion) .* (max(1,zs.n0a) + spe) .* option.eioniz .* 1.602176462e-19;
    else
        pioniz_e  = (1 - fcxion) .* (max(1,zs.n0a) + spe).* z0eioniz_div(zs.tebord,zs.nebord) .* 1.602176462e-19;
    end
    zs.pioniz_i  = fcxion .* max(1,zs.n0a) .* ticx .* 1.602176462e-19;
    zs.pioniz = pioniz_e + zs.pioniz_i;
    %disp('here')
else
    fpe     = min(1,zs.taup ./ max(zs.taup,zs.taue)) .* zs.frac_pellet;
    switch option.configuration
    case {0,1}
	fpe    =  fpe ./ option.fn0a;
    case {2,3}
	fpe    =  fpe ./ (zs.xpoint .* option.fn0a_div + (~zs.xpoint) .* option.fn0a);
    otherwise
	fpe    =  fpe ./ option.fn0a_div;
    end
    if option.eioniz > 0
	    pioniz_e = (1 - fcxion) .* max(1,zs.n0a) .* option.eioniz .* 1.602176462e-19 ./  max(0.1,1 - fpe);
    else
	    pioniz_e = (1 - fcxion) .* max(1,zs.n0a) .* z0eioniz_div(zs.tebord,zs.nebord) .* 1.602176462e-19 ./  max(0.1,1 - fpe);
    end 
    zs.pioniz_i  = fcxion .* max(1,zs.n0a) .* ticx .* 1.602176462e-19;
    zs.pioniz = pioniz_e + zs.pioniz_i;
end
%figure(21);plot(cons.temps,zs.pioniz,'k',cons.temps,pioniz_e,'r',cons.temps,zs.pioniz_i,'b');drawnow

% calcul de wdia
% le (2/3) vient desobservation sur JET , la pente dans le tanh doit etre ajustee plus finement
switch option.icrh_model
case 'Dumont-Vu'
    % already computed
otherwise
    frpar_icrh  = (1/3) + (1/3) .*  max(0,tanh(zs.einj_icrh ./ zs.ecrit_icrh));
end
frpar_nbi   = (1/3) + (2/3) .* abs(real(zs.mu0_nbi)) .*  max(0,tanh(option.einj ./ max(eps,real(zs.ecrit_nbi)))) + ...
              sqrt(-1) .* ((1/3) + (2/3) .* abs(imag(zs.mu0_nbi)) .*  max(0,tanh(option.einj2 ./ max(eps,imag(zs.ecrit_nbi)))));
zs.wdia = zs.wth + zs.esup_fus + (3/2) .* zs.esup_icrh .* (1 -  frpar_icrh)  + ...
          (3/2) .* real(zs.esup_nbi) .* (1 - real(frpar_nbi)) +  (3/2) .* imag(zs.esup_nbi) .* (1 - imag(frpar_nbi)) + ...
          (1/3) .* zs.esup_lh;
% le (1/3) pour LH provient d'un calcul DKE fait par J.Decker pour TS.	  
%zs.w      = zs.wth + (zs.esup_fus + zs.esup_icrh + zs.esup_nbi + zs.esup_lh);
% fast electron due to runaway electrons addedt to LH terms
if option.runaway ~= 0
    % assuming v// ~ c
    zs.wdia = zs.wdia  + zs.irun .* phys.me .* phys.c ./ phys.e;
end


% calcul du boostrap et de l'ohmique (calcul complet)
if mode == 0
   li = cons.li;
elseif option.limode ~= 0
   li = zs.li;
else
   li = cons.li;
end
if mode == 0
   ipin = max(1,cons.ip);
elseif option.vloop ==0
   ipin = max(1,cons.ip);
else
   if any(zs.ip > (2 .* cons.ip))
	%disp('Metis : Ip >  2 * Ip_refrence in mode vloop');
	fprintf('I');
   end
   if option.evolution == 0
  	 ipin = min(2 .* cons.ip,max(1,zs.ip));
   else
   	 ipin = max(1,cons.ip);	 
   end
end
if option.berror ~= 0
        %cons.temps(cons.ip < 1)
	ipin(cons.ip < 1) = cons.ip(cons.ip < 1);
end
ipin(1) = cons.ip(1);

swlh = abs(zs.plh ./ max(option.pth_min,zs.pohm)) > 0.2;


switch option.icrh_model
case 'Dumont-Vu'
    pfweh = puissance_fw;
    picrh = zs.picrh_th - pfweh;
    tune_frac = 0;
otherwise
    pfweh = max(0,min(1,option.fwcd ~=0)) .* zs.picrh_th + (option.fwcd ==0) .* zs.pel_icrh./2;
    picrh = zs.picrh_th - pfweh;
    tune_frac = option.icrh_width;
end
if option.zeff > 1
	zeffin = zs.zeffsc;
else
	% la consigne est le zeff lineique, il faut travailer avec le volumique moyen
	zeffin = real(cons.zeff) .* zs.zmszl;
end

% nan et imag
z0dnanimag
%save('sat7','option','cons','geo','zs','amorti','profli');
% le probleme de precision numerique arrive apres cette ligne


% directivite
switch option.lhmode
case 0
   directivity = abs(option.etalh);
case 3
   directivity = abs(option.etalh);
case 4
   directivity = abs(option.etalh);
otherwise
   directivity = 0.75;
end
directivity = max(0,min(1,directivity));
xlhin =  option.xlh * ones(size(cons.temps));
dlhin =  option.dlh * ones(size(cons.temps));

% test des champs optionnels pour le couplage avec FREEBIE
if isfield(cons,'flux_ip')
	flux_ip = cons.flux_ip;
else
	flux_ip = NaN * cons.ip;
end
if isfield(cons,'L_ext')
	L_ext = cons.L_ext;
else
	L_ext = NaN * cons.ip;
end
if isfield(cons,'K_min')
	K_min = cons.K_min;
else
	K_min = NaN * cons.ip;
end
if isfield(cons,'d0_ext')
	d0_ext = cons.d0_ext;
else
	d0_ext = NaN * cons.ip;
end
if fwr == 1; disp('zboot0');end
evolution_boot = option.evolution;
if option.evolution == 1
	switch option.dwdt_method
	case 'freebie'
		evolution_boot  = evolution_boot + sqrt(-1);
	end
end

% utilisation de la pression // si anysotrope 
% ref: Equilibrium Analysis of Tokamak Discharges with Anisotropic Pressure, W Zwingmann, L-G Eriksson and P Stubberfield
% W Zwingmann et al 2001 Plasma Phys. Control. Fusion 43 1441 doi:10.1088/0741-3335/43/11/302
switch option.equi_ppar
case 1
    wptot = max(0,zs.w - zs.wdia) + zs.wth; 
case 2
    wptot = zs.w; 
case 3
    wptot = max(0,zs.w - zs.wdia) + zs.wth + zs.wrot; 
otherwise
    wptot = zs.wdia;
end

if (option.refined_ptot == 1) && isfield(profli,'ptot')
   % compute more accurate suprathermal pressure
   ve = ones(size(profli.xli));

   % NBI contribution
   psupra = eps + ((real(zs.esup_nbi) ./ max(1,trapz(profli.xli,real(profli.nbinesource) .* profli.vpr,2))) * ve) .* real(profli.nbinesource);
   psupra = psupra + ((imag(zs.esup_nbi)./ max(1,trapz(profli.xli,imag(profli.nbinesource) .* profli.vpr,2))) * ve) .* imag(profli.nbinesource);

   % ICRH
   psupra = psupra + ((zs.esup_icrh ./ max(1,trapz(profli.xli,profli.picrh .* profli.vpr,2))) * ve) .* real(profli.picrh);

   % fusion
   psupra = psupra + ((zs.esup_fus ./ max(1,trapz(profli.xli,profli.pfus .* profli.vpr,2))) * ve) .* real(profli.pfus);
   
   
   % LH
   psupra = psupra + ((zs.esup_lh ./ max(1,trapz(profli.xli,profli.plh .* profli.vpr,2))) * ve) .* real(profli.plh);

   % rotation part
   switch option.equi_ppar
   case 3
      psupra = psupra + ((zs.wrot ./ max(1,trapz(profli.xli,real(profli.omega) .^ 2 .* profli.vpr,2))) * ve) .* real(profli.omega) .^ 2;   
   end

   % normalisation
   wsupra = max(0,wptot - zs.wth);
   %figure(21);clf;plot(cons.temps,any(psupra > 1 ,2) .* (wsupra ./ max(1,trapz(profli.xli,psupra .* profli.vpr,2))));drawnow
   psupra = psupra .* ((wsupra ./ max(1,trapz(profli.xli,psupra .* profli.vpr,2))) * ve);
   %figure(21);plot(profli.xli,psupra);drawnow;
   profli.ptot = profli.ptot + sqrt(-1) .* psupra;
end

% anti aliasing for ST
qdds_inboot = option.qdds + sqrt(-1) .* option.s1crit;
% if (option.evolution == 1) && (option.qdds < 0 )
%     keyboard
%     % minimal time between to crash
%     if isappdata(0,'MINIMAL_TAU_ST')
%         %%%zs.indice_inv(1:2)
%         if any(zs.indice_inv(1:2) > 1)
%             indlst = find(zs.indice_inv(1:2) > 0,1);
%             %cons.temps(indlst)
%             setappdata(0,'LAST_TIME_ST',cons.temps(indlst));
%             qdds_inboot = 0;
%         end
%         if isappdata(0,'LAST_TIME_ST')
%             %disp([cons.temps(3) - getappdata(0,'LAST_TIME_ST'),getappdata(0,'MINIMAL_TAU_ST')])
%             if getappdata(0,'MINIMAL_TAU_ST') > (cons.temps(3) - getappdata(0,'LAST_TIME_ST'))
%                 qdds_inboot = 0;
%             end
%         end
%     else
%         % must have two time slices between crash
%         if any(zs.indice_inv(1:2) > 0)
%             qdds_inboot = 0;
%         end
%     end
% end

if isfield(zs,'kidds_evol') && (length(zs.kidds_evol) == length(cons.temps))
    % backward compatibility for unitialised call of evolution mode
    if any(~isfinite(zs.kidds_evol))
        zs.kidds_evol(~isfinite(zs.kidds_evol)) = 1;
    end
    kidds_inboot = option.kidds + sqrt(-1) .* zs.kidds_evol;
else
    kidds_inboot = option.kidds;
end

% Hydrogen minoritary
switch option.mino
    case 'H'
        cmin_loc = option.cmin;
    otherwise
        cmin_loc = 0;
end

%save('loc4');
% transport coefficient for turbulent resistivity
if isfield(profli,'dn')
    D0_eta = profli.dn ./ max(eps,max(profli.dn,[],2) * ones(1,size(profli.dn,2)));
else
    D0_eta = 0;
end
[zs.ip,zs.iboot,zs.iohm,zs.pohm,zs.RR,zs.vloop,zs.qa,zs.q95,zs.qmin,zs.q0,zs.betap,zs.piqj, zs.wbp,zs.dwbpdt, ...
 zs.asser,zs.tauj,zs.li,zs.tauip,zs.hitb,zs.xitb,voidate,zs.aitb,zs.hmhd,voidte0,zs.fwcorr,zs.xlh,zs.dlh, ...
 zs.zeff,zs.zmszl,zs.qeff,zs.efficiency,zs.vmes,zs.difcurconv,zs.ipar,void,profli,zs.phiplasma,zs.indice_inv, ...
 zs.poynting,zs.kidds_evol] = ...
        zboot0diff(option.vloop,zs.xpoint,option.vref,cons.temps,zs.nem,zs.tem,zeffin, ...
                   geo.R,geo.a,geo.K,geo.d,geo.b0, ...
                   ipin,zs.icd,li,zs.wth,zs.ane,zs.ate,zs.vp,zs.sp,zs.sext,zs.hitb,option.transitoire,zs.taue,zs.tauip, ...
                   zs.peri,zs.inbicd,zs.xnbi,zs.piqnbi,zs.ieccd,zs.xeccd,zs.ifwcd,zs.ilh,xlhin,dlhin,zs.ifus, ...
                   zs.xfus,zs.jxfus,zs.j0fus,option.runaway,zs.tebord,zs.nebord,swlh,zs.pped,pfweh,picrh, ...
                   zs.plh_th,zs.pecrh,zs.pfus_th,zs.pnbi_th,zs.tite,wptot,zs.nim,option.wlh,option.gaz,option.npar0,option.freqlh,...
		           zs.pcyclo,zs.pbrem,zs.prad,zs.vloop,zu1,zu2,option.zimp,option.zmax,option.rimp,option.zeff,zs.nhem,option.frhe0, ...
		           zs.pion_icrh,zs.pion_nbi,zs.pion_fus,option.kishape + sqrt(-1) .* option.ki_expo ,qdds_inboot, ...
		           kidds_inboot,option.modeboot,Rsepa,Zsepa,amorti, ...
		           directivity,abs(zs.plhrip)./max(1,max(abs(zs.plhrip),zs.plh_th)),zs.drmdt,zs.fracmino, ....
		           zs.xres,zs.pioniz,zs.irun,zs.difcurconv,option.laochange,evolution_boot,cons.flux,option.breakdown,option.sitb, ...
                   zs.indice_inv,zs.meff,option.mode_expo_inte,option.cronos_regul,option.tswitch, ...
                   option.bootmul .* fact_confinement,option.ffit_ped,option.upshiftmode,option.fupshift,option.ddsmode, ....
                   option.w1,option.epsq,option.npar_neg, ...
                   option.itb_sensitivity,option.itb_slope_max,option.faccu,option.xieorkie,option.berror,profli, ...
                   flux_ip,L_ext,K_min,d0_ext,tune_frac,option.width_ecrh,option.hollow, option.moments_mode,option.q0_dds_trig, ...
                   option.betap1crit,option.force_spitzer,option.neutral_friction,option.f_eta_turb .* D0_eta, ...
                   zs.pioniz_i,option.Sn_fraction,option.cor_rel_spitzer,cons.iso,zs.nTm,cmin_loc,option.natural_nD_o_nH,option.collapse);
   
               
% prise en compte du burn-through
%zs.betap = zs.betap .* fact_confinement;
zs.hitb = 1 + (zs.hitb - 1) .* fact_confinement;
%figure(21);clf;plot(cons.temps,zs.hitb,cons.temps,fact_confinement);drawnow

% remplissage de rm et drmdt
zs.rm       = profli.rmx(:,end);
zs.drmdt    = z0dxdt(zs.rm,cons.temps);
if option.evolution == 1
	switch option.dwdt_method
	case 'freebie'
		zs.drmdt    = z0dxdt_freebie(zs.rm,cons.temps);
		%zs.drmdt(:) = polyder(polyfit(cons.temps,zs.rm,1));
	end
end

if option.limode == 0
   zs.li = cons.li;
end
if option.sitb == 0
   zs.hitb = ones(size(cons.temps));
   zs.aitb = ones(size(cons.temps));
   zs.xitb = zeros(size(cons.temps));
%  elseif option.sitb >= 2
%     % ajout du  mode ion chaud
%     % ajout d'une securite pour eviter le declenchement pendant des choc oh
%     hotion  = tanh(max(1,zs.pion) ./ max(1,3 .* zs.pohm)) .*...
%               tanh(zs.tite .* ((zs.pion - zs.pel + zs.pei) ./ max(1,zs.pion + zs.pel)));
%     zs.hitb = zs.hitb + hotion .* (hotion > 0);
end
if option.sitb ~= 3
   zs.hmhd = ones(size(cons.temps));
end


% securite anti depassement de limite sur la temperature (itb trop forte)
zs.hitb = 1 + (zs.hitb - 1) .* min(1,max(0,(1 - tanh((zs.te0 - 1e5) ./ 0.66e4))));

% effet de la mhd
if option.smhd == 0
   smhd = 4 .* zs.li;
elseif option.smhd < 0
   smhd =  abs(option.smhd) .* zs.li;
else
   smhd = option.smhd;
end
fmhd    = (betan0  - smhd/100) ./ betan0.*3;
fmhd    = tanh(fmhd.* (fmhd > 0));
zs.hmhd = min(zs.hmhd,max(0.1,1 - fmhd .* (cons.temps >= option.tmhd)));

% perte de la bariere a proximite de la disruption
zs.hitb   = zs.hitb .* (zs.disrup == 0) + (zs.disrup ~= 0);

% securite 1er temps
if option.evolution == 0
	zs.hitb(1) = 1;
	zs.hmhd(1) = 1;
end

zs.betaptot  = zs.betap .* zs.w ./ zs.wth;
indbad = find((3 .* geo.R ./ geo.a) < zs.betaptot);
if ~isempty(indbad)
    zs.betaptot(indbad) = zs.betap(indbad);
end
zs.ini  = zs.icd + zs.iboot;


% cas grad_ped = 2 
if (option.grad_ped == 3) 
	% nouveau calcul des sources
	profli.source_ion = profli.pnbi_ion + profli.picrh_ion + profli.pfus_ion - profli.pioniz_i;
	profli.source_el  = profli.plh  + profli.pnbi + profli.pecrh + profli.pfweh + profli.picrh + profli.pfus + profli.pohm - ...
                    	    profli.pbrem - profli.prad - (profli.pioniz - profli.pioniz_i)  - profli.pcyclo - profli.source_ion;	
end


% changing  xieshape during breakdown to mimic smaller plasma
if (option.berror > 0) && (mode ~= 0) && isfield(profli,'tep') && (option.kishape ~= 0)
   mask = max(1 / 30,double((x_ioniz * ones(size(profli.xli)) - ones(size(cons.temps)) * profli.xli) < 0));
   mask = mask ./ (mean(mask,2) * ones(size(profli.xli)));
   %figure(21);clf;plot(profli.xli,mask);set(gca,'xlim',[0 3]);set(gca,'ylim',[-0.5 1.5]);drawnow
   profli.xieshape = mask .* profli.xieshape;
   profli.xieshape_itb = mask .* profli.xieshape_itb;
   zs.hitb(x_ioniz < 1) =  max(1 + sqrt(eps),zs.hitb(x_ioniz < 1));
end




%runaway
if (option.runaway > 0) &  (option.transitoire == 1)
	[zs.irun,vl1] = z0irun(cons,geo,zs,profli,option);
	zs.vloop(1)   = vl1(1);
else
	zs.irun = zeros(size(cons.ip));
end

% shafranov shift (formule plasma circulaire + correction Lao + Miller)
if isfield(profli,'Raxe')
	zs.d0  =  profli.Raxe(:,1) - profli.Raxe(:,end);
else
	zs.d0   = geo.a .^ 2  ./ 2 ./ geo.R  .* (2.*(geo.K .^ 2 + 1)./ (3 .* geo.K .^ 2 + 1) .* ...
           ( zs.betaptot + zs.li ./ 2) +  ...
           0.5 .* (geo.K .^ 2 - 1)./ (3 .* geo.K .^ 2 + 1)) .* ( 1 - geo.a ./ geo.R);
end
% nan et imag
z0dnanimag
%save('sat8','option','cons','geo','zs','amorti','profli');
% le probleme de precision numerique arrive avant cette ligne

% calcul de la puissance de fusion (calcul avec les effet de profil et la section efficace complete)
if fwr == 1; disp('zfus0');end
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    ftnbi = 0 .* cons.ftnbi;
else
    switch option.gaz
        case {1,2,4}
            ftnbi = 0 .* cons.ftnbi;
        otherwise
            ftnbi = cons.ftnbi;
    end
end


% here we add new function for pB11 and DHe3
switch option.gaz
    case 5
%         here zs.nhem means density of helium-3,
%         ftnbi means fraction of heium-3 in beam.
        if strcmp(option.mino ,'He3')&& any(cons.iso > 0)
            [pfus0,salpha0,zs.pfus,zs.salpha,ifus,zs.xfus,jxfus,j0fus,zs.taus_he,ecrit_he_void,zs.pfus_nbi,zs.pfus_loss,profli.jfusshape,~,palf0,profli.salf,profli.palf] = ...
                zfus0tae_dhe(zs.nDm,zs.nhem,zs.nTm,zs.tem,zs.nem,zs.zeff,zs.tite,geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,zs.vp,zs.sp, ...
                zs.pnbi_th,zs.taus_nbi,zs.ecrit_nbi,option.einj,option.einj2 ,ftnbi, ...
                zs.pion_icrh,zs.taus_icrh,zs.ecrit_nbi,zs.einj_icrh,cons.temps,cons.pnbi,zs.d0,zs.qa,zs.qmin,zs.te0,zs.nebord,zs.tebord,zs.pped,zs.nim, ...
                zs.wth,option.tae,option.nb_nbi,option.fspot,option.e_shielding,profli,option.fpolarized,option.forced_H_NBI,option.te_max);
        else
            [pfus0,salpha0,zs.pfus,zs.salpha,ifus,zs.xfus,jxfus,j0fus,zs.taus_he,ecrit_he_void,zs.pfus_nbi,zs.pfus_loss,profli.jfusshape,~,palf0,profli.salf,profli.palf] = ...
                zfus0tae_dhe(zs.nDm,zs.nhem,zs.nTm,zs.tem,zs.nem,zs.zeff,zs.tite,geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,zs.vp,zs.sp, ...
                zs.pnbi_th,zs.taus_nbi,zs.ecrit_nbi,option.einj,option.einj2,ftnbi, ...
                zeros(size(zs.pion_icrh)),zs.taus_icrh,zs.ecrit_nbi,zs.einj_icrh,cons.temps,cons.pnbi,zs.d0,zs.qa,zs.qmin,zs.te0,zs.nebord,zs.tebord,zs.pped,zs.nim,...
                zs.wth,option.tae,option.nb_nbi,option.fspot,option.e_shielding,profli,option.fpolarized,option.forced_H_NBI,option.te_max);
        end
        
        % neutron dd
        [zs.ndd,zs.ndd_th,zs.ndd_nbi_th,zs.ndd_nbi_nbi,zs.pddfus,proton_dd,zs.pnbi_icrh,zs.einj_nbi_icrh] = z0neutron_dd(option,cons,zs,profli); 
        Salf_ddNBI.he3 = zs.ndd_nbi_th + zs.ndd_nbi_nbi;
        Salf_ddNBI.p   = proton_dd.nbi_th + proton_dd.nbi_nbi;
        Salf_ddNBI.t   = proton_dd.nbi_th + proton_dd.nbi_nbi;
        pfus_ddNBI.he3 = Salf_ddNBI.he3 .* fus.ddn_he3 .* 1.602176462e-19;
        pfus_ddNBI.p   = Salf_ddNBI.p   .* fus.ddp_p  .* 1.602176462e-19;
        pfus_ddNBI.t   = Salf_ddNBI.t   .* fus.ddp_t  .* 1.602176462e-19;
        
        zs.pfus_he4_DHe3 = pfus0.he4_DHe3;
        zs.pfus_p_DHe3   = pfus0.p_DHe3;
        zs.pfus_he3_DDn  = pfus0.he3_DDn + pfus_ddNBI.he3;
        zs.pfus_t_DDp    = pfus0.t_DDp   + pfus_ddNBI.t;
        zs.pfus_p_DDp    = pfus0.p_DDp   + pfus_ddNBI.p;
        zs.pfus_he4_DT   = pfus0.he4_DT;
        zs.salpha_he4    = salpha0.he4;
        zs.salpha_p      = salpha0.p     + Salf_ddNBI.p;
        zs.salpha_t      = salpha0.t     + Salf_ddNBI.t;
        zs.salpha_he3    = salpha0.he3   + Salf_ddNBI.he3;
        zs.salpha_n      = salpha0.n;
        profli.palf0_he4_DHe3 = palf0.he4_DHe3;
        profli.palf0_p_DHe3   = palf0.p_DHe3;
        profli.palf0_he3_DDn  = palf0.he3_DDn;
        profli.palf0_t_DDp    = palf0.t_DDp;
        profli.palf0_p_DDp    = palf0.p_DDp; 
        profli.palf0_he4_DT   = palf0.he4_DT;
        profli.palf0_tot      = palf0.tot;

    case 11
         [zs.pfus,zs.salpha,ifus,zs.xfus,jxfus,j0fus,zs.taus_he,ecrit_he_void,zs.pfus_nbi,zs.pfus_loss,profli.jfusshape,profli.salf,profli.palf] = ...
 	     zfus0tae_pB11(zs.nDm,zs.nTm,zs.tem,zs.nem,zs.zeff,zs.tite,geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,zs.vp,zs.sp, ...
         zs.pnbi_th,zs.taus_nbi,zs.ecrit_nbi,option.einj,option.einj2 ,ftnbi, ...,
         zs.pion_icrh,zs.taus_icrh,zs.ecrit_nbi,zs.einj_icrh,cons.temps,cons.pnbi,zs.d0,zs.qa,zs.qmin, ...
         zs.te0,zs.nebord,zs.tebord,zs.pped,zs.nim,zs.wth,option.tae,option.nb_nbi,option.fspot,option.e_shielding,profli,option.fpolarized,option.forced_H_NBI,option.mino);
    otherwise
        if strcmp(option.mino ,'T')&& any(cons.iso > 0)
            [zs.pfus,zs.salpha,ifus,zs.xfus,jxfus,j0fus,zs.taus_he,ecrit_he_void,zs.pfus_nbi,zs.pfus_loss,profli.jfusshape,profli.salf,profli.palf] = ...
                zfus0tae(zs.nDm,zs.nTm,zs.tem,zs.nem,zs.zeff,zs.tite,geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,zs.vp,zs.sp, ...
                zs.pnbi_th,zs.taus_nbi,zs.ecrit_nbi,option.einj,option.einj2 ,ftnbi, ...
                zs.pion_icrh,zs.taus_icrh,zs.ecrit_nbi,zs.einj_icrh,cons.temps,cons.pnbi,zs.d0,zs.qa,zs.qmin, ...
                zs.te0,zs.nebord,zs.tebord,zs.pped,zs.nim,zs.wth,option.tae,option.nb_nbi,option.fspot,option.e_shielding,profli,option.fpolarized,option.forced_H_NBI);
        else
            [zs.pfus,zs.salpha,ifus,zs.xfus,jxfus,j0fus,zs.taus_he,ecrit_he_void,zs.pfus_nbi,zs.pfus_loss,profli.jfusshape,profli.salf,profli.palf] = ...
                zfus0tae(zs.nDm,zs.nTm,zs.tem,zs.nem,zs.zeff,zs.tite,geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,zs.vp,zs.sp, ...
                zs.pnbi_th,zs.taus_nbi,zs.ecrit_nbi,option.einj,option.einj2,ftnbi, ...
                zeros(size(zs.pion_icrh)),zs.taus_icrh,zs.ecrit_nbi,zs.einj_icrh,cons.temps,cons.pnbi,zs.d0,zs.qa, ...
                zs.qmin,zs.te0,zs.nebord,zs.tebord,zs.pped,zs.nim,zs.wth,option.tae,option.nb_nbi,option.fspot,option.e_shielding,profli,option.fpolarized,option.forced_H_NBI);
            
        end
end

% calcul du courant induit par les alpha de fusion 
taufus   = min(zs.taus_he,zs.esup_fus ./ max(1,zs.pfus)) .* (zs.pfus > (1e-3 .* max(1,max(zs.ploss,zs.pin))));
zs.ifus  = ifus  .* taufus;
zs.jxfus = jxfus .* taufus;
zs.j0fus = j0fus .* taufus;
zs.ifus  = zs.ifus .* (zs.ifus < zs.ip);

% addding contribution from ICRH fast ions
switch option.icrh_model
    case 'Dumont-Vu'
        if isfield(profli,'jfusshape')
            jfus    = profli.jfusshape;
            jfus    = jfus .* (jfus > 0);
            jfus    = jfus .* ((zs.ifus./  max(1,trapz(profli.xli,profli.spr .* jfus,2))) * ones(size(profli.xli)));
            profli.jfusshape = j_fast_mino + jfus;
            zs.ifus     = zs.ifus + i_fast_mino;
            zs.ifus  = zs.ifus .* (zs.ifus < zs.ip);
        end
end

% correction de la puissance de fusion  
if option.transitoire == 1
	zs.pfus  = zs.pfus + zs.pddfus + zs.pttfus;
end

% neutron dd
switch option.gaz
    case 5
	% already computed
    otherwise
	[zs.ndd,zs.ndd_th,zs.ndd_nbi_th,zs.ndd_nbi_nbi,zs.pddfus,void_proton_dd,zs.pnbi_icrh,zs.einj_nbi_icrh] = z0neutron_dd(option,cons,zs,profli); 
end
% neutron from TT
switch option.gaz
  case {3,5}
    [~,~,~,~,zs.pttfus,~,picrh_nbi_tt,~] = ...
	             z0neutron_tt(option,cons,zs,profli);
    zs.pnbi_icrh = zs.pnbi_icrh +  picrh_nbi_tt;           
otherwise
    zs.pttfus = zeros(size(zs.pfus));
end

% calcul du rayonnement
% if option.gaz == 4
%    zgaz = 2;
% else
%    zgaz = 1;
% end
switch option.gaz
    case 1
        zgaz = 1;
    case 2
        zgaz = 1;
    case 3
        zgaz = 1;
    case 4
        zgaz = 2;
    case 5
        zgaz = mean((1  + 4 .* cons.iso)   ./  (1 + 2 .* cons.iso));
    case 11
        zgaz = mean((1  + 25 .* cons.iso)   ./ (1 + 5.* cons.iso));
end

% selon la configuration
switch option.configuration
case {2,3,4}
	telim       = zs.xpoint .* zs.telim;
	nesol_ave   = (1 - option.fnesol) .* zs.nebord + option.fnesol .* zs.nelim;	
otherwise
	telim       = 0 .* zs.telim;
	nesol_ave   = zs.nebord;	
end
%  switch option.gaz
%  case   1
%  	ya0_w = z0yield('H','W',zs.telim);
%  case 2
%  	ya0_w = z0yield('D','W',zs.telim);
%  case 3
%  	ya0_w = z0yield('D','W',zs.telim) .* (1 -cons.iso) + z0yield('T','W',zs.telim) .* cons.iso;
%  case 4
%  	ya0_w = z0yield('He','W',zs.telim);
%  end
%  ya0_w = abs(option.faccu .* ya0_w);

switch option.gaz
    case 11
        nboronm = zs.nTm;
    otherwise
        nboronm = zeros(size(zs.nem));
end

if fwr == 1; disp('zrad0');end

[zs.prad,zs.pbrem,zs.pcyclo,zs.pradsol,profli] = zrad0(zs.nem,zs.tem,zs.zeff ./ zs.zmszl,option.zmax,option.zimp,option.rimp,zgaz, ...
                                            geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,option.rw,zs.vp,zs.sext, ...
					    (option.modeh > 0),option.frad,zs.te0,zs.nbar,option.matthews, ...
                                            telim,nesol_ave,zs.taue + sqrt(-1) .* zs.taup,zs.dsol,0,option.z_prad,option.sol_rad, ...
                                            option.gaunt,option.noncoronal,cons.temps,profli,option.Sn_fraction,option.te_max, ...
                                            option.cor_rel_brem,nboronm);


% adding gas puff ionization 
if any(gas_puff_mem > 0)
    if option.eioniz > 0
	     zs.pradsol = zs.pradsol + max(0,1 - option.eta_gas_puff) .* gas_puff_mem .*  option.eioniz .* 1.602176462e-19;
    else
	     zs.pradsol = zs.pradsol + max(0,1 - option.eta_gas_puff) .* gas_puff_mem .* z0eioniz_div(zs.telim,zs.nelim) .* 1.602176462e-19;
    end 
end


%  [void.prad,void.pbrem,void.pcyclo,void.pradsol] = zrad0(zs.nem,zs.tem,zs.zeff ./ zs.zmszl,option.zmax,option.zimp,option.rimp,zgaz, ...
%                                              geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,option.rw,zs.vp,zs.sext, ...
%  					    (option.modeh > 0),option.frad,zs.te0,zs.nbar,option.matthews, ...
%                                              telim,nesol_ave,zs.taue,zs.dsol,0,option.z_prad,option.sol_rad,option.gaunt,0,cons.temps,profli,option.Sn_fraction);
%  
%  figure(21);clf;
%  subplot(3,1,1);
%  plot(cons.temps,zs.prad,'r',cons.temps,void.prad,'b');
%  subplot(3,1,2);
%  plot(cons.temps,zs.pbrem,'r',cons.temps,void.pbrem,'b');
%  subplot(3,1,3);
%  plot(cons.temps,zs.pradsol,'r',cons.temps,void.pradsol,'b');
%  drawnow

% calcul de la disruption radiative
%pradtot   = option.fprad .* zs.prad + zs.pcyclo + zs.pbrem;
pradtot   = zs.prad + zs.pcyclo + zs.pbrem;
zs.disrup = (pradtot > (zs.pin + zs.dwdt .* (zs.dwdt >0))) | ((zs.pel-min(0,zs.pei)) <= (0.01 .* zs.pin));
zs.disrup = double(zs.disrup);

% en mode limiteur le rayonnement dans la sol est augmente par la distribution de densite	
% securite
switch option.ploss_exp
case 'no_prad'
	zs.pradsol = max(1,min(0.99 .* zs.pin,zs.pradsol));
	zs.prad = max(1,min(pi .* zs.pin,zs.prad));
case 'max_power'
	zs.pradsol = max(1,min(0.99 .* zs.ploss,zs.pradsol));
	zs.prad = max(1,min(pi .* zs.pin,zs.prad));
case 'max(pel)+max(pion)'
 	zs.pradsol = max(1,min(0.99 .* zs.ploss,zs.pradsol));
	zs.prad = max(1,min(pi .* zs.pin,zs.prad));
otherwise
	zs.prad = max(1,min(0.99 .* zs.pin,zs.prad));
	zs.pradsol = max(1,min(0.99 .* zs.ploss,zs.pradsol));
end

% pour maintenir la coherence
if isfield(profli,'vpr') && isfield(profli,'fprad')
	prad_in       = trapz(profli.xli,profli.fprad .* profli.vpr,2);
	profli.prad   = profli.fprad .* ((zs.prad ./max(eps,prad_in)) * ones(1,length(profli.xli)));
	profli.pcyclo = profli.pcyclo .* ((zs.pcyclo ./ max(1,trapz(profli.xli,profli.pcyclo .* profli.vpr,2))) * ones(1,length(profli.xli)));
end

% calcul de la qualite de l'alignement des courants
% ref : M.R. Wade NF 43 (2003) p 634-
if isfield(profli,'jni') & isfield(profli,'nep') & isfield(profli,'tep')
	zs.ialign = 1 - trapz(profli.xli,profli.nep ./ max(1,profli.tep) .* abs(profli.jeff - profli.jni),2) ./ ...
      			max(eps,trapz(profli.xli,profli.nep ./ max(1,profli.tep) .* abs(profli.jeff),2));
	zs.ialign = min(1,max(0,zs.ialign));
end

% donnees 0D extraites des profiles pour compatibilite ITM
% (donnees redondnates) 
zs.tibord     = profli.tip(:,end);
zs.nibord     = profli.nip(:,end);
zs.xped       = 0.95 .* zs.modeh + (~zs.modeh);
zs.teped      = profli.tep(:,end - 1) .* zs.modeh + profli.tep(:,end ) .* (~zs.modeh);
zs.tiped      = profli.tip(:,end - 1) .* zs.modeh + profli.tip(:,end ) .* (~zs.modeh);
zs.neped      = profli.nep(:,end - 1) .* zs.modeh + profli.nep(:,end ) .* (~zs.modeh);
zs.niped      = profli.nip(:,end - 1) .* zs.modeh + profli.nip(:,end ) .* (~zs.modeh);
zs.edgeflux   = 2 .* pi .* profli.psi(:,end);

if (mode == 0) | (option.transitoire == 0) 
	[zs.wrad,zs.snbi,void_sicrh,void_sfus,void_sripth, ...
         void_sriplh,void_sripicrh,void_sturb,void_fact,zs.wrot,zs.slh,profli] = ...
                                        z0rot3(zs,option,cons,geo,profli);
end


% taue definition donnees dans :
% Collisional transport in magnetized plasma , P. Helander, and D.J. Sigmar, Cambridge University Press, 2002 (p 164)
paux   = (zs.pfus + zs.picrh + zs.plh + real(zs.pnbi) + imag(zs.pnbi) + zs.pecrh) - (zs.pbrem + zs.pcyclo + zs.prad);
iv     = zs.ip .* profli.vpr(:,end) ./ profli.spr(:,end) .* profli.epar(:,end);
bpol   = sqrt( profli.grho2r2) .* abs(pdederive(profli.xli,profli.psi,0,2,2,1)) ./ (profli.rmx(:,end) * ones(size(profli.xli)));
xcor   = profli.rmx(:,end) .* trapz(profli.xli, ...
          bpol .* pdederive(cons.temps,bpol,2,2,1,1) + ...
         (profli.fdia - (profli.fdia(:,end) *  ones(size(profli.xli))))  .* profli.ri.* ...
         pdederive(cons.temps,profli.fdia .* profli.ri,2,2,1,1) ,2) ./ (4 .* pi .* 1e-7);
%
dt      = diff(cons.temps);
dt      = cat(1,dt(1),dt);
% ? tt      = (min(cons.temps):dt:max(cons.temps))';
u       = paux + iv - xcor;
dw      = zs.dwdt .* dt;
span    = cat(2,0.1 - 9.9 / 1000,linspace(0.1,10,1001),10 + 9.9 / 1000);
taueref = zs.taue;
tauetab = mean(cat(2,taueref,cat(1,taueref(2:end),taueref(end)),cat(1,taueref(1),taueref(1:end-1))),2) * span;
vs = ones(1,size(tauetab,2));
dd = abs(dw * vs  - (zs.w * vs - tauetab .* (u *vs)) .* (exp(- (dt * vs) ./ tauetab) -1));
ddmin = min(dd,[],2);
taue_alt = sum(tauetab .* (dd == (ddmin *vs)),2) ./ max(1,sum(dd == (ddmin *vs),2));
taue_alt(taue_alt == max(tauetab,[],2)) = 0;
taue_alt(taue_alt == min(tauetab,[],2)) = 0;
zs.taue_alt =  taue_alt;

% nan et imag
z0dnanimag
%save('sat9','option','cons','geo','zs','amorti','profli');
%keyboard