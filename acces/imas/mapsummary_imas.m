%% --------------------------------------------------------
%% CONVERT DATA FROM METIS TO SUMMARY IDS (MAPPING)
%% --------------------------------------------------------
function summary = mapsummary_imas ...
    (z0dstruct,data_zerod,profil0d,texte_diary ...
     ,error_flag,summary,run,occurrence,sigma_B0_eff,sigma_bvac_r)
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

% Sn
if ~isfield(z0dstruct.z0dinput.option,'Sn_fraction')
    z0dstruct.z0dinput.option.Sn_fraction = 0;
end

%% FOR METIS TIME-INDEPENDENT VARIABLES WHICH BECOME TIME-DEPENDENT IN THE IDS
ntime = length(data_zerod.temps);
vt = ones(size(data_zerod.temps));

%% FUELLING CALCULATION
fuelling_data = metis_gaz_imas(z0dstruct,data_zerod,profil0d);

% precomputation
amin                   = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.a,profil0d.temps,'pchip','extrap');
magnetic_shear         = pdederive(profil0d.xli,profil0d.qjli,0,2,2,1) ./ profil0d.qjli .* (ones(size(profil0d.qjli,1),1) * profil0d.xli);
%
% densite ioniques
%
ve    = ones(size(profil0d.xli));
nDm   = interp1_imas(data_zerod.temps,data_zerod.nDm,profil0d.temps,'pchip','extrap');
nTm   = interp1_imas(data_zerod.temps,data_zerod.nTm,profil0d.temps,'pchip','extrap');
nDp   = max(1,profil0d.n1p .* ((nDm./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
nTp   = max(1,profil0d.n1p .* ((nTm./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
switch z0dstruct.z0dinput.option.gaz
    case 11
        nbp   = nTp;
        nTp   = zeros(size(nTp));
        nHp   = max(1e13,profil0d.n1p - nDp);
    otherwise
        nbp   = zeros(size(nTp));
        nHp   = max(1e13,profil0d.n1p - nTp - nDp);
end
switch z0dstruct.z0dinput.option.gaz
    case 5
        nhep3  = profil0d.nhep;
        nhep   = option.frhe0 .* profil0d.nep;
    otherwise
        nhep  = max(1e13,profil0d.nhep);
        nhep3 = zeros(size(profil0d.nhep));
end
%nHp   = max(1,profil0d.n1p - nTp - nDp);
%nhep  = max(1,profil0d.nhep);
nz1p  = max(1,profil0d.nzp);
nz2p  = max(1,profil0d.nzp .* z0dstruct.z0dinput.option.rimp);   
if  z0dstruct.z0dinput.option.Sn_fraction > 0
    nwp   = (1 - z0dstruct.z0dinput.option.Sn_fraction) .* max(1,profil0d.nwp);
    nsnp  = z0dstruct.z0dinput.option.Sn_fraction .* max(1,profil0d.nwp);
else
    nwp   = max(1,profil0d.nwp);
    nsnp  = zeros(size(nwp)); 
end
%
switch z0dstruct.z0dinput.option.mino
    case 'He3'
        switch z0dstruct.z0dinput.option.gaz
            case 4
                nHe3m = z0dstruct.z0dinput.option.cmin .* data_zerod.nhem;
                nHem  = max(0,data_zerod.nhem - nHe3m);
            case 5
                nHe3m = data_zerod.nhem;
                nHem  = z0dstruct.z0dinput.option.frhe0 .* data_zerod.nem;
                
            otherwise
                nHe3m = z0dstruct.z0dinput.option.cmin .* data_zerod.n1m;
                nHem  = max(0,data_zerod.nhem - nHe3m);
        end
    otherwise
        nHem  = data_zerod.nhem;
        nHe3m = 0 .* nHem;
end
frhe3  = nHe3m ./ max(1e11,nHe3m + nHem);
frhe3  = interp1_imas(data_zerod.temps,frhe3,profil0d.temps,'pchip','extrap') * ve;
nions  = NaN * ones(length(profil0d.temps),length(profil0d.xli),10);
% switch z0dstruct.z0dinput.option.mino
% case 'He3'
% 	switch z0dstruct.z0dinput.option.gaz
% 	case 4
% 		nHe3m = z0dstruct.z0dinput.option.cmin .* data_zerod.nhem;
% 		nHem  = max(0,data_zerod.nhem - nHe3m);
% 	otherwise
% 		nHe3m = z0dstruct.z0dinput.option.cmin .* data_zerod.n1m;
% 		nHem  = max(0,data_zerod.nhem - nHe3m);
% 	end
% otherwise
% 	nHem  = data_zerod.nhem;
% 	nHe3m = 0 .* nHem;
% end
% frhe3  = nHe3m ./ max(1e11,nHe3m + nHem);
% frhe3  = interp1_imas(data_zerod.temps,frhe3,profil0d.temps,'pchip','extrap') * ve;
% if  z0dstruct.z0dinput.option.Sn_fraction > 0
%     nions  = NaN * ones(length(profil0d.temps),length(profil0d.xli),9);
% else
%     nions  = NaN * ones(length(profil0d.temps),length(profil0d.xli),8);
% end
nions(:,:,1) = nHp;
nions(:,:,2) = nDp;
nions(:,:,3) = nTp;
switch z0dstruct.z0dinput.option.gaz
   case 5
        nions(:,:,4) = nhep3;
        nions(:,:,5) = nhep;
   otherwise
       nions(:,:,4) = nhep .* frhe3;
       nions(:,:,5) = nhep .* (1 - frhe3);
       %nhep3        = nhep .* (1 - frhe3);       
       %nhep         = nhep .* frhe3;
end
nions(:,:,6) = nz1p;
nions(:,:,7) = nz2p;
nions(:,:,8) = nwp;
nions(:,:,9) = nsnp;  
nions(:,:,10) = nbp;

for k = 1:size(nions,3)
  nions_vol(:,k) = trapz(profil0d.xli,profil0d.vpr .* squeeze(nions(:,:,k)),2) ./ trapz(profil0d.xli,profil0d.vpr,2);
end

% ion masse and charge table
% H
A(1) = 1; 
Z(1) = 1;
label{1} = 'H';
% D
A(2) = 2; 
Z(2) = 1;
label{2} = 'D';
% T
A(3) = 3; 
Z(3) = 1;
label{3} = 'T';
% He3
A(4) = 3; 
Z(4) = 2;
label{4} = 'He3';
% He4
A(5) = 4; 
Z(5) = 2;
label{5} = 'He4';
% imp 1
Z(6) = z0dstruct.z0dinput.option.zimp;
[A(6),label{6}] = getafromz_imas(Z(6));
% imp 2
Z(7) = z0dstruct.z0dinput.option.zmax;
[A(7),label{7}] = getafromz_imas(Z(7));
% W
Z(8) = 74;
[A(8),label{8}] = getafromz_imas(Z(8));
% Sn
Z(9) = 50;
[A(9),label{9}] = getafromz_imas(Z(9));
% boron
Z(10) = 5;
A(11) = 11;
label{10} = 'B';

% rotation 
[rtor,vtor,vpol,omega,mtor] = z0rot_imas(data_zerod,profil0d,z0dstruct.z0dinput.option,frhe3,z0dstruct.z0dinput.geo,z0dstruct.z0dinput.cons);

% fast particles densities
[psupra,psupra_para,psupra_perp,nfast] = nfast4imas(data_zerod,profil0d,z0dstruct.z0dinput.option,frhe3,z0dstruct.z0dinput.geo,z0dstruct.z0dinput.cons);

% neutrals
% cold first and hot second
switch z0dstruct.z0dinput.option.gaz
    case 1
        amn = [1,1];
        zn  = [1,1];
        ion_index = [1,1];
        neutlab   = {'H2','H2'};
    case 2
        amn = [2,2];
        zn  = [1,1];
        ion_index = [2,2];
        neutlab   = {'D2','D2'};
    case 3
        amn = [2,3,2,3];
        zn  = [1,1,1,1];
        ion_index = [2,3,2,3];
        neutlab   = {'D2','T2','D2','T2'};
    case 4
        amn = [2,2];
        zn  = [4,4];
        ion_index = [5,5];
        neutlab   = {'He4','He4'};
    case 5
        amn = [2,3,2,3];
        zn  = [1,2,1,2];
        ion_index = [2,4,2,4];
        neutlab   = {'D2','He3','D2','He3'};
        
    case 11
        amn = [1,11,1,11];
        zn  = [1,5,1,5];
        ion_index = [1,10,1,10];
        neutlab   = {'H2','B','H2','B'};
        
    otherwise
        error('plasma composition not yet implemented in METIS4IMAS');
end
%
n_neutrals  = NaN * ones(length(profil0d.temps),length(profil0d.xli),length(amn));
t_neutrals  = NaN * ones(length(profil0d.temps),length(profil0d.xli),length(amn));
% utilisation de la vitesse du son (plus stable numeriquement, generalement peu differente)
t_hot  = (profil0d.tip + profil0d.tep .* profil0d.zeff)./ 2;
t_cold = interp1_imas(data_zerod.temps,data_zerod.telim,profil0d.temps,'pchip','extrap') * ones(size(profil0d.xli));

%
switch z0dstruct.z0dinput.option.gaz
    case {3,5}
        iso = interp1_imas(z0dstruct.z0dinput.cons.temps,real(z0dstruct.z0dinput.cons.iso),profil0d.temps,'pchip','extrap');
        iso = iso * ones(size(profil0d.xli));
        n_neutrals(:,:,1)   = profil0d.n0m ./ (1+iso);
        t_neutrals(:,:,1)   = t_cold;
        n_neutrals(:,:,3)   = profil0d.n0  ./ (1+iso);
        t_neutrals(:,:,3)   = t_hot;
        n_neutrals(:,:,2)   = profil0d.n0m ./ (1+iso) .* iso;
        t_neutrals(:,:,2)   = t_cold;
        n_neutrals(:,:,4)   = profil0d.n0  ./ (1+iso) .* iso;
        t_neutrals(:,:,4)   = t_hot;
        % take into account molecule (for H and D) and number of charge for He4.
        n_neutrals = n_neutrals ./ 2;
        if z0dstruct.z0dinput.option.gaz == 5
            warning('nHe3onD & nTonD not yet implemented !');
        end
    case 11
        iso = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.cons.iso,profil0d.temps,'pchip','extrap');
        iso = iso * ones(size(profil0d.xli));
        n_neutrals(:,:,1)   = profil0d.n0m ./ (1+iso) / 2;
        t_neutrals(:,:,1)   = t_cold;
        n_neutrals(:,:,3)   = profil0d.n0  ./ (1+iso) / 2;
        t_neutrals(:,:,3)   = t_hot;
        n_neutrals(:,:,2)   = profil0d.n0m ./ (1+iso) .* iso / 5;
        t_neutrals(:,:,2)   = t_cold;
        n_neutrals(:,:,4)   = profil0d.n0  ./ (1+iso) .* iso / 5;
        t_neutrals(:,:,4)   = t_hot;
    otherwise
        n_neutrals(:,:,1)   = profil0d.n0m;
        t_neutrals(:,:,1)   = t_cold;
        n_neutrals(:,:,2)   = profil0d.n0;
        t_neutrals(:,:,2)   = t_hot;
        % take into account molecule (for H and D) and number of charge for He4.
        n_neutrals = n_neutrals ./ 2;
end
% take into account molecule (for H and D) and number of charge for He4.
%%n_neutrals = n_neutrals ./ 2;

% compute neutron rates
[neutron.time,neutron.ntot,neutron.ndd,neutron.ndd_th,neutron.ndd_nbi_th,neutron.ndd_nbi_nbi, ...
 neutron.ndt,neutron.ndt_th,neutron.ndt_nbi_th,neutron.ndt_nbi_nbi,neutron.ndt_nbi_icrh, ...
 neutron.ntt,neutron.ntt_th,neutron.ntt_nbi_th,neutron.ntt_nbi_nbi] = imas_neutron(z0dstruct);

% compute fusion power from neutrons
[fusion.pfusion_total,fusion.pfusion_DT,fusion.pfusion_DD,fusion.pfusion_TT,fusion.ntt, ...
          fusion.frac_dt_thermal,fusion.frac_dt_beam_plasma,fusion.frac_dt_beam_beam,fusion.frac_dt_icrh, ...
          fusion.frac_dd_thermal,fusion.frac_dd_beam_plasma,fusion.frac_dd_beam_beam, ...
          fusion.frac_tt_thermal,fusion.frac_tt_beam_plasma,fusion.frac_tt_beam_beam]  =  ...
          compute_total_fusion_power(z0dstruct);

% done
configs.config.value = 0;
switch z0dstruct.z0dinput.option.configuration
    case 0
        configs.config.value = 2;
    case 1
        configs.config.value = 4;
    case 2
        configs.config.value = 2 + 5 .* double(any(data_zerod.modeh>0));
    case 3
        configs.config.value =  4 + 3 .* double(any(data_zerod.modeh>0));
    case 4
        configs.config.value = 7;
end
      
configs.ecrh_freq.value     = [];
configs.ecrh_loc.value      = mean(z0dstruct.z0dinput.cons.xece) ;
configs.ecrh_mode.value     = []; 
configs.ecrh_tor_ang.value  = (z0dstruct.z0dinput.option.sens .* (pi / 8) ) * vt; 
configs.ecrh_pol_ang.value  = (z0dstruct.z0dinput.option.angle_ece ./180 .* pi) * vt;
configs.ecrh_harm.value     = [];
configs.enbi.value          = z0dstruct.z0dinput.option.einj * vt ;
configs.r_nbi.value         = z0dstruct.z0dinput.option.rtang * vt ;
configs.grad_b_drift.value  = 1;
configs.icrh_freq.value     = z0dstruct.z0dinput.option.freq *vt ;
configs.icrh_phase.value    = [];

if z0dstruct.z0dinput.option.fwcd == 2	
  configs.icrh_scheme = 'FW';	
elseif z0dstruct.z0dinput.option.fwcd == -1
  configs.icrh_scheme = 'FW_CCD';	
elseif z0dstruct.z0dinput.option.fwcd == 1
  configs.icrh_scheme = 'FW_CD';	
else
  configs.icrh_scheme = sprintf('%s_min_%d',z0dstruct.z0dinput.option.mino,ceil(mean(data_zerod.harm)));
end

configs.LH_freq.value    = (z0dstruct.z0dinput.option.freqlh .* 1e9) * vt ;
configs.LH_npar.value	  = z0dstruct.z0dinput.option.npar0 *vt; 
configs.pellet_ang.value = [];
configs.pellet_v.value   = [];
configs.pellet_nba.value = [];

% champs manquant
switch z0dstruct.z0dinput.option.scaling
    case 0
        configs.lmode_sc 	= 'ITERH-96P(th)';
        configs.hmode_sc 	= 'ITERH-98P(y,2)';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
        configs.pedestal_sc  	= '(W from ITERH-98P(y,2) - W from ITERH-96P(th)) / 2';
    case 1
        configs.lmode_sc 	= 'OH scaling law (G. Bracco and K. Thomsen, NF, 37 ,1997)';
        configs.hmode_sc 	= 'ITERH-98P(y,2)';
        configs.core_sc  	=  'W from L-mode or H-mode - Wpedestal';
        configs.pedestal_sc  	='(W from ITERH-98P(y,2) - W from ITERH-96P(th)) / 2';
    case 2
        configs.lmode_sc 	= 'W from core scaling from McDonald NF 47 (2007) p147';
        configs.hmode_sc 	= 'W from core scaling  +W from  pedestal scaling from McDonald NF 47 (2007) p147';
        configs.core_sc  	= 'core scaling from McDonald NF 47 (2007) p147';
        configs.pedestal_sc  	= 'pedestal scaling from McDonald NF 47 (2007) p147';
    case 3
        configs.lmode_sc 	= 'W form GS03 - W from pedestal scaling';
        configs.hmode_sc 	= 'GS03 (ITPA @AIEA 2004 from McDonald  PPCF 46 ,2004,A215- , #11';
        configs.core_sc  	= 'W form GS03 - W from pedestal scaling';
        configs.pedestal_sc  	= 'T. Hoang for ITPA 2004';
    case 4
        configs.lmode_sc 	= 'fit of experimental measurement of Wdia';
        configs.hmode_sc 	= 'fit of experimental measurement of Wdia';
        configs.core_sc  	= 'fit of experimental measurement of Wdia - Wpedestal';
        configs.pedestal_sc  	= '(W from ITERH-98P(y,2) - W from ITERH-96P(th)) / 2';
    case 5
        configs.lmode_sc 	= 'ITERH-96P(th)';
        configs.hmode_sc 	= 'ITERH-EIV(y,2)  (Mc Donald, NF 47 (2007) 147-174)';
        configs.core_sc  	= 'difference between Wtot and Wpedestal';
        configs.pedestal_sc  	= 'W from ITERH-EIV(y,2) - W from ITERH-96P(th)';
    case 6
        configs.lmode_sc 	= 'OH scaling law from Wesson Tokamak with transition to ITERH-96P(th) when additional heating is applied';
        configs.hmode_sc 	= 'ITERH-98P(y,2)';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
        configs.pedestal_sc  	='(W from ITERH-98P(y,2) - W from ITERH-96P(th)) / 2';
    case 7
        configs.lmode_sc 	= 'ITERH-98P(y,2)/2';
        configs.hmode_sc 	= 'ITERH-98P(y,2)';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
        configs.pedestal_sc  	= 'ITERH-98P(y,2)/2';
    case 8
        configs.lmode_sc 	= 'user defined scaling';
        configs.hmode_sc 	= 'user defined scaling';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
        configs.pedestal_sc  	= 'user defined scaling';
    case 9
        configs.lmode_sc 	= 'ITERH-96P(th)';
        configs.hmode_sc 	= 'J_pol = 0 (Hybrid or advance summary)';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
        configs.pedestal_sc  	= 'difference between H mode and L mode';
    case 10
        configs.lmode_sc 	= 'Sauter & Martin';
        configs.hmode_sc 	= 'quasi analytical scaling law';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
        configs.pedestal_sc  	= 'difference between H mode and L mode';
    case 11
        configs.lmode_sc 	= 'Sauter & Martin';
        configs.hmode_sc 	= 'Elbeze EIV 2005';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
        configs.pedestal_sc  	= 'difference between H mode and L mode';
    case {12,13}
        configs.lmode_sc 	= 'Sauter & Martin';
        configs.hmode_sc 	= 'Elbeze EIV 2005';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal with limitation on beta_N';
        configs.pedestal_sc  	= 'difference between H mode and L mode';
    case 14
        configs.lmode_sc 	= 'Robust scaling L1 (A.Murari, NF 57,2017,120617)';
        configs.hmode_sc 	= 'Robust scaling H1  (A.Murari, NF 57,2017,120617)';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
        configs.pedestal_sc = 'difference between H mode and L mode';
    case 15
        configs.lmode_sc 	= 'Robust scaling L1 (A.Murari, NF 57,2017,120617)';
        configs.hmode_sc 	= 'Robust scaling H1  (A.Murari, NF 57,2017,120617)';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal with limitation on beta_N';
        configs.pedestal_sc = 'difference between H mode and L mode';
    case 16
        configs.lmode_sc 	= 'ITER89-P';
        configs.hmode_sc 	= 'ITERH-98P(y,2)';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal with limitation on beta_N';
        configs.pedestal_sc = 'difference between H mode and L mode';
    case 17
        configs.lmode_sc 	= 'ITER89-P';
        configs.hmode_sc 	= 'ITER89-P + Cordey tau_{pedestal}';
        configs.core_sc  	= 'ITER89-P';
        configs.pedestal_sc = 'Cordey tau_{pedestal}';
    case 18
        configs.lmode_sc 	= 'ITER89-P';
        configs.hmode_sc 	= 'Petty 2008 gyroBohm scaling in H-mode)';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal with limitation on beta_N';
        configs.pedestal_sc = 'difference between H mode and L mode';
    case 19
        configs.lmode_sc 	= 'ITER89-P';
        configs.hmode_sc 	= 'ITER89-P + W_ped scaling incorporating experimental data and prediction from MHD code for ITER and DEMO';
        configs.core_sc  	= 'ITER89-P';
        configs.pedestal_sc = ' W_ped scaling incorporating experimental data and prediction from MHD code for ITER and DEMO';
    case {20,21}
        configs.lmode_sc 	= 'Elbeze 32nd EPS 2005 P1.033, Tarragona (EIV) Bohm';
        configs.hmode_sc 	= 'ITPA20 scaling for ITER 2020 with triangularity dependance';
        configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal with limitation on beta_N';
        configs.pedestal_sc = 'difference between H mode and L mode';
        
end

if z0dstruct.z0dinput.option.tauhemul  == 0
    configs.helium_sc 	= 'ITER,physics basis, Nuclear Fusion 39, 1999, p 2226';
else
    configs.helium_sc 	= sprintf('%g * tau_E',z0dstruct.z0dinput.option.tauhemul);
end

switch z0dstruct.z0dinput.option.zeff
    case 0
        configs.impurity_sc = 'None, line averaged Zeff given';
    case 1
        configs.impurity_sc = 'None, line averaged Zeff given + impurities neoclassical accumulation';
    case 2
        configs.impurity_sc = 'undefined scaling ; must be not used (reserved)';
    case 3
        configs.impurity_sc = 'undefined scaling ; must be not used (reserved)';
    case 4
        configs.impurity_sc = 'undefined scaling ; must be not used (reserved)';
    case 5
        configs.impurity_sc = 'Tore Supra scaling law';
    case 6
        configs.impurity_sc = 'Matthews scaling law';
    otherwise
        configs.impurity_sc = 'other scaling law';
end

switch z0dstruct.z0dinput.option.l2hscaling
    case 0
        configs.l2h_sc		= 'LH99(1)';
    case 1
        configs.l2h_sc		= 'LH2002 without Zeff effect';
    case 2
        configs.l2h_sc		= 'LH2002 with Zeff effect';
    case 3
        configs.l2h_sc		= 'YR Martin 2008';
    case 4
        configs.l2h_sc		= 'NLM-7 Murari 2012';
    case 5
        configs.l2h_sc		= 'NLM-11 Murari 2012';
    case 6
        configs.l2h_sc		= 'Jpol change of signe in edge region (E. R. Solano rule)';
    case 10
        configs.l2h_sc		= 'multimachine scaling law from Murami 2013 (BUEMS)';
    case 28
        configs.l2h_sc		= 'Low density case - Ryter et al, Nucl. Fusion 54 (2014) 083003 (9pp), equation 4';
    otherwise
        configs.l2h_sc		= 'criterion base on plasma rotation ( Gamma_ExB / Gamma_ITG)';
end
if z0dstruct.z0dinput.option.taurotmul == 0
  configs.tor_rot_sc	= 'ion heat confinement time';
else
  configs.tor_rot_sc	= sprintf(' somptaneous rotation given by Rice NF  47, (2007) p 1618-1624 + NBI torque with %g * tau_E confinement time',z0dstruct.z0dinput.option.taurotmul);
end
%
configs.coordinate	= 'sqrt(PSI_TOR/(pi * B0))';
% matieres 
%  configs.wall_mat 	= '';
%  configs.evap_mat 	= '';
%  configs.lim_mat 	= '';
%  configs.div_mat 	= '';

%% create an empty IDS
% already created
%summary = ids_gen('summary');


%% DEFINE EMPTY CODE AND IDS_PROPERTIES STRUCTURES
if ~isfield(summary,'code')
  summary.code = code_empty_imas;
end
if ~isfield(summary,'ids_properties')
  summary.ids_properties = ids_properties_empty_imas;
end

%% DATE (MSUAL: ADD DATE IN IDSS)
%if ~defined_imas(summary.ids_properties.pulse_type)
%  summary.ids_properties.pulse_type = sprintf('%f',clock2julday);
%end

%% COMMENT
if ~defined_imas(summary.ids_properties.comment)
  summary.ids_properties.comment = 'IMAS implementation of METIS';
end

%% HOMOGENEOUS_TIME
if ~defined_imas(summary.ids_properties.homogeneous_time)
  summary.ids_properties.homogeneous_time = 1;
end

if isfield(summary,'data_entry')

  %% USER
  if ~defined_imas(summary.data_entry.user)
    summary.data_entry.user = summary.data_entry.machine;
  end

  %% MACHINE
  if ~defined_imas(summary.data_entry.machine)
    summary.data_entry.machine =  z0dstruct.z0dinput.machine;
  end

  %% SHOT NUMBER
  if ~defined_imas(summary.data_entry.pulse)
    summary.data_entry.pulse = real(z0dstruct.z0dinput.shot(1));
  end

  %% RUN NUMBER
  if ~defined_imas(summary.data_entry.run)
    summary.data_entry.run = real(run);
  end
end

%
% temps
%
summary.time = data_zerod.temps;

% creation des donnees pour codeparam
z0dstruct.z0dinput.option.configs=configs;
data = z0dstruct.z0dinput.option;
tpn = tempname;
xml_write(tpn,data);
% lecture du fichier
fid = fopen(tpn,'r');
if fid > 0
  parameters = char(fread(fid,Inf,'char')');
  fclose(fid);
else
  parameters = 'unable to read parameters xml file';
end
delete(tpn);
%
summary.code.name       = 'METIS4IMAS';
summary.code.version    = num2str(zinebversion);
summary.code.parameters = parameters;
%
data             = [];
data.diboot      = data_zerod.diboot;
data.dw          = data_zerod.dw;
data.dpfus       = data_zerod.dpfus;
data.dini        = data_zerod.dini;
data.difcurconv  = data_zerod.difcurconv;
data.stf         = data_zerod.stf;
data.nb          = data_zerod.nb;
data.texte_diary = texte_diary;
tpn = tempname;
xml_write(tpn,data);
% lecture du fichier
fid = fopen(tpn,'r');
if fid > 0
  parameters = char(fread(fid,Inf,'char')');
  fclose(fid);
else
  parameters = 'unable to read parameters xml file';
end
delete(tpn);
summary.code.output_diag  = parameters;
%
summary.code.output_flag = error_flag * ones(size(data_zerod.temps));


%% global_quantities
summary.global_quantities.ip.value = data_zerod.ip;
summary.global_quantities.current_non_inductive.value = data_zerod.ini;
summary.global_quantities.current_bootstrap.value = data_zerod.iboot; % + data_zerod.ifus;
summary.global_quantities.current_ohm.value = data_zerod.iohm;
summary.global_quantities.current_alignment.value = data_zerod.ialign;
summary.global_quantities.v_loop.value = data_zerod.vloop;
summary.global_quantities.li.value = data_zerod.li;
if length(data_zerod.temps) == 1
  summary.global_quantities.beta_tor.value = (2/3) .* (data_zerod.w ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0(end) .^ 2 ./2 ./ phys.mu0); 
else
  summary.global_quantities.beta_tor.value = (2/3) .* (data_zerod.w ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0 .^ 2 ./2 ./ phys.mu0); 
end
summary.global_quantities.beta_tor_norm.value = 100 .* data_zerod.betan;
summary.global_quantities.beta_pol.value = data_zerod.betaptot;
summary.global_quantities.energy_diamagnetic.value = data_zerod.wdia;
summary.global_quantities.energy_total.value = data_zerod.w;
summary.global_quantities.energy_thermal.value = data_zerod.wth;
summary.global_quantities.energy_b_field_pol.value = data_zerod.wbp;
summary.global_quantities.volume.value = data_zerod.vp;
summary.global_quantities.h_mode.value = data_zerod.modeh;
%  if length(profil0d.temps) == 1
%      r0 = profil0d.Raxe(end,end);
%  else
%      r0              = interp1_imas(profil0d.temps,profil0d.Raxe(:,end),data_zerod.temps,'pchip','extrap');
%      if length(data_zerod.temps) > 1
%  	  r0              = trapz(data_zerod.temps,data_zerod.ip .* r0) ./ trapz(data_zerod.temps,data_zerod.ip);
%      end
%  end
%  summary.global_quantities.r0.value = r0;
%  if length(data_zerod.temps) == 1
%      summary.global_quantities.b0.value = z0dstruct.z0dinput.geo.b0(end) .* z0dstruct.z0dinput.geo.R(end) ./ r0;
%  else
%      summary.global_quantities.b0.value = z0dstruct.z0dinput.geo.b0 .* z0dstruct.z0dinput.geo.R ./ r0;
%  end
rb0 =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,data_zerod.temps,'pchip','extrap');
r0  =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,data_zerod.temps,'pchip','extrap');
summary.global_quantities.r0.value = mean(z0dstruct.z0dinput.geo.R);
summary.global_quantities.b0.value = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ summary.global_quantities.r0.value;



[h98,h96] = compute_H98(z0dstruct.zerod,z0dstruct.z0dinput.geo,z0dstruct.z0dinput.cons);
if length(data_zerod.temps) == 1
    summary.global_quantities.h_98.value = h98(end);
else
    summary.global_quantities.h_98.value = h98;
end
summary.global_quantities.tau_energy.value = data_zerod.taue;
summary.global_quantities.tau_helium.value = data_zerod.tauhe;
summary.global_quantities.tau_resistive.value = data_zerod.ip;
summary.global_quantities.resistance.value = data_zerod.RR;


%% LINE-AVERAGED QUANTITIES
summary.line_average.n_e.value 		= data_zerod.nbar;
summary.line_average.zeff.value 	= data_zerod.zmszl ./ max(eps,data_zerod.zeff);
summary.line_average.t_e.value          = interp1_imas(profil0d.temps,trapz(profil0d.xli,profil0d.tep,2),data_zerod.temps,'pchip','extrap');
summary.line_average.t_i_average.value  = interp1_imas(profil0d.temps,trapz(profil0d.xli,profil0d.tip,2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i_total.value    = interp1_imas(profil0d.temps,trapz(profil0d.xli,profil0d.nip,2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i.hydrogen.value = interp1_imas(profil0d.temps,trapz(profil0d.xli,squeeze(nions(:,:,1)),2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i.deuterium.value = interp1_imas(profil0d.temps,trapz(profil0d.xli,squeeze(nions(:,:,2)),2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i.tritium.value   = interp1_imas(profil0d.temps,trapz(profil0d.xli,squeeze(nions(:,:,3)),2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i.helium_3.value = interp1_imas(profil0d.temps,trapz(profil0d.xli,squeeze(nions(:,:,4)),2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i.helium_4.value = interp1_imas(profil0d.temps,trapz(profil0d.xli,squeeze(nions(:,:,5)),2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i.tungsten.value = interp1_imas(profil0d.temps,trapz(profil0d.xli,squeeze(nions(:,:,8)),2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i.tin.value = interp1_imas(profil0d.temps,trapz(profil0d.xli,squeeze(nions(:,:,9)),2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i.boron.value = interp1_imas(profil0d.temps,trapz(profil0d.xli,squeeze(nions(:,:,10)),2),data_zerod.temps,'pchip','extrap');
summary.line_average.n_i.berylium.value = zeros(size(data_zerod.temps));
summary.line_average.n_i.lithium.value  = zeros(size(data_zerod.temps));
summary.line_average.n_i.carbon.value   = zeros(size(data_zerod.temps));
summary.line_average.n_i.nitrogen.value = zeros(size(data_zerod.temps));
summary.line_average.n_i.neon.value     = zeros(size(data_zerod.temps));
summary.line_average.n_i.argon.value    = zeros(size(data_zerod.temps));
summary.line_average.n_i.xenon.value    = zeros(size(data_zerod.temps));
summary.line_average.n_i.oxygen.value   = zeros(size(data_zerod.temps));
clear as labels
for ispec=6:7
    [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
    nions_ispec = interp1_imas(profil0d.temps,trapz(profil0d.xli,squeeze(nions(:,:,ispec)),2),data_zerod.temps,'pchip','extrap');
    switch labels{ispec}
        case 'Be'
            summary.line_average.n_i.berylium.value = nions_ispec;
        case 'Li'
            summary.line_average.n_i.lithium.value  = nions_ispec;
        case 'C'
            summary.line_average.n_i.carbon.value   = nions_ispec;
        case 'N'
            summary.line_average.n_i.nitrogen.value = nions_ispec;
        case 'Ne'
            summary.line_average.n_i.neon.value     = nions_ispec;
        case 'Ar'
            summary.line_average.n_i.argon.value    = nions_ispec;
        case 'Xe'
            summary.line_average.n_i.xenon.value    = nions_ispec;
        case 'O'
            summary.line_average.n_i.oxygen.value   = nions_ispec;
        otherwise
            warning('plasma composition not yet implemented in summary IDS');
    end
end
%% VOLUME AVERAGED QUANTITIES
summary.volume_average.t_e.value 		= data_zerod.tem;
summary.volume_average.t_i.value 		= data_zerod.tem .* data_zerod.tite;
summary.volume_average.n_e.value 		= data_zerod.nem;
summary.volume_average.n_i_total.value 		= data_zerod.nim;
summary.volume_average.zeff.value 	        = data_zerod.zeff; 
summary.volume_average.n_i.hydrogen.value       = interp1_imas(profil0d.temps,nions_vol(:,1),data_zerod.temps,'pchip','extrap');
summary.volume_average.n_i.deuterium.value      = interp1_imas(profil0d.temps,nions_vol(:,2),data_zerod.temps,'pchip','extrap');
summary.volume_average.n_i.tritium.value        = interp1_imas(profil0d.temps,nions_vol(:,3),data_zerod.temps,'pchip','extrap');
summary.volume_average.n_i.helium_3.value       = interp1_imas(profil0d.temps,nions_vol(:,4),data_zerod.temps,'pchip','extrap');
summary.volume_average.n_i.helium_4.value       = interp1_imas(profil0d.temps,nions_vol(:,5),data_zerod.temps,'pchip','extrap');
summary.volume_average.n_i.tungsten.value       = interp1_imas(profil0d.temps,nions_vol(:,8),data_zerod.temps,'pchip','extrap');
summary.volume_average.n_i.tin.value       = interp1_imas(profil0d.temps,nions_vol(:,9),data_zerod.temps,'pchip','extrap');
summary.volume_average.n_i.boron.value       = interp1_imas(profil0d.temps,nions_vol(:,10),data_zerod.temps,'pchip','extrap');
summary.volume_average.n_i.berylium.value = zeros(size(data_zerod.temps));
summary.volume_average.n_i.lithium.value  = zeros(size(data_zerod.temps));
summary.volume_average.n_i.carbon.value   = zeros(size(data_zerod.temps));
summary.volume_average.n_i.nitrogen.value = zeros(size(data_zerod.temps));
summary.volume_average.n_i.neon.value     = zeros(size(data_zerod.temps));
summary.volume_average.n_i.argon.value    = zeros(size(data_zerod.temps));
summary.volume_average.n_i.xenon.value    = zeros(size(data_zerod.temps));
summary.volume_average.n_i.oxygen.value   = zeros(size(data_zerod.temps));
clear as labels
for ispec=6:7
  [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
  nions_ispec = interp1_imas(profil0d.temps,nions_vol(:,ispec),data_zerod.temps,'pchip','extrap');
  switch labels{ispec}
   case 'Be'
    summary.volume_average.n_i.berylium.value = nions_ispec;
   case 'Li'
    summary.volume_average.n_i.lithium.value  = nions_ispec;
   case 'C'
    summary.volume_average.n_i.carbon.value   = nions_ispec; 
   case 'N'
    summary.volume_average.n_i.nitrogen.value = nions_ispec;
   case 'Ne'
    summary.volume_average.n_i.neon.value     = nions_ispec; 
   case 'Ar'
    summary.volume_average.n_i.argon.value    = nions_ispec; 
   case 'Xe'
    summary.volume_average.n_i.xenon.value    = nions_ispec; 
   case 'O'
    summary.volume_average.n_i.oxygen.value   = nions_ispec;
   otherwise
    warning('plasma composition not yet implemented in summary IDS');
  end
end

% fusion data
summary.fusion.power.value                             = data_zerod.pfus;
summary.fusion.current.value                           = data_zerod.ifus;
summary.fusion.neutron_fluxes.total.value              = interp1_imas(neutron.time,neutron.ntot,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.thermal.value            = interp1_imas(neutron.time,neutron.ndd_th + neutron.ndt_th,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.dd.total.value           = interp1_imas(neutron.time,neutron.ndd,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.dd.thermal.value         = interp1_imas(neutron.time,neutron.ndd_th,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.dd.beam_thermal.value    = interp1_imas(neutron.time,neutron.ndd_nbi_th,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.dd.beam_beam.value       = interp1_imas(neutron.time,neutron.ndd_nbi_nbi,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.dt.total.value           = interp1_imas(neutron.time,neutron.ndt,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.dt.thermal.value         = interp1_imas(neutron.time,neutron.ndt_th,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.dt.beam_thermal.value    = interp1_imas(neutron.time,neutron.ndt_nbi_th,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.dt.beam_beam.value       = interp1_imas(neutron.time,neutron.ndt_nbi_nbi,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.tt.total.value           = interp1_imas(neutron.time,neutron.ntt,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.tt.thermal.value         = interp1_imas(neutron.time,neutron.ntt_th,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.tt.beam_thermal.value    = interp1_imas(neutron.time,neutron.ntt_nbi_th,data_zerod.temps,'pchip','extrap');
summary.fusion.neutron_fluxes.tt.beam_beam.value       = interp1_imas(neutron.time,neutron.ntt_nbi_nbi,data_zerod.temps,'pchip','extrap');
% new definition 
summary.fusion.neutron_rates = summary.fusion.neutron_fluxes;
% power from neutrons
if isfield(z0dstruct,'profil0d')
    summary.fusion.neutron_power_total.value = max(0,interp1_imas(z0dstruct.profil0d.temps,fusion.pfusion_total,data_zerod.temps,'pchip','extrap') - summary.fusion.power.value);
else
    summary.fusion.neutron_power_total.value = max(0,interp1_imas(z0dstruct.profil.temps,fusion.pfusion_total,data_zerod.temps,'pchip','extrap') - summary.fusion.power.value);
end
% runaways
summary.runaways.particles.value                = data_zerod.irun ./ phys.e ./ phys.c;
summary.runaways.current.value                  = data_zerod.irun;

%SOL
summary.scrape_off_layer.t_e_decay_length.value                 = data_zerod.dsol;
summary.scrape_off_layer.t_i_average_decay_length.value         = data_zerod.dsol;
summary.scrape_off_layer.n_e_decay_length.value                 = data_zerod.dsol;
summary.scrape_off_layer.n_i_total_decay_length.value           = data_zerod.dsol;
summary.scrape_off_layer.heat_flux_e_decay_length.value         = data_zerod.dsol;
summary.scrape_off_layer.heat_flux_i_decay_length.value         = data_zerod.dsol;
summary.scrape_off_layer.power_radiated.value                   = data_zerod.pradsol;
pn0 = z0dstruct.z0dinput.option.temp_vac .* phys.k .* profil0d.n0m(:,end) + phys.e .* profil0d.tep(:,end) .* profil0d.n0(:,end);
summary.scrape_off_layer.pressure_neutral.value                 = interp1_imas(profil0d.temps,pn0,data_zerod.temps,'pchip','extrap');

%% CENTRAL VALUES
summary.local.magnetic_axis.position.rho_tor_norm = zeros(length(data_zerod.te0),1);
summary.local.magnetic_axis.position.rho_tor      = zeros(length(data_zerod.te0),1);     
summary.local.magnetic_axis.position.psi          = interp1_imas(profil0d.temps,profil0d.psi(:,1),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.position.r            = interp1_imas(profil0d.temps,profil0d.Raxe(:,1),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.position.z            = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.z0,data_zerod.temps,'pchip','extrap');              
summary.local.magnetic_axis.zeff.value            = interp1_imas(profil0d.temps,profil0d.zeff(:,1),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.momentum_tor.value    = interp1_imas(profil0d.temps,profil0d.rtor(:,1),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.q.value               = data_zerod.q0;
summary.local.magnetic_axis.magnetic_shear.value  = interp1_imas(profil0d.temps,magnetic_shear(:,1),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.t_e.value             = data_zerod.te0;
summary.local.magnetic_axis.t_i_average.value     = interp1_imas(profil0d.temps,profil0d.tip(:,1),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.n_e.value             = data_zerod.ne0;
summary.local.magnetic_axis.n_i_total.value       = data_zerod.ni0;
summary.local.magnetic_axis.position.psi          = interp1_imas(profil0d.temps,profil0d.psi(:,1),data_zerod.temps,'pchip','extrap');
%% FILL VTOR FOR KNOWN SPECIES
summary.local.magnetic_axis.velocity_tor.hydrogen.value = interp1_imas(profil0d.temps,vtor(:,1,1),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.velocity_tor.deuterium.value = interp1_imas(profil0d.temps,vtor(:,1,2),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.velocity_tor.tritium.value = interp1_imas(profil0d.temps,vtor(:,1,3),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.velocity_tor.helium_3.value = interp1_imas(profil0d.temps,vtor(:,1,4),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.velocity_tor.helium_4.value = interp1_imas(profil0d.temps,vtor(:,1,5),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.velocity_tor.tungsten.value = interp1_imas(profil0d.temps,vtor(:,1,8),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.velocity_tor.tin.value = interp1_imas(profil0d.temps,vtor(:,1,9),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.velocity_tor.boron.value = interp1_imas(profil0d.temps,vtor(:,1,10),data_zerod.temps,'pchip','extrap'); 
%% FILL VTOR TESTING ION SPECIES FOR POSITIONS 6 AND 7
summary.local.magnetic_axis.velocity_tor.berylium.value = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.velocity_tor.lithium.value  = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.velocity_tor.carbon.value   = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.velocity_tor.nitrogen.value = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.velocity_tor.neon.value     = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.velocity_tor.argon.value    = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.velocity_tor.xenon.value    = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.velocity_tor.oxygen.value   = zeros(size(data_zerod.temps));
clear as labels
for ispec=6:7
  [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
  switch labels{ispec}
   case 'Be'
    summary.local.magnetic_axis.velocity_tor.berylium.value = interp1_imas(profil0d.temps,vtor(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Li'
    summary.local.magnetic_axis.velocity_tor.lithium.value  = interp1_imas(profil0d.temps,vtor(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'C'
    summary.local.magnetic_axis.velocity_tor.carbon.value   = interp1_imas(profil0d.temps,vtor(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'N'
    summary.local.magnetic_axis.velocity_tor.nitrogen.value = interp1_imas(profil0d.temps,vtor(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ne'
    summary.local.magnetic_axis.velocity_tor.neon.value     = interp1_imas(profil0d.temps,vtor(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ar'
    summary.local.magnetic_axis.velocity_tor.argon.value    = interp1_imas(profil0d.temps,vtor(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Xe'
    summary.local.magnetic_axis.velocity_tor.xenon.value    = interp1_imas(profil0d.temps,vtor(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'O'
    summary.local.magnetic_axis.velocity_tor.oxygen.value   = interp1_imas(profil0d.temps,vtor(:,1,ispec),data_zerod.temps,'pchip','extrap');
   otherwise
    warning('plasma composition not yet implemented in summary IDS');
  end
end
summary.local.magnetic_axis.n_i.hydrogen.value   = interp1_imas(profil0d.temps,nions(:,1,1),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.n_i.deuterium.value  = interp1_imas(profil0d.temps,nions(:,1,2),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.n_i.tritium.value    = interp1_imas(profil0d.temps,nions(:,1,3),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.n_i.helium_3.value   = interp1_imas(profil0d.temps,nions(:,1,4),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.n_i.helium_4.value   = interp1_imas(profil0d.temps,nions(:,1,5),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.n_i.berylium.value   = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.n_i.lithium.value    = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.n_i.carbon.value     = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.n_i.nitrogen.value   = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.n_i.neon.value       = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.n_i.argon.value      = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.n_i.xenon.value      = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.n_i.oxygen.value      = zeros(size(data_zerod.temps));
summary.local.magnetic_axis.n_i.tungsten.value   = interp1_imas(profil0d.temps,nions(:,1,8),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.n_i.tin.value   = interp1_imas(profil0d.temps,nions(:,1,9),data_zerod.temps,'pchip','extrap');
summary.local.magnetic_axis.n_i.boron.value   = interp1_imas(profil0d.temps,nions(:,1,10),data_zerod.temps,'pchip','extrap');
clear as labels
for ispec=6:7
  [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
  switch labels{ispec}
   case 'Be'
    summary.local.magnetic_axis.n_i.berylium.value = interp1_imas(profil0d.temps,nions(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Li'
    summary.local.magnetic_axis.n_i.lithium.value  = interp1_imas(profil0d.temps,nions(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'C'
    summary.local.magnetic_axis.n_i.carbon.value   = interp1_imas(profil0d.temps,nions(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'N'
    summary.local.magnetic_axis.n_i.nitrogen.value = interp1_imas(profil0d.temps,nions(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ne'
    summary.local.magnetic_axis.n_i.neon.value     = interp1_imas(profil0d.temps,nions(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ar'
    summary.local.magnetic_axis.n_i.argon.value    = interp1_imas(profil0d.temps,nions(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Xe'
    summary.local.magnetic_axis.n_i.xenon.value    = interp1_imas(profil0d.temps,nions(:,1,ispec),data_zerod.temps,'pchip','extrap');
   case 'O'
    summary.local.magnetic_axis.n_i.oxygen.value   = interp1_imas(profil0d.temps,nions(:,1,ispec),data_zerod.temps,'pchip','extrap');
   otherwise
    warning('plasma composition not yet implemented in summary IDS');
  end
end

%% EDGE QUANTITIES
summary.local.edge.t_e.value             = data_zerod.tebord; 
summary.local.edge.t_i_average.value     = data_zerod.tibord;  
summary.local.edge.n_e.value             = data_zerod.nebord; 
summary.local.edge.n_i_total.value       = data_zerod.nibord;
summary.local.edge.position.psi          = interp1_imas(profil0d.temps,profil0d.psi(:,end),data_zerod.temps,'pchip','extrap') ;
summary.local.edge.position.rho_tor      = data_zerod.rm; 
summary.local.edge.position.rho_tor_norm = ones(size(data_zerod.rm));
summary.local.edge.q.value               = data_zerod.qa;
summary.local.edge.zeff.value            = interp1_imas(profil0d.temps,profil0d.zeff(:,1),data_zerod.temps,'pchip','extrap');
summary.local.edge.momentum_tor.value    = interp1_imas(profil0d.temps,profil0d.rtor(:,end),data_zerod.temps,'pchip','extrap');
summary.local.edge.magnetic_shear.value  = interp1_imas(profil0d.temps,magnetic_shear(:,end),data_zerod.temps,'pchip','extrap');
%% FILL VTOR FOR KNOWN SPECIES
summary.local.edge.velocity_tor.hydrogen.value = interp1_imas(profil0d.temps,vtor(:,end,1),data_zerod.temps,'pchip','extrap');
summary.local.edge.velocity_tor.deuterium.value = interp1_imas(profil0d.temps,vtor(:,end,2),data_zerod.temps,'pchip','extrap');
summary.local.edge.velocity_tor.tritium.value = interp1_imas(profil0d.temps,vtor(:,end,3),data_zerod.temps,'pchip','extrap');
% temporary, until spelling fixed in the summary IDS
summary.local.edge.velocity_tor.titrium.value = interp1_imas(profil0d.temps,vtor(:,end,3),data_zerod.temps,'pchip','extrap');
summary.local.edge.velocity_tor.helium_3.value = interp1_imas(profil0d.temps,vtor(:,end,4),data_zerod.temps,'pchip','extrap');
summary.local.edge.velocity_tor.helium_4.value = interp1_imas(profil0d.temps,vtor(:,end,5),data_zerod.temps,'pchip','extrap');
summary.local.edge.velocity_tor.tungsten.value = interp1_imas(profil0d.temps,vtor(:,end,8),data_zerod.temps,'pchip','extrap');
summary.local.edge.velocity_tor.tin.value = interp1_imas(profil0d.temps,vtor(:,end,9),data_zerod.temps,'pchip','extrap');
summary.local.edge.velocity_tor.boron.value = interp1_imas(profil0d.temps,vtor(:,end,10),data_zerod.temps,'pchip','extrap');

%% FILL VTOR TESTING ION SPECIES FOR POSITIONS 6 AND 7
summary.local.edge.velocity_tor.berylium.value = zeros(size(data_zerod.temps));
summary.local.edge.velocity_tor.lithium.value  = zeros(size(data_zerod.temps));
summary.local.edge.velocity_tor.carbon.value   = zeros(size(data_zerod.temps));
summary.local.edge.velocity_tor.nitrogen.value = zeros(size(data_zerod.temps));
summary.local.edge.velocity_tor.neon.value     = zeros(size(data_zerod.temps));
summary.local.edge.velocity_tor.argon.value    = zeros(size(data_zerod.temps));
summary.local.edge.velocity_tor.xenon.value    = zeros(size(data_zerod.temps));
summary.local.edge.velocity_tor.oxygen.value   = zeros(size(data_zerod.temps));
clear as labels
for ispec=6:7
  [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
  switch labels{ispec}
   case 'Be'
    summary.local.edge.velocity_tor.berylium.value = interp1_imas(profil0d.temps,vtor(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'Li'
    summary.local.edge.velocity_tor.lithium.value  = interp1_imas(profil0d.temps,vtor(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'C'
    summary.local.edge.velocity_tor.carbon.value   = interp1_imas(profil0d.temps,vtor(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'N'
    summary.local.edge.velocity_tor.nitrogen.value = interp1_imas(profil0d.temps,vtor(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ne'
    summary.local.edge.velocity_tor.neon.value     = interp1_imas(profil0d.temps,vtor(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ar'
    summary.local.edge.velocity_tor.argon.value    = interp1_imas(profil0d.temps,vtor(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'Xe'
    summary.local.edge.velocity_tor.xenon.value    = interp1_imas(profil0d.temps,vtor(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'O'
    summary.local.edge.velocity_tor.oxygen.value   = interp1_imas(profil0d.temps,vtor(:,end,ispec),data_zerod.temps,'pchip','extrap');
   otherwise
    warning('plasma composition not yet implemented in summary IDS');
  end
end

summary.local.edge.n_i.hydrogen.value   = interp1_imas(profil0d.temps,nions(:,end,1),data_zerod.temps,'pchip','extrap');
summary.local.edge.n_i.deuterium.value  = interp1_imas(profil0d.temps,nions(:,end,2),data_zerod.temps,'pchip','extrap');
summary.local.edge.n_i.tritium.value    = interp1_imas(profil0d.temps,nions(:,end,3),data_zerod.temps,'pchip','extrap');
summary.local.edge.n_i.helium_3.value   = interp1_imas(profil0d.temps,nions(:,end,4),data_zerod.temps,'pchip','extrap');
summary.local.edge.n_i.helium_4.value   = interp1_imas(profil0d.temps,nions(:,end,5),data_zerod.temps,'pchip','extrap');
summary.local.edge.n_i.berylium.value   = zeros(size(data_zerod.temps));
summary.local.edge.n_i.lithium.value    = zeros(size(data_zerod.temps));
summary.local.edge.n_i.carbon.value     = zeros(size(data_zerod.temps));
summary.local.edge.n_i.nitrogen.value   = zeros(size(data_zerod.temps));
summary.local.edge.n_i.neon.value       = zeros(size(data_zerod.temps));
summary.local.edge.n_i.argon.value      = zeros(size(data_zerod.temps));
summary.local.edge.n_i.xenon.value      = zeros(size(data_zerod.temps));
summary.local.edge.n_i.oxygen.value      = zeros(size(data_zerod.temps));
summary.local.edge.n_i.tungsten.value   = interp1_imas(profil0d.temps,nions(:,1,8),data_zerod.temps,'pchip','extrap');
summary.local.edge.n_i.tin.value   = interp1_imas(profil0d.temps,nions(:,1,9),data_zerod.temps,'pchip','extrap');
summary.local.edge.n_i.boron.value   = interp1_imas(profil0d.temps,nions(:,1,10),data_zerod.temps,'pchip','extrap');

clear as labels
for ispec=6:7
  [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
  switch labels{ispec}
   case 'Be'
    summary.local.edge.n_i.berylium.value = interp1_imas(profil0d.temps,nions(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'Li'
    summary.local.edge.n_i.lithium.value  = interp1_imas(profil0d.temps,nions(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'C'
    summary.local.edge.n_i.carbon.value   = interp1_imas(profil0d.temps,nions(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'N'
    summary.local.edge.n_i.nitrogen.value = interp1_imas(profil0d.temps,nions(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ne'
    summary.local.edge.n_i.neon.value     = interp1_imas(profil0d.temps,nions(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ar'
    summary.local.edge.n_i.argon.value    = interp1_imas(profil0d.temps,nions(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'Xe'
    summary.local.edge.n_i.xenon.value    = interp1_imas(profil0d.temps,nions(:,end,ispec),data_zerod.temps,'pchip','extrap');
   case 'O'
    summary.local.edge.n_i.oxygen.value   = interp1_imas(profil0d.temps,nions(:,end,ispec),data_zerod.temps,'pchip','extrap');
   otherwise
    warning('plasma composition not yet implemented in summary IDS');
  end
end

%summary.local.edge.position.phi       = interp1_imas(profil0d.temps,profil0d.phi(:,end),data_zerod.temps,'pchip','extrap'); 
%summary.local.edge.drho_edge_dt.value = data_zerod.drmdt;
%summary.local.edge.neutral_flux.value = data_zerod.n0a;
%summary.local.edge.phi_plasma.value   = data_zerod.phiplasma;
% Miss n_i, zeff, momentum_tor, magnetic_shear
%

% Not adapted: wait for variables to exist in the summary IDS
%  summary.global_quantities.energy.w_tot.value 		= data_zerod.w;
%  summary.global_quantities.energy.w_b_pol.value 		= data_zerod.wbp;
%  if length(data_zerod.temps) == 1
%    wdia  = z0dstruct.zerod.wdia;
%    t     = z0dstruct.zerod.temps;
%    summary.global_quantities.energy.dwdia_dt.value 	= (wdia(end) - wdia(end-1)) ./ (t(end) - t(end -1));
%  else
%    summary.global_quantities.energy.dwdia_dt.value 	= z0dxdt(data_zerod.wdia,data_zerod.temps);
%  end
%  summary.global_quantities.energy.w_b_tor_pla.value 	= data_zerod.wmagtor;
%  summary.global_quantities.energy.w_th.value 		= data_zerod.wth;
%  summary.global_quantities.energy.dwtot_dt.value 	= data_zerod.dwdt;
%  summary.global_quantities.energy.dwbpol_dt.value 	= data_zerod.dwbpdt;
%  summary.global_quantities.energy.dwbtorpla_dt.value 	= data_zerod.dwmagtordt;
%  summary.global_quantities.energy.dwth_dt.value 		= data_zerod.dwthdt;
%  summary.global_quantities.energy.esup_icrhtot.value 	= data_zerod.esup_icrh;
%  summary.global_quantities.energy.esup_icrhper.value 	= ((1/3)+(1/3).*max(0,tanh(data_zerod.einj_icrh./data_zerod.ecrit_icrh))).*data_zerod.esup_icrh;
%  summary.global_quantities.energy.esup_nbitot.value 	= data_zerod.esup_nbi;
%  summary.global_quantities.energy.esup_nbiperp.value 	= ((1/3) + (2/3) .*  sin(z0dstruct.z0dinput.option.angle_nbi./180*pi) .*  ...
%  						          max(0,tanh(z0dstruct.z0dinput.option.einj ./ data_zerod.ecrit_nbi))) .* data_zerod.esup_nbi;
%  summary.global_quantities.energy.esup_lhcd.value 	= data_zerod.esup_lh;
%  summary.global_quantities.energy.esup_alpha.value 	= data_zerod.esup_fus;
%  
%  summary.global_quantities.energy_diamagnetic.value = data_zerod.wdia;

%% EQGOMETRY NOT IN SUMMARY IDS, KEPT HERE FOR NOT LOSING THE INFORMATION
%  summary.eqgeometry.source              = 'METIS';
%  %   for k=1:length(data_zerod.temps)
%  %   	if data_zerod.xpoint(k) == 1
%  %   		summary.eqgeometry.boundarytype{k}  = 'separatrix';
%  %   	else 
%  %   		summary.eqgeometry.boundarytype{k}  = 'limiter';
%  %   	end
%  %   end
%  %summary.eqgeometry.boundarytype = sprintf('%f',data_zerod.xpoint);
%  summary.eqgeometry.boundarytype = data_zerod.xpoint;
%  if isfield(profil0d,'Rsepa') && isfield(profil0d,'Zsepa')
%    summary.eqgeometry.boundary(1).r 		= interp1_imas(profil0d.temps,profil0d.Rsepa,data_zerod.temps,'nearest','extrap');
%    summary.eqgeometry.boundary(1).z 		= interp1_imas(profil0d.temps,profil0d.Zsepa,data_zerod.temps,'nearest','extrap');
%  else
%    summary.eqgeometry.boundary.r 		= [];
%    summary.eqgeometry.boundary.z 		= [];
%    
%  end
%  if length(data_zerod.temps) == 1
%    summary.eqgeometry.geom_axis.r 		= z0dstruct.z0dinput.geo.R(end);
%    summary.eqgeometry.geom_axis.z 		= z0dstruct.z0dinput.geo.z0(end);
%    summary.eqgeometry.a_minor 	 		= z0dstruct.z0dinput.geo.a(end);
%    summary.eqgeometry.elongation 	 		= z0dstruct.z0dinput.geo.K(end);
%    summary.eqgeometry.tria_upper 	 		= z0dstruct.z0dinput.geo.d(end);
%    summary.eqgeometry.tria_lower 	 		= z0dstruct.z0dinput.geo.d(end);
%  else
%    summary.eqgeometry.geom_axis.r 		= z0dstruct.z0dinput.geo.R;
%    summary.eqgeometry.geom_axis.z 		= z0dstruct.z0dinput.geo.z0;
%    summary.eqgeometry.a_minor 	 		= z0dstruct.z0dinput.geo.a;
%    summary.eqgeometry.elongation 	 		= z0dstruct.z0dinput.geo.K;
%    summary.eqgeometry.tria_upper 	 		= z0dstruct.z0dinput.geo.d;
%    summary.eqgeometry.tria_lower 	 		= z0dstruct.z0dinput.geo.d;
%  end
%  summary.eqgeometry.xpts.r 			= [];
%  summary.eqgeometry.xpts.z 			= [];
%  summary.eqgeometry.left_low_st.r 		= [];
%  summary.eqgeometry.left_low_st.z 		= [];
%  summary.eqgeometry.right_low_st.r 		= [];
%  summary.eqgeometry.right_low_st.z 		= [];
%  summary.eqgeometry.left_up_st.r 		= []; 
%  summary.eqgeometry.left_up_st.z 		= []; 
%  summary.eqgeometry.right_up_st.r 		= [];
%  summary.eqgeometry.right_up_st.z 		= [];
%  summary.eqgeometry.active_limit.r 		= [];  
%  summary.eqgeometry.active_limit.z 		= [];
%  %
%  
%  summary.global_quantities.ip.value = data_zerod.ip;
%  if length(data_zerod.temps) == 1
%    ip    = z0dstruct.zerod.ip;
%    t     = z0dstruct.zerod.temps;
%    %summary.global_param.dip_dt.value 	= (ip(end) - ip(end-1)) ./ (t(end) - t(end -1));
%  else
%    %summary.global_param.dip_dt.value 		= z0dxdt(data_zerod.ip,data_zerod.temps); 
%  end
%  summary.global_quantities.beta_pol.value = data_zerod.betaptot;
%  mu0 =  4 .* pi .* 1e-7;
%  if length(data_zerod.temps) == 1
%    summary.global_quantities.beta_tor.value = (2/3) .* (data_zerod.w ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0(end) .^ 2 ./2 ./ mu0); 
%  else
%    summary.global_quantities.beta_tor.value = (2/3) .* (data_zerod.w ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0 .^ 2 ./2 ./ mu0); 
%  end	
%  summary.global_quantities.beta_tor_norm.value = data_zerod.betan;
%  summary.global_quantities.li.value = data_zerod.li; 
%  summary.global_quantities.volume.value = data_zerod.vp;
%  %summary.global_param.area_pol.value 	= data_zerod.sp;
%  %summary.global_param.area_ext.value 	= data_zerod.sext; 
%  %summary.global_param.len_sepa.value 	= data_zerod.peri; 
%  %summary.global_param.beta_pol_th.value 	= data_zerod.betap;
%  %if length(data_zerod.temps) == 1
%  %  summary.global_param.beta_tor_th.value	=  (2/3) .* (data_zerod.wth ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0(end) .^ 2 ./2 ./ mu0); 
%  %else
%  %  summary.global_param.beta_tor_th.value	=  (2/3) .* (data_zerod.wth ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0 .^ 2 ./2 ./ mu0); 
%  %end
%  %summary.global_param.beta_n_th.value 	= summary.energy.w_th.value ./ max(eps,summary.energy.w_tot.value ) .* data_zerod.betan;
%  %summary.global_param.disruption.value 	= data_zerod.disrup;
%  %summary.global_param.mode_h.value 		= data_zerod.modeh;
%  %summary.global_param.s_alpha.value 		= data_zerod.salpha;

%% EC, FIRST LAUNCHER
summary.heating_current_drive.ec{1}.power.value = data_zerod.pecrh;
%summary.heating_current_drive.ec{1}.frequency.value
if length(profil0d.temps) == 1
  summary.heating_current_drive.ec{1}.position.value =interp1_imas(profil0d.xli,profil0d.rmx ./ profil0d.rmx(end,end),z0dstruct.z0dinput.cons.xece(end),'pchip','extrap');
else
  [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
  rho_ece = griddata_imas(tt,xx,profil0d.rmx ./ (profil0d.rmx(:,end) * ones(size(profil0d.xli))),z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.cons.xece,'cubic');
  summary.heating_current_drive.ec{1}.position.value =  interp1_imas(z0dstruct.z0dinput.cons.temps,rho_ece,data_zerod.temps,'pchip','extrap');
end
%summary.heating_current_drive.ec{1}.polarisation.value
%summary.heating_current_drive.ec{1}.harmonic.value
%summary.heating_current_drive.ec{1}.angle_tor.value
%ssummary.heating_current_drive.ec{1}.angle_pol.value
summary.heating_current_drive.ec{1}.current.value  = data_zerod.ieccd;
%summary.heating_current_drive.ec{1}.energy_fast.value
%% REMINDER: SECOND EC LAUNCHER STORED IN LH WHEN LHMODE=5
if z0dstruct.z0dinput.option.lhmode == 5
  summary.heating_current_drive.ec{2} = summary.heating_current_drive.ec{1};
  summary.heating_current_drive.ec{2}.power.value = data_zerod.plh;
  if length(profil0d.temps) == 1
    summary.heating_current_drive.ec{2}.position.value =interp1_imas(profil0d.xli,profil0d.rmx ./ profil0d.rmx(end,end),z0dstruct.z0dinput.option.xlh,'pchip','extrap');
  else
    [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
    rho_ece = griddata_imas(tt,xx,profil0d.rmx ./ (profil0d.rmx(:,end) * ones(size(profil0d.xli))),data_zerod.temps,z0dstruct.z0dinput.option.xlh * ones(size(data_zerod.temps)),'cubic');
    summary.heating_current_drive.ec{1}.position.value =  rho_ece;
  end
  summary.heating_current_drive.ec{1}.current.value  = data_zerod.ilh;
else
  summary.heating_current_drive.lh{1}.power.value       = data_zerod.plh;
  summary.heating_current_drive.lh{1}.frequency.value   = z0dstruct.z0dinput.option.freqlh * ones(size(data_zerod.temps)) .* 1e9;
  dd     = abs(profil0d.plh - max(profil0d.plh,[],2) * ones(size(profil0d.xli)));
  mask   = double(dd == 0);
  rho_lh = sum(profil0d.rmx .* mask,2) ./ max(1,sum(mask,2));
  if length(rho_lh) == 1
      summary.heating_current_drive.lh{1}.position.value  =  rho_lh; 
  else 
      summary.heating_current_drive.lh{1}.position.value  = interp1_imas(profil0d.temps,rho_lh,data_zerod.temps,'pchip','extrap');
  end
  summary.heating_current_drive.lh{1}.n_parallel.value  = z0dstruct.z0dinput.option.npar0 * ones(size(data_zerod.temps));
  summary.heating_current_drive.lh{1}.power.value       = data_zerod.plh;
  summary.heating_current_drive.lh{1}.current.value     = data_zerod.ilh;
  summary.heating_current_drive.lh{1}.energy_fast.value = data_zerod.esup_lh;
end
%% ICRH
summary.heating_current_drive.ic{1}.power.value = data_zerod.picrh;
summary.heating_current_drive.ic{1}.frequency.value = z0dstruct.z0dinput.option.freq * ones(size(data_zerod.temps)) .* 1e6;
if length(profil0d.temps) == 1
    summary.heating_current_drive.ic{1}.position.value = interp1_imas(profil0d.xli,profil0d.rmx ./ profil0d.rmx(end,end),data_zerod.xres,'pchip','extrap');
else
    [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
    rho_icrh = griddata_imas(tt,xx,profil0d.rmx ./ (profil0d.rmx(:,end) * ones(size(profil0d.xli))),data_zerod.temps,data_zerod.xres,'cubic');
    summary.heating_current_drive.ic{1}.position.value = rho_icrh;
end
summary.heating_current_drive.ic{1}.n_tor.value = z0dstruct.z0dinput.option.nphi * ones(size(data_zerod.temps));
%summary.heating_current_drive.ic{1}.k_perpendicular.value = 
%summary.heating_current_drive.ic{1}.e_field_plus_minus_ratio.value
summary.heating_current_drive.ic{1}.harmonic.value = data_zerod.harm;
%summary.heating_current_drive.ic{1}.phase.value
summary.heating_current_drive.ic{1}.current.value = data_zerod.ifwcd;
summary.heating_current_drive.ic{1}.energy_fast.value = data_zerod.esup_icrh;

%% NBI
if isfield(z0dstruct.z0dinput.option,'forced_H_NBI') && (z0dstruct.z0dinput.option.forced_H_NBI ~= 0)
    ftnbi = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.cons.ftnbi,data_zerod.temps,'pchip','extrap');
    summary.heating_current_drive.nbi{1}.power.value            = real(data_zerod.pnbi) .* 0;
    summary.heating_current_drive.nbi{1}.current.value          = real(data_zerod.inbicd) .* 0;
    summary.heating_current_drive.nbi{1}.tangency_radius.value  = z0dstruct.z0dinput.option.rtang;
    rext = abs(profil0d.Raxe(:,end) .* (1 + profil0d.epsi(:,end)) - z0dstruct.z0dinput.option.rtang * ones(size(profil0d.temps)));
    zext = z0dstruct.z0dinput.option.zext .* ones(size(profil0d.temps));
    if length(profil0d.temps) == 1
        kext = interp1_imas(profil0d.xli,profil0d.kx,zext,'pchip','extrap');
    else
        [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
        kext = griddata_imas(tt,xx,profil0d.kx,profil0d.temps,zext,'cubic');
    end
    zext = zext .* kext .*  profil0d.Raxe(:,end) .* profil0d.epsi(:,end);
    angle_ext = atan(zext ./ rext);
    summary.heating_current_drive.nbi{1}.angle.value            = mean(angle_ext);
    summary.heating_current_drive.nbi{1}.direction.value        = sign(z0dstruct.z0dinput.option.angle_nbi);
    summary.heating_current_drive.nbi{1}.energy.value           = z0dstruct.z0dinput.option.einj .* ones(size(summary.time));
    summary.heating_current_drive.nbi{1}.species.z_n.value      = 1;
    summary.heating_current_drive.nbi{1}.species.a.value        = 2;
    %
    summary.heating_current_drive.nbi{2} = summary.heating_current_drive.nbi{1};
    summary.heating_current_drive.nbi{2}.power.value            = real(data_zerod.pnbi);
    summary.heating_current_drive.nbi{2}.current.value          = real(data_zerod.inbicd);
    summary.heating_current_drive.nbi{2}.tangency_radius.value  = z0dstruct.z0dinput.option.rtang;
    rext = abs(profil0d.Raxe(:,end) .* (1 + profil0d.epsi(:,end)) - z0dstruct.z0dinput.option.rtang * ones(size(profil0d.temps)));
    zext = z0dstruct.z0dinput.option.zext .* ones(size(profil0d.temps));
    if length(profil0d.temps) == 1
        kext = interp1_imas(profil0d.xli,profil0d.kx,zext,'pchip','extrap');
    else
        [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
        kext = griddata_imas(tt,xx,profil0d.kx,profil0d.temps,zext,'cubic');
    end
    zext = zext .* kext .*  profil0d.Raxe(:,end) .* profil0d.epsi(:,end);
    angle_ext = atan(zext ./ rext);
    summary.heating_current_drive.nbi{2}.angle.value            = mean(angle_ext);
    summary.heating_current_drive.nbi{2}.direction.value        = sign(z0dstruct.z0dinput.option.angle_nbi);
    summary.heating_current_drive.nbi{2}.energy.value           = z0dstruct.z0dinput.option.einj .* ones(size(summary.time));
    % summary.heating_current_drive.nbi{2}.beam_current_fraction.value
    % summary.heating_current_drive.nbi{2}.beam_power_fraction.value
    summary.heating_current_drive.nbi{2}.species.z_n.value      = 1;
    summary.heating_current_drive.nbi{2}.species.a.value = 1;
    
    if z0dstruct.z0dinput.option.nb_nbi > 1
        summary.heating_current_drive.nbi{3} = summary.heating_current_drive.nbi{1};
        summary.heating_current_drive.nbi{3}.power.value            = imag(data_zerod.pnbi)   .* 0;
        summary.heating_current_drive.nbi{3}.current.value          = imag(data_zerod.inbicd) .* 0;
        summary.heating_current_drive.nbi{3}.tangency_radius.value  = z0dstruct.z0dinput.option.rtang2;
        rext = abs(profil0d.Raxe(:,end) .* (1 + profil0d.epsi(:,end)) - z0dstruct.z0dinput.option.rtang2 * ones(size(profil0d.temps)));
        zext = z0dstruct.z0dinput.option.zext2 .* ones(size(profil0d.temps));
        if length(profil0d.temps) == 1
            kext = interp1_imas(profil0d.xli,profil0d.kx,zext,'pchip','extrap');
        else
            [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
            kext = griddata_imas(tt,xx,profil0d.kx,profil0d.temps,zext,'cubic');
        end
        zext = zext .* kext .*  profil0d.Raxe(:,end) .* profil0d.epsi(:,end);
        angle_ext = atan(zext ./ rext);
        summary.heating_current_drive.nbi{3}.angle.value            = mean(angle_ext);
        summary.heating_current_drive.nbi{3}.direction.value        = sign(z0dstruct.z0dinput.option.angle_nbi2);
        summary.heating_current_drive.nbi{3}.energy.value           = z0dstruct.z0dinput.option.einj2 .* ones(size(summary.time));
        summary.heating_current_drive.nbi{3}.species.z_n.value      = 1;
        summary.heating_current_drive.nbi{3}.species.a.value        = 2;
        % summary.heating_current_drive.nbi{2}.beam_current_fraction.value
        % summary.heating_current_drive.nbi{2}.beam_power_fraction.value
        summary.heating_current_drive.nbi{4} = summary.heating_current_drive.nbi{1};
        summary.heating_current_drive.nbi{4}.power.value            = imag(data_zerod.pnbi) ;
        summary.heating_current_drive.nbi{4}.current.value          = imag(data_zerod.inbicd);
        summary.heating_current_drive.nbi{4}.tangency_radius.value  = z0dstruct.z0dinput.option.rtang2;
        rext = abs(profil0d.Raxe(:,end) .* (1 + profil0d.epsi(:,end)) - z0dstruct.z0dinput.option.rtang2 * ones(size(profil0d.temps)));
        zext = z0dstruct.z0dinput.option.zext2 .* ones(size(profil0d.temps));
        if length(profil0d.temps) == 1
            kext = interp1_imas(profil0d.xli,profil0d.kx,zext,'pchip','extrap');
        else
            [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
            kext = griddata_imas(tt,xx,profil0d.kx,profil0d.temps,zext,'cubic');
        end
        zext = zext .* kext .*  profil0d.Raxe(:,end) .* profil0d.epsi(:,end);
        angle_ext = atan(zext ./ rext);
        summary.heating_current_drive.nbi{4}.angle.value            = mean(angle_ext);
        summary.heating_current_drive.nbi{4}.direction.value        = sign(z0dstruct.z0dinput.option.angle_nbi2);
        summary.heating_current_drive.nbi{4}.energy.value           = z0dstruct.z0dinput.option.einj2 .* ones(size(summary.time));
        % summary.heating_current_drive.nbi{2}.beam_current_fraction.value
        % summary.heating_current_drive.nbi{2}.beam_power_fraction.value
        summary.heating_current_drive.nbi{4}.species.z_n.value      = 1;
        summary.heating_current_drive.nbi{4}.species.a.value = 1;
    end
else
    ftnbi = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.cons.ftnbi,data_zerod.temps,'pchip','extrap');
    summary.heating_current_drive.nbi{1}.power.value            = real(data_zerod.pnbi) .* (1 - real(ftnbi));
    summary.heating_current_drive.nbi{1}.current.value          = real(data_zerod.inbicd) .* (1 - real(ftnbi));
    summary.heating_current_drive.nbi{1}.tangency_radius.value  = z0dstruct.z0dinput.option.rtang;
    rext = abs(profil0d.Raxe(:,end) .* (1 + profil0d.epsi(:,end)) - z0dstruct.z0dinput.option.rtang * ones(size(profil0d.temps)));
    zext = z0dstruct.z0dinput.option.zext .* ones(size(profil0d.temps));
    if length(profil0d.temps) == 1
        kext = interp1_imas(profil0d.xli,profil0d.kx,zext,'pchip','extrap');
    else
        [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
        kext = griddata_imas(tt,xx,profil0d.kx,profil0d.temps,zext,'cubic');
    end
    zext = zext .* kext .*  profil0d.Raxe(:,end) .* profil0d.epsi(:,end);
    angle_ext = atan(zext ./ rext);
    summary.heating_current_drive.nbi{1}.angle.value            = mean(angle_ext);
    summary.heating_current_drive.nbi{1}.direction.value        = sign(z0dstruct.z0dinput.option.angle_nbi);
    summary.heating_current_drive.nbi{1}.energy.value           = z0dstruct.z0dinput.option.einj .* ones(size(summary.time));
    % summary.heating_current_drive.nbi{1}.beam_current_fraction.value
    % summary.heating_current_drive.nbi{1}.beam_power_fraction.value
    switch z0dstruct.z0dinput.option.gaz
        case 11
            summary.heating_current_drive.nbi{1}.species.z_n.value      = 1;
            summary.heating_current_drive.nbi{1}.species.a.value        = 1;
        otherwise
            
            summary.heating_current_drive.nbi{1}.species.z_n.value      = 1;
            summary.heating_current_drive.nbi{1}.species.a.value        = 2;
    end
    %
    summary.heating_current_drive.nbi{2} = summary.heating_current_drive.nbi{1};
    summary.heating_current_drive.nbi{2}.power.value            = real(data_zerod.pnbi)   .* real(ftnbi);
    summary.heating_current_drive.nbi{2}.current.value          = real(data_zerod.inbicd) .* real(ftnbi);
    summary.heating_current_drive.nbi{2}.tangency_radius.value  = z0dstruct.z0dinput.option.rtang;
    rext = abs(profil0d.Raxe(:,end) .* (1 + profil0d.epsi(:,end)) - z0dstruct.z0dinput.option.rtang * ones(size(profil0d.temps)));
    zext = z0dstruct.z0dinput.option.zext .* ones(size(profil0d.temps));
    if length(profil0d.temps) == 1
        kext = interp1_imas(profil0d.xli,profil0d.kx,zext,'pchip','extrap');
    else
        [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
        kext = griddata_imas(tt,xx,profil0d.kx,profil0d.temps,zext,'cubic');
    end
    zext = zext .* kext .*  profil0d.Raxe(:,end) .* profil0d.epsi(:,end);
    angle_ext = atan(zext ./ rext);
    summary.heating_current_drive.nbi{2}.angle.value            = mean(angle_ext);
    summary.heating_current_drive.nbi{2}.direction.value        = sign(z0dstruct.z0dinput.option.angle_nbi);
    summary.heating_current_drive.nbi{2}.energy.value           = z0dstruct.z0dinput.option.einj .* ones(size(summary.time));
    % summary.heating_current_drive.nbi{2}.beam_current_fraction.value
    % summary.heating_current_drive.nbi{2}.beam_power_fraction.value
    summary.heating_current_drive.nbi{2}.species.z_n.value      = 1;
    switch z0dstruct.z0dinput.option.gaz
        case 3
            summary.heating_current_drive.nbi{2}.species.z_n.value      = 1;
            summary.heating_current_drive.nbi{2}.species.a.value = 3;
        case 5
            summary.heating_current_drive.nbi{2}.species.z_n.value      = 2;
            summary.heating_current_drive.nbi{2}.species.a.value = 3;
        case 11
            summary.heating_current_drive.nbi{2}.species.z_n.value      = 5;
            summary.heating_current_drive.nbi{2}.species.a.value = 11;
        otherwise
            summary.heating_current_drive.nbi{2}.species.z_n.value      = 1;
            summary.heating_current_drive.nbi{2}.species.a.value = 1;
    end
    
    if z0dstruct.z0dinput.option.nb_nbi > 1
        summary.heating_current_drive.nbi{3} = summary.heating_current_drive.nbi{1};
        summary.heating_current_drive.nbi{3}.power.value            = imag(data_zerod.pnbi)   .* (1 - imag(ftnbi));
        summary.heating_current_drive.nbi{3}.current.value          = imag(data_zerod.inbicd) .* (1 - imag(ftnbi));
        summary.heating_current_drive.nbi{3}.tangency_radius.value  = z0dstruct.z0dinput.option.rtang2;
        rext = abs(profil0d.Raxe(:,end) .* (1 + profil0d.epsi(:,end)) - z0dstruct.z0dinput.option.rtang2 * ones(size(profil0d.temps)));
        zext = z0dstruct.z0dinput.option.zext2 .* ones(size(profil0d.temps));
        if length(profil0d.temps) == 1
            kext = interp1_imas(profil0d.xli,profil0d.kx,zext,'pchip','extrap');
        else
            [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
            kext = griddata_imas(tt,xx,profil0d.kx,profil0d.temps,zext,'cubic');
        end
        zext = zext .* kext .*  profil0d.Raxe(:,end) .* profil0d.epsi(:,end);
        angle_ext = atan(zext ./ rext);
        summary.heating_current_drive.nbi{3}.angle.value            = mean(angle_ext);
        summary.heating_current_drive.nbi{3}.direction.value        = sign(z0dstruct.z0dinput.option.angle_nbi2);
        summary.heating_current_drive.nbi{3}.energy.value           = z0dstruct.z0dinput.option.einj2 .* ones(size(summary.time));
        switch z0dstruct.z0dinput.option.gaz
            case 11
                summary.heating_current_drive.nbi{3}.species.z_n.value      = 1;
                summary.heating_current_drive.nbi{3}.species.a.value        = 1;
            otherwise
                
                summary.heating_current_drive.nbi{3}.species.z_n.value      = 1;
                summary.heating_current_drive.nbi{3}.species.a.value        = 2;
        end
        % summary.heating_current_drive.nbi{2}.beam_current_fraction.value
        % summary.heating_current_drive.nbi{2}.beam_power_fraction.value
        summary.heating_current_drive.nbi{4} = summary.heating_current_drive.nbi{1};
        summary.heating_current_drive.nbi{4}.power.value            = imag(data_zerod.pnbi)   .* imag(ftnbi);
        summary.heating_current_drive.nbi{4}.current.value          = imag(data_zerod.inbicd) .* imag(ftnbi);
        summary.heating_current_drive.nbi{4}.tangency_radius.value  = z0dstruct.z0dinput.option.rtang2;
        rext = abs(profil0d.Raxe(:,end) .* (1 + profil0d.epsi(:,end)) - z0dstruct.z0dinput.option.rtang2 * ones(size(profil0d.temps)));
        zext = z0dstruct.z0dinput.option.zext2 .* ones(size(profil0d.temps));
        if length(profil0d.temps) == 1
            kext = interp1_imas(profil0d.xli,profil0d.kx,zext,'pchip','extrap');
        else
            [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
            kext = griddata_imas(tt,xx,profil0d.kx,profil0d.temps,zext,'cubic');
        end
        zext = zext .* kext .*  profil0d.Raxe(:,end) .* profil0d.epsi(:,end);
        angle_ext = atan(zext ./ rext);
        summary.heating_current_drive.nbi{4}.angle.value            = mean(angle_ext);
        summary.heating_current_drive.nbi{4}.direction.value        = sign(z0dstruct.z0dinput.option.angle_nbi2);
        summary.heating_current_drive.nbi{4}.energy.value           = z0dstruct.z0dinput.option.einj2 .* ones(size(summary.time));
        % summary.heating_current_drive.nbi{2}.beam_current_fraction.value
        % summary.heating_current_drive.nbi{2}.beam_power_fraction.value
        switch z0dstruct.z0dinput.option.gaz
            case 3
                summary.heating_current_drive.nbi{4}.species.z_n.value      = 1;
                summary.heating_current_drive.nbi{4}.species.a.value = 3;
            case 5
                summary.heating_current_drive.nbi{4}.species.z_n.value      = 2;
                summary.heating_current_drive.nbi{4}.species.a.value = 3;
            case 11
                summary.heating_current_drive.nbi{4}.species.z_n.value      = 5;
                summary.heating_current_drive.nbi{4}.species.a.value = 11;
            otherwise
                summary.heating_current_drive.nbi{4}.species.z_n.value      = 1;
                summary.heating_current_drive.nbi{4}.species.a.value = 1;
        end
        
%         summary.heating_current_drive.nbi{4}.species.z_n.value      = 1;
%         if z0dstruct.z0dinput.option.gaz == 3
%             summary.heating_current_drive.nbi{4}.species.a.value = 3;
%         else
%             summary.heating_current_drive.nbi{4}.species.a.value = 1;
%         end
    end
end
%% ITB LOCAL QUANTITIES
if length(profil0d.temps) == 1
  xx = profil0d.xli;
  tt = profil0d.temps;
else
  [xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
end
summary.local.itb.q.value = griddata_imas(tt,xx,profil0d.qjli,data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.t_e.value = griddata_imas(tt,xx,profil0d.tep,data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.t_i_average.value = griddata_imas(tt,xx,profil0d.tip,data_zerod.temps,data_zerod.xitb,'cubic'); 
summary.local.itb.n_e.value = griddata_imas(tt,xx,profil0d.nep,data_zerod.temps,data_zerod.xitb,'cubic'); 
summary.local.itb.n_i_total.value = griddata_imas(tt,xx,profil0d.nip,data_zerod.temps,data_zerod.xitb,'cubic'); 
summary.local.itb.position.psi = griddata_imas(tt,xx,profil0d.psi,data_zerod.temps,data_zerod.xitb,'cubic');  
summary.local.itb.position.rho_tor = griddata_imas(tt,xx,profil0d.rmx,data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.position.rho_tor_norm = summary.local.itb.position.rho_tor ./ summary.local.edge.position.rho_tor;
summary.local.itb.zeff.value            = griddata_imas(tt,xx,profil0d.zeff,data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.momentum_tor.value    = griddata_imas(tt,xx,profil0d.rtor,data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.magnetic_shear.value  = griddata_imas(tt,xx,magnetic_shear,data_zerod.temps,data_zerod.xitb,'cubic');
%width = interp1_imas(profil0d.temps,mean(double(profil0d.xieshape_itb < profil0d.xieshape),2),data_zerod.temps,'pchip','extrap') .* data_zerod.rm;
%summary.local.itb.width_itb.value 		= width;
%% FILL VTOR FOR KNOWN SPECIES
summary.local.itb.velocity_tor.hydrogen.value = griddata_imas(tt,xx,vtor(:,:,1),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.velocity_tor.deuterium.value = griddata_imas(tt,xx,vtor(:,:,2),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.velocity_tor.tritium.value = griddata_imas(tt,xx,vtor(:,:,3),data_zerod.temps,data_zerod.xitb,'cubic');
% temporary, until spelling fixed in the summary IDS
summary.local.itb.velocity_tor.helium_3.value = griddata_imas(tt,xx,vtor(:,:,4),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.velocity_tor.helium_4.value = griddata_imas(tt,xx,vtor(:,:,5),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.velocity_tor.tungsten.value = griddata_imas(tt,xx,vtor(:,:,8),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.velocity_tor.tin.value = griddata_imas(tt,xx,vtor(:,:,9),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.velocity_tor.boron.value = griddata_imas(tt,xx,vtor(:,:,10),data_zerod.temps,data_zerod.xitb,'cubic');

%% FILL VTOR TESTING ION SPECIES FOR POSITIONS 6 AND 7
summary.local.itb.velocity_tor.berylium.value = zeros(size(data_zerod.temps));
summary.local.itb.velocity_tor.lithium.value  = zeros(size(data_zerod.temps));
summary.local.itb.velocity_tor.carbon.value   = zeros(size(data_zerod.temps));
summary.local.itb.velocity_tor.nitrogen.value = zeros(size(data_zerod.temps));
summary.local.itb.velocity_tor.neon.value     = zeros(size(data_zerod.temps));
summary.local.itb.velocity_tor.argon.value    = zeros(size(data_zerod.temps));
summary.local.itb.velocity_tor.xenon.value    = zeros(size(data_zerod.temps));
summary.local.itb.velocity_tor.oxygen.value   = zeros(size(data_zerod.temps));
clear as labels
for ispec=6:7
  [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
  switch labels{ispec}
   case 'Be'
    summary.local.itb.velocity_tor.berylium.value = griddata_imas(tt,xx,vtor(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'Li'
    summary.local.itb.velocity_tor.lithium.value  = griddata_imas(tt,xx,vtor(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'C'
    summary.local.itb.velocity_tor.carbon.value   = griddata_imas(tt,xx,vtor(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'N'
    summary.local.itb.velocity_tor.nitrogen.value = griddata_imas(tt,xx,vtor(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'Ne'
    summary.local.itb.velocity_tor.neon.value     = griddata_imas(tt,xx,vtor(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'Ar'
    summary.local.itb.velocity_tor.argon.value    = griddata_imas(tt,xx,vtor(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'Xe'
    summary.local.itb.velocity_tor.xenon.value    = griddata_imas(tt,xx,vtor(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'O'
    summary.local.itb.velocity_tor.oxygen.value   = griddata_imas(tt,xx,vtor(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   otherwise
    warning('plasma composition not yet implemented in summary IDS');
  end
end

summary.local.itb.n_i.hydrogen.value = griddata_imas(tt,xx,nions(:,:,1),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.n_i.deuterium.value = griddata_imas(tt,xx,nions(:,:,2),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.n_i.tritium.value = griddata_imas(tt,xx,nions(:,:,3),data_zerod.temps,data_zerod.xitb,'cubic');
% temporary, until spelling fixed in the summary IDS
summary.local.itb.n_i.helium_3.value = griddata_imas(tt,xx,nions(:,:,4),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.n_i.helium_4.value = griddata_imas(tt,xx,nions(:,:,5),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.n_i.tungsten.value = griddata_imas(tt,xx,nions(:,:,8),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.n_i.tin.value = griddata_imas(tt,xx,nions(:,:,9),data_zerod.temps,data_zerod.xitb,'cubic');
summary.local.itb.n_i.boron.value = griddata_imas(tt,xx,nions(:,:,10),data_zerod.temps,data_zerod.xitb,'cubic');

%% FILL VTOR TESTING ION SPECIES FOR POSITIONS 6 AND 7
summary.local.itb.n_i.berylium.value = zeros(size(data_zerod.temps));
summary.local.itb.n_i.lithium.value  = zeros(size(data_zerod.temps));
summary.local.itb.n_i.carbon.value   = zeros(size(data_zerod.temps));
summary.local.itb.n_i.nitrogen.value = zeros(size(data_zerod.temps));
summary.local.itb.n_i.neon.value     = zeros(size(data_zerod.temps));
summary.local.itb.n_i.argon.value    = zeros(size(data_zerod.temps));
summary.local.itb.n_i.xenon.value    = zeros(size(data_zerod.temps));
summary.local.itb.n_i.oxygen.value   = zeros(size(data_zerod.temps));
clear as labels
for ispec=6:7
  [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
  switch labels{ispec}
   case 'Be'
    summary.local.itb.n_i.berylium.value = griddata_imas(tt,xx,nions(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'Li'
    summary.local.itb.n_i.lithium.value  = griddata_imas(tt,xx,nions(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'C'
    summary.local.itb.n_i.carbon.value   = griddata_imas(tt,xx,nions(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'N'
    summary.local.itb.n_i.nitrogen.value = griddata_imas(tt,xx,nions(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'Ne'
    summary.local.itb.n_i.neon.value     = griddata_imas(tt,xx,nions(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'Ar'
    summary.local.itb.n_i.argon.value    = griddata_imas(tt,xx,nions(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'Xe'
    summary.local.itb.n_i.xenon.value    = griddata_imas(tt,xx,nions(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   case 'O'
    summary.local.itb.n_i.oxygen.value   = griddata_imas(tt,xx,nions(:,:,ispec),data_zerod.temps,data_zerod.xitb,'cubic');
   otherwise
    warning('plasma composition not yet implemented in summary IDS');
  end
end



%summary.local.itb.itb_type.value 		= ((data_zerod.hitb .* data_zerod.hmhd) > 0) .* 8;  % a completer ulterieument
%

%  %% LIMITER / DIVERTOR / WALL
%  summary.local.lim_div_wall.t_e.value = data_zerod.telim;
%  summary.local.lim_div_wall.t_i_average.value = [];
%  summary.local.lim_div_wall.n_e.value = data_zerod.nelim;
%  summary.local.lim_div_wall.n_i_total.value = []; 
%  %summary.local.lim_div_wall.p_peak_div.value 	= data_zerod.peakdiv;
%  %summary.local.lim_div_wall.surf_temp.value 	= []; 
%  %summary.local.lim_div_wall.p_lim_div.value 	= data_zerod.plim;
%  %summary.local.lim_div_wall.p_rad_div.value 	= [];
%  %summary.local.lim_div_wall.wall_temp.value 	= []; 
%  %summary.local.lim_div_wall.wall_state.value 	= [];
%  %summary.local.lim_div_wall.detach_state.value 	= double(((data_zerod.prad + data_zerod.pbrem + data_zerod.pcyclo + data_zerod.pioniz + data_zerod.pradsol) >= data_zerod.ploss) | (data_zerod.telim < 1));
%  %summary.local.lim_div_wall.pump_flux.value 	= fuelling_data.pumping;

%zeffne = trapz(profil0d.xli,profil0d.nep .* profil0d.zeff,2) ./ max(eps,trapz(profil0d.xli,profil0d.nep,2));
%summary.line_average.ne_zeff_line.value 	= interp1_imas(profil0d.temps,zeffne,data_zerod.temps,'pchip','extrap'); 
%if length(data_zerod.temps) == 1
%  nbar = z0dstruct.zerod.nbar;
%  t     = z0dstruct.zerod.temps;
%  summary.line_average.dne_line_dt.value 	= (nbar(end) - nbar(end-1)) ./ (t(end) - t(end -1));
%else
%  summary.line_average.dne_line_dt.value 	= z0dxdt(data_zerod.nbar,data_zerod.temps);
%end
%

%  %% NEUTRONS (TO UPDATE, IF MY PROPOSAL FOR THE MODIFICATION OF THE SUMMARY IDS IS ACCEPTED
%  summary.neutron_fluxes.dd.total.value = data_zerod.ndd;
%  summary.neutron_fluxes.dd.thermal.value	= data_zerod.ndd_th; 
%  %summary.neutron_fluxes.dd.nbi_th.value	= data_zerod.ndd_nbi_th; 
%  %summary.neutron_fluxes.dd.nbi_nbi.value = data_zerod.ndd_nbi_nbi;
%  summary.neutron_fluxes.dt.total.value = data_zerod.pfus ./(3.56e6 .* 1.602176462e-19);
%  summary.neutron_fluxes.dt.thermal.value = max(0,data_zerod.pfus - data_zerod.pfus_nbi) ./ (3.56e6 .* 1.602176462e-19); 

%  %% AT RHO = 0.95 (NOT ADAPTED HERE, KEPT ONLY FOR NOT LOSING THE INFORMATION)
%  tt = profil0d.temps * ones(size(profil0d.xli));
%  psin = (profil0d.psi - profil0d.psi(:,end) *  ones(size(profil0d.xli))) ./ ((profil0d.psi(:,1) - profil0d.psi(:,end)) *  ones(size(profil0d.xli)));
%  v95  = 0.95 .* ones(size(data_zerod.temps));
%  summary.local.ninety_five.q_95.value 		= griddata_imas(tt,psin,profil0d.qjli,data_zerod.temps,v95,'cubic');
%  summary.local.ninety_five.elong_95.value 	= griddata_imas(tt,psin,profil0d.kx,data_zerod.temps,v95,'cubic');
%  summary.local.ninety_five.tria_95.value 	= griddata_imas(tt,psin,profil0d.dx,data_zerod.temps,v95,'cubic');
%  summary.local.ninety_five.tria_up_95.value 	= [];
%  summary.local.ninety_five.tria_lo_95.value 	= []; 
%  summary.local.ninety_five.te_95.value 		= griddata_imas(tt,psin,profil0d.tep,data_zerod.temps,v95,'cubic'); 
%  summary.local.ninety_five.ti_95.value 		= griddata_imas(tt,psin,profil0d.tip,data_zerod.temps,v95,'cubic'); 
%  summary.local.ninety_five.ne_95.value 		= griddata_imas(tt,psin,profil0d.nep,data_zerod.temps,v95,'cubic'); 
%  summary.local.ninety_five.ni_95.value 		= griddata_imas(tt,psin,profil0d.nip,data_zerod.temps,v95,'cubic'); 
%  summary.local.ninety_five.phi_95.value 		= griddata_imas(tt,psin,profil0d.phi,data_zerod.temps,v95,'cubic');   
%  summary.local.ninety_five.rho_95.value 		= griddata_imas(tt,psin,profil0d.rmx,data_zerod.temps,v95,'cubic');  
%  summary.local.ninety_five.vtor_95.value 	= griddata_imas(tt,psin,profil0d.vtor,data_zerod.temps,v95,'cubic');

%% PEDESTAL VALUES
summary.local.pedestal.t_e_ped.value 	= data_zerod.teped;
summary.local.pedestal.t_i_average.value= data_zerod.tiped;
summary.local.pedestal.n_e.value 	= data_zerod.neped; 
summary.local.pedestal.n_i_total.value 	= data_zerod.niped; 
summary.local.pedestal.position.psi 	= interp1_imas(profil0d.temps,profil0d.psi(:,end-1),data_zerod.temps,'pchip','extrap'); 
%					  (~data_zerod.modeh) .* summary.local.edge.position.phi; 
summary.local.pedestal.position.rho_tor	= interp1_imas(profil0d.temps,profil0d.rmx(:,end-1),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.position.rho_tor_norm = summary.local.pedestal.position.rho_tor ./ summary.local.edge.position.rho_tor; 
summary.local.pedestal.q.value 	        = interp1_imas(profil0d.temps,profil0d.qjli(:,end-1),data_zerod.temps,'pchip','extrap');  
summary.local.pedestal.zeff.value       = interp1_imas(profil0d.temps,profil0d.zeff(:,end-1),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.momentum_tor.value = interp1_imas(profil0d.temps,profil0d.rtor(:,end-1),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.magnetic_shear.value = interp1_imas(profil0d.temps,magnetic_shear(:,end-1),data_zerod.temps,'pchip','extrap');

%summary.local.pedestal.pressure_ped.value = data_zerod.modeh .* data_zerod.pped;
%% FILL VTOR FOR KNOWN SPECIES
summary.local.pedestal.velocity_tor.hydrogen.value  = interp1_imas(profil0d.temps,vtor(:,end-1,1),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.velocity_tor.deuterium.value = interp1_imas(profil0d.temps,vtor(:,end-1,2),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.velocity_tor.tritium.value   = interp1_imas(profil0d.temps,vtor(:,end-1,3),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.velocity_tor.helium_3.value  = interp1_imas(profil0d.temps,vtor(:,end-1,4),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.velocity_tor.helium_4.value  = interp1_imas(profil0d.temps,vtor(:,end-1,5),data_zerod.temps,'pchip','extrap'); 
summary.local.pedestal.velocity_tor.tungsten.value  = interp1_imas(profil0d.temps,vtor(:,end-1,8),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.velocity_tor.tin.value  = interp1_imas(profil0d.temps,vtor(:,end-1,9),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.velocity_tor.boron.value  = interp1_imas(profil0d.temps,vtor(:,end-1,10),data_zerod.temps,'pchip','extrap');

%% FILL VTOR TESTING ION SPECIES FOR POSITIONS 6 AND 7
summary.local.pedestal.velocity_tor.berylium.value = zeros(size(data_zerod.temps));
summary.local.pedestal.velocity_tor.lithium.value  = zeros(size(data_zerod.temps));
summary.local.pedestal.velocity_tor.carbon.value   = zeros(size(data_zerod.temps));
summary.local.pedestal.velocity_tor.nitrogen.value = zeros(size(data_zerod.temps));
summary.local.pedestal.velocity_tor.neon.value     = zeros(size(data_zerod.temps));
summary.local.pedestal.velocity_tor.argon.value    = zeros(size(data_zerod.temps));
summary.local.pedestal.velocity_tor.xenon.value    = zeros(size(data_zerod.temps));
summary.local.pedestal.velocity_tor.oxygen.value   = zeros(size(data_zerod.temps));
clear as labels
for ispec=6:7
  [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
  switch labels{ispec}
   case 'Be'
    summary.local.pedestal.velocity_tor.berylium.value = interp1_imas(profil0d.temps,vtor(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Li'
    summary.local.pedestal.velocity_tor.lithium.value  = interp1_imas(profil0d.temps,vtor(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'C'
    summary.local.pedestal.velocity_tor.carbon.value   = interp1_imas(profil0d.temps,vtor(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'N'
    summary.local.pedestal.velocity_tor.nitrogen.value = interp1_imas(profil0d.temps,vtor(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ne'
    summary.local.pedestal.velocity_tor.neon.value     = interp1_imas(profil0d.temps,vtor(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ar'
    summary.local.pedestal.velocity_tor.argon.value    = interp1_imas(profil0d.temps,vtor(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Xe'
    summary.local.pedestal.velocity_tor.xenon.value    = interp1_imas(profil0d.temps,vtor(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'O'
    summary.local.pedestal.velocity_tor.oxygen.value   = interp1_imas(profil0d.temps,vtor(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   otherwise
    warning('plasma composition not yet implemented in summary IDS');
  end
end

summary.local.pedestal.n_i.hydrogen.value   = interp1_imas(profil0d.temps,nions(:,end-1,1),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.n_i.deuterium.value  = interp1_imas(profil0d.temps,nions(:,end-1,2),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.n_i.tritium.value    = interp1_imas(profil0d.temps,nions(:,end-1,3),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.n_i.helium_3.value   = interp1_imas(profil0d.temps,nions(:,end-1,4),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.n_i.helium_4.value   = interp1_imas(profil0d.temps,nions(:,end-1,5),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.n_i.berylium.value   = zeros(size(data_zerod.temps));
summary.local.pedestal.n_i.lithium.value    = zeros(size(data_zerod.temps));
summary.local.pedestal.n_i.carbon.value     = zeros(size(data_zerod.temps));
summary.local.pedestal.n_i.nitrogen.value   = zeros(size(data_zerod.temps));
summary.local.pedestal.n_i.neon.value       = zeros(size(data_zerod.temps));
summary.local.pedestal.n_i.argon.value      = zeros(size(data_zerod.temps));
summary.local.pedestal.n_i.xenon.value      = zeros(size(data_zerod.temps));
summary.local.pedestal.n_i.oxygen.value      = zeros(size(data_zerod.temps));
summary.local.pedestal.n_i.tungsten.value   = interp1_imas(profil0d.temps,nions(:,1,8),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.n_i.tin.value   = interp1_imas(profil0d.temps,nions(:,1,9),data_zerod.temps,'pchip','extrap');
summary.local.pedestal.n_i.boron.value   = interp1_imas(profil0d.temps,nions(:,1,10),data_zerod.temps,'pchip','extrap');

clear as labels
for ispec=6:7
  [as(ispec),labels{ispec}]=getafromz_imas(Z(ispec));
  switch labels{ispec}
   case 'Be'
    summary.local.pedestal.n_i.berylium.value = interp1_imas(profil0d.temps,nions(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Li'
    summary.local.pedestal.n_i.lithium.value  = interp1_imas(profil0d.temps,nions(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'C'
    summary.local.pedestal.n_i.carbon.value   = interp1_imas(profil0d.temps,nions(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'N'
    summary.local.pedestal.n_i.nitrogen.value = interp1_imas(profil0d.temps,nions(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ne'
    summary.local.pedestal.n_i.neon.value     = interp1_imas(profil0d.temps,nions(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Ar'
    summary.local.pedestal.n_i.argon.value    = interp1_imas(profil0d.temps,nions(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'Xe'
    summary.local.pedestal.n_i.xenon.value    = interp1_imas(profil0d.temps,nions(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   case 'O'
    summary.local.pedestal.n_i.oxygen.value   = interp1_imas(profil0d.temps,nions(:,end-1,ispec),data_zerod.temps,'pchip','extrap');
   otherwise
    warning('plasma composition not yet implemented in summary IDS');
  end
end

%% NOT IN THE SUMMARY IDS BUT DON'T WANT TO LOSE THE INFORMATION
%  if length(data_zerod.temps) == 1
%    summary.references.plh.value 		= z0dstruct.z0dinput.cons.plh(end);
%    summary.references.picrh.value 	= z0dstruct.z0dinput.cons.picrh(end);
%    summary.references.pecrh.value 	= z0dstruct.z0dinput.cons.pecrh(end); 
%    summary.references.pnbi.value 	= z0dstruct.z0dinput.cons.pnbi(end);  
%    summary.references.ip.value		= z0dstruct.z0dinput.cons.ip(end); 
%  else
%    summary.references.plh.value 		= z0dstruct.z0dinput.cons.plh;
%    summary.references.picrh.value 	= z0dstruct.z0dinput.cons.picrh;
%    summary.references.pecrh.value 	= z0dstruct.z0dinput.cons.pecrh; 
%    summary.references.pnbi.value 	= z0dstruct.z0dinput.cons.pnbi;  
%    summary.references.ip.value		= z0dstruct.z0dinput.cons.ip; 
%  end
%  if length(data_zerod.temps) == 1
%    summary.references.bvac_r.value 	= sigma_bvac_r .* z0dstruct.z0dinput.geo.R(end) .* z0dstruct.z0dinput.geo.b0(end); 
%  else
%    summary.references.bvac_r.value 	= sigma_bvac_r .* z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0; 
%  end 
%  if length(data_zerod.temps) == 1
%    summary.references.zeffl.value 	= z0dstruct.z0dinput.cons.zeff(end);
%    summary.references.nbar.value		= real(z0dstruct.z0dinput.cons.nbar(end));
%    if imag(z0dstruct.z0dinput.cons.nbar(end))
%      summary.references.gas_puff.value   = imag(z0dstruct.z0dinput.cons.nbar(end));
%    end
%    summary.references.xecrh.value 	= z0dstruct.z0dinput.cons.xece(end);
%    summary.references.pol_flux.value 	= z0dstruct.z0dinput.cons.flux(end) .* 2 .* pi; 
%    summary.references.enhancement.value 	= z0dstruct.z0dinput.cons.hmore(end); 
%    summary.references.isotopic.value 	= z0dstruct.z0dinput.cons.iso(end); 
%    summary.references.nbi_td_ratio.value	= z0dstruct.z0dinput.cons.ftnbi(end);
%  else
%    summary.references.zeffl.value 	= z0dstruct.z0dinput.cons.zeff;
%    summary.references.nbar.value		= real(z0dstruct.z0dinput.cons.nbar);
%    if any(imag(z0dstruct.z0dinput.cons.nbar(end)))
%      summary.references.gas_puff.value   = imag(z0dstruct.z0dinput.cons.nbar);
%    end
%    summary.references.xecrh.value 	= z0dstruct.z0dinput.cons.xece;
%    summary.references.pol_flux.value 	= z0dstruct.z0dinput.cons.flux .* 2 .* pi; 
%    summary.references.enhancement.value 	= z0dstruct.z0dinput.cons.hmore; 
%    summary.references.isotopic.value 	= z0dstruct.z0dinput.cons.iso; 
%    summary.references.nbi_td_ratio.value	= z0dstruct.z0dinput.cons.ftnbi;
%  end

%% SOL QUANTITIES
%  if isfield(data_zerod,'dsol')
%    summary.local.sol.n_e.value 		= data_zerod.dsol;
%  elseif length(data_zerod.temps) == 1
%    summary.local.sol.n_e.value 		= z0dstruct.z0dinput.geo.a(end) ./ 100;
%  else
%    summary.local.sol.n_e.value 		= z0dstruct.z0dinput.geo.a ./ 100;
%  end
%  summary.local.sol.n_i_total.value 	= summary.local.sol.n_e.value;
%  summary.local.sol.qe.value 		= summary.local.sol.n_e.value ./ 0.62;
%  summary.local.sol.qi.value 		= summary.local.sol.qe.value;
%  summary.local.sol.t_e.value 		= 3 .* summary.local.sol.qe.value;  
%  summary.local.sol.t_i_average.value 	= 3 .* summary.local.sol.qe.value;
%  summary.local.sol.p_rad.value		= data_zerod.pradsol;
%  summary.local.sol.gaz_puff       	= fuelling_data.gas_puff;



% 
% composition
% 
%% KEEP IT BUT NOT AVAILABLE IN THE SUMMARY IDS
%  summary.composition.amn 		= [1,2,3,3,4,ceil(7/3 .* z0dstruct.z0dinput.option.zimp),ceil(7/3 .* z0dstruct.z0dinput.option.zmax),183.84];
%  summary.composition.zn	 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax,74];
%  summary.composition.zion 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax,74];
%  summary.composition.imp_flag 		= zeros(size(summary.composition.zion));
%  summary.composition.rot_imp_flag 	= [0,0,0,0,0,1,0,0];
%  switch z0dstruct.z0dinput.option.gaz
%   case 1
%    summary.composition.pellet_amn = 1;
%    summary.composition.pellet_zn 	= 1;
%    summary.composition.nbi_amn 	= [1,2];
%    summary.composition.nbi_zn 	= [1,1];
%  case 2
%   summary.composition.pellet_amn = 2;
%   summary.composition.pellet_zn 	= 1;
%   summary.composition.nbi_amn 	= [1,2];
%   summary.composition.nbi_zn 	= [1,1];
%   case 3
%    summary.composition.pellet_amn = [2,3];
%    summary.composition.pellet_zn 	= [1,1];
%    summary.composition.nbi_amn 	= [2,3];
%    summary.composition.nbi_zn 	= [1,1];
%   case 4
%    summary.composition.pellet_amn = [1,2];
%    summary.composition.pellet_zn 	= [1,1];
%    summary.composition.nbi_amn 	= [1,2];
%    summary.composition.nbi_zn 	= [1,1];
%   otherwise
%    error(sprintf('unknown gaz option %d',z0dstruct.z0dinput.option.gaz));
%  end 
%  summary.compositions = copy_composition_to_compositionstype(summary.composition);
summary.local.limiter.t_e.value  = data_zerod.telim .* double(data_zerod.xpoint == 0) - 9.0e40 .* double(data_zerod.xpoint ~= 0); 
summary.local.limiter.n_e.value  = data_zerod.nelim .* double(data_zerod.xpoint == 0) - 9.0e40 .* double(data_zerod.xpoint ~= 0); 
summary.local.divertor_plate{1}.t_e.value  = data_zerod.telim .* double(data_zerod.xpoint ~= 0) - 9.0e40 .* double(data_zerod.xpoint == 0); 
summary.local.divertor_plate{1}.n_e.value  = data_zerod.nelim .* double(data_zerod.xpoint ~= 0) - 9.0e40 .* double(data_zerod.xpoint == 0); 
switch z0dstruct.z0dinput.option.sol_model
case '2_points'
      zeff_div = summary_2points(z0dstruct,data_zerod,profil0d);
      summary.local.divertor_plate{1}.zeff.value = zeff_div .* double(data_zerod.xpoint ~= 0) - 9.0e40 .* double(data_zerod.xpoint == 0); 
end
%  summary.local.limiter.t_i_average.value
%  summary.local.limiter.n_i_total.value
%  summary.local.limiter.zeff.value
%  summary.local.limiter.flux_expansion.value
%  summary.local.limiter.n_i.hydrogen.value
%  summary.local.limiter.n_i.deuterium.value
%  summary.local.limiter.n_i.tritium.value
%  summary.local.limiter.n_i.helium_3.value
%  summary.local.limiter.n_i.helium_4.value
%  summary.local.limiter.n_i.berylium.value
%  summary.local.limiter.n_i.lithium.value
%  summary.local.limiter.n_i.carbon.value
%  summary.local.limiter.n_i.nitrogen.value
%  summary.local.limiter.n_i.neon.value
%  summary.local.limiter.n_i.argon.value
%  summary.local.limiter.n_i.xenon.value
%  summary.local.limiter.n_i.oxygen.value
%  summary.local.limiter.n_i.tungsten.value

%  summary.local.divertor_plate{1}.t_i_average.value
%  summary.local.divertor_plate{1}.n_i_total.value
%  summary.local.divertor_plate{1}.flux_expansion.value
%  summary.local.divertor_plate{1}.n_i.hydrogen.value
%  summary.local.divertor_plate{1}.n_i.deuterium.value
%  summary.local.divertor_plate{1}.n_i.tritium.value
%  summary.local.divertor_plate{1}.n_i.helium_3.value
%  summary.local.divertor_plate{1}.n_i.helium_4.value
%  summary.local.divertor_plate{1}.n_i.berylium.value
%  summary.local.divertor_plate{1}.n_i.lithium.value
%  summary.local.divertor_plate{1}.n_i.carbon.value
%  summary.local.divertor_plate{1}.n_i.nitrogen.value
%  summary.local.divertor_plate{1}.n_i.neon.value
%  summary.local.divertor_plate{1}.n_i.argon.value
%  summary.local.divertor_plate{1}.n_i.xenon.value
%  summary.local.divertor_plate{1}.n_i.oxygen.value
%  summary.local.divertor_plate{1}.n_i.tungsten.value



%
% end wrapping IDS summary
% 



function [h98,h96] = compute_H98(zs,geo,cons)

% rapport pour le facteur d'amplification (Ealpha +En)/Ealpha
rfan = (3.56e6 + 14.03e6) ./ 3.56e6 ;
pin  = zs.pin ./ 1e6;
ploss = max(1e-6,pin - (zs.pbrem + zs.pcyclo + (1/3) .* zs.prad + zs.pioniz) ./ 1e6);
ip = zs.ip ./ 1e6;
Bt = geo.b0;
ne = cons.nbar./ 1e19;
R  =  geo.R;
a  =  geo.a;
ep  = a ./ R;
K   = geo.K; 
Vp = zs.vp;
Ka  = Vp ./ (2*pi^2.*R.*a.^2);
meff = zs.meff;


% loi standard ITER
% ITERH-96P(th)        
tauthl  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
      R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s       

% ITERH-98P(y,2)        
tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
	R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s    
	
% output
h98    = zs.taue ./ tauh;
h96    = zs.taue ./ tauthl;

%% IMAS part start here
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


function zeff_div = summary_2points(z0dstruct,data_zerod,profil0d)

% precomputation
% donnees du modele a 2 points
zs      = data_zerod;
option  = z0dstruct.z0dinput.option; 
geo     = z0dstruct.z0dinput.geo; 
cons    = z0dstruct.z0dinput.cons; 
profli  = profil0d;
% compatibilite
if ~isfield(option,'yield_model')
 	option.yield_model      = 'Javev';
end

% puissance conduite a la separatrice
pl        = max(zs.pin .* sqrt(eps),zs.pin - zs.prad - zs.pbrem - zs.pcyclo - zs.pioniz);
% fraction perdue en volume dans la sol par rayonnement:
%fesol   = max(0,min(1, (zs.pradsol + max(0,1 - option.fprad) .* zs.prad) ./ max(1,pl)));
switch option.sol_rad
case 'coupled'
	fesol   = max(0,min(1, (zs.pradsol + max(0,1 - option.fprad) .* zs.prad) ./ max(1,pl)));
        pradsol = zs.pradsol + max(0,1 - option.fprad) .* zs.prad;
otherwise
	fesol   = max(0,min(1, zs.pradsol ./ max(1,pl)));
        pradsol = zs.pradsol;
end

% these E. Tsitrone
lclim = pi .* geo.R .* zs.qa;
lcpol = pi .* geo.R;
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

if isfield(zs,'dsol')
	dsol = zs.dsol;
elseif option.sol_lscale  == 0
    dsol        = geo.a ./ 100;
elseif option.sol_lscale  > 0
    dsol        = geo.a .* option.sol_lscale;
else
    dsol        = - geo.R .* option.sol_lscale;
end


% flux //  (formula 5.64 et 5.75 Stangeby)
x      = profli.xli;
ve     =  ones(size(x));
Raxea  = interp1(profli.temps,profli.Raxe,zs.temps,'pchip','extrap');
Fa     = interp1(profli.temps,profli.fdia,zs.temps,'pchip','extrap');
psi    = interp1(profli.temps,profli.psi,zs.temps,'pchip','extrap');
rmx    = interp1(profli.temps,profli.rmx,zs.temps,'pchip','extrap');
zeffp  = interp1(profli.temps,profli.zeff,zs.temps,'pchip','extrap');
nzp    = interp1(profli.temps,profli.nzp,zs.temps,'pchip','extrap');
nwp    = interp1(profli.temps,profli.nwp,zs.temps,'pchip','extrap');
nhep   = interp1(profli.temps,profli.nhep,zs.temps,'pchip','extrap');
n1p    = interp1(profli.temps,profli.n1p,zs.temps,'pchip','extrap');
nep    = interp1(profli.temps,profli.nep,zs.temps,'pchip','extrap');
Qe     = interp1(profli.temps,profli.qe(:,end),zs.temps,'pchip','extrap');
Qi     = interp1(profli.temps,profli.qi(:,end),zs.temps,'pchip','extrap');
rext         = Raxea + geo.a * x;
btor         = Fa ./ rext;
grho         = abs((rmx(:,end) * ve)./ max(eps,pdederive(x,rext,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,psi,0,2,2,1) ./ rext .* grho ./ (rmx(:,end) * ve);
ut           = atan(abs(bpol(:,end) ./  btor(:,end)));
% flux //  (formula 5.64 et 5.75 Stangeby)
%bpola = interp1(profli.temps,profli.bpol(:,end),zs.temps,'pchip','extrap');
%Fa = interp1(profli.temps,profli.fdia(:,end),zs.temps,'pchip','extrap');
%Raxea = interp1(profli.temps,profli.Raxe(:,end),zs.temps,'pchip','extrap');

%ut = atan(bpola ./  Fa .* Raxea);
%qpl_tot     = pl  ./ (4 .* pi  .* Raxea(:,end) .* dsol .* sin(ut));
Asol_para = 4 .* pi  .* Raxea(:,end) .* dsol .* sin(ut);
switch option.sol_model
case '2_points'
	% rien
otherwise
	warndlg('the 2 points model is not used in this simulation','2 points model');
	option.sol_model = '2_points';
end
option_mem = option;
option.plot2points = 'Yes';
[tebord,nelim,telim,qpl_target,err,nb,indbad,fmom,qpl_rad_div,qpl_neutral_div,qpl_tot_in,pl_in,zeff_div,gamma,mach_target,prad_loc,pradsol_loc,fcond,profli] = ...
    z0convergence_2points_dic(option,cons,geo,zs,profli);
option = option_mem;
qpl_rad_sol = option.fpower .* pradsol ./ (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut));
qpl_tot     = option.fpower .* pl ./ (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut));
qpl_in_max  = option.fpower .* zs.pin ./ (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut));
%(4*pi*rp*sol_width/q95)
%qpar=p_sep*1e6/(4*pi*rp*sol_width/q95); % W/m2 estimate of parallel q
%L=pi*q95*rp;
qpl_target  = qpl_tot -  qpl_rad_sol - qpl_rad_div - qpl_neutral_div;

fie = 1 + zs.tibord ./ zs.tebord;
tite_loc = zs.tibord ./ zs.tebord;
%fpe = min(1,max(0.1,zs.pel ./ (zs.pel + zs.pion)));
fpe = min(1,max(0.1,Qe ./ (Qe + Qi)));



if option.fmom == 0
	% longueur de ionisation pres des plaques
	[svi1s,svi2s,svcx,svrec,sii,sss,Ass] = z0sectionh(zs.telim,zs.telim .* (zs.tibord ./ zs.tebord));  % attention ici le rapport ti/te est calculer a l'exterieur, tebord ne doit pas etre mis a jour
	%% equilibre entre 1s et 2s pour le neutres de centre 
	alphas = nelim .* sss ./ (nelim .* sss + Ass);
	% etat d'equilibre 1s/2s
	sviss  = svi2s .* alphas + svi1s .* (1 - alphas);
	alpha  = (sviss  + sii)  ./ (sviss  + sii+ svcx);
	%flimx  = min(1,exp(telim - option.eioniz));
	fmom_corr = max(0.1,min(1,0.17 + 2 .* (alpha ./ (alpha + 1)) .^ ((alpha + 1) ./ 2)));
else
	fmom_corr = option.fmom .* ones(size(zs.telim));		
end


% flux de matiere dans la sol
switch option.gaz
    case 1
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 1);
    case 2
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 2);
    case {3,5}
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
            (2 + 3.* real(cons.iso)) .* (1 + real(cons.iso)));
        if option.gaz == 5
                warning('nHe3onD & nTonD not yet implemented !');
        end
    case 4
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 4);
    case 11
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
            (1 + 11 .* cons.iso) .* (1 + cons.iso));
end
flux_target = mach_target .* cs0 .* sqrt(zs.telim) .* zs.nelim;
neu =  (4 .* zs.telim .* zs.nelim) ./ (fmom .* fie) ./ zs.tebord;

if option.mach_corr == 1
	ftm = (gamma .* (1 + mach_target .^ 2) ./ (2 .* abs(mach_target) .* (gamma - 1 + mach_target .^ 2))) .^ 2;
	fnm = 2 ./ (ftm .* (1 + mach_target .^ 2));
	neu = neu ./ ftm ./ fnm;
        flux_target = flux_target .* mach_target;
end
flux_output = interp1(profli.temps,profli.ge(:,end) .* profli.grho2(:,end) .* profli.vpr_tor(:,end),zs.temps,'pchip','extrap');

maskx = ones(size(zs.temps));
maskx(zs.xpoint == 0) = NaN;
