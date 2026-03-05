function [yneo_out] = neo_std(neoparams,option,z0dinput,profil,itime,twp);
currentpath=pwd;

ee = 1.6e-19;   

% Initialisations
Chii = zeros(1,length(profil.xli));

yneo_out.Gion1 = zeros(1,length(profil.xli));
yneo_out.Gion2 = zeros(1,length(profil.xli));
yneo_out.Gion3 = zeros(1,length(profil.xli));
yneo_out.Gion4 = zeros(1,length(profil.xli));
yneo_out.Gion5 = zeros(1,length(profil.xli));

yneo_out.Gion1_nclass = zeros(1,length(profil.xli));
yneo_out.Gion2_nclass = zeros(1,length(profil.xli));
yneo_out.Gion3_nclass = zeros(1,length(profil.xli));
yneo_out.Gion4_nclass = zeros(1,length(profil.xli));
yneo_out.Gion5_nclass = zeros(1,length(profil.xli));

yneo_out.Qion1 = zeros(1,length(profil.xli));
yneo_out.Qion1_nclass = zeros(1,length(profil.xli));

% Path to NEO routines, need to provide something for a TGLF install function
path_to_METIS_coef = fullfile(fileparts(which('install_QLKNN10D')));
path_to_METIS_NEO = [path_to_METIS_coef '/NEO_TOOLS'];
addpath(path_to_METIS_NEO);

%neo_input = [path_to_METIS_NEO '/input_templates'];
neoin = read_neoinput('input_templates',[path_to_METIS_NEO '/']);

%neoparams = declare_NEO_parameters(); % Call the function which consists of the input parameters for TGLF.

neoin.SIM_MODEL = neoparams.SIM_MODEL;
neoin.EQUILIBRIUM_MODEL = neoparams.EQUILIBRIUM_MODEL;
neoin.COLLISION_MODEL = neoparams.COLLISION_MODEL;
neoin.N_ENERGY = neoparams.N_ENERGY;
neoin.N_XI = neoparams.N_XI;
neoin.N_THETA = neoparams.N_THETA;

irok = find(profil.xli>0);
irok = 1:length(profil.xli);
Q = interp1(profil.temps,profil.qjli(:,irok),profil.temps(itime));
insize = size(Q);
Dw = zeros(neoin.N_SPECIES,length(profil.xli(:,irok)));
Vw = zeros(neoin.N_SPECIES,length(profil.xli(:,irok)));
Flux = zeros(neoin.N_SPECIES,length(profil.xli(:,irok)));


%Two runs per species, one normal profile and another shifted profile
if neoparams.DV == 2  
  %Loop for each positions
  for ii=2:(length(profil.xli(:,irok)))

    %Loop for considered impurities
    for tt = [0, 1, 2, 3, 4, 5]  %for normal run, z1 et z2
      
      if tt == 2 
        calc_DV = 2;
      else  calc_DV = tt;
      end
      disp(calc_DV);
      [neoin_mod.(sprintf('s%d', tt)), ...
       ynorm.(sprintf('s%d', tt)), ...
       yspec.(sprintf('s%d', tt))] = build_local_neo_inputs(option,z0dinput,profil,neoin,profil.xli(ii),twp,1,calc_DV);

      if isempty(getappdata(0,'local'))
         runpath = [path_to_METIS_NEO '/runs_neo/'];
      else
         runpath = [getappdata(0,'local') '/runs_neo/'];
      end
      
      neo_input_scan = sprintf('neo_tmp_radial_s%d_%d', tt, ii);
  
      make_neoinput(neoin_mod.(sprintf('s%d', tt)), neo_input_scan, runpath);
      disp(['==============     RUNS NEO (' sprintf('tt=%d', tt) ') at r/a=' num2str(round(profil.xli(ii) * 100) / 100) '   =================='])
      
      [sout, tout] = unix(sprintf(['module purge; module load gcc/13.2.0 lapack/3.9.1-gcc-13.2.0 openmpi/4.1.1-gcc-13.2.0; export GACODE_PLATFORM=IRFM; export GACODE_ROOT=/ZONE_TRAVAIL/PM258588/codes/gacode; . ${GACODE_ROOT}/shared/bin/gacode_setup; ' ...
                                  'cd %s; nice -n 9 neo -e ''./'' -n 16'], [runpath neo_input_scan]));
      
      disp('==============     ENDS OF NEO RUN        ==================')
  
      neoout.(sprintf('s%d', tt)) = read_neo_transport([runpath neo_input_scan], 0);
      nclass_out.(sprintf('s%d', tt)) = read_neo_out_nclass([runpath neo_input_scan]);

    end

    %Computation of the transport coefficients
    for tt = [1, 2, 3, 4, 5]
      %run without shifted profils
      flux_SI_0=neoout.(sprintf('s%d', 0)).(sprintf('pflux%d', tt+1))*ynorm.(sprintf('s%d', 0)).nnorm*ynorm.(sprintf('s%d', 0)).vnorm;
      n_0=neoin_mod.(sprintf('s%d', 0)).(sprintf('DENS_%d', tt+1))*ynorm.(sprintf('s%d', 0)).nnorm;
      ALn_0=neoin_mod.(sprintf('s%d', 0)).(sprintf('DLNNDR_%d', tt+1));
      
      if tt == 4
        %for H
        flux_SI_pp_s=neoout.(sprintf('s%d', tt)).(sprintf('pflux%d', tt+1))*ynorm.(sprintf('s%d', tt)).nnorm*ynorm.(sprintf('s%d', tt)).vnorm;
        ALn_pp_s=neoin_mod.(sprintf('s%d', tt)).(sprintf('DLNNDR_%d', tt+1));
        Dw(tt+1,ii)=(1/n_0)*(ynorm.(sprintf('s%d', tt)).anorm*flux_SI_pp_s-ynorm.(sprintf('s%d', 0)).anorm*flux_SI_0)/(ALn_pp_s-ALn_0);
        Vw(tt+1,ii)=ynorm.(sprintf('s%d', 0)).anorm*flux_SI_0/n_0-ALn_0*Dw(tt,ii);
        Flux(tt+1,ii)=neoout.(sprintf('s%d', 0)).(sprintf('pflux%d', tt+1))*ynorm.(sprintf('s%d', 0)).nnorm*ynorm.(sprintf('s%d', 0)).vnorm;

        %for the main ions also
        flux_SI_0_i=neoout.(sprintf('s%d', 0)).(sprintf('pflux%d', 1))*ynorm.(sprintf('s%d', 0)).nnorm*ynorm.(sprintf('s%d', 0)).vnorm;
        n_0_i=neoin_mod.(sprintf('s%d', 0)).(sprintf('DENS_%d', 1))*ynorm.(sprintf('s%d', 0)).nnorm;
        ALn_0_i=neoin_mod.(sprintf('s%d', 0)).(sprintf('DLNNDR_%d', 1));
        flux_SI_pp_s=neoout.(sprintf('s%d', tt)).(sprintf('pflux%d', 1))*ynorm.(sprintf('s%d', tt)).nnorm*ynorm.(sprintf('s%d', tt)).vnorm;
        ALn_pp_s=neoin_mod.(sprintf('s%d', tt)).(sprintf('DLNNDR_%d', 1));
        Dw(1,ii)=(1/n_0_i)*(ynorm.(sprintf('s%d', tt)).anorm*flux_SI_pp_s-ynorm.(sprintf('s%d', 0)).anorm*flux_SI_0_i)/(ALn_pp_s-ALn_0_i);
        Vw(1,ii)=ynorm.(sprintf('s%d', 0)).anorm*flux_SI_0_i/n_0_i-ALn_0_i*Dw(1,ii);
        Flux(1,ii)=neoout.(sprintf('s%d', 0)).(sprintf('pflux%d', 1))*ynorm.(sprintf('s%d', 0)).nnorm*ynorm.(sprintf('s%d', 0)).vnorm;
      else
        %shifted components (a=anorm)
        flux_SI_pp_s=neoout.(sprintf('s%d', tt)).(sprintf('pflux%d', tt+1))*ynorm.(sprintf('s%d', tt)).nnorm*ynorm.(sprintf('s%d', tt)).vnorm;
        ALn_pp_s=neoin_mod.(sprintf('s%d', tt)).(sprintf('DLNNDR_%d', tt+1));
        Dw(tt+1,ii)=(1/n_0)*(ynorm.(sprintf('s%d', tt)).anorm*flux_SI_pp_s-ynorm.(sprintf('s%d', 0)).anorm*flux_SI_0)/(ALn_pp_s-ALn_0);
        Vw(tt+1,ii)=ynorm.(sprintf('s%d', 0)).anorm*flux_SI_0/n_0-ALn_0*Dw(tt,ii);
        Flux(tt+1,ii)=neoout.(sprintf('s%d', 0)).(sprintf('pflux%d', tt+1))*ynorm.(sprintf('s%d', 0)).nnorm*ynorm.(sprintf('s%d', 0)).vnorm;
      end
      
    end
  end

  %For others species
  for pp = 1:6
    D(:,:,pp) = reshape(Dw(pp,:),insize);
    V(:,:,pp)   = reshape(Vw(pp,:),insize);
    G(:,:,pp) = reshape(Flux(pp,:),insize);
  end

  % D(:,:,6) = zeros(size(D(:,:,2)));
  % V(:,:,6)   = zeros(size(D(:,:,2)));
  % G(:,:,6) = zeros(size(D(:,:,2)));

  yneo_out.G=G;
  yneo_out.D=D;
  yneo_out.V=V;
  
%One run for shifted profil and one normal (only for traces impurities)
elseif neoparams.DV == 1  
  %Loop for each positions
  for ii=2:(length(profil.xli(:,irok)))

    %Loop for considered impurities
    for tt = [0,1]  %for normal run, z1 et z3
      
      if tt == 1 
        calc_DV = 6;
      else  calc_DV = tt;
      end
      disp(calc_DV);
      [neoin_mod.(sprintf('s%d', tt)), ...
      ynorm.(sprintf('s%d', tt)), ...                     
      yspec.(sprintf('s%d', tt))] = build_local_neo_inputs(option,z0dinput,profil,neoin,profil.xli(ii),twp,1,calc_DV);
       
      if isempty(getappdata(0,'local'))
         runpath = [path_to_METIS_NEO '/runs_neo/'];
      else
         runpath = [getappdata(0,'local') '/runs_neo/'];
      end
      
      neo_input_scan = sprintf('neo_tmp_radial_s%d_%d', tt, ii);
  
      make_neoinput(neoin_mod.(sprintf('s%d', tt)), neo_input_scan, runpath);
      disp(['==============     RUNS NEO (' sprintf('tt=%d', tt) ') at r/a=' num2str(round(profil.xli(ii) * 100) / 100) '   =================='])
            [sout, tout] = unix(sprintf(['module purge; module load gcc/13.2.0 lapack/3.9.1-gcc-13.2.0 openmpi/4.1.1-gcc-13.2.0; export GACODE_PLATFORM=IRFM; export GACODE_ROOT=/ZONE_TRAVAIL/PM258588/codes/gacode; . ${GACODE_ROOT}/shared/bin/gacode_setup; ' ...
					   'cd %s; nice -n 9 neo -e ''./'' -n 16'], [runpath neo_input_scan]));	    
      
      disp('==============     ENDS OF NEO RUN        ==================')
  
      neoout.(sprintf('s%d', tt)) = read_neo_transport([runpath neo_input_scan], 0);
      nclass_out.(sprintf('s%d', tt)) = read_neo_out_nclass([runpath neo_input_scan]);
      
    end

    for tt = 1:neoin.N_SPECIES
      %run without shifted profils
      flux_SI_0=neoout.(sprintf('s%d', 0)).(sprintf('pflux%d', tt))*ynorm.(sprintf('s%d', 0)).nnorm*ynorm.(sprintf('s%d', 0)).vnorm;
      n_0=neoin_mod.(sprintf('s%d', 0)).(sprintf('DENS_%d', tt))*ynorm.(sprintf('s%d', 0)).nnorm;
      ALn_0=neoin_mod.(sprintf('s%d', 0)).(sprintf('DLNNDR_%d', tt));
      
      
      %shifted components (a=anorm)
      flux_SI_pp_s=neoout.(sprintf('s%d', 1)).(sprintf('pflux%d', tt))*ynorm.(sprintf('s%d', 1)).nnorm*ynorm.(sprintf('s%d', 1)).vnorm;
      ALn_pp_s=neoin_mod.(sprintf('s%d', 1)).(sprintf('DLNNDR_%d', tt));
      Dw(tt,ii)=(1/n_0)*(ynorm.(sprintf('s%d', 1)).anorm*flux_SI_pp_s-ynorm.(sprintf('s%d', 0)).anorm*flux_SI_0)/(ALn_pp_s-ALn_0);
      Vw(tt,ii)=ynorm.(sprintf('s%d', 0)).anorm*flux_SI_0/n_0-ALn_0*Dw(tt,ii);
      Flux(tt,ii) = flux_SI_0;

    end
  end

  %Computation of the transport coefficients
  for pp = 1:neoin.N_SPECIES 
    NNout(:,1) = Dw(pp,:);
    NNout(:,2) = Vw(pp,:);
    NNout(:,3) = Flux(pp,:);
    D(:,:,pp) = reshape(NNout(:,1),insize);
    V(:,:,pp)   = reshape(NNout(:,2),insize);
    G(:,:,pp) = reshape(NNout(:,3),insize);
  end
  yneo_out.G=G;
  yneo_out.D=D;
  yneo_out.V=V;

else 
  for ii=2:length(profil.xli(:,irok))

    [neoin_mod ynorm yspec] = build_local_neo_inputs(option,z0dinput,profil,neoin,profil.xli(ii),twp,1,0);

    if isempty(getappdata(0,'local'))
      runpath = [path_to_METIS_NEO '/runs_neo/'];
    else
      runpath = [getappdata(0,'local') '/runs_neo/'];
    end
      
    neo_input_scan = ['neo_tmp_radial' num2str(ii)];

    make_neoinput(neoin_mod,neo_input_scan,runpath);
				% cd([runpath neo_input_scan])
    disp(['==============     RUNS NEO at r/a=' num2str(round(profil.xli(ii)*100)/100) '   =================='])
    [sout, tout] = unix(sprintf(['module purge; module load gcc/13.2.0 lapack/3.9.1-gcc-13.2.0 openmpi/4.1.1-gcc-13.2.0; export GACODE_PLATFORM=IRFM; export GACODE_ROOT=/ZONE_TRAVAIL/PM258588/codes/gacode; . ${GACODE_ROOT}/shared/bin/gacode_setup; cd %s; nice -n 9 neo -e "./" -n 16'], [runpath neo_input_scan]));	

    disp('==============     ENDS OF NEO RUN        ==================')
    
    neoout = read_neo_transport([runpath neo_input_scan],0);
    nclass_out=read_neo_out_nclass([runpath neo_input_scan]);

    for pp = 1:neoin.N_SPECIES 
      flux_SI_pp=neoout.(sprintf('pflux%d', pp))*ynorm.nnorm*ynorm.vnorm;
      Flux(pp,ii)=flux_SI_pp;
    end

    yneo_out.Qion1(ii) = neoout.qflux1.*ynorm.nnorm.*ynorm.tnorm.*ynorm.vnorm*ee;
    yneo_out.Qion1_nclass(ii) = nclass_out.qflux(1).*ynorm.nnorm.*ynorm.tnorm.*ynorm.vnorm*ee;

    for kk=1:neoin.N_SPECIES-1
      str_ev = ['yneo_out.Gion' num2str(kk) '(ii) = neoout.pflux' num2str(kk) '.*ynorm.nnorm.*ynorm.vnorm;'];
      eval(str_ev);
      str_ev = ['yneo_out.Gion' num2str(kk) '_nclass(ii) = nclass_out.pflux(' num2str(kk) ').*ynorm.nnorm.*ynorm.vnorm;'];
      eval(str_ev);   
    end

  end
  for pp = 1:neoin.N_SPECIES 
    NNout(:,1) = Flux(pp,:);
    G(:,:,pp) = reshape(NNout(:,1),insize);
  end
  yneo_out.G=G;
  yneo_out.D=NaN.*G;
  yneo_out.V=NaN.*G;  
end

cd(currentpath)


