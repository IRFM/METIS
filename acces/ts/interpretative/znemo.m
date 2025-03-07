function [sortie,particules,memoire] = znemo(param,proto,cons,geo,equi,injection, ...
                                          prof,neo,impur,phys,compo,gene,memoire)
%% ----------------------------------------------------------------------------------------
%% PURPOSE:
%% CALCULATION OF NBI SOURCES USING SINBAD-RELOADED AND SPOT
%% ----------------------------------------------------------------------------------------
%% CALL:
%% -----
%%   [SORTIE,MATIERE,MEMOIRE] = ZNEMO(PARAM,PROTO,CONS,GEO,EQUI,INJECTION, ...
%%                              PROF,NEO,IMPUR,PHY,COMPOSITION,GENE,MEMOIRE);
%% ----------------------------------------------------------------------------------------
%% INPUT:
%% ------
%%  PARAM     = FUNCTION PARAMETERS TUNED VIA THE INTERFACE
%%  PROTO     = PROTOTYPE FOR THE SOURCE STRUCTURE (FILLED WITH ZEROS)
%%  CONS      = POWER FOR EACH INJECTOR (W) (DATA.CONS.IDN)
%%  GEO       = PLASMA GEOMETRY (DATA.GEO)
%%  EQUI      = PLASMA EQUILIBRIUM (DATA.EQUI)
%%  INJECTION = NOT USED HERE, VARIABLE FOR GAS INJECTION, (DATA.CONS.C)
%%  PROF      = PLASMA KINETIC PROFILES (DATA.PROF) 
%%  NEO       = NOT USED HERE, NEOCLASSIC DATA (DATA.NEO)
%%  IMPUR     = STRUCTURE FOR IMPURITY CONCENTRATION (DATA.IMPUR)
%%  PHYS      = PHYSICAL CONSTANTS (PARAM.PHYS)
%%  COMPO     = PLASMA COMPOSITION (SPECIES CHARGE AND MASS) (PARAM.COMPO)
%%  GENE      = MISCELLANEOUS GLOBAL PAREMETERS (PARAM.GENE)
%%  MEMOIRE   = STRUCTURE CONTAINING THE LATEST CALCULATION OF THE MODULE
%%              (DATAK.MEMOIRE.IDN)
%% ----------------------------------------------------------------------------------------
%% OUTPUT:
%% -------
%%  SORTIE    = SOURCE STRUCTURE (SAME FORMAT AS THE INITIALIZATION STRUCTURE "PROTO")
%%  MATIERE   = MATTER SOURCE FOR EACH SPECIE
%%  MEMOIRE   = STRUCTURE CONTAINING THE LATEST CALCULATION OF THE MODULE
%%              (DATAK.MEMOIRE.IDN)
%% ----------------------------------------------------------------------------------------

%% INITIALIZATION MODE
if nargin <=1

  if nargin == 0
    pini_number = 1;
  else
    pini_number = param;
  end

  valeur.machine_nom           = 'ITER';
  type.machine_nom             = 'string';
  borne.machine_nom            = {'TS','JET','ITER','DIIID','AUG','COMPASS','JT-60U','JT-60SA','DEMO','FAST','WEST'};
  defaut.machine_nom           = 'ITER';
  info.machine_nom             = 'Tokamak'; 
  
  valeur.fokker                = 'FAST';
  type.fokker                  = 'string';
  borne.fokker                 = {'FAST','SPOT','RISK'};
  defaut.fokker                = 'FAST';
  info.fokker                  = 'Fokker-Planck module associated with NEMO'; 
  
  valeur.korbit                = 0;
  type.korbit                  = 'integer';
  borne.korbit                 = [0,1];
  defaut.korbit                = 0;
  info.korbit                  = 'Orbit width correction in NEMO (1) or not (0)'; 
  
  valeur.directivity           = ones(1,pini_number);
  type.directivity             = 'integer';
  borne.directivity            = {-1,1};
  defaut.directivity           = 1;
  info.directivity             = '1 -> co-current, -1 -> counter current';

  valeur.energie               = 80e3 .* ones(1,pini_number);
  type.energie                 = 'float';
  borne.energie                = [1,2e6];
  defaut.energie               = 80;
  info.energie                 = 'Beam energy (eV)';

  valeur.fraction1             = 0.8  .* ones(1,pini_number);
  type.fraction1               = 'float';
  borne.fraction1              = [0,1];
  defaut.fraction1             = 0.8;
  info.fraction1               = 'Fraction of full energy ions (-)';

  valeur.fraction2             = 0.15 .* ones(1,pini_number);
  type.fraction2               = 'float';
  borne.fraction2              = [0,1];
  defaut.fraction2             = 0.15;
  info.fraction2               = 'Fraction of half energy ions (-)';

  valeur.fraction3             = 0.05 .* ones(1,pini_number);
  type.fraction3               = 'float';
  borne.fraction3              = [0,1];
  defaut.fraction3             = 0.05;
  info.fraction3               = 'Fraction of third energy ions (-)';

  % a changer (cf. align)
  valeur.injector_config       = 'On-axis';
  type.injector_config         = 'string';
  borne.injector_config        = {'On-axis','Off-axis','On-off-axis','On-axis_2box'};
  defaut.injector_config       = 'On-axis';
  info.injector_config         = 'Configuration of the injector (not used if only one config. is possible)';

  valeur.align                 = zeros(1,pini_number);
  type.align                   = 'integer';
  borne.align                  = {0,1};
  defaut.align                 = 0;
  info.align                   = 'Injector alignment: 0=STANDARD, 1=UPSHIFTED (not used if only one config possible)';

  % VARIABLE NOT USED FOR NOW
  valeur.type                  = 3;
  type.type                    = 'integer';
  borne.type                   = {1,2,3,4,5};
  defaut.type                  = 3;
  info.type                    = '(NOT USED) interaction:1)DT=>4HE+N 2)DD=>T+P 3)DD=>3HE+N 4)TT=>4HE+2N 5)D3HE=>4HE+P';
  
  % VARIABLE NOT USED FOR NOW (=> NO CHOICE: H ISOTOPES ONLY)
  valeur.charge                = 1 .* ones(1,pini_number);
  type.charge                  = 'integer';
  borne.charge                 = {1};
  defaut.charge                = 1;
  info.charge                  = '(NOT USED) Z charge of the injected ion (not a large choice for now) (-)';

  valeur.masse                 = 2 .* ones(1,pini_number);
  type.masse                   = 'integer';
  borne.masse                  = {1,2,3,4};
  defaut.masse                 = 2;
  info.masse                   = 'A atomic mass of the injected ion (U)';

  valeur.n_out_profiles        = 20;
  type.n_out_profiles          = 'integer';
  borne.n_out_profiles         = [5 50];
  defaut.n_out_profiles        = [20];
  info.n_out_profiles          = 'Resolution of ouptut 1D-profiles';
  
  valeur.n_output_2d           = 41;
  type.n_output_2d             = 'integer';
  borne.n_output_2d            = [20 200];
  defaut.n_output_2d           = 41;
  info.n_output_2d             = 'Resolution of ouptut 2D (R,Z) grids';

  valeur.n_pitch_resol         = 101;
  type.n_pitch_resol           = 'integer';
  borne.n_pitch_resol          = [21 101];
  defaut.n_pitch_resol         = 101;
  info.n_pitch_resol           = 'Resolution of ouptut pitch angle profiles';

  valeur.spot_rlong            = 21;
  type.spot_rlong              = 'real';
  borne.spot_rlong             = [10 101];
  defaut.spot_rlong            = 21;
  info.spot_rlong              = 'Size of radial output array in SPOT (-)'; 

  valeur.spot_nproc            = 1;
  type.spot_nproc              = 'integer';
  borne.spot_nproc             = [1,100];
  defaut.spot_nproc            = 1;
  info.spot_nproc              = 'Number of processors used in SPOT (if MPI is present) (-)';

  valeur.spot_ncreated         = 10;
  type.spot_ncreated           = 'real';
  borne.spot_ncreated          = [1,1000];
  defaut.spot_ncreated         = 10;
  info.spot_ncreated           = 'Number of MC particles generated at each SPOT time step (-)';

%%%

  valeur.spot_pitch_scat         = 1;
  type.spot_pitch_scat           = 'real';
  borne.spot_pitch_scat          = {0,1};
  defaut.spot_pitch_scat         = 1;
  info.spot_pitch_scat           = 'Activate/disable pitch angle scattering';

  valeur.spot_collision         = 1;
  type.spot_collision           = 'real';
  borne.spot_collision          = {0,1};
  defaut.spot_collision         = 1;
  info.spot_collision           = 'Activate/disable collisions (slowing down + energy diffusion)';

  valeur.spot_engy_diff         = 1;
  type.spot_engy_diff           = 'real';
  borne.spot_engy_diff          = {0,1};
  defaut.spot_engy_diff         = 1;
  info.spot_engy_diff           = 'Activate/disable energy diffusion';

%%%
  
  valeur.spot_verbose          = 'No';
  type.spot_verbose            = 'string';
  borne.spot_verbose           = {'Yes','No'};
  defaut.spot_verbose          = 'Yes';
  info.spot_verbose            = 'Flag = 0 for silent SPOT, =1 for talkative SPOT';

  valeur.spot_init             = 'Zero';
  type.spot_init               = 'string'; 
  borne.spot_init              = {'Steady-state','Zero'};
  defaut.spot_init             = 'Zero';
  info.spot_init               = 'Status of fast ion particle distribution function for the first call to the module';

  valeur.spot_save             = 'No';
  type.spot_save               = 'string';
  borne.spot_save              = {'Yes','No'};
  defaut.spot_save             = 'No';
  info.spot_save               = 'Save context of SPOT for tests (WARNING: it fills the disk very quickly....)';

  valeur.spot_anomalous        = 'No';
  type.spot_anomalous          = 'string';
  borne.spot_anomalous         = {'Yes','No'};
  defaut.spot_anomalous        = 'No';
  info.spot_anomalous          = 'YES to account for the fast ion anomalous transport in SPOT';

  valeur.spot_ano_dcoef1       = 0.;
  type.spot_ano_dcoef1         = 'real';
  borne.spot_ano_dcoef1        = [0,10];
  defaut.spot_ano_dcoef1       = 0.;
  info.spot_ano_dcoef1         = 'Diffusion coefficient at the centre for fast ion anomalous transport in SPOT';

  valeur.spot_ano_dcoef2       = 0.;
  type.spot_ano_dcoef2         = 'real';
  borne.spot_ano_dcoef2        = [0,10];
  defaut.spot_ano_dcoef2       = 0.;
  info.spot_ano_dcoef2         = 'Diffusion coefficient at the edge for fast ion anomalous transport in SPOT';

  valeur.spot_ano_rhocut       = 0.5;
  type.spot_ano_rhocut         = 'real';
  borne.spot_ano_rhocut        = [0,1];
  defaut.spot_ano_rhocut       = 0.5;
  info.spot_ano_rhocut         = 'Radial limit to separate diffusion coefficients at the centre and at the edge';

  valeur.spot_ano_vconv        = 0.;
  type.spot_ano_vconv          = 'real';
  borne.spot_ano_vconv         = [0,10];
  defaut.spot_ano_vconv        = 0.;
  info.spot_ano_vconv          = 'Maximum convection velocity for fast ion anomalous transport in SPOT';

  valeur.rtang        = 5;
  type.rtang          = 'real';
  borne.rtang         = [0,10];
  defaut.rtang        = 0;
  info.rtang          = '(DEMO only) Tangency radius (m)';

  valeur.zshift        = 5;
  type.zshift          = 'real';
  borne.zshift         = [-10,10];
  defaut.zshift        = 0;
  info.zshift          = '(DEMO only) Vertical shift of injector (m)';

  valeur.tilt        = 0;
  type.tilt          = 'real';
  borne.tilt         = [-10,10];
  defaut.tilt        = 0;
  info.tilt          = '(DEMO only) Tilt angle (degrees)';

  valeur.save                  = 'No';
  type.save                    = 'string';
  borne.save                   = {'Yes','No'};
  defaut.save                  = 'No';
  info.save                    = 'save context for test';
  
  interface.ts  = '';
  interface.jet = '';

  sortie.valeur      = valeur;
  sortie.type        = type;
  sortie.borne       = borne;
  sortie.defaut      = defaut;
  sortie.info        = info;
  sortie.interface   = interface;
  sortie.description = 'Calculation of NBI sources with NEMO and SPOT';
  sortie.help        = ''; % HELP FILE, IF IT EXISTS
  sortie.gui         = ''; % NAME OF THE SPECIFIC GRAPHIC INTERFACE, IF IT EXISTS
  sortie.controle    = ''; % NAME OF THE CONTROL FUNCTION, IF IT EXISTS

  return
end

if isfield(param,'machine_nom')
  machine_nom = param.machine_nom;
else
  machine_nom = [];
end

%% SAVE INITAL STATE OF PARAM STRUCTURE TO RESTAURE IT AT THE END
param_save = param;

%% INITIALIZATIONS
sortie     = proto;
particules = zeros(1,gene.nbrho,gene.nbg);

param.pini_number = length(cons);

%% GO FURTHER ONLY IF THERE IS SOME NBI POWER
if sum(cons(:)) > 50e3
  
  %% CALL NEMO FOR CALCULATION OF NEUTRAL BEAM SOURCE (INJECTION, DEPOSITION, ATTENUATION)
  [output_nemo] = nemocall(equi,gene,prof,impur,compo,geo,param,cons,memoire.data);
	 %output_nemo.rout
	 %output_nemo.torque
	 %gene.x

  %% ROTATION SOURCE (COMING FROM NEMO) (KG.M-1.S-2 = N.M-2)
  sortie.wb = interp1(output_nemo.rout/output_nemo.rout(end),sum(output_nemo.torque,1),gene.x);
  
  %% SAVE THE WORKSPACE IN THE DEBUG MODE
  if strcmp(param.save,'Yes')
    [chemin,void] = fileparts(gene.rapsauve);
    [dum,filename,dum]=fileparts(gene.origine);
    sprintf('save %s/last_nemo_%s_%s_%s.mat',chemin,getenv('USER'),int2str(fix(gene.t*1000)),filename)
    eval(sprintf('save %s/last_nemo_%s_%s_%s.mat',chemin,getenv('USER'),int2str(fix(gene.t*1000)),filename))
  end

  %% CALL SPOT FOR FOKKER-PLANCK CALCULATION (FAST ION PROPAGATION AND RELAXATION)
  if strcmp(param.fokker,'SPOT')

    [output_spot] = spotcall(output_nemo,geo,equi,gene,param,memoire.data,cons,compo,prof,impur,phys);

    %% ------------------------------------------------------------------------------------------
    %% EVALUATE THE FRACTION OF POWER BELOW THE THRESHOLD THAT IS TRANSFERRED TO THE PLASMA IONS
    %% ------------------------------------------------------------------------------------------
    
    %% THRESHOLD ENERGY (EV)
    eb0 = mean(param.energie);
    
    %% CRITICAL ENERGY (EV)
    amass     = mean(param.masse);
    zcharge   = mean(param.charge);
    xhout     = output_spot.inprhoout ./ equi.rhomax;              % NORMALISED RHO USED FOR PROFILES
    inpneout  = interp1(gene.x,prof.ne,xhout);                     % ELECTRON DENSITY PROFILE
    inpteout  = interp1(gene.x,prof.te,xhout);                     % ELECTRON TEMPERATURE PROFILE
    inptiout  = interp1(gene.x,prof.ti,xhout);                     % ION TEMPERATURE PROFILE
    inpn1out  = interp1(gene.x,squeeze(impur.impur(:,:,1)),xhout); % DENSITY PROFILE OF SPECIES 1
    inpn2out  = interp1(gene.x,squeeze(impur.impur(:,:,2)),xhout); % DENSITY PROFILE OF SPECIES 2
    inpn3out  = interp1(gene.x,squeeze(impur.impur(:,:,3)),xhout); % DENSITY PROFILE OF SPECIES 3
    inpn4out  = interp1(gene.x,squeeze(impur.impur(:,:,4)),xhout); % DENSITY PROFILE OF SPECIES 4
    inpn5out  = interp1(gene.x,squeeze(impur.impur(:,:,5)),xhout); % DENSITY PROFILE OF SPECIES 5
    logl_prof = 24-log((inpneout*1e-6).^0.5./inpteout);            % (-)
    tslow_ana = 6.27e8*amass*(inpteout.^1.5)./(inpneout*1e-6*zcharge^2.*logl_prof);  % SLOWING-DOWN TIME
    inpzi = compo.z; % Z OF ALL SPECIES
    inpai = compo.a; % A OF ALL SPECIES
    gzbar = (inpn1out*inpzi(1).^2/inpai(1)+inpn2out*inpzi(2).^2/ ... % ZBAR OF THE PLASMA
	     inpai(2)+inpn3out*inpzi(3).^2/inpai(3)+inpn4out* ...
	     inpzi(4).^2/inpai(4)+inpn5out*inpzi(5).^2/inpai(5))./inpneout;
    estixal  = 14.8 .* inpteout*1.e-3 .* amass .*gzbar.^(2./3.); % CRITICAL ENERGY
    vcri_ana = sqrt(2000*1.6e-19*estixal/(amass*1.6726e-27));    % CRITICAL VELOCITY

    %% CORRECTION OF VCRI DUE TO NON-CONSTANT COULOMB LOGARITHM
    logle = 15.2-0.5*log(inpneout/1e20)+log(inpteout)+log(zcharge);
    logli = 17.3-0.5*log(inpneout/1e20)+1.5*log(inptiout)+log(zcharge);
    vcri_ana = vcri_ana.*(logli./logle).^(1./3.);
    ecri = 0.5*amass*1.67e-27.*vcri_ana.^2./1.6e-19;

    %% FRACTION OF POWER TO THERMAL IONS
    if eb0 > 0.
      frac = phi_th_frac(eb0./ecri);
    else
      frac = 1.;
    end
    elpower1dp  = -output_spot.elpower1d'-(1-frac').*output_spot.thpower1d';
    ionpower1dp = -output_spot.ionpower1d'-frac'.*output_spot.thpower1d';    
    
    %% OUTPUT MANAGEMENT
    rhovec        = output_spot.inprhoout / max(output_spot.inprhoout);
    sortie.el     = zbornes(interp1(rhovec,elpower1dp,gene.x),0,inf,0);
    sortie.ion    = zbornes(interp1(rhovec,ionpower1dp,gene.x),0,inf,0);
    sortie.j      = zbornes(interp1(rhovec,output_spot.icurpara1d(:,1'),gene.x),0,inf,0);
    %sortie.psupra = zbornes(interp1(rhovec,output_spot.pressperp1d,gene.x),0,inf,0);
    %% HERE I USE A WRONG DEFINITION ON PURPOSE TO CORRECT THE EQUATION
    %% IN COEUR/ZPROFDIVERS FOR THE CALCULATION OF PTOT: THE PARALLEL
    %% PRESSURE HAS TO BE USED, INSTEAD OF THE ISOTROPIC (I.E
    %% PERPENDICULAR) COMPONENT, FOR PSUPRA
    sortie.psupra = zbornes(interp1(rhovec,output_spot.presspara1d',gene.x),0,inf,0);

    sortie.paniso = 0.5*(zbornes(interp1(rhovec,output_spot.presspara1d', ...
		    gene.x),0,inf,0)-zbornes(interp1(rhovec,output_spot.pressperp1d',gene.x),0,inf,0));
    sortie.ne     = zbornes(interp1(rhovec,output_spot.initrate1d',gene.x),0,inf,0);
    sortie.nfast  = zbornes(interp1(rhovec,output_spot.nalpha1d',gene.x),0,inf,0);
    kthermion     = zbornes(interp1(rhovec,output_spot.kthermalpha1d',gene.x),0,inf,0);

    %% ROTATION SOURCE (COMING FROM SPOT) (KG.M-1.S-2 = N.M-2)
    sortie.w = interp1(rhovec,output_spot.torque1d',gene.x);

    %% ERR = SHINETHROUG LOSSES (W)
    sortie.err = sum(output_nemo.pshine_tot);
    sortie.neutron.dd = zeros(size(sortie.j));

    %% FILL VARIABLES FOR SOURCE OF MATTER WITH THERMALIZED IONS
    for iinj=1:param.pini_number
      zcharge = param.charge(iinj);
      amass   = param.masse(iinj);
      index   = find((compo.z == zcharge) &(compo.a == amass));
      particules(1,:,index) = particules(1,:,index) + kthermion/param.pini_number;
    end

    %% RUN THE FAST FOKKER-PLANCK CALCULATION
  elseif strcmp(param.fokker,'FAST')
  
    disp('----------------------------')
    disp('RUN FAST FOKKER CALCULATION')
    disp('----------------------------')

    nidn               = 0 .* impur.impur;
    taus               = 0 .* gene.x;

    param.orbit = 'No';

    for iinj=1:param.pini_number
      if cons(iinj)~=0
	
	% charge et masse
	ag = param.masse(iinj);
	zg = param.charge(iinj);
	
	for iene=1:3
	  
	  % pitch angle
	  pitch_int  = squeeze(output_nemo.pitchvec(iene,:,iinj));
		%output_nemo.pitchvec
		%output_nemo.rout
		%pitch_int
		%gene.x
	  pitch_in   = zbornes(interp1(output_nemo.rout/output_nemo.rout(end),pitch_int,gene.x),0,inf,0);
	  
	  if iene==1
	    partic_ene_frac = param.fraction1;
	  elseif iene==2
	    partic_ene_frac = param.fraction2;
	  else
	    partic_ene_frac = param.fraction3;
	  end
	  
	  einj     = output_nemo.input.injected_energy(iene,iinj);
	  pdep_int = squeeze(output_nemo.hvec_pini(iinj,iene,:)*einj.*phys.e) ...
	      .*partic_ene_frac(iinj);
	  pdep_in  = zbornes(interp1(output_nemo.rout/output_nemo.rout(end),pdep_int,gene.x),0,inf,0);
	  
	  [pdep,ploss,pdep_ion,pdep_el,jnbi,wnbi,snbi,taus_nbi,tauseff,psupra,paniso,pin]= ...
	      pseudo_FP(gene.x,pdep_in,abs(pitch_in),ag,zg,einj,param.orbit, ...
	      equi.a(end),equi.F(end) .* equi.ri(end), ...
	      prof.ne,prof.te,prof.ni,prof.zeff,impur.impur,equi.vpr.*equi.rhomax, ...
	      equi.raxe,equi.q,equi.ftrap,compo,phys);

	  sortie.el  = sortie.el  + pdep_el;
	  sortie.ion = sortie.ion + pdep_ion;
	  sortie.j   = sortie.j   + jnbi .* sign(pitch_in);
	  % MS 12.01.09: rotation calculated by NEMO (filled later)
	  % sortie.w   = sortie.w   + wnbi;
	  sortie.ne  = sortie.ne  + pdep_in ./ einj ./ phys.e;
	  sortie.psupra  = sortie.psupra  + psupra;
	  sortie.paniso  = sortie.paniso  + paniso;

	  taus = taus + tauseff .* pin;
			
	  indg =  find((ag == compo.a) & (zg == compo.z),1,'first');
	  particules(1,:,indg) = particules(1,:,indg) + snbi;     
	  sortie.err = sortie.err + ploss;

	end % loop over energy fractions iene
      end   % test cons~=0
    end     % loop over pinis iinj

    %% ROTATION SOURCE (COMING FROM NEMO) (KG.M-1.S-2 = N.M-2)
    sortie.w = interp1(output_nemo.rout/output_nemo.rout(end),sum(output_nemo.torque,1),gene.x);
    
    %% start cooking
    if 2==2
    pitchmoyt =  squeeze(mean(mean(output_nemo.pitchvec,1),3));
    pitchmoy   = zbornes(interp1(output_nemo.rout/output_nemo.rout(end),pitchmoyt,gene.x),0,inf,0);

    % normalisation de taus
    taus = taus ./ max(1,sum(abs(cons)));

    % gestion du temps et des memoires
    if ~isfield(memoire,'data')
      memoire.data = [];
    end
    if isempty(memoire.data) || ~isfield(memoire.data,'sinbadfast')
      if strcmp(param.spot_init,'Zero')
	param.debut = 'Before';
      elseif strcmp(param.spot_init,'Steady-state')
	param.debut = 'During';
      end
      switch param.debut
       case 'Before'
	% initialisation des memoires
	memoire.t = gene.t;
	memoire.data.input.cons(1,:)  = 0 .* cons(:);
	memoire.data.input.temps = gene.t - gene.dt;
	memoire.data.input.dt    = gene.dt;
	memoire.data.input.taus  = taus;
	memoire.data.output      = sortie;
	memoire.data.particules  = particules;
	sortie = proto;	
	particules  = zeros(size(particules));
       otherwise
	% initialisation des memoires
	memoire.t = gene.t;
	memoire.data.input.cons(1,:)  = cons(:);
	memoire.data.input.temps = gene.t - gene.dt;
	memoire.data.input.dt    = gene.dt;
	memoire.data.input.taus  = taus;
	memoire.data.output      = sortie;	
	memoire.data.particules  = particules;
      end
      memoire.data.sinbadfast = 1;
    else
      % remplissage des memoires
      memoire.t = gene.t;
      memoire.data.input.cons(end+1,:)    = cons(:);
      memoire.data.input.temps(end+1,1)   = gene.t;
      memoire.data.input.dt(end+1,1)      = gene.dt;
      memoire.data.input.taus(end+1,:)    = taus;
      memoire.data.particules(end+1,:,:)  = particules;
      noms = fieldnames(sortie);
      for k =1:length(noms)
	if ~isstruct(sortie.(noms{k}))
	  memoire.data.output.(noms{k})(end+1,:) = sortie.(noms{k});
	end
      end

      % diffusion anormale 
      switch param.spot_anomalous
       case 'Yes'
	memoire.data.input.Dan = param.spot_ano_dcoef1 .* (param.spot_ano_rhocut > gene.x) + ...
	    param.spot_ano_dcoef2 .* (param.spot_ano_rhocut <= gene.x);

	indc    = find(param.spot_ano_rhocut<= gene.x,1);
	fxp     = linspace(0,1,indc);	
	fxm     = linspace(1,0,length(gene.x) - indc + 1);							
	fx      = cat(2,fxp,fxm(2:end)); 
	memoire.data.input.Van = - param.spot_ano_vconv .*  fx;
	
	% calcul de la diffusion 
	difu = memoire.data.input.Dan;
	conv = memoire.data.input.Van;
		
	% diffusion des electron
	rep = zdiffusedensite(memoire.data.output.el(end-1,:) .* memoire.data.input.taus(end-1) ./ 2,0,0 ...
            ,equi.rhomax,equi.drhomaxdt,equi.grho2,equi.vpr,equi.dvprdt,difu,difu,conv,conv ...
            ,memoire.data.output.el(end-1,:),memoire.data.output.el(end,:),gene.x,memoire.data.input.dt(end) ...
            ,max(1e-6,memoire.data.input.taus(end) ./ 20),0.5);		    
	memoire.data.output.el(end,:)   = rep(end,:) ./ taus(end,:) .* 2;
   
	% diffusion des ions
	rep = zdiffusedensite(memoire.data.output.ion(end-1,:) .* memoire.data.input.taus(end-1) ./ 2,0,0 ...
            ,equi.rhomax,equi.drhomaxdt,equi.grho2,equi.vpr,equi.dvprdt,difu,difu,conv,conv ...
            ,memoire.data.output.ion(end-1,:),memoire.data.output.ion(end,:),gene.x,memoire.data.input.dt(end) ...
            ,max(1e-6,memoire.data.input.taus(end) ./ 20),0);		    
	memoire.data.output.ion(end,:)   = rep(end,:) ./ taus(end,:) .* 2;
   
	% diffusion du courant
	rep = zdiffusedensite(memoire.data.output.j(end-1,:) .* memoire.data.input.taus(end-1),0,0 ...
            ,equi.rhomax,equi.drhomaxdt,equi.grho2,equi.vpr,equi.dvprdt,difu,difu,conv,conv ...
            ,memoire.data.output.j(end-1,:),memoire.data.output.j(end,:),gene.x,memoire.data.input.dt(end) ...
            ,max(1e-6,memoire.data.input.taus(end) ./ 10),0);		    
	memoire.data.output.j(end,:)   = rep(end,:) ./ taus(end,:);
   
	% diffusion du moment
	rep = zdiffusedensite(memoire.data.output.w(end-1,:) .* memoire.data.input.taus(end-1),0,0 ...
            ,equi.rhomax,equi.drhomaxdt,equi.grho2,equi.vpr,equi.dvprdt,difu,difu,conv,conv ...
            ,memoire.data.output.w(end-1,:),memoire.data.output.w(end,:),gene.x,memoire.data.input.dt(end) ...
            ,max(1e-6,memoire.data.input.taus(end) ./ 10),0);		    

	frac = equi.ftrap .* pitchmoy;
	memoire.data.output.w(end,:) = rep(end,:) ./ taus(end,:).* (1 - frac) + frac .*  memoire.data.output.w(end,:);

	for k=1:size(particules,3)
	  snbi = squeeze(memoire.data.particules(:,:,k));
	  if any(snbi(:))
	    rep = zdiffusedensite(memoire.data.output.ion(end-1,:) .* memoire.data.input.taus(end-1) ./ 2,0,0 ...
               ,equi.rhomax,equi.drhomaxdt,equi.grho2,equi.vpr,equi.dvprdt,difu,difu,conv,conv ...
               ,memoire.data.output.ion(end-1,:),memoire.data.output.ion(end,:),gene.x,memoire.data.input.dt(end) ...
               ,max(1e-6,memoire.data.input.taus(end) ./ 20),0);
	    memoire.data.particules(end,:,k)   = rep(end,:) ./ taus(end,:) .* 2;
	  end
	end

	noms = fieldnames(sortie);
	for k =1:length(noms)
	  if ~isstruct(sortie.(noms{k}))
	    sortie.(noms{k}) = memoire.data.output.(noms{k})(end,:);
	  end
	end

	particules = memoire.data.particules(end,:,:);

      end % (test anomalous diffusion)

      % champs utiles
      el    = memoire.data.output.el;
      ion   = memoire.data.output.ion;
      w     = memoire.data.output.w;
      j     = memoire.data.output.j;
      taus  = max(1e-6,memoire.data.input.taus);
      temps = memoire.data.input.temps;
      
      % calcul des ode
      % P -> taus/2
      % J -> taus
      % W passantes -> taus
      % W piegees -> immediat
      [trep,rep] = fast_ode(temps,el,taus ./ 2 ,0 .* gene.x);
      sortie.el  = rep(end,:) ./ taus(end,:) .* 2;
      pfast      = max(1,3 ./ 2 .* sortie.psupra + sortie.paniso);
      a2         = 1 - sortie.psupra ./ pfast;
      sortie.psupra = rep(end,:) .* (1 - a2);
      sortie.paniso = rep(end,:) .* a2 - 0.5 .* sortie.psupra;
      sortie.psupra = rep(end,:);
      [trep,rep] = fast_ode(temps,ion,taus ./ 2 ,0 .* gene.x);
      sortie.ion  = rep(end,:) ./ taus(end,:) .* 2;
      [trep,rep] = fast_ode(temps,j,taus ,0 .* gene.x);
      sortie.j   = rep(end,:) ./ taus(end,:);
      
      frac       = equi.ftrap .* pitchmoy;
      [trep,rep] = fast_ode(temps,w,taus ,0 .* gene.x);
      sortie.w   = rep(end,:) ./ taus(end,:) .* (1 - frac) + w(end,:) .* frac;
      
      % normalisation sur la source d'ions de particules
      for k=1:size(particules,3)
	snbi = squeeze(memoire.data.particules(:,:,k));
	[trep,rep] = fast_ode(temps,snbi,taus ./ 2 ,0 .* gene.x);
	particules(:,:,k)  = rep(end,:) ./ taus(end,:) .* 2;
      end

    end % (test memoire structure filled)
    end % test 2==3
    %% end cooking

    %% SAVE THE WORKSPACE IN THE DEBUG MODE
    if strcmp(param.save,'Yes')
      [chemin,void] = fileparts(gene.rapsauve);
      [dum,filename,dum]=fileparts(gene.origine);
      sprintf('save %s/last_fastfok_%s_%s_%s.mat',chemin,getenv('USER'),int2str(fix(gene.t*1000)),filename)
      eval(sprintf('save %s/last_fastfok_%s_%s_%s.mat',chemin,getenv('USER'),int2str(fix(gene.t*1000)),filename))
    end

  %% RUN THE RISK FOKKER-PLANCK CALCULATION
  elseif strcmp(param.fokker,'RISK')
  
    disp('-----------------------------')
    disp('RUN RISK FOKKER CALCULATION')
    disp('-----------------------------')

    [output_risk,cons] = riskcall(output_nemo,geo,equi,gene,param,memoire.data,cons,compo,prof,impur,phys);
    
    %% ------------------------------------------------------------------------------------------
    %% EVALUATE THE FRACTION OF POWER BELOW THE THRESHOLD THAT IS TRANSFERRED TO THE PLASMA IONS
    %% ------------------------------------------------------------------------------------------
    
    %% THRESHOLD ENERGY (EV)
    %% MS205274 PROVISOIRE LE TEMPS QUE CETTE QUANTITE SOIT EXPORTEE DEPUIS RISK
    eb0 = prof.ti(1);
    
    %% CRITICAL ENERGY (EV)
    amass     = mean(param.masse);
    zcharge   = mean(param.charge);
    xhout     = output_risk.rhonorm;                              % NORMALISED RHO USED FOR PROFILES
    inpneout  = interp1(gene.x,prof.ne,xhout);                     % ELECTRON DENSITY PROFILE
    inpteout  = interp1(gene.x,prof.te,xhout);                     % ELECTRON TEMPERATURE PROFILE
    inptiout  = interp1(gene.x,prof.ti,xhout);                     % ION TEMPERATURE PROFILE
    inpn1out  = interp1(gene.x,squeeze(impur.impur(:,:,1)),xhout); % DENSITY PROFILE OF SPECIES 1
    inpn2out  = interp1(gene.x,squeeze(impur.impur(:,:,2)),xhout); % DENSITY PROFILE OF SPECIES 2
    inpn3out  = interp1(gene.x,squeeze(impur.impur(:,:,3)),xhout); % DENSITY PROFILE OF SPECIES 3
    inpn4out  = interp1(gene.x,squeeze(impur.impur(:,:,4)),xhout); % DENSITY PROFILE OF SPECIES 4
    inpn5out  = interp1(gene.x,squeeze(impur.impur(:,:,5)),xhout); % DENSITY PROFILE OF SPECIES 5
    logl_prof = 24-log((inpneout*1e-6).^0.5./inpteout);            % (-)
    tslow_ana = 6.27e8*amass*(inpteout.^1.5)./(inpneout*1e-6*zcharge^2.*logl_prof);  % SLOWING-DOWN TIME
    inpzi = compo.z; % Z OF ALL SPECIES
    inpai = compo.a; % A OF ALL SPECIES
    gzbar = (inpn1out*inpzi(1).^2/inpai(1)+inpn2out*inpzi(2).^2/ ... % ZBAR OF THE PLASMA
	     inpai(2)+inpn3out*inpzi(3).^2/inpai(3)+inpn4out* ...
	     inpzi(4).^2/inpai(4)+inpn5out*inpzi(5).^2/inpai(5))./inpneout;
    estixal  = 14.8 .* inpteout*1.e-3 .* amass .*gzbar.^(2./3.); % CRITICAL ENERGY
    vcri_ana = sqrt(2000*1.6e-19*estixal/(amass*1.6726e-27));    % CRITICAL VELOCITY

    %% CORRECTION OF VCRI DUE TO NON-CONSTANT COULOMB LOGARITHM
    logle = 15.2-0.5*log(inpneout/1e20)+log(inpteout)+log(zcharge);
    logli = 17.3-0.5*log(inpneout/1e20)+1.5*log(inptiout)+log(zcharge);
    vcri_ana = vcri_ana.*(logli./logle).^(1./3.);
    ecri = 0.5*amass*1.67e-27.*vcri_ana.^2./1.6e-19;

    %% FRACTION OF POWER TO THERMAL IONS
    %frac = phi_th_frac(eb0./ecri);
    %elpower1dp  = output_risk.elpower1d+(1-frac').*output_risk.thpower1d;
    %ionpower1dp = output_risk.ionpower1d+frac'.*output_risk.thpower1d;    
    elpower1dp  = output_risk.elpower1d;
    ionpower1dp = output_risk.ionpower1d+output_risk.thpower1d;
    
    %% OUTPUT MANAGEMENT
    sortie.el     = zbornes(interp1(xhout,elpower1dp,gene.x),0,inf,0);
    sortie.ion    = zbornes(interp1(xhout,ionpower1dp,gene.x),0,inf,0);
    sortie.j      = zbornes(interp1(xhout,output_risk.jnbcd,gene.x),0,inf,0);

    %sortie.psupra = zbornes(interp1(xhout,output_risk.pressperp1d,gene.x),0,inf,0);
    %% HERE I USE A WRONG DEFINITION ON PURPOSE TO CORRECT THE EQUATION
    %% IN COEUR/ZPROFDIVERS FOR THE CALCULATION OF PTOT: THE PARALLEL
    %% PRESSURE HAS TO BE USED, INSTEAD OF THE ISOTROPIC (I.E
    %% PERPENDICULAR) COMPONENT, FOR PSUPRA
    sortie.psupra = zbornes(interp1(xhout,output_risk.presspara1d,gene.x),0,inf,0);
    sortie.paniso = 0.5*(zbornes(interp1(xhout,output_risk.presspara1d, ...
		    gene.x),0,inf,0)-zbornes(interp1(xhout,output_risk.pressperp1d,gene.x),0,inf,0));
    %% MS205274: IT MAY HAPPEN THAT PANISO IS SLIGHTLY <0, PUT EQUAL TO ZERO THEN
    sortie.paniso(find(sortie.paniso<0)) = 0;
    
    sortie.ne     = zbornes(interp1(output_nemo.rout./output_nemo.rout(end),sum(output_nemo.hvec,1),gene.x),0,inf,0);
    kthermion     = zbornes(interp1(xhout,output_risk.thnumber1d,gene.x),0,inf,0);

    %% ROTATION SOURCE (COMING FROM RISK) (KG.M-1.S-2 = N.M-2)
    sortie.w = interp1(xhout,output_risk.torque1d,gene.x);

    %% ERR = SHINETHROUGH LOSSES (W)
    sortie.err = sum(output_nemo.pshine_tot);
    sortie.neutron.dd = zeros(size(sortie.j));

    %% FILL VARIABLES FOR SOURCE OF MATTER WITH THERMALIZED IONS
    for iinj=1:param.pini_number
      zcharge = param.charge(iinj);
      amass   = param.masse(iinj);
      index   = find((compo.z == zcharge) &(compo.a == amass));
      particules(1,:,index) = particules(1,:,index) + kthermion/param.pini_number;
    end
    
    %% SAVE THE WORKSPACE IN THE DEBUG MODE
    if strcmp(param.save,'Yes')
      [chemin,void] = fileparts(gene.rapsauve);
      [dum,filename,dum]=fileparts(gene.origine);
      sprintf('save %s/last_risk_%s_%s_%s.mat',chemin,getenv('USER'),int2str(fix(gene.t*1000)),filename)
      eval(sprintf('save %s/last_risk_%s_%s_%s.mat',chemin,getenv('USER'),int2str(fix(gene.t*1000)),filename))
    end

  else
    disp('-------------------------------------------------------------')
    disp('NO CALCULATION IMPLEMENTED FOR THIS FOKKER-PLANCK CALCULATION')
    disp('-------------------------------------------------------------')
    keyboard
    return
  end % test SPOT or FAST FOKKER or RISK

else
  disp('-----------------------------------------------------------')
  disp('NO SIGNIFICANT NBI POWER YET => DO NOT BOTHER CAPTAIN NEMO!')
  disp('-----------------------------------------------------------')
  sortie.el     = 0.*ones(size(gene.x));
  sortie.ion    = 0.*ones(size(gene.x));
  sortie.j      = 0.*ones(size(gene.x));
  sortie.psupra = 0.*ones(size(gene.x));
  sortie.paniso = 0.*ones(size(gene.x));
  sortie.ne     = 0.*ones(size(gene.x));
  sortie.w      = 0.*ones(size(gene.x));
  kthermion     = 0.*ones(size(gene.x));
  sortie.err    = 0.;
  sortie.neutron.dd = 0.*ones(size(gene.x));

end % test cons~=0

%% SAVE EVERYTHING INTO THE MEMORY FOR NEMO
if exist('output_nemo')
  memoire.data.output_nemo = output_nemo;
end

%% SAVE EVERYTHING INTO THE MEMORY FOR SPOT
if strcmp(param.fokker,'SPOT')
  memoire.t               = gene.t;
  memoire.data.sortie     = sortie;
  memoire.data.particules = particules;
  memoire.data.cons       = cons;
  if exist('output_spot')
    memoire.data.memoryspot = output_spot;
  end
end

%% SAVE EVERYTHING INTO THE MEMORY FOR RISK
if strcmp(param.fokker,'RISK')
  memoire.t               = gene.t;
  memoire.data.sortie     = sortie;
  memoire.data.particules = particules;
  memoire.data.cons       = cons;
  if exist('output_risk')
    memoire.data.output_risk = output_risk;
  end
end

%% SAVE DISK SPACE FOR THE TEST FILE STORED IN CVS
if strcmp(param.spot_init,'Zero_testmode')
  memoire.data.memoryspot = [];
  memoire.t = 0.;
  memoire.data.input=[];
end

%% RESTAURE INITIAL STATE OF PARAM STRUCTURE
param=param_save;


