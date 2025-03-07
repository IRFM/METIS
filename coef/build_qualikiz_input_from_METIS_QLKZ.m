% Builds QuaLiKiz standalone from QuaLiKiz in METIS
% Need to change all path for generic use
function [] = build_qualikiz_input_from_METIS_QLKZ(NNin,runname_qlk,metis_qlkz_path,qlkparams,Sn_fraction)

runname=runname_qlk;
currentpath = pwd;

path_to_METIS_QLKZ = metis_qlkz_path;

% Creates runs directory if it doesn't already exist 
tmp_dir = tempdir;
runpath=[tmp_dir runname];

if isdir('runpath') == 0
  eval(['!mkdir ',runpath]);
end

eval(['!cp ',path_to_METIS_QLKZ,'/parameters_metis_qlkz.m ',runpath,'/parameters.m']);
eval(['!cp ',path_to_METIS_QLKZ,'/subjob ',runpath,'/subjob']);

% Use latest compiled QuaLiKiz version
eval(['!cp ',path_to_METIS_QLKZ,'/QLKZ/QuaLiKiz/QuaLiKiz ',runpath,'/','QuaLiKiz']);

% Load standard parameters
parameters_metis_qlkz;

phys_meth = qlkparams.phys_meth;
coll_flag = qlkparams.coll_flag;
rot_flag = qlkparams.rot_flag;
numsols   = qlkparams.numsols;
relacc1   = qlkparams.relacc1;
relacc2   = qlkparams.relacc2;
ETGmult = qlkparams.ETG_mult;
collmult = qlkparams.coll_mult;
separateflux = qlkparams.separateflux;

scann = length(NNin.ne);

clear Tix;
clear Zi;
clear Ai;
clear ninorm;
clear Ati;
clear Ani;
clear ion_type;
clear anis;
clear danisdr;

x=NNin.rova;
rho=x;

% Small and big radii at LCFS
Rmin = NNin.amin;
R0 = NNin.R0;

qx = NNin.q;
Ro = NNin.Raxe;
Bo = abs(NNin.B0);

smag = NNin.shearr;
Tex = NNin.te;
Nex = NNin.ne;
Ate = NNin.rlte;
Ane = NNin.rlne;
anise = 1.*ones(scann,1); 
danisedr = 0.*ones(scann,1);

Tix(:,1) = NNin.te.*NNin.tite;
% Hard limit in QuaLiKiz standalone on Temperature
Tex(find(Tex<=0.0015)) = 0.002;

Tix(find(Tix(:,1)<=0.0015)) = 0.002;
Tix(:,2) = Tix(:,1);
Tix(:,3) = Tix(:,1);
Tix(:,4) = Tix(:,1);
Tix(:,5) = Tix(:,1);

Ati(:,1) = NNin.rlti;
Ati(:,2) = Ati(:,1);
Ati(:,3) = Ati(:,1);
Ati(:,4) = Ati(:,1);
Ati(:,5) = Ati(:,1);

ninorm(:,1) = NNin.nip./Nex; % Temporary, changed with quasineutrality
ninorm(:,2) = NNin.nz1p./Nex;
ninorm(:,3) = NNin.nz2p./Nex;
ninorm(:,4) = NNin.nwp./Nex;
ninorm(:,5) = NNin.nminp./Nex;

Zi(:,1) = 1.*ones(scann,1);
Zi(:,2)=NNin.z1;
Zi(:,3)=NNin.z2;
Zi(:,4)=NNin.z3;
Zi(:,5)=1.*ones(scann,1);

% Temporary need to get it from METIS
Ai(:,1) = 2.*ones(scann,1);
Ai(:,2)=2*Zi(:,2);
Ai(:,3)=2*Zi(:,3);
Ai(:,4)= round((1 - Sn_fraction) .* 183.84  +  Sn_fraction .* 118.71) .*ones(scann,1);
Ai(:,5)=1.*ones(scann,1);

Ani(:,1) = NNin.rlne;
Ani(:,2) = NNin.rlnz1;
Ani(:,3) = NNin.rlnz2;
Ani(:,4) = NNin.rlnz3;
Ani(:,1) = (Ane'-sum(Zi(:,2:4).*Ani(:,2:4).*ninorm(:,2:4),2))./(1+NNin.cmin');
Ani(:,5) = Ani(:,1);

% Normalisation with sound speed Tref = Te
Machtor = 0.*Ane;
Machpar = Machtor;
Autor = 0.*Ane;
Aupar = Autor;
gammaE = 0.*Ane;

mu0 = 4*pi*1e-7;
alphai = sum(ninorm.*Nex'.*Tix.*(Ati+Ani),2);
%alphax = 0.*ones(scann,1);
alphax = qx.^2.*(Nex.*Tex.*(Ate+Ane)+alphai')./(Bo.^2./(2*mu0))*1e19*1.6*1e-16;

% Could change it for fast minority
anis(1:scann,1:5)=1.*ones(scann,size(Ai,2)); %Tperp/Tpar at LFS. Tix defined as Tpar
danisdr(1:scann,1:5)=0.*ones(scann,size(Ai,2)); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar

ion_type(1:scann,1:5)=1.*ones(scann,size(Ai,2)); %Tperp/Tpar at LFS. Tix defined as Tpar

%%%%%%%%%%%%%%%%%% Write binary input files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(sprintf('%s',['!mkdir ' runpath '/input']))
eval(sprintf('%s',['!mkdir ' runpath '/output']))
eval(sprintf('%s',['!mkdir ' runpath '/output/primitive']))
eval(sprintf('%s',['!mkdir ' runpath '/debug']))

%qualikiz_names;
kc  =  1;
name{kc} = 'dimx'; 	p{kc}  =  length(x);                  kc = kc+1;%p{1} (-) Number of radial or scan points
name{kc} = 'dimn';	p{kc}  =  length(kthetarhos);         kc = kc+1;%p{2} (-) Number of wavenumbers
name{kc} = 'nions';	p{kc}  =  nions;                      kc = kc+1;%p{3} (-) Number of ions in system

%Flag input and metadata
name{kc} = 'phys_meth';		p{kc}  =  phys_meth;	            kc = kc+1;%p{4} (-) Flag for additional calculation (default 0.0)
name{kc} = 'coll_flag';		p{kc}  =  coll_flag;	            kc = kc+1;%p{5} (-) Flag for collisionality (default 0.0)
name{kc} = 'write_primi';	p{kc}  =  write_primi;	            kc = kc+1;% Flag for writing primitive outputs (default 1)
name{kc} = 'rot_flag';		p{kc}  =  rot_flag;	            kc = kc+1;%p{6} (-) Flag  for rotation (default 0.0)
name{kc} = 'verbose';		p{kc}  =  verbose;	            kc = kc+1;%p{7} (-) Flag  for setting level of output verbosity
name{kc} = 'numsols';		p{kc}  =  numsols;	            kc = kc+1;%p{8} (-) Number of solutions requested
name{kc} = 'relacc1';		p{kc}  =  relacc1;	 	    kc = kc+1;%p{9} (-) Relative accuracy in 1D integrals
name{kc} = 'relacc2';		p{kc}  =  relacc2;	 	    kc = kc+1;%p{10} (-) Relative accuracy in 2D integrals
name{kc} = 'absacc1';		p{kc}  =  absacc1;                    kc = kc+1;%p{15}(-) collisionality multiplier (for testing)
name{kc} = 'absacc2';		p{kc}  =  absacc2;                    kc = kc+1;%p{15}(-) collisionality multiplier (for testing)
name{kc} = 'maxruns';		p{kc}  =  maxruns; 	            kc = kc+1;%p{11} (-) Number of runs jumping directly to Newton between contour checks
name{kc} = 'maxpts';		p{kc}  =  maxpts;                     kc = kc+1;%p{12}(-) Number of integrand evaluations done in 2D integral
name{kc} = 'timeout';		p{kc}  =  timeout;                    kc = kc+1;%p{13}(-) Upper time limit (s) beyond which solutions are not sought after at a given wavenumber and radius
name{kc} = 'ETGmult';		p{kc}  =  ETGmult;                    kc = kc+1;%p{14}(-) ETG multiplier (for testing)
name{kc} = 'collmult';		p{kc}  =  collmult;                   kc = kc+1;%p{15}(-) collisionality multiplier (for testing)

%Geometry input
name{kc} = 'kthetarhos';	p{kc}  =  kthetarhos;                 kc = kc+1;%p{16} (-) Wave spectrum input: Vector (dimn)
name{kc} = 'x';			p{kc}  =  x;                          kc = kc+1;%p{17} (-) radial normalised coordinate (midplane average)
name{kc} = 'rho';		p{kc}  =  rho;                        kc = kc+1;%p{18} (-) normalized toroidal flux coordinate
name{kc} = 'rhomin';		p{kc}  =  rhomin;                     kc = kc+1;%p{18} (-) normalized toroidal flux coordinate
name{kc} = 'rhomax';		p{kc}  =  rhomax;                     kc = kc+1;%p{18} (-) normalized toroidal flux coordinate
name{kc} = 'Ro';		p{kc}  =  Ro;		            kc = kc+1;%p{19} (m) Major radius. Radial profile due to Shafranov shift
name{kc} = 'Rmin';		p{kc}  =  Rmin;	                    kc = kc+1;%p{20} (m) Geometric minor radius. Assumed to be a midplane average at LCFS. Currently a profile but should probably be shifted to a scalar
name{kc} = 'Bo';		p{kc}  =  Bo;		            kc = kc+1;%p{21} (T) Likely not very rigorous to use this sqrt(<B��>) for calculating the Larmor radius % quite close to <Bphi> in practice however 
name{kc} = 'R0';		p{kc}  =  R0;		            kc = kc+1;%p{22} (m) Geometric major radius used for normalizations
name{kc} = 'q';			p{kc}  =  qx;               	    kc = kc+1;%p{23} (-) Vector (radial grid x(aa))
name{kc} = 'smag';		p{kc}  =  smag;            	    kc = kc+1;%p{24} (-) Vector (radial grid x(aa))  q is a flux surface quantity --> makes sense to consider s  =  rho/q dq/drho
name{kc} = 'alpha';		p{kc}  =  alphax;            	    kc = kc+1;%p{25} (-) Vector (radial grid x(aa)) 

%Rotation input
name{kc} = 'Machtor';		p{kc}  =  Machtor;            	    kc = kc+1;%p{26} (-) Vector (radial grid x(aa)) 
name{kc} = 'Autor';		p{kc}  =  Autor;            	    kc = kc+1;%p{27} (-) Vector (radial grid x(aa)) 
name{kc} = 'Machpar';		p{kc}  =  Machpar;            	    kc = kc+1;%p{28} (-) Vector (radial grid x(aa)) 
name{kc} = 'Aupar';		p{kc}  =  Aupar;            	    kc = kc+1;%p{29} (-) Vector (radial grid x(aa)) 
name{kc} = 'gammaE';		p{kc}  =  gammaE;            	    kc = kc+1;%p{30} (-) Vector (radial grid x(aa)) 

%Electron input
name{kc} = 'Te';		p{kc}  =  Tex;       		    kc = kc+1;%p{31} (keV) Vector (radial grid x(aa))
name{kc} = 'ne';		p{kc}  =  Nex;    		    kc = kc+1;%p{32} (10^19 m^-3) Vector (radial grid x(aa))
name{kc} = 'Ate';		p{kc}  =  Ate;                        kc = kc+1;%p{33} (-) Vector (radial grid x(aa))
name{kc} = 'Ane';		p{kc}  =  Ane;                	    kc = kc+1;%p{34} (-) Vector (radial grid x(aa))
name{kc} = 'typee';		p{kc}  =  el_type;                    kc = kc+1;%p{35} Kinetic or adiabatic
name{kc} = 'anise';		p{kc}  =  anise;                      kc = kc+1;%p{36}  Tperp/Tpar at LFS
name{kc} = 'danisdre';		p{kc}  =  danisedr;                   kc = kc+1;%p{37}  d/dr(Tperp/Tpar) at LFS

%Ion inputs (can be for multiple species)
name{kc} = 'Ai';		p{kc}  =  Ai;	            kc = kc+1;%p{38} (-) Ion mass
name{kc} = 'Zi';		p{kc}  =  Zi;     	    kc = kc+1;%p{39} (-) Ion charge
name{kc} = 'Ti';		p{kc}  =  Tix;                kc = kc+1;%p{40} (keV) Vector (radial grid x(aa))
name{kc} = 'normni';		p{kc}  =  ninorm;             kc = kc+1;%p{41} ni/ne Vector (radial grid x(aa))
name{kc} = 'Ati';		p{kc}  =  Ati;      	    kc = kc+1;%p{42}  (-) Vector (radial grid x(aa))
name{kc} = 'Ani';		p{kc}  =  Ani;                kc = kc+1;%p{43}  (-) Vector (radial grid x(aa))  check calculation w.r.t. Qualikiz electroneutrality assumption
name{kc} = 'typei'; 		p{kc}  =  ion_type;           kc = kc+1;%p{44}  Kinetic, adiabatic, tracer
name{kc} = 'anisi';		p{kc}  =  anis;               kc = kc+1;%p{45}  Tperp/Tpar at LFS
name{kc} = 'danisdri'; 		p{kc}  =  danisdr;            kc = kc+1;%p{46}  d/dr(Tperp/Tpar) at LFS
name{kc} = 'separateflux';	p{kc}  =  separateflux;       kc = kc+1;%p{47}  Output seperate fluxes
name{kc} = 'simple_mpi_only'; 	p{kc}  =  simple_mpi_only;    kc = kc+1;%p{48}  Output seperate fluxes
name{kc} = 'integration_routine';	p{kc}  =  integration_routine; kc = kc+1;%p{15}(-) nag

nargu = kc-1;
stringind=[]; %no string indexes. Kept if adding any future string inputs

cd(runpath)

% Write binary files
for ii=1:nargu
  if any(ii == stringind)
    eval(sprintf('fid=fopen(\047input/p%d.txt\047,\047wb\047);',ii));
    eval(sprintf('fwrite(fid,p{%d});',ii));
    fclose(fid);
  else  
    eval(sprintf('fid=fopen(\047input/%s.bin\047,\047wb\047);',name{ii}));
    eval(sprintf('fwrite(fid,p{%d},\047double\047);',ii));
    fclose(fid);
  end
end

% Launch on interactive node (can be made neater with batch jobs)
tic
disp('============== Starts QuaLiKiz run ===============')
[sout,tout] = unix(sprintf('module purge \n module load intel/2018 mpi/2018 \n mpiexec -np 8 ./QuaLiKiz'));
toc
disp('============== End of QuaLiKiz run ===============')

cd(currentpath)

