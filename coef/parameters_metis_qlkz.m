phys_meth = 2.0;% Method on additional integrals on particle flux if not 0
coll_flag = 1.0;% Flag for collisionless simulations if 0
verbose = 1;
rot_flag = 0.0;
write_primi = 1.0;% Do not write primitive outputs if 0
numsols   = 3;  % Number of solutions in the output
nprocs    = 1;  % Number of processors used in the calculation
maxruns   = 1; %number of runs between contour checks
maxpts    = 5e5; %number of integrand evaluations done in 2D integral
relacc1   = 1e-3; %relative accuracy in 1D integrals
relacc2   = 2e-2; %relative accuracy in 2D integrals
timeout = 600;
ETGmult = 1.0;
collmult = 1.0;
separateflux=1.0;
simple_mpi_only=0.0; %for debugging; if 1 then no openMP used (F. Casson QLK development)
write_primi = 1.0;% Do not write primitive outputs if 0
integration_routine = 1;
absacc1 = 0;
absacc2 = 0;
rhomin = 0;
rhomax = 1;

ntheta=64; %resolution of parallel direction
numecoefs = 13; %number of outputs in the ecoefs matrix

set_ninorm1=1;      %flag for automating main ion concentration for maintaining QN
set_Ani1=1; %flag for automating main ion gradient for maintaining QN
set_QN_grad=1; %flag for maintaining quasineutrality of gradients.

%Set the number and range of wave number points
%kthetarhos = [linspace(0.1,1,10) 1.5 2 3 4 linspace(6,36,8)];
kthetarhos = linspace(0.1,2,10);

numn = length(kthetarhos); %number of ktheta in the system
scann=21;  %Number of points in parameter scan (radial points)

%NOTE: for general scans, any of the below can be changed to a vector of size (ones(1,scann))
%e.g. for a q-profile scan: qx=linspace(1,4,scann)

Bo = 3.7.*ones(scann,1); %magnetic field
Ro = 1.71.*ones(scann,1); %major radius
R0 = 1.71;
Rmin = 1.*ones(scann,1); %minor radius

x=linspace(0.,1,scann);	% x is the r/a throughout the scan: impacts the fraction of trapped particles	
rho=x;

% Following depends on exp profiles
qx = 2.00*ones(scann,1); %set q-profile
smag=1*ones(scann,1);
alphax=0.*ones(scann,1); %MHD alpha

%Electrons
Tex=1.*ones(scann,1); % kev
Nex=5.*ones(scann,1); % ne is in units of 10^19
Ate=9.*ones(scann,1);
Ane = 3.*ones(scann,1); 
el_type = 1; % 1 for active, 2 for adiabatic, 3 for adiabatic passing at ion scales (kthetarhosmax<2)
anise = 1.*ones(scann,1); 
danisedr = 0.*ones(scann,1); 
iind=1;

% Ions
ion_name{iind}='Main D';
Ai(1:scann,iind)=2.*ones(scann,1);
Zi(1:scann,iind)=1.*ones(scann,1);
Tix(1:scann,iind)=1.*ones(scann,1); % keV
ninorm(1:scann,iind)=1.*ones(scann,1);  % ni/ne arbitrary for main ions (will be rewritten to ensure quasineutrality)
Ati(1:scann,iind)=9.*ones(scann,1); %logarithmic gradient normalized to R
Ani(1:scann,iind) = 3.*ones(scann,1); %arbitrary for main ions, will be rewritten to ensure quasineutrality of gradients
ion_type(1:scann,iind)=1.*ones(scann,1); %1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
anis(1:scann,iind)=1.*ones(scann,1); %Tperp/Tpar at LFS. Tix defined as Tpar
danisdr(1:scann,iind)=0.*ones(scann,1); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
iind=iind+1;

% Tungsten
ion_name{iind}='W1';
Ai(1:scann,iind)=184.*ones(scann,1);
Zi(1:scann,iind)=45.*ones(scann,1);
Tix(1:scann,iind)=1.*ones(scann,1); % keV
ninorm(1:scann,iind)=1e-10*ones(scann,1);  ; % ni/ne arbitrary for main ions (will be rewritten to ensure quasineutrality)
Ati(1:scann,iind)= 9.*ones(scann,1); %logarithmic gradient normalized to R
Ani(1:scann,iind) = 1.*ones(scann,1); %arbitrary for main ions, will be rewritten to ensure quasineutrality of gradients
ion_type(1:scann,iind)=1.*ones(scann,1); %1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
anis(1:scann,iind)=1.*ones(scann,1); %Tperp/Tpar at LFS. Tix defined as Tpar
danisdr(1:scann,iind)=0.*ones(scann,1); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
iind=iind+1;

% Nitrogen
ion_name{iind}='N';
Ai(1:scann,iind)=14.*ones(scann,1);
Zi(1:scann,iind)=7.*ones(scann,1);
Tix(1:scann,iind)=1.*ones(scann,1); % keV
ninorm(1:scann,iind)=1e-10*ones(scann,1);  ; % ni/ne arbitrary for main ions (will be rewritten to ensure quasineutrality)
Ati(1:scann,iind)= 9.*ones(scann,1); %logarithmic gradient normalized to R
Ani(1:scann,iind) = 1.*ones(scann,1); %arbitrary for main ions, will be rewritten to ensure quasineutrality of gradients
ion_type(1:scann,iind)=1.*ones(scann,1); %1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
anis(1:scann,iind)=1.*ones(scann,1); %Tperp/Tpar at LFS. Tix defined as Tpar
danisdr(1:scann,iind)=0.*ones(scann,1); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
iind=iind+1;

% Second impurity in METIS (with rimp)
ion_name{iind}='Cu';
Ai(1:scann,iind)=1.*ones(scann,1);
Zi(1:scann,iind)=1.*ones(scann,1);
Tix(1:scann,iind)=1.*ones(scann,1); % keV
ninorm(1:scann,iind)=1*ones(scann,1);  ; % ni/ne arbitrary for main ions (will be rewritten to ensure quasineutrality)
Ati(1:scann,iind)= 9.*ones(scann,1); %logarithmic gradient normalized to R
Ani(1:scann,iind) = 1.*ones(scann,1); %arbitrary for main ions, will be rewritten to ensure quasineutrality of gradients
ion_type(1:scann,iind)=1.*ones(scann,1); %1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
anis(1:scann,iind)=1.*ones(scann,1); %Tperp/Tpar at LFS. Tix defined as Tpar
danisdr(1:scann,iind)=0.*ones(scann,1); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
iind=iind+1;

% Minority
ion_name{iind}='H';
Ai(1:scann,iind)=1.*ones(scann,1);
Zi(1:scann,iind)=1.*ones(scann,1);
Tix(1:scann,iind)=1.*ones(scann,1); % keV
ninorm(1:scann,iind)=1e-10*ones(scann,1);  ; % ni/ne arbitrary for main ions (will be rewritten to ensure quasineutrality)
Ati(1:scann,iind)= 9.*ones(scann,1); %logarithmic gradient normalized to R
Ani(1:scann,iind) = 1.*ones(scann,1); %arbitrary for main ions, will be rewritten to ensure quasineutrality of gradients
ion_type(1:scann,iind)=1.*ones(scann,1); %1 for active, 2 for adiabatic, 3 for tracer (also won't be included in quasineutrality checks), 4 for tracer but included in Zeff
anis(1:scann,iind)=1.*ones(scann,1); %Tperp/Tpar at LFS. Tix defined as Tpar
danisdr(1:scann,iind)=0.*ones(scann,1); %d/dr(Tperp/Tpar) at LFS. Tix defined as Tpar
iind=iind+1;

nions=iind-1; %number of ions in system

%rotation inputs. Intrinsic assumption that all species rotate together
Machtor = 0.0.*ones(scann,1); %Main ion mach number of the toroidal rotation. All others are scaled to this.
Machpar = 0.0.*ones(scann,1); %Main ion mach number of the toroidal rotation. All others are scaled to this.
gammaE = 0.0.*ones(scann,1); %ExB shear. Defined with the sqrt(2)!
Aupar = 0.*ones(scann,1); %Parallel velocity shear. -1111 means that the consistent value with toroidal rotation is calculated.   
Autor = 0.*ones(scann,1); %Parallel velocity shear. -1111 means that the consistent value with toroidal rotation is calculated.   









