function [Chie,Chii,De,Ve] = run_tglf_3sp(tglfparams,post,itime,iplt);

if ~exist('iplt')
  iplt=0;
end

% To do:
% Force consistency with different switches like sat2 and xnumodel=3

%%%% Initialisations of the transport flux and its coefficients
Chie = zeros(1,length(post.profil0d.xli));   % electron heat diffusivity
Chii = zeros(1,length(post.profil0d.xli));   % ion heat diffusivity

Deff_elec = zeros(1,length(post.profil0d.xli));   % effective diffusion coefficient
Deff_imp = zeros(1,length(post.profil0d.xli));   % effective diffusion coefficient

De   = zeros(1,length(post.profil0d.xli));   % electron particle diffusion coefficient
Di   = zeros(1,length(post.profil0d.xli));   % ion particle diffusion coefficient
Dimp = zeros(1,length(post.profil0d.xli));

Ve   = zeros(1,length(post.profil0d.xli));   % electron convective coefficient
Vimp   = zeros(1,length(post.profil0d.xli));

Gelec = zeros(1,length(post.profil0d.xli)); % collective elctron particle flux
Gion1 = zeros(1,length(post.profil0d.xli)); % collective ion (main) particle flux
Gion2 = zeros(1,length(post.profil0d.xli)); % collective ion (here for impurity) particle flux 

Gelec_diffusive = zeros(1,length(post.profil0d.xli));
Gelec_convective = zeros(1,length(post.profil0d.xli));
Gimp_diffusive = zeros(1,length(post.profil0d.xli));
Gimp_convective= zeros(1,length(post.profil0d.xli));


Qelec = zeros(1,length(post.profil0d.xli)); % elctron heat flux
Qion1 = zeros(1,length(post.profil0d.xli)); % ion (main) heat flux 
Qion2 = zeros(1,length(post.profil0d.xli)); % ion (here for impurity) heat flux

% Path to TGLF routines, need to provide something for a TGLF install function
path_to_METIS_coef = fullfile(fileparts(which('install_QLKNN10D')));
path_to_METIS_TGLF = [path_to_METIS_coef '/TGLF_TOOLS'];
addpath(path_to_METIS_TGLF);

%disp(path_to_METIS_coef);
%disp(path_to_METIS_TGLF);

% To store the output folders
if isempty(getappdata(0,'local'))
  path_runs = [path_to_METIS_TGLF '/runs_tglf_std'];
  path_plots = [path_to_METIS_TGLF '/plots_tglf_std'];  
else
  path_runs = [getappdata(0,'local') '/runs_tglf_std'];
  path_plots = [getappdata(0,'local') '/plots_tglf_std'];  
end

unix(sprintf(['mkdir -p ' path_runs]));
unix(sprintf(['mkdir -p ' path_plots]));

%%%%% Read input TGLF template and make the changes needed >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tglf_input = [path_to_METIS_TGLF '/input_templates'];
tglfin = read_tglfinput_3sp(tglf_input);

tglfin.UNITS 		    = tglfparams.UNITS;
tglfin.FILTER 		    = tglfparams.FILTER;
tglfin.NS 		    = 3; %tglfparams.NS;  		
tglfin.USE_TRANSPORT_MODEL  = tglfparams.USE_TRANSPORT_MODEL;
tglfin.GEOMETRY_FLAG        = tglfparams.GEOMETRY_FLAG;
tglfin.USE_BPER             = tglfparams.USE_BPER;
tglfin.USE_BPAR 	    = tglfparams.USE_BPAR;
tglfin.USE_MHD_RULE 	    = tglfparams.USE_MHD_RULE;
tglfin.USE_BISECTION	    = tglfparams.USE_BISECTION;
tglfin.USE_INBOARD_DETRAPPED= tglfparams.USE_INBOARD_DETRAPPED;
tglfin.SAT_RULE 	    = tglfparams.SAT_RULE;
tglfin.KYGRID_MODEL 	    = tglfparams.KYGRID_MODEL; 
tglfin.XNU_MODEL 	    = tglfparams.XNU_MODEL;
tglfin.VPAR_MODEL           = tglfparams.VPAR_MODEL;
tglfin.SIGN_BT              = tglfparams.SIGN_BT;
tglfin.SIGN_IT              = tglfparams.SIGN_IT;
tglfin.KY                   = tglfparams.KY;
tglfin.NEW_EIKONAL          = tglfparams.NEW_EIKONAL;
tglfin.VEXB_SHEAR 	    = tglfparams.VEXB_SHEAR;
tglfin.DEBYE 		    = tglfparams.DEBYE;
tglfin.IFLUX 		    = tglfparams.IFLUX;
tglfin.IBRANCH		    = tglfparams.IBRANCH;
tglfin.NMODES 		    = tglfparams.NMODES;
tglfin.NBASIS_MAX 	    = tglfparams.NBASIS_MAX;
tglfin.NBASIS_MIN 	    = tglfparams.NBASIS_MIN;
tglfin.NXGRID 		    = tglfparams.NXGRID;
tglfin.NKY 		    = tglfparams.NKY;
tglfin.ADIABATIC_ELEC       = tglfparams.ADIABATIC_ELEC; 
tglfin.ALPHA_MACH 	    = tglfparams.ALPHA_MACH;
tglfin.ALPHA_E 		    = tglfparams.ALPHA_E;
tglfin.ALPHA_P 		    = tglfparams.ALPHA_P;
tglfin.ALPHA_QUENCH	    = tglfparams.ALPHA_QUENCH;
tglfin.XNU_FACTOR	    = tglfparams.XNU_FACTOR;
tglfin.DEBYE_FACTOR	    = tglfparams.DEBYE_FACTOR;
tglfin.ETG_FACTOR 	    = tglfparams.ETG_FACTOR;
tglfin.WRITE_WAVEFUNCTION_FLAG= tglfparams.WRITE_WAVEFUNCTION_FLAG;

%%%%% PHYSICS CONSTANTS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
phys.e   =   1.602176462e-19;
phys.mD  =  3.345e-27;                 % mass of deuterium
phys.mu0 =   4*pi*1e-7;                % permeability in vacuum (H/m) (definition)
phys.kB  =   1.3806503e-23;            % Boltzmann Constant (J/K)
phys.ua  =   1.66053873e-27;           % in kgs 


%%%%% Call the profiles from METIS:>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% from profil0d in post (METIS outputs):
xli	     	 = post.profil0d.xli;   	% normalised radius (1,nr=21)
time_profil0d= post.profil0d.temps; 	% (nt1,1)

%from 'geo' and 'option' in z0dinput in post, Geometrical parameters (METIS inputs):  
minor_radius = post.z0dinput.geo.a;  	  % (nt0,1) minor radius of the plasma waveform (m).
major_radius = post.z0dinput.geo.R;  	  % (nt0,1) major radius of the plasma waveform (m).
vertical_pos = post.z0dinput.geo.R*0; 	  % (nt0,1) vertical position of the plasma waveform (m) set to 0 here.
Btoroidal    = post.z0dinput.geo.b0;	  % (nt0,1) vacuum magnetic toroidal field measured at geo.R(T)
rimp         = post.z0dinput.option.rimp; % ratio of density b/w main impurity and another impurity
main_ion     = post.z0dinput.option.gaz;  % main ion 1>H, 2>D, 3>T, 4>He
zlightimp    = post.z0dinput.option.zimp; % charge of other mpurity

% from profil0d in post
kappa 	   = post.profil0d.kx; 	                     % elongation (nt1,nr)
delta      = post.profil0d.dx;	                     % triangularity (nt1,nr)
epsilon    = post.profil0d.epsi;
R_maj	   = post.profil0d.Raxe;	             % (nt1,nr) major radius of centre of each flux surface  (m)
r_min      = post.profil0d.epsi.*post.profil0d.Raxe; % (nt1,nr) minor radius of centre of each flux surface  (m) 
q          = post.profil0d.qjli; 	             % Safety factor (nt1,nr)
TeeV       = post.profil0d.tep;			     % electron temperature (eV)
TieV       = post.profil0d.tip;			     % ion temperature (eV)
Te         = post.profil0d.tep.*phys.e; 	     % (nt1,nr) (J)
Ti         = post.profil0d.tip.*phys.e;              % (nt1,nr) (J)
Tz 	   = post.profil0d.tip.*phys.e; 	     % (nt1,nr) (J)
n_e	   = post.profil0d.nep; 	  	     % (nt1,nr) electron density (m^-3) 
n_lightimp = post.profil0d.nzp;%.*1e-13;     	     % (nt1,nr) density of light impurity (m^-3)
vtor       = post.profil0d.vtor;      	             % (nt1,nr) Toroidal rotation speed (m/s)
ptotal     = post.profil0d.ptot; 		     % (nt1,nr) Total pressure profile (Pa)
Phi 	   = post.profil0d.phi;     	             % (nt1,nr) Toroidal magnetic flux (Wb)
Psi 	   = post.profil0d.psi;    	             % (nt1,nr) Poloidal  magnetic flux (Wb)
Bpol 	   = post.profil0d.bpol;
xie	   = post.profil0d.xie;			     % electron heat diffusivity estimation using power balance	
xii	   = post.profil0d.xii;			     % ion heat diffusivity estimation using power balance	

% from zerod in post (Time dependent Scaler Data given in output file):
time_zerod   = post.zerod.temps; % (nt0,1)


%%%%% Interpolation and Reshaping on profiles 1D time grid.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
a0_interp   = interp1(time_zerod,minor_radius,time_profil0d,'spline','extrap');  
R0_interp   = interp1(time_zerod,major_radius,time_profil0d,'spline','extrap'); 
z0_interp   = interp1(time_zerod,vertical_pos,time_profil0d,'spline','extrap');
Btor_interp = interp1(time_zerod,Btoroidal,time_profil0d,'spline','extrap');
mp = 1.6726e-27;
try
  mass_dt     = interp1(time_zerod,(2 .* (1 - post.z0dinput.cons.iso) +  3 .* post.z0dinput.cons.iso)*mp,time_profil0d);
catch
  fiso = post.zerod.nTm./(post.zerod.nDm+post.zerod.nTm);
  mass_dt     = interp1(time_zerod,(2 .* (1 - fiso) +  3 .* fiso)*mp,time_profil0d);  
end

xli_rep = repmat(xli,length(post.profil0d.temps),1);
a0      = repmat(a0_interp,1,length(xli)); 	%(nt1,nr) (m)
R0      = repmat(R0_interp,1,length(xli)); 	%(nt1,nr) (m)
z0      = repmat(z0_interp,1,length(xli)); 	%(nt1,nr) (m)
Btor    = repmat(Btor_interp,1,length(xli));
R       = R_maj + (xli_rep .* a0); 		% (nt1,nr) R = R0 + r Cos(theta=0) (m).


%%%%% 
tabmat = load('Lz_zave.mat');
tabmat = tabmat.tabmat;

impurities = fieldnames(tabmat);
for i = 1:length(impurities)
	impurity = impurities{i};
	current_impurity = getfield(tabmat, impurity); % Access each impurity structure

	if current_impurity.Z == zlightimp
		mass_lightimp = current_impurity.A./2;
            for ll=1:length(time_profil0d)
  		  		charge_lightimp(ll,:) = round(interp1(current_impurity.data(:,1)*1e3,current_impurity.data(:,3),TeeV(ll,:)));
            end
    end
	
end

%%%%% imposing quasi neutrality: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
n_i = (n_e - n_lightimp.*charge_lightimp)/(tglfin.ZS_2);

%%%%%% Gradients d(quantity)/dRMIN_LOC:>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dtedxli  	= pdederive(xli,Te,0,2,2,1);  		% (nt1,nr) electron tmperature gradient
dtidxli  	= pdederive(xli,Ti,0,2,2,1);  		% (nt1,nr) ion tmperature gradient
dnedxli  	= pdederive(xli,n_e,0,2,2,1); 		% (nt1,nr) electron density gradient
dnidxli  	= pdederive(xli,n_i,0,2,2,1); 		% (nt1,nr) ion density gradient
dnlightimpdxli  = pdederive(xli,n_lightimp,0,2,2,1); 	% (nt1,nr) main impurity density gradient
dvtordxli 	= pdederive(xli,(vtor./R),0,2,2,1); 	% (nt1,nr) d(toroidal angular rotation frequency)/dxli, xli = r*a
dqdxli    	= pdederive(xli,q,0,2,2,1);   		% (nt1,nr) d(safety factor)/dxli
dkdxli 	  	= pdederive(xli,kappa,0,2,2,1); 	% (nt1,nr) d(kappa)/dxli
dddxli    	= pdederive(xli,delta,0,2,2,1); 	% (nt1,nr) d(delta)/dxli
dpdxli    	= pdederive(xli,ptotal,0,2,2,1); 	% (nt1,nr)
dPsidxli  	= pdederive(xli,Psi,0,2,2,1); 

%%%%% Some physics quantities>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
SIGN_Itorr = tglfin.SIGN_IT; 		                 % Itorr is the plasma current, positive for counter-clockwise 
Cs	   = sqrt(Te./phys.mD); 			 % (nt1,nr)deuterium sound speed
B_unit	   = abs((q ./ (xli_rep.* a0.*a0)) .* dPsidxli); % (nt1,nr) psi is poloidal magnetic field.
Omegai	   = phys.e * B_unit ./phys.mD;  		 % (nt1,nr)
rhol_ref   = Cs./Omegai;  				 % (nt1,nr)

%%%%%% Collisionality and Betae (Reference: Wesson)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
lnlamda = 24.0 - (0.5 .* log(n_e ./ 1e6)) + log(TeeV);       
taue	= ((3.44e5) .* ((TeeV).^1.5))./((n_e ./ 1e6) .* lnlamda);  	%  sec
xnuei   = 0.75 .* sqrt(pi)./taue ;         				% 1/sec 
xnuei   = xnuei .* a0 ./ Cs ;            				% normalized electron-ion collision rate
betae   = (8.0 * pi .* (n_e./ 1e6) .* (Te*1e7)) ./(B_unit.^2 .*1e8); 		% electron beta unit

%%%%% SPECIES VECTORS (Cant be changed by the user.)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rlte	    = -dtedxli ./ Te;  					% (nt1,nr) Normalised electron temperature gradient calculation
rlti  	    = -dtidxli ./ Ti;  					% (nt1,nr) Normalised ion temperature gradient calculation
rlne  	    = -dnedxli ./ n_e; 					% (nt1,nr) Normalised electron density gradient calculation
rlni        = -dnidxli ./ n_i; 					% (nt1,nr) Normalised ion density gradient calculation
rlnlightimp = -dnlightimpdxli ./ n_lightimp; 			% (nt1,nr) Normalised main inpurity density gradient
TeTe  	    = Te./Te;          					% (nt1,nr) Normalised temperature calculation TAUS_S: 
TiTe  	    = Ti./Te;        	 				% (nt1,nr)
TzTe        = Tz./Te;		 				% (nt1,nr)
nene        = n_e./n_e;        					% Normalised density calculation AS_S: 
nine        = n_i./n_e;	 					% (nt1,nr) 	
nlightimpne = n_lightimp./n_e;	 				% (nt1,nr)
Vpar_e	    =  (SIGN_Itorr .* R_maj .* vtor) ./ (R .* Cs);	% VPAR_1 Species array of parallel velocities:
Vpar_shear_e= (- SIGN_Itorr .* R_maj .* dvtordxli) ./ Cs; 	% VPAR_SHEAR_1: Normalised parallel velocity gradient. 

%%%%% Calculate zeff
zeff        = nine.*tglfin.ZS_2.^2 + nlightimpne.*charge_lightimp.^2 ;

%%%%% MILLER GEOMETRY PARAMETERS: changes with space r and time t>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rmin_loc    = xli_rep;                    					% (nt1,nr) Flux surface centroid minor radius r/a 
rmaj_loc    = R_maj./a0;                 					%(nt1,nr) Flux surface centroid major radius R0/a0.
zmaj_loc    = z0./a0;                     					%(nt1,nr) Flux surface centroid elevation Zmaj/a0
q_loc	    = abs(q);                       					% (nt1,nr) Absolute value of the safety factor
q_prime_loc = (q ./rmin_loc).* dqdxli; 						%(nt1,nr) 
p_prime_loc = (phys.mu0/4*pi) .* (q./ (rmin_loc .* (B_unit.^2))) .* dpdxli;
drmindx_loc = 1;                       						%  x=r_min
drmajdx_loc = pdederive(xli,R_maj,0,2,2,1)./a0;
dzmajdx_loc = 0;
kappa_loc   = kappa;
s_kappa_loc = (rmin_loc ./ kappa_loc) .* dkdxli; 
delta_loc   = delta;

s_delta_loc = rmin_loc .*  dddxli;
zeta_loc    = 0;
s_zeta_loc  = 0;
kx0_loc     = 0;  

gammas_max_all = zeros(1,length(xli)); 
rlne_all = zeros(1,length(xli));
rlte_all = zeros(1,length(xli));  
rlnimp_all = zeros(1,length(xli));
rltimp_all = zeros(1,length(xli));  


jj = itime; 

% Main Ions:
switch main_ion
		case 1 
			tglfin.ZS_2 = 1;
			tglfin.MASS_2 = 0.5; % H
			
		case 2
			tglfin.ZS_2 = 1;
			tglfin.MASS_2 = 1;   % D
			
		case 3
		        tglfin.ZS_2 = 1;
   		        tglfin.MASS_2 = mass_dt(jj)./(2*mp); % DT
			
		case 4
			tglfin.ZS_2 = 2;
			tglfin.MASS_2 = 2;   % He

end

%%%%% Loop over the radial points for the diffusion coefficients vs r/a
for ii = 2:length(xli)  
disp([' charge_lightimp = ', num2str(charge_lightimp(jj,ii))]);
disp([' mass_lightimp = ', num2str(mass_lightimp)]);

	disp(['================================ RUNS TGLF r/a=' num2str(xli(ii)) '================================']);

	% Create the directory for this radial point 
	tglf_input_scan_radial = [path_runs '/tglf_tmp_radial_' num2str(ii)];
        unix(sprintf(['mkdir -p ' tglf_input_scan_radial]));  	  

	% Parameters Initialization:
	rln_offset = 0.01;
	R_squared_threshold = 0.9899;
	max_refinement_iterations = 5;
    	refinement_iteration = 0;
	%disp(['rlne_variations:',num2str(rlne_variations)]);
	%disp(class(rlne_variations));
	%disp(rlne_variations(2));
	
	while refinement_iteration < max_refinement_iterations
		 
		rln_variations = [0, rln_offset];
		rlne_values = []; 
		%rlnimp_values = [];
		elec_flux_values = [];
		%imp_flux_values = [];
	

		% Loop over rln variations
		for offset = rln_variations
    			modified_rlne = rlne(jj, ii) + offset;
			%modified_rlnimp = rlnlightimp(jj,ii) + offset;

			tglf_input_scan_rln = [tglf_input_scan_radial '/rlne_' num2str(modified_rlne)];
                        unix(sprintf(['mkdir -p ' tglf_input_scan_rln]));  	  			  

			% Setup for the inpus of tglf inputs for this rln
			% ELECTRONS
			tglfin.ZS_1         = -1;    	 
			tglfin.MASS_1       = 2.723100e-04;
			tglfin.RLNS_1       = modified_rlne; %rlne(jj,ii);
			tglfin.RLTS_1       = rlte(jj,ii);
			tglfin.TAUS_1       = TeTe(jj,ii); 
			tglfin.AS_1	    = nene(jj,ii);
			tglfin.VPAR_1 	    = Vpar_e(jj,ii);
			tglfin.VPAR_SHEAR_1 = Vpar_shear_e(jj,ii);

			% #2 MAIN IONS:
			tglfin.RLNS_2       = rlni(jj,ii);
			tglfin.RLTS_2       = rlti(jj,ii);
			tglfin.TAUS_2       = TiTe(jj,ii);
			tglfin.AS_2 	    = nine(jj,ii);
			tglfin.VPAR_2 	    = Vpar_e(jj,ii);
			tglfin.VPAR_SHEAR_2 = Vpar_shear_e(jj,ii);

			% #3 Some light impurity like C, N, Be
    			tglfin.ZS_3         = charge_lightimp(jj,ii); 
			tglfin.MASS_3	    = mass_lightimp;	
			tglfin.RLNS_3 	    = rlnlightimp(jj,ii);   % modified_rlnimp;
 			tglfin.RLTS_3	    = rlti(jj,ii);
			tglfin.TAUS_3	    = TzTe(jj,ii);
			tglfin.AS_3         = nlightimpne(jj,ii);
			tglfin.VPAR_3	    = Vpar_e(jj,ii);	
			tglfin.VPAR_SHEAR_3 = Vpar_shear_e(jj,ii);

			% General inputs
			tglfin.BETAE 	   = betae(jj,ii);
			tglfin.XNUE	   = xnuei(jj,ii);
			tglfin.ZEFF 	   = zeff(jj,ii);
			tglfin.RMIN_LOC    = rmin_loc(jj,ii);
			tglfin.RMAJ_LOC    = rmaj_loc(jj,ii);
			tglfin.ZMAJ_LOC    = zmaj_loc(jj,ii); 
			tglfin.Q_LOC       = q_loc(jj,ii);
			tglfin.Q_PRIME_LOC = q_prime_loc(jj,ii);
			tglfin.P_PRIME_LOC = p_prime_loc(jj,ii);
			tglfin.DRMINDX_LOC = drmindx_loc; 
			tglfin.DRMAJDX_LOC = drmajdx_loc(jj,ii);
			tglfin.DZMAJDX_LOC = dzmajdx_loc;
			tglfin.KAPPA_LOC   = kappa_loc(jj,ii);
			tglfin.S_KAPPA_LOC = s_kappa_loc(jj,ii); 
			tglfin.DELTA_LOC   = delta_loc(jj,ii);
			tglfin.S_DELTA_LOC = s_delta_loc(jj,ii);
			tglfin.ZETA_LOC    = zeta_loc;
			tglfin.S_ZETA_LOC  = s_zeta_loc; 
			tglfin.KX0_LOC     = kx0_loc; 
		

			make_tglfinput(tglf_input_scan_rln,tglfin); %make_tglfinput is a function to generate the inputs for tglf 


			% Run TGLF
			cd(tglf_input_scan_rln)
			!rm -f out.tglf*
        		disp([ '===TGLF for rlne=' num2str(modified_rlne)])
			
			[sout,tout] = unix(sprintf(['module purge; module load gcc/13.2.0 lapack/3.9.1-gcc-13.2.0 openmpi/4.1.1-gcc-13.2.0; export GACODE_PLATFORM=IRFM; export GACODE_ROOT=/ZONE_TRAVAIL/PM258588/codes/gacode; . ${GACODE_ROOT}/shared/bin/gacode_setup; nice -n 9 tglf -e ./ -n ' num2str(round(tglfparams.numproc))]));
       	 		disp('===  ENDS OF TGLF RUN  ===')
                        disp(['=== RUN DATA STORED IN ' tglf_input_scan_rln])

                 	try
                          tglfout = read_tglf_outrun_3sp(tglf_input_scan_rln);
                	catch
                 	  tglfout = read_tglf_outrun_gbflux_3sp(tglf_input_scan_rln);
                  	end
			
			Gelec_per_rlne = tglfout.Ge_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Cs(jj,ii); % Flux in SI units
			Gimp_per_rlnimp = tglfout.G2_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Cs(jj,ii); % Flux in SI units

			rlne_values = [rlne_values,modified_rlne];
			%rlnimp_values = [rlnimp_values,modified_rlnimp];

			elec_flux_values = [elec_flux_values, Gelec_per_rlne];
			%imp_flux_values = [imp_flux_values, Gimp_per_rlnimp];
		 

			% store flux corresponding to original rln
        		if offset == 0
          			Gelec(ii) = Gelec_per_rlne;
				Gion1(ii) = tglfout.G1_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Cs(jj,ii); 
				Gion2(ii) = Gimp_per_rlnimp;
				
				Qelec(ii) = tglfout.Qe_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Te(jj,ii) * Cs(jj,ii);
				Qion1(ii) = tglfout.Q1_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Te(jj,ii) * Cs(jj,ii);
				Qion2(ii) = tglfout.Q2_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Te(jj,ii) * Cs(jj,ii);
	
			
				tglfout_eigens = read_eigenvalues(tglf_input_scan_rln);
				tglfout_ky = read_ky_spectrum(tglf_input_scan_rln);
				gamma1(ii,:) =  tglfout_eigens.gamma1;
				frequency1(ii,:) = tglfout_eigens.freq1; 
				gamma2(ii,:) =  tglfout_eigens.gamma2;
				frequency2(ii,:) = tglfout_eigens.freq2; 
				ky_vals(ii,:) = tglfout_ky.ky_values;

				%disp([' ky_values = ', num2str(ky_vals(ii,:))]);
				%disp([' gamma1_values = ', num2str(gamma1(ii,:))]);


			end % end of if condition for storing the flux and eigenvalues
				
	
		end % end for loop for rlne variation loop 

	
	% Perform the linear fitting for G(SI) = (D*n_e/a)*rlne + V*n_e 
	elec_coeffs = polyfit(rlne_values, elec_flux_values, 1); 	% linear fit
	%imp_coeffs = polyfit(rlnimp_values, imp_flux_values, 1);
	
	rlne_fit = linspace(min(rlne_values),max(rlne_values),100); % defined to make a smooth plot
	%rlnimp_fit = linspace(min(rlnimp_values),max(rlnimp_values),100); % defined to make a smooth plot

	elec_flux_fit = polyval(elec_coeffs,rlne_fit); 		% gives fitted values of fluxes
	%imp_flux_fit = polyval(imp_coeffs,rlnimp_fit);	

	% Calculate goodness of fit (R-squared)
	elec_residuals = elec_flux_values - polyval(elec_coeffs,rlne_values);
        elec_SS_residual = sum(elec_residuals.^2);
        elec_SS_total = sum((elec_flux_values - mean(elec_flux_values)).^2);
        elec_R_squared = 1 - (elec_SS_residual / elec_SS_total);
	%disp([' R^2 for elec = ', num2str(elec_R_squared)]);
	disp(['R^2 = ', num2str(elec_R_squared)]);


		%imp_residuals   = imp_flux_values - polyval(imp_coeffs,rlnimp_values);
    		%imp_SS_residual = sum(imp_residuals.^2);
   		%imp_SS_total    = sum((imp_flux_values - mean(imp_flux_values)).^2);
    		%imp_R_squared   = 1 - (imp_SS_residual / imp_SS_total);
		%disp([' R^2 for imp = ', num2str(imp_R_squared)]);




	% Check if the fit is consistent
        if elec_R_squared >= R_squared_threshold %elec_R_squared >= R_squared_threshold  &&
       		disp(['Good fit']);
		
			Gelec_diffusive(ii) = elec_coeffs(1).*rlne(jj,ii);
			Gelec_convective(ii) = elec_coeffs(2);

			%Gimp_diffusive(ii) = imp_coeffs(1).*rlnlightimp(jj,ii);
			%Gimp_convective(ii) = imp_coeffs(2);


			%Calculate transport coefficients
			Delec = elec_coeffs(1) .*a0(jj,ii) ./ n_e(jj, ii);  % Slope gives D
			Velec = elec_coeffs(2) ./ n_e(jj, ii);  % Intercept divided by density is V

			%Dion2 = imp_coeffs(1) .*a0(jj,ii) ./ n_lightimp(jj, ii);  % Slope gives D
			%Vion2 = imp_coeffs(2) ./ n_lightimp(jj, ii);  % Intercept divided by density is V
			
			De(ii) = Delec;
			Ve(ii) = Velec;

			%Dimp(ii) = Dion2;
			%Vimp(ii) = Vion2;

			


			%figure;
			%plot(rlnimp_values,imp_flux_values,'-bo', 'Linewidth',2,'DisplayName','Data');
			%hold on;
			%plot(rlnimp_fit,imp_flux_fit,'r-', 'Linewidth',2,'DisplayName',['Fit, R^2 =',num2str(imp_R_squared)]);
			%xlabel('Normalised impurity density gradient');
			%ylabel('G (SI)');
			%legend('location', 'best');
			%title(['D = ' num2str(Dion2) ' and V = ' num2str(Vion2) 'at radial point' num2str(xli(ii))]);
			%set(gca,'FontSize',30);
			%grid on;
			%saveas(gcf, fullfile(path_to_plots, ['Flux_vs_rlnimp_ii_' num2str(ii) '.jpg']));
			%saveas(gcf, fullfile(path_to_plots, ['Flux_vs_rlnimp_ii_' num2str(ii) '.fig']));
			%close;
			disp(['R^2 = ', num2str(elec_R_squared)]);

		break;
	else
		disp('Poor Fit')
		refinement_iteration = refinement_iteration + 1;
            	rlne_offset = rlne_offset/2;
        
	end % if else for the R^2 > R^2_treshold

    end % while end

ky_cutoff = 1;
ky_mask = ky_vals(ii, :) < ky_cutoff;
%disp([' ky_mask = ', num2str(ky_mask)]);
ky_valid = ky_vals(ii, ky_mask);
%disp([' ky_valid = ', num2str(ky_valid)]);
gamma_valid = gamma1(ii, ky_mask);
%disp([' gamma_valid = ', num2str(gamma_valid)]);

if ~isempty(gamma_valid)
  [gamma_max, idx_max] = max(gamma_valid);
  ky_at_gamma_max = ky_valid(idx_max);
else
  gamma_max=0;
  ky_at_gamma_max = NaN;
end

gammas_max_all(ii) = gamma_max;

% Store or display
%fprintf('r/a = %.3f | gamma_max = %.4f at ky = %.3f\n', xli(ii), gamma_max, ky_at_gamma_max);

%rlne_all(ii) = rlne(jj,ii);
%rlnimp_all(ii) = rlnlightimp(jj,ii);
%rlte_all(ii) = rlte(jj,ii);
%rltimp_all(ii) = rlti(jj,ii);
	
end % end of radial loop

% Compute the effective diffusion coefficients
Chie = - Qelec./n_e(jj,:)./(dtedxli(jj,:)./a0(jj,:));
Chii = - Qion1./n_i(jj,:)./(dtidxli(jj,:)./a0(jj,:));

% Magnetic axis
Chie(1)=0;
Chii(1)=0;

Deff_elec = - Gelec./(dnedxli(jj,:)./a0(jj,:)); % from TGLF flux
Deff_ion = - Gion1./(dnidxli(jj,:)./a0(jj,:));
Deff_imp = - Gion2./(dnlightimpdxli(jj,:)./a0(jj,:));

Deff_elec_rln =  Gelec.*a0(jj,:)./n_e(jj,:)./rlne(jj,:);
Deff_imp_rln =   Gion2.*a0(jj,:)./n_lightimp(jj,:)./rlnlightimp(jj,:);

% For my own peace of mind.
gradne = zeros(1,length(post.profil0d.xli));
gradne  = - n_e(jj,:) .* rlne(jj,:);
gradnimp = zeros(1,length(post.profil0d.xli));
gradnimp  = - n_lightimp(jj,:) .* rlnlightimp(jj,:);

e_flux_tot_calculated_rln = zeros(1,length(post.profil0d.xli)); %G(SI) = (D*n_e/a)*rlne + V*n_e
e_flux_tot_calculated_rln = (De .* n_e(jj,:) .* rlne(jj,:)./a0(jj,:)) + (Ve .* n_e(jj,:));

e_flux_tot_calculated_grad = zeros(1,length(post.profil0d.xli)); 
e_flux_tot_calculated_grad = (- De .* (dnedxli(jj,:)./a0(jj,:))) + (Ve .* n_e(jj,:));

%%%%%%%%%%%%%%%%%%%%%%%%Plot Fluxes

figure;
subplot(2,2,1)
plot(xli, Qelec, '-ro', 'LineWidth', 2,'DisplayName','el.');
hold on;
plot(xli, Qion1, '-ko', 'LineWidth', 2,'DisplayName','main ion');
hold on;
plot(xli, Qion2, '-bo', 'LineWidth', 2,'DisplayName','imp.');
xlabel('r/a');
ylabel('Heat flux (W/m^2)');
legend('location', 'best');
set(gca,'FontSize',14);
grid on;

subplot(2,2,2)
plot(xli, Gelec, '-ro', 'LineWidth', 2,'DisplayName','el.');
hold on;
plot(xli, Gion1, '-ko', 'LineWidth', 2,'DisplayName','main ion');
hold on;
plot(xli, Gion2, '-bo', 'LineWidth', 2, 'DisplayName','imp.');
xlabel('r/a');
ylabel('Particle Flux (m^{-2} s^{-1})');
legend('location', 'best');
set(gca,'FontSize',14);
grid on;

subplot(2,2,3)
plot(xli, Gelec_diffusive, '-ro', 'LineWidth', 2);
xlabel('r/a');
ylabel('Diffusive flux (m^{-2} s^{-1})');
%xlim([0 0.9]);
legend('location', 'best');
set(gca,'FontSize',14);
grid on;


subplot(2,2,4)
plot(xli, Gelec_convective, '-ro', 'LineWidth', 2);
xlabel('r/a');
ylabel('Convective flux (m^{-2} s^{-1})');
%xlim([0 0.9]);
legend('location', 'best');
set(gca,'FontSize',14);
grid on;

%%%%%%%%%%%%%%%%%%%%Plot transport coefficients 
figure;
subplot(2,2,1)
plot(xli,Chie,'rs-','linewidth',2);
hold on;
plot(xli, Chii, 'ko-', 'LineWidth', 2);
hold on; 
plot(xli, post.profil0d.xie(jj,:), 'r--', 'LineWidth', 2);
hold on;
plot(xli, post.profil0d.xii(jj,:), 'b--', 'LineWidth', 2);
xlabel('r/a', 'FontSize', 14)
ylabel('\chi (m^2/s)', 'FontSize', 14)
set(gca,'FontSize',14)
legend('elec(TGLF)', 'main ion(TGLF)','elec(PB)', 'main ion(PB)')
grid on;

subplot(2,2,2)
plot(xli, Deff_elec, '-ro', 'LineWidth', 2, 'DisplayName','el.'); hold on;
plot(xli, Deff_imp, '-bo', 'LineWidth', 2, 'DisplayName','imp.');   hold on; 
plot(xli, Deff_elec_rln, '--ro', 'LineWidth', 2, 'DisplayName','Deff elec rln'); hold on;
plot(xli, Deff_imp_rln, '--bo', 'LineWidth', 2, 'DisplayName','Deff imp rln');
xlabel('r/a');
ylabel('D_{eff} (m^2/s)');
legend('location', 'best');
%xlim([0 0.9]);
set(gca,'FontSize',14);
grid on;

subplot(2,2,3)
plot(xli, De, '-ro', 'LineWidth', 2);
xlabel('r/a');
ylabel('D_e (m^2/s)');
legend('location', 'best');
%xlim([0 0.9]);
set(gca,'FontSize',14);
grid on;

subplot(2,2,4)
plot(xli, Ve, '-ro', 'LineWidth', 2);
xlabel('r/a');
ylabel('V_e (m/s)');
legend('location', 'best');
%xlim([0 0.9]);
set(gca,'FontSize',14);
grid on;

%%%%%%%%%%%%%%%%%%%%% Plot profiles
figure;
subplot(2,2,1) 
plot(xli,n_e(jj,:),'-ro', 'Linewidth',2,'DisplayName','n_e');
hold on;
plot(xli,n_i(jj,:),'-ko', 'Linewidth',2,'DisplayName','n_i');
hold on;
plot(xli,n_lightimp(jj,:),'-bo', 'Linewidth',2,'DisplayName','n_{imp}');
xlabel('r/a');
ylabel('n (m^{-3})');
%xlim([0 0.9]);
legend('location', 'best');
set(gca,'FontSize',14);
grid on;

subplot(2,2,2) 
plot(xli,TeeV(jj,:),'-ro', 'Linewidth',2,'DisplayName','T_e');
hold on;
plot(xli,TieV(jj,:),'-ko', 'Linewidth',2,'DisplayName','T_i and T_z');
xlabel('r/a');
ylabel('T (eV)');
%xlim([0 0.9]);
legend('location', 'best');
set(gca,'FontSize',14);
grid on;

subplot(2,2,3) 
plot(xli,rlne(jj,:),'-ro', 'Linewidth',2,'DisplayName','el.');
hold on;
plot(xli,rlni(jj,:),'-ko', 'Linewidth',2,'DisplayName','main ion');
hold on;
plot(xli,rlnlightimp(jj,:),'-bo', 'Linewidth',2,'DisplayName','imp.');
xlabel('r/a');
ylabel('a/L_n');
%xlim([0 0.9]);
legend('location', 'best');
set(gca,'FontSize',14);
grid on;

subplot(2,2,4) 
plot(xli,rlte(jj,:),'-ro', 'Linewidth',2,'DisplayName','el.');
hold on;
plot(xli,rlti(jj,:),'-ko', 'Linewidth',2,'DisplayName','ions');
hold on;
xlabel('r/a');
ylabel('a/L_T');
%xlim([0 0.9]);
legend('location', 'best');
set(gca,'FontSize',14);
grid on;


% Create slider plot for gamma1 (growth rate for the mode 1)
hFig_gamma = figure;
hAx_gamma = axes('Parent', hFig_gamma);
figure(hFig_gamma)
plot(ky_vals(2, :), gamma1(2, :), 'rs-', 'linewidth', 2);
%xlim([0 30]);
%ylim([-0.1 1.5]);
grid on;
xlabel('ky \rho_s', 'FontSize', 14);
ylabel('\gamma (c_s/a)', 'FontSize', 14);
set(gca, 'FontSize', 14);
%title('Gamma vs ky at radial point 2 ');
title(sprintf('Gamma vs ky at radial point 2 (r/a = %.3f)', xli(2)));

% Add slider
hSlider_gamma = uicontrol('Style', 'slider', 'Min', 2, 'Max', length(xli), 'Value', 2, ...
    'Units', 'normalized', 'Position', [0.2 0.02 0.6 0.05], ...
    'SliderStep', [1/(length(xli)-1), 1/(length(xli)-1)], ...
    'Callback', @(src, event) updatePlot_gamma(src, hAx_gamma, ky_vals, gamma1, xli));

%function to update the plot
    function updatePlot_gamma(src, hAx_gamma, ky_vals, gamma1, xli)
        idx = round(get(src, 'Value'));
        if idx < 2
            idx = 2;
        elseif idx > length(xli)
            idx = length(xli);
        end
        
				% Update the plot
        plot(hAx_gamma, ky_vals(idx, :), gamma1(idx, :), 'rs-', 'linewidth', 2);
        %xlim(hAx, [0 30]);
        %ylim(hAx, [-0.1 1.5]);
        grid(hAx_gamma, 'on');
        xlabel(hAx_gamma, 'ky \rho_s', 'FontSize', 14);
        ylabel(hAx_gamma, '\gamma (c_s/a)', 'FontSize', 14);
        set(hAx_gamma, 'FontSize', 14);
	title(hAx_gamma, sprintf('Gamma vs ky at radial point %d (r/a = %.3f)', idx, xli(idx)));
    end


% Create slider plot for frequency1 (for the mode 1)
hFig_frequency = figure;
hAx_frequency = axes('Parent', hFig_frequency);
figure(hFig_frequency)
plot(ky_vals(2, :), frequency1(2, :), 'rs-', 'linewidth', 2);
%xlim([0 30]);
%ylim([-0.1 1.5]);
grid on;
xlabel('ky \rho_s', 'FontSize', 14);
ylabel('\omega (c_s/a)', 'FontSize', 14);
set(gca, 'FontSize', 14);
%title('frequency vs ky at radial point 2');
title(sprintf('Frequency vs ky at radial point 2 (r/a = %.3f)', xli(2)));

% Add slider
hSlider_frequency = uicontrol('Style', 'slider', 'Min', 2, 'Max', length(xli), 'Value', 2, ...
    'Units', 'normalized', 'Position', [0.2 0.02 0.6 0.05], ...
    'SliderStep', [1/(length(xli)-1), 1/(length(xli)-1)], ...
    'Callback', @(src, event) updatePlot_frequency(src, hAx_frequency, ky_vals, frequency1, xli));


    function updatePlot_frequency(src, hAx_frequency, ky_vals, frequency1, xli)
        idx = round(get(src, 'Value'));
        if idx < 2
            idx = 2;
        elseif idx > length(xli)
            idx = length(xli);
        end

        plot(hAx_frequency, ky_vals(idx, :), frequency1(idx, :), 'bo-', 'linewidth', 2);
        %xlim(hAx, [0 30]);
        %ylim(hAx, [-0.1 1.5]);
        grid(hAx_frequency, 'on');
        xlabel(hAx_frequency, 'ky \rho_s', 'FontSize', 14);
        ylabel(hAx_frequency, '\omega (c_s/a)', 'FontSize', 14);
        set(hAx_frequency, 'FontSize', 14);
	title(hAx_frequency, sprintf('Frequency vs ky at radial point %d (r/a = %.3f)', idx, xli(idx)));
    end

% Create GIF for gamma1 vs ky across the radial positions
gif_filename = fullfile(path_plots, 'Gamma1_vs_ky.gif');

gamma_dir = fullfile(path_plots, 'Gamma1_Plots_FromGIF');
unix(sprintf(['mkdir -p ' gamma_dir]));  	  			    

delay_time = 1.2; % Time between frames in seconds
fig_gif = figure('Visible','off'); % Optional: set 'on' to watch as it's created
				 
for idx = 2:length(xli)
  plot(ky_vals(idx, :), gamma1(idx, :), 'rs-', 'LineWidth', 2);
  
    xlabel('ky \rho_s', 'FontSize', 14);
    ylabel('gamma (c_s/a)', 'FontSize', 14);
    title(sprintf('Gamma vs ky at radial point %d (r/a = %.3f)', idx, xli(idx)));
    grid on;
    set(gca, 'FontSize', 14);

	
    % Capture frame as image
    frame = getframe(fig_gif);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write to the GIF File
    if idx == 2
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end

close(fig_gif); % Close the figure used for GIF
disp(['GIF saved at: ' gif_filename]);

% Create GIF for the frequency1 vs ky accross the radial positions
gif_filename_freq = fullfile(path_plots, 'Frequency1_vs_ky.gif');

freq_dir = fullfile(path_plots, 'Frequency1_Plots_FromGIF');
unix(sprintf(['mkdir -p ' freq_dir]));  	  			      

delay_time = 1.2; % seconds between frames

fig_gif_freq = figure('Visible','off'); % Set to 'on' if you want to watch it live
for idx = 2:length(xli)
    plot(ky_vals(idx, :), frequency1(idx, :), 'bo-', 'LineWidth', 2);
    xlabel('ky \rho_s', 'FontSize', 14);
    ylabel('\omega (c_s/a)', 'FontSize', 14);
    title(sprintf('Frequency vs ky at radial point %d (r/a = %.3f)', idx, xli(idx)));
    grid on;
    set(gca, 'FontSize', 14);

    % Capture frame
    frame = getframe(fig_gif_freq);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write to GIF
    if idx == 2
        imwrite(imind, cm, gif_filename_freq, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
    else
        imwrite(imind, cm, gif_filename_freq, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end

close(fig_gif_freq); % Close the figure used for GIF
disp(['GIF saved at: ' gif_filename_freq]);


end


