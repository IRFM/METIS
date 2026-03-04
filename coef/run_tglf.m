function [Chie,Chii,De,Ve] = run_tglf(tglfparams,post,itime,iplt);

if ~exist('iplt')
  iplt=0;
end

% To do:
% Force consistency with different switches like sat2 and xnumodel=3

%%%% Initialisations
Chie = zeros(1,length(post.profil0d.xli));
Chii = zeros(1,length(post.profil0d.xli));
De = zeros(1,length(post.profil0d.xli));
Di = zeros(1,length(post.profil0d.xli));
Ve = zeros(1,length(post.profil0d.xli));

Gelec = zeros(1,length(post.profil0d.xli));
Gion1 = zeros(1,length(post.profil0d.xli));
Gion2 = zeros(1,length(post.profil0d.xli));

Qelec = zeros(1,length(post.profil0d.xli));
Qion1 = zeros(1,length(post.profil0d.xli));
Qion2 = zeros(1,length(post.profil0d.xli));

% Path to TGLF routines, need to provide something for a TGLF install function
path_to_METIS_coef = fullfile(fileparts(which('install_QLKNN10D')))
path_to_METIS_TGLF = [path_to_METIS_coef '/TGLF_TOOLS'];
addpath(path_to_METIS_TGLF);

persistent path_runs_root

if isempty(path_runs_root)
    % Base "runs" directory (shared root)
    if isempty(getappdata(0,'local'))
      base_runs = fullfile(path_to_METIS_TGLF, 'runs_tglf');    
    else
      base_runs = fullfile(getappdata(0,'local'), 'runs_tglf');    
    end  

    unix(sprintf(['mkdir -p ' base_runs]));

    % Unique ID based on date & time for this MATLAB session
    timestamp = datestr(now,'yyyymmdd_HHMMSS');  % e.g. 20251118_231530123   HHMMSSFFF')
    path_runs_root = fullfile(base_runs, ['run_' timestamp]);
    
    unix(sprintf(['mkdir -p ' path_runs_root]));    

    fprintf('====> TGLF stored in run directory %s\n', path_runs_root);
end

%param = declare_TGLF_op_parameters;
%tglfparams  = param.valeur;

%%%%% Read input TGLF template  using>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
tglfin.ADIABATIC_ELEC 	    = tglfparams.ADIABATIC_ELEC; 
tglfin.ALPHA_MACH 	    = tglfparams.ALPHA_MACH;
tglfin.ALPHA_E 		    = tglfparams.ALPHA_E;
tglfin.ALPHA_P 		    = tglfparams.ALPHA_P;
tglfin.ALPHA_QUENCH	    = tglfparams.ALPHA_QUENCH;
tglfin.XNU_FACTOR	    = tglfparams.XNU_FACTOR;
tglfin.DEBYE_FACTOR	    = tglfparams.DEBYE_FACTOR;
tglfin.ETG_FACTOR 	    = tglfparams.ETG_FACTOR;
tglfin.WRITE_WAVEFUNCTION_FLAG= tglfparams.WRITE_WAVEFUNCTION_FLAG;

% Add some consistency checks e.g. XNU_MODEL and SAT_RULE


%%%%% PHYSICS CONSTANTS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
phys.e   =   1.602176462e-19;
phys.mD  =  3.345e-27;                 % mass of deuterium
phys.mu0 =   4*pi*1e-7;                % permeability in vacuum (H/m) (definition)
phys.kB  =   1.3806503e-23;            % Boltzmann Constant (J/K)
phys.ua  =   1.66053873e-27;           % in kgs 


%%%%% Call the profiles from METIS:>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% from profil0d in post (output of METIS):
xli	     = post.profil0d.xli;   	% normalised radius (1,nr=21)
time_profil0d= post.profil0d.temps; 	% (nt1,1)

%from z0dinput in post, Geometrical parameters (Inputs for METIS):  
minor_radius = post.z0dinput.geo.a;  	  % (nt0,1) minor radius of the plasma waveform (m).
major_radius = post.z0dinput.geo.R;  	  % (nt0,1) major radius of the plasma waveform (m).
vertical_pos = post.z0dinput.geo.R*0; 	  % (nt0,1) vertical position of the plasma waveform (m) set to 0 here.
Btoroidal    = post.z0dinput.geo.b0;	  % (nt0,1) vacuum magnetic toroidal field measured at geo.R(T)
rimp         = post.z0dinput.option.rimp; % ratio of density b/w main impurity and another impurity
main_ion     = post.z0dinput.option.gaz;  % main ion 1>H, 2>D, 3>T, 4>He
zlightimp    = post.z0dinput.option.zimp; % charge of other impurity

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
n_lightimp = post.profil0d.nzp;     		     % (nt1,nr) density of light impurity (m^-3)
vtor       = post.profil0d.vtor;      	             % (nt1,nr) Toroidal rotation speed (m/s)
ptotal     = post.profil0d.ptot; 		     % (nt1,nr) Total pressure profile (Pa)
Phi 	   = post.profil0d.phi;     	             % (nt1,nr) Toroidal magnetic flux (Wb)
Psi 	   = post.profil0d.psi;    	             % (nt1,nr) Poloidal  magnetic flux (Wb)
Bpol 	   = post.profil0d.bpol;
xie	   = post.profil0d.xie;			     % electron diffusivity estimation using power balance	
xii	   = post.profil0d.xii;			     % ion diffusivity estimation using power balance	

% from zerod in post (Time dependent Scaler Data given in output file):
time_zerod   = post.zerod.temps; % (nt0,1)


%%%%% Interpolation and Reshaping on profiles 1D time grid.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
dvtordxli 	= pdederive(xli,(vtor./R),0,2,2,1); 	%(nt1,nr) d(toroidal angular rotation frequency)/dxli, xli = r*a
dqdxli    	= pdederive(xli,q,0,2,2,1);   		% (nt1,nr) d(safety factor)/dxli
dkdxli 	  	= pdederive(xli,kappa,0,2,2,1); 	% (nt1,nr) d(kappa)/dxli
dddxli    	= pdederive(xli,delta,0,2,2,1); 	% (nt1,nr) d(delta)/dxli
dpdxli    	= pdederive(xli,ptotal,0,2,2,1); 	% (nt1,nr)
dPsidxli  	= pdederive(xli,Psi,0,2,2,1); 

%%%%% Some physics quantities>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
SIGN_Itorr = tglfin.SIGN_IT;		%1 ;			 % Itorr is the plasma current, positive for counter-clockwise 
Cs	   = sqrt(Te./phys.mD); 			 % (nt1,nr)deuterium sound speed
B_unit	   = abs((q ./ (xli_rep.* a0.*a0)) .* dPsidxli); % (nt1,nr) psi is poloidal magnetic field.
Omegai	   = phys.e * B_unit ./phys.mD;  		 % (nt1,nr)
rhol_ref   = Cs./Omegai;  				 % (nt1,nr)

%%%%%% Collisionality and Betae: 
lnlamda = 24.0 - (0.5 .* log(n_e ./ 1e6)) + log(TeeV);       
taue	= ((3.44e5) .* ((TeeV).^1.5))./((n_e ./ 1e6) .* lnlamda);  	%  sec
xnuei   = 0.75 .* sqrt(pi)./taue ;         				% 1/sec 
xnuei   = xnuei .* a0 ./ Cs ;            				% normalized electron-ion collision rate
betae   = (8.0 * pi .* (n_e./ 1e6) .* (Te*1e7)) ./(B_unit.^2 .*1e8); 		% electron beta unit

%%%%% SPECIES VECTORS (Cant be changed by the user.)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
Vpar_shear_e= (- SIGN_Itorr .* R_maj .* (dvtordxli./a0)) ./ Cs; 	% VPAR_SHEAR_1: Normalised parallel velocity gradient. 


%%%%% Calculate zeff
zeff        = nine.*tglfin.ZS_2.^2 + nlightimpne.*charge_lightimp.^2 ;

%%%%% MILLER GEOMETRY PARAMETERS: changes with space r and time t>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rmin_loc    = xli_rep;                    					% (nt1,nr) Flux surface centroid minor radius r/a 
rmaj_loc    = R_maj./a0;                 					%(nt1,nr) Flux surface centroid major radius R0/a0.
zmaj_loc    = z0./a0;                     					%(nt1,nr) Flux surface centroid elevation Zmaj/a0
q_loc	    =  abs(q);                       					% (nt1,nr) Absolute value of the safety factor
q_prime_loc = (q ./rmin_loc).* dqdxli; 						%(nt1,nr) 
p_prime_loc = (phys.mu0/(4*pi)) .* (q./ (rmin_loc .* (B_unit.^2))) .* dpdxli;
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

%%%%% Passing the value of time chose by the user
jj = itime;

%%%%% Species with constant charges but types are varying from METIS.
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

%%%%% Loop over the radial domain

for ii = 2:length(xli)   % doing a loop over radial domain

    disp(['========================= RUNS TGLF r/a=' num2str(xli(ii)) '=========================']);

    % Create the directory for this radial point     
    tglf_input_scan_radial = fullfile(path_runs_root,['tglf_tmp_radial_' num2str(ii)])
    unix(sprintf(['mkdir -p ' tglf_input_scan_radial]));        

    % >>>>>>>>>>>>>>>>>>>>>>> rlne scan parameters >>>>>>>>>>>>>>>>>>>>>>>>
    rln_offset = 0.01;  		% small perturbation to rlne
    rlne_variations = [0, rln_offset];  

    rlne_values = [];
    elec_flux_values = [];

    % Loop over the rlne variations
    rln_count=0;
    for offset = rlne_variations
        rln_count=rln_count+1;
        modified_rlne = rlne(jj, ii) + offset;
	tglf_input_scan_rln = [tglf_input_scan_radial '/rlne_' num2str(rln_count)];
        unix(sprintf(['mkdir -p ' tglf_input_scan_rln]));        	

        % Update TGLF input for this rlne
        tglfin.ZS_1 = -1;    	 
        tglfin.MASS_1 = 2.723100e-04;
        tglfin.RLNS_1 = modified_rlne;
        tglfin.RLTS_1 = rlte(jj,ii);
        tglfin.TAUS_1 = TeTe(jj,ii);
        tglfin.AS_1 = nene(jj,ii);
        tglfin.VPAR_1 = Vpar_e(jj,ii);
        tglfin.VPAR_SHEAR_1 = Vpar_shear_e(jj,ii);

        % MAIN IONS
        tglfin.RLNS_2 = rlni(jj,ii);
        tglfin.RLTS_2 = rlti(jj,ii);
        tglfin.TAUS_2 = TiTe(jj,ii);
        tglfin.AS_2 = nine(jj,ii);
        tglfin.VPAR_2 = Vpar_e(jj,ii);
        tglfin.VPAR_SHEAR_2 = Vpar_shear_e(jj,ii);

        % LIGHT IMPURITY
        tglfin.ZS_3 = charge_lightimp(jj,ii); 
        tglfin.MASS_3 = mass_lightimp;	
        tglfin.RLNS_3 = rlnlightimp(jj,ii);
        tglfin.RLTS_3 = rlti(jj,ii);
        tglfin.TAUS_3 = TzTe(jj,ii);
        tglfin.AS_3 = nlightimpne(jj,ii);
        tglfin.VPAR_3 = Vpar_e(jj,ii);	
        tglfin.VPAR_SHEAR_3 = Vpar_shear_e(jj,ii);

        % Geometry and other plasma parameters 
        tglfin.BETAE = betae(jj,ii);
        tglfin.XNUE = xnuei(jj,ii);
        tglfin.ZEFF = zeff(jj,ii);
        tglfin.RMIN_LOC = rmin_loc(jj,ii);
        tglfin.RMAJ_LOC = rmaj_loc(jj,ii);
        tglfin.ZMAJ_LOC = zmaj_loc(jj,ii); 
        tglfin.Q_LOC = q_loc(jj,ii);
        tglfin.Q_PRIME_LOC = q_prime_loc(jj,ii);
        tglfin.P_PRIME_LOC = p_prime_loc(jj,ii);
        tglfin.DRMINDX_LOC = drmindx_loc; 
        tglfin.DRMAJDX_LOC = drmajdx_loc(jj,ii);
        tglfin.DZMAJDX_LOC = dzmajdx_loc;
        tglfin.KAPPA_LOC = kappa_loc(jj,ii);
        tglfin.S_KAPPA_LOC = s_kappa_loc(jj,ii); 
        tglfin.DELTA_LOC = delta_loc(jj,ii);
        tglfin.S_DELTA_LOC = s_delta_loc(jj,ii);
        tglfin.ZETA_LOC = zeta_loc;
        tglfin.S_ZETA_LOC = s_zeta_loc; 
        tglfin.KX0_LOC = kx0_loc; 

        make_tglfinput(tglf_input_scan_rln, tglfin);
	
        cd(tglf_input_scan_rln);
	!rm -f out.tglf*
        disp([ '===TGLF for rlne=' num2str(modified_rlne)]);

	[sout,tout] = unix(sprintf(['module purge; module load gcc/13.2.0 lapack/3.9.1-gcc-13.2.0 openmpi/4.1.1-gcc-13.2.0; export GACODE_PLATFORM=IRFM; export GACODE_ROOT=/ZONE_TRAVAIL/PM258588/codes/gacode; . ${GACODE_ROOT}/shared/bin/gacode_setup; nice -n 9 tglf -e ./ -n ' num2str(round(tglfparams.numproc))]))

        disp('=== TGLF RUN COMPLETE ===');
	try
          tglfout = read_tglf_outrun_3sp(tglf_input_scan_rln);
	catch
	  tglfout = read_tglf_outrun_gbflux_3sp(tglf_input_scan_rln);
	end
	
        % Compute particle flux (SI units) at this rlne
        Gelec_current = tglfout.Ge_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Cs(jj,ii);

        rlne_values = [rlne_values, modified_rlne];
        elec_flux_values = [elec_flux_values, Gelec_current];

        % Store heat fluxes for chie, chii (only for original rlne)
        if offset == 0
            Gelec(ii) = Gelec_current;

            Qelec(ii) = tglfout.Qe_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Te(jj,ii) * Cs(jj,ii);
            Qion1(ii) = tglfout.Q1_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Te(jj,ii) * Cs(jj,ii);
	    Qion2(ii) = tglfout.Q2_gB * (rhol_ref(jj,ii)/a0(jj,ii))^2 * n_e(jj,ii) * Te(jj,ii) * Cs(jj,ii);
        end

    end % rlne scan loop

    % >>>>>>>>>>>>>>>>>>>>>>> Fit Gelec vs rlne >>>>>>>>>>>>>>>>>>>>>>>>
    elec_coeffs = polyfit(rlne_values, elec_flux_values, 1); % [slope, intercept]
    De(ii) = elec_coeffs(1) * a0(jj,ii) / n_e(jj,ii);   % slope ----> diffusion
    Ve(ii) = elec_coeffs(2) / n_e(jj,ii);               % intercept ----> pinch

    % Heat diffusivities (for reference)
    %Chie(ii) = -Qelec(ii) / n_e(jj,ii) / (dtedxli(jj,ii) / a0(jj,ii));
    %Chii(ii) = -Qion1(ii) / n_i(jj,ii) / (dtidxli(jj,ii) / a0(jj,ii));

end % radial loop

    	%delete(fullfile(tglf_input_scan, '*')); % Delete all files in this radial position directory
	%rmdir(tglf_input_scan); 		% Remove the subdirectory itself
	
Chie = -Qelec./n_e(jj,:)./(dtedxli(jj,:)./a0(jj,:));
Chii = -Qion1./n_i(jj,:)./(dtidxli(jj,:)./a0(jj,:));

% Magnetic axis
Chie(1)=0;
Chii(1)=0;

% Put the total flux in diffusive part
%De = - Gelec./dnedxli(jj,:)./a0(jj,:); 
%Di = - Gion1./dnidxli(jj,:)./a0(jj,:); 

if ~iplt
  return
end

%{
% Create slider plot for gamma1 (growth rate for the mode 1)
%hFig_gamma = figure;
%hAx_gamma = axes('Parent', hFig_gamma);
%plot(ky_vals(2, :), gamma1(2, :), 'rs-', 'linewidth', 2);
%xlim([0 30]);
%ylim([-0.1 1.5]);
%grid on;
%xlabel('ky', 'FontSize', 20);
%ylabel('gamma', 'FontSize', 20);
%set(gca, 'FontSize', 15);
%title('Gamma vs ky at radial point 2');

% Add slider
hSlider_gamma = uicontrol('Style', 'slider', 'Min', 2, 'Max', length(xli), 'Value', 2, ...
    'Units', 'normalized', 'Position', [0.2 0.02 0.6 0.05], ...
    'SliderStep', [1/(length(xli)-1), 1/(length(xli)-1)], ...
    'Callback', @(src, event) updatePlot_gamma(src, hAx_gamma, ky_vals, gamma1, xli));

    function updatePlot_gamma(src, hAx_gamma, ky_vals, gamma1, xli)
        idx = round(get(src, 'Value'));
        if idx < 2
            idx = 2;
        elseif idx > length(xli)
            idx = length(xli);
        end
        
        plot(hAx_gamma, ky_vals(idx, :), gamma1(idx, :), 'rs-', 'linewidth', 2);
        %xlim(hAx, [0 30]);
        %ylim(hAx, [-0.1 1.5]);
        grid(hAx_gamma, 'on');
        xlabel(hAx_gamma, 'ky', 'FontSize', 20);
        ylabel(hAx_gamma, 'gamma1', 'FontSize', 20);
        set(hAx_gamma, 'FontSize', 15);
        title(hAx_gamma, sprintf('Gamma1 vs ky at radial point %d', idx));
    end


% Create slider plot for frequency1 (for the mode 1)
hFig_frequency = figure;
hAx_frequency = axes('Parent', hFig_frequency);
plot(ky_vals(2, :), frequency1(2, :), 'rs-', 'linewidth', 2);
%xlim([0 30]);
%ylim([-0.1 1.5]);
grid on;
xlabel('ky', 'FontSize', 20);
ylabel('frequency', 'FontSize', 20);
set(gca, 'FontSize', 15);
title('frequency vs ky at radial point 2');

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
        
        plot(hAx_frequency, ky_vals(idx, :), frequency1(idx, :), 'rs-', 'linewidth', 2);
        %xlim(hAx, [0 30]);
        %ylim(hAx, [-0.1 1.5]);
        grid(hAx_frequency, 'on');
        xlabel(hAx_frequency, 'ky', 'FontSize', 20);
        ylabel(hAx_frequency, 'frequency1', 'FontSize', 20);
        set(hAx_frequency, 'FontSize', 15);
        title(hAx_frequency, sprintf('frequency1 vs ky at radial point %d', idx));
    end



figure
subplot(2,2,1)
hold on;
plot(xli, Qelec,'rs-','linewidth',2.5);
plot(xli, Qion1, 'bo-', 'LineWidth', 2.5) ;
plot(xli, Qion2, 'g*-', 'LineWidth', 2.5) ;
grid on;
hold off;
xlabel('r/a', 'FontSize', 20);
ylabel('Q','FontSize', 20);
set(gca,'FontSize',15);
legend('elec', 'main ion','light imp');

subplot(2,2,2)
hold on;
plot(xli,Gelec,'rs-','linewidth',2);
plot(xli, Gion1, 'bo-', 'LineWidth', 2.5) ;
plot(xli, Gion2, 'g*-', 'LineWidth', 2.5) ;
grid on;
hold off;
xlabel('r/a', 'FontSize', 20);
ylabel('G', 'FontSize', 20);
set(gca,'FontSize',15);
legend('elec', 'main ion','light imp');


figure;
subplot(2,2,1)
hold on;
plot(xli,Chie,'rs-','linewidth',2.5)
plot(xli, Chii, 'bo-', 'LineWidth', 2.5) 
plot(xli, post.profil0d.xie(itime,:), 'r--', 'LineWidth', 2.5)
plot(xli, post.profil0d.xii(itime,:), 'b--', 'LineWidth', 2.5) 
grid on
hold off;
xlabel('r/a', 'FontSize', 20)
ylabel('Chi', 'FontSize', 20)
set(gca,'FontSize',15)
legend('elec(TGLF)', 'main ion(TGLF)','elec(PB)', 'main ion(PB)')

subplot(2,2,2)                      
hold on;
plot(xli,De,'rs-','linewidth',2.5)
%plot(xli, Di, 'bo-', 'LineWidth', 2.5) 
grid on
hold off;
xlabel('r/a', 'FontSize', 20)
ylabel('D', 'FontSize', 20)
set(gca,'FontSize',15)
legend('elec')

subplot(2,2,3)                      
hold on;
plot(xli,Ve,'rs-','linewidth',2.5)
%plot(xli, Di, 'bo-', 'LineWidth', 2.5) 
grid on
hold off;
xlabel('r/a', 'FontSize', 20)
ylabel('V', 'FontSize', 20)
set(gca,'FontSize',15)
legend('elec')
%}






