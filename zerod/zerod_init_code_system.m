% Read data from code system results.
function z0dinput = zerod_init_code_system(mode_exp,shot,gaz,temps,z0dinput);

[file,path]=uigetfile('*.*','Code System : choose a file ?');
drawnow
if ~ischar(file)
	disp('Operation cancelled');
	z0dinput = [];
	return
end
filename = strcat(path,file);

% test type of code system
text = fileread(filename);
if ~isempty(strfind(text,'* PROCESS *'))
    % case PROCESS filename
    disp('Reading Process results file');
    % template for METIS parameters (to tune the composition)
    z0dinput = sycomore2metis;
    parameters_filename = z0dinput.option
    % name of Sycomore result file
    process_filename =filename;
    % computation mode (1 = fast for testing; 0 = full)
    computation_mode = 0;
    % reading sycomore results
    data = read_process_data(process_filename);

    % impurities
    liste = {
    {'Ag',47,107.8682 },
    {'Al',13,26.98154 },
    {'Au',79,196.96654},
    {'B',  5,10.81   },
    {'Be', 4,9.01    },
    {'C', 6,12.011   },
    {'Cd',48,112.41},
    {'Co',27,58.93},
    {'Cr',24,51.996},
    {'Cu',29,63.546},
    {'Fe',26,55.845},
    {'Ge',32,72.61},
    {'Hf',72,178.49},
    {'Ir',77,192.217},
    {'Mn',25,54.94 },
    {'Mo',42,95.94 },
    {'Nb',41,92.90638},
    {'Ni',28,58.6934},
    {'Os',76,190.2  },
    {'Pb',82,207.2},
    {'Pd',46,106.42},
    {'Pt',78,195.08},
    {'Re',75,186.21},
    {'Rh',45,102.9055},
    {'Ru',44,101.07},
    {'Si',14,28.0855},
    {'Sn',50,118.69},
    {'Th',90,232.04},
    {'Ti',22,47.9},
    {'Ta',73,180.9479},
    {'U', 92,238.03 },
    {'V', 23,50.9415},
    {'W', 74,183.84},
    {'Zn',30,65.38},
    {'Zr',40,91.22},
    {'Ar', 18,39.95},
    {'D',   1,2},
    {'T',   1,3},
    {'Ga', 31,69.75},
    {'H',   1,1},
    {'He',  2,4},
    {'He3', 2,3},
    {'Kr', 36,83.8},
    {'Li',  3,6.94},
    {'Na', 11,22.99},
    {'Ne', 10,20.18},
    {'O',   8,16.00},
    {'Xe', 54,131.3},
    {'N',7,14},
    };
    for k = 1:length(liste)
	name_liste{k} = liste{k}{1};
	num_liste(k)  = liste{k}{2};
	a_liste(k)  = liste{k}{3};
    end

    A = [];
    Z = [];
    label= {};
    fraction = [];
    for k=1:20
        name = sprintf('fimp_%2.2d',k);
        if isfield(data.plasma,name)
            ind = strmatch(data.plasma.(name).unit,name_liste,'exact');
            if isempty(ind)
                if data.plasma.(name).comment(3) == '_'
                    data.plasma.(name).unit = data.plasma.(name).comment(2);
                else
                    data.plasma.(name).unit = data.plasma.(name).comment(2:3);
                end
                ind = strmatch(data.plasma.(name).unit,name_liste,'exact');
            end
            if isempty(ind)
                if data.plasma.(name).line(3) == '_'
                    data.plasma.(name).unit = data.plasma.(name).line(2);
                else
                    data.plasma.(name).unit = data.plasma.(name).line(2:3);
                end
                ind = strmatch(data.plasma.(name).unit,name_liste,'exact');
            end
            if data.plasma.(name).value ~= 0
                label{end+1} = name_liste{ind};
                Z(end+1) = num_liste(ind);
                A(end+1) = round(a_liste(ind));
                fraction(end+1) = data.plasma.(name).value;
            end
        end
    end

    % mapping the data
    sycomore_data.ip = data.plasma.plascuro1D6.value .* 1e6;                % flat top plasma current (A).
    sycomore_data.q95 = data.plasma.q95.value;                % flat top safety factor
    sycomore_data.rb0 = data.plasma.bt.value .* data.plasma.rmajor.value;         % magnetic rigidity (flat top)
    sycomore_data.available_flux = abs(data.volt_second.vstot.value);     % Poloidal available flux (Wb)
    sycomore_data.device = sprintf('DEMO-R%d-a%d-RBt%d-Ip%d',ceil(data.plasma.rmajor.value*100),ceil(data.plasma.rminor.value*100),ceil(sycomore_data.rb0),ceil(sycomore_data.ip/1e6)); % device name (used to generate file name and comments)
    sycomore_data.scaling  = 0;             % code for energy plasma content scaling law (same as METIS)
    sycomore_data.H_H  = data.plasma.hfact.value;                 % enhancement factor for energy content on flat top
    sycomore_data.f_Greenwald = data.plasma.dnla_gw.value;       % Greenwald density fraction on flat top
    if isfield(sycomore_data,'shot')
      sycomore_data.shot = sycomore_data.shot;                 % shot number, will be used to write data in UAL
      sycomore_data.run  = sycomore_data.run + 1;                 % run  number, will be used to write data in UAL
    else
      sycomore_data.shot = ceil(rand*1e4);                 % shot number, will be used to write data in UAL
      sycomore_data.run  = 1;                 % run  number, will be used to write data in UAL
    end
    if isfield(data.plasma,'faccd_')
        sycomore_data.f_ni = data.plasma.faccd_.value + data.plasma.bootipf_.value;                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
    else
        sycomore_data.f_ni = data.hcd.faccd.value + data.hcd.bootipf.value;                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
    end
    if sycomore_data.f_ni  > 0.99
	sycomore_data.f_ni = 1;
    end
    sycomore_data.zeff      = data.plasma.zeff.value - 4 .* sum(fraction(Z == 2));          % line averaged Zeff without He ashes contribution
    ind_impur = find(Z>2);
    ind_imp   = ind_impur(fraction(ind_impur) == max(fraction(ind_impur)));
    ind_imp   = ind_imp(1);
    ind_impur(ind_impur == ind_imp) = [];
    ind_max   = ind_impur(fraction(ind_impur) == max(fraction(ind_impur)));
    ind_max   = ind_max(1);
    if Z(ind_max) == 74
      % computation of rimp
      %f4 = ((sycomore_data.zeff - 1) - fraction(ind_imp) .* (Z(ind_imp).^ 2 - Z(ind_imp)) - fraction(ind_max) .* (Z(ind_max).^ 2 - Z(ind_max))) ./ 12;
      parameters_filename.zimp = Z(ind_imp);
      parameters_filename.zmax = Z(ind_imp);
      sycomore_data.rimp = 0.1;               % ratio  between berillium and argon in core plasma
      sycomore_data.cW               = fraction(ind_max);  % tungsten concentration in core plasma   , must be provide by Sycomore 
      fzmax = fraction(ind_imp);
    else
      parameters_filename.zimp = Z(ind_imp);
      parameters_filename.zmax = Z(ind_max);
      sycomore_data.rimp = fraction(ind_max) ./ fraction(ind_imp);               % ratio  between berillium and argon in core plasma
      fzmax = fraction(ind_max);  
      sycomore_data.cW               = 1e-5;  % tungsten concentration in core plasma   , must be provide by Sycomore 
    end
    sycomore_data.Recycling = 0.99;         % recycling  @ divertor 
    sycomore_data.tau_He_o_tau_E_core = 3;  % ratio between core confinement time of He over energy confinement time.
    sycomore_data.rw = 0.7;                 % cyclotron radiation reflection coefficient
    if isfield(data.hcd,'pnbeam')
        sycomore_data.PNBI = data.hcd.pnbeam.value .* 1e6;              % NBI power used in Sycomore simulation (W, for ITER like default  case contains also ICRH power, if PNBI <  0, scale on ITER power).
    else
        sycomore_data.PNBI = 0;              % NBI power used in Sycomore simulation (W, for ITER like default  case contains also ICRH power, if PNBI <  0, scale on ITER power).        
    end
    if isfield(data.hcd,'echpwr')
        sycomore_data.PECRH = data.hcd.echpwr.value .* 1e6;   % ECRH power 
    else
        sycomore_data.PECRH = 0;             
    end
    sycomore_data.fR_target = 1;            % position of outer target in unit of R0.
    if isfield(data.hcd,'enbeam')
        sycomore_data.nbi_e_inj = data.hcd.enbeam.value .* 1e3;          % neutral beam energy injection.
    else
        sycomore_data.nbi_e_inj = 2e6;          % neutral beam energy injection.
    end
    sycomore_data.signe = 1;                % sign of Bt . J_phi: must be provided by the machine description
    sycomore_data.flux_expansion = 3;       % divertor flux expansion 
    sycomore_data.target_orientation = 90;  % poloidal orientation of outer divertor target  (degrees)
    sycomore_data.S_factor = 3e-3;          % divertor flux spreading factor (m) 
    % There no enhancement of impurities concentration id divertor. have been checked  in soldiv module of Sycomore.
    if isfield(data,'divertor')
         sycomore_data.fzmax_div        = 100 .* max(0,(data.divertor.zeffso.value -  data.plasma.zeff.value) ./ parameters_filename.zmax .^ 2 - fzmax);     % Argon concentration enhancement in divertor (%)  , must be provide by Sycomore
         sycomore_data.f_DSOL             = data.divertor.delw.value./3e-3;   % multiplicator applied to the SOL width when it is defined by a scaling law, for H mode only (Dsol_Hmode = factor_scale * Goldston scaling).
    else
         sycomore_data.fzmax_div        = 0;    
        sycomore_data.f_DSOL            = 0;   % multiplicator applied to the SOL width when it is defined by a scaling law, for H mode only (Dsol_Hmode = factor_scale * Goldston scaling).
    end
    sycomore_data.f_ne_LCFS          = 2;   % factor applied to edge scaling law for density:\nif > 0, ne_edge =  nea_factor * LCFS_denstity_scaling_law;\nif < 0,  ne_edge =  abs(nea_factor) * n_bar'electron density at LCFS (m^-3)  , must be provide by Sycomore
    sycomore_data.couplage         = 0.15   % coupling coefficient given the ratio between poloidal flux consumption in the plasam and in the central selenoid (CS) during breakdown to take into account dissipated flux in passive structure.
    sycomore_data.duration         = 7200;  % estimated duration of the shot in s
    % use constant P/R for allowed shine througth (and first orbit losses)
    % take < 5% of 33 MW for ITER @ 6.2m -> 2e5 W/m
    sycomore_data.shine_through_limit =  1 .*  sycomore_data.PNBI/33;  % limit of power lost in shine through (W)

    %  optionnal fields 
    sycomore_data.tokamak  = '';           % new UAL tokamak database name
    sycomore_data.user     = '';           % new UAL user database name
    sycomore_data.dataversion = '';        % selected a differente version of data (not the last one)
    sycomore_data.occurrence = '';         % input cpo scenario occurrence in Kepler (default = [])
    sycomore_data.path4metis_files = '';   % path for files save as a postprocessing of METIS (METIS data, figures, ...); if left empty, then no wrinting       
 
    %  structure details for sepa_option (with data example for ITER):
    K = data.plasma.kappa.value;
    d = data.plasma.triang.value;
    Kref = (1.687 + 2.001) / 2;
    dref = (0.466 + 0.568) / 2;
    sepa_option.rxup      = 0.466 .* d ./ dref;     % upper triangularity (minor radius unit)
    sepa_option.zxup      = 1.687 .* K / Kref;    % upper altitude X point (minor radius unit)
    sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
    sepa_option.ra        = data.plasma.rmajor.value;       % major radius R0 (m) [6.2]
    sepa_option.za        = 0;       % altitude of the magnetic axis (m) [0.9]
    sepa_option.a         = data.plasma.rminor.value;         % minor radius (m) [2]
    sepa_option.rxdo      = 0.568 .* d ./ dref;     % lower triangularity (minor radius unit)
    sepa_option.zxdo      = 2.001 .* K / Kref;       % lower altitude X point (minor radius unit)
    sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
    sepa_option.b0        = data.tf_coils.bmaxtf.value;      % magnetic field at R0
    if isfield(data.radial_build,'TFcoilinboardleg')
        sepa_option.delta     = data.plasma.rmajor.value - data.plasma.rminor.value - data.radial_build.TFcoilinboardleg.value(2,1);      % magnetic field at R0
    else
        sepa_option.delta     = data.plasma.rmajor.value - data.plasma.rminor.value - data.radial_build.tfcth.value(2);      % magnetic field at R0
    end
    % metis call
    z0dinput = sycomore2metis(sepa_option,parameters_filename,sycomore_data,computation_mode);
    
else
    % case SYCOMORE
    disp('Reading Sycomore results file');
    % name of Sycomore result file
    sycomore_filename = filename;
    % parameters
    parameters_filename =[]; % to use default tuning
    computation_mode = 1;
    % reading sycomore results
    data = read_sycomore_output_file(sycomore_filename);
    % mapping the data
    sycomore_data.ip = data.Plasmacurrent .* 1e6;                % flat top plasma current (A).
    sycomore_data.q95 = data.q95;                % flat top safety factor
    sycomore_data.rb0 = data.Btonaxis .* data.Majorradius;         % magnetic rigidity (flat top)
    mu0 = 4 .* pi .* 1e-7;
    if ~isfield(data,'flux_CS') || ~isfinite(data.flux_CS) || (data.flux_CS == 0) 
    disp('=========> available_flux not provided by SYCOMORE (CS contribution)');
    % from Jean-Luc Paper
    acs = (data.OuterradiusofCS - data.InnerradiusofCS) ./ data.Numberoflayers
    Ics = data.Coilcurrent
    jcs = Ics  ./ acs .^ 2 
    bmax_cs = mu0 .* (data.OuterradiusofCS -  data.InnerradiusofCS) .* jcs
    flux_cs = 2 .* pi ./ 3 .* bmax_cs .* (data.OuterradiusofCS .^ 2 + data.InnerradiusofCS .^ 2 + data.OuterradiusofCS .* data.InnerradiusofCS)
    data.flux_CS = flux_cs;
    end
    if ~isfield(data,'flux_bvert') || ~isfinite(data.flux_bvert) || (data.flux_bvert == 0) 
      disp('=========> available_flux not provided by SYCOMORE (Bvert contribution)');
      %flux du champ vertical
      if ~isfield(data,'Fractionnon_inductivecurrent');
 	  li = 0.85;     
      elseif isnumeric(data.Fractionnon_inductivecurrent)
	  li = 1 - 0.4 .* data.Fractionnon_inductivecurrent;
      else
	  li = 0.85;
      end
      % calcul Wesson p 120-123.
      bv =  mu0 .* data.Plasmacurrent .* 1e6 ./ (4 .* pi .* data.Majorradius) .* (log(8 .* data.Majorradius ./((data.Poloidalcrosssection/pi) .^ 0.5)) +  ... 
	    data.NormalizedBeta ./ data.NormalizedthermalBeta .* data.PoloidalthermalBeta + li ./ 2 - 3/2)
      % from Jean-Luc Paper
      flux_bv =  pi .* data.Majorradius .^ 2 .* bv
      data.flux_Bvert = flux_bv;
    end 
    sycomore_data.available_flux = data.flux_CS + data.flux_Bvert;
    sycomore_data.device = sprintf('DEMO-R%d-a%d-RBt%d-Ip%d',ceil(data.Majorradius*100),ceil(data.Minorradius*100),ceil(sycomore_data.rb0),ceil(sycomore_data.ip/1e6)); % device name (used to generate file name and comments)
    sycomore_data.scaling  = 0;             % code for energy plasma content scaling law (same as METIS)
    sycomore_data.H_H  = data.H98;                 % enhancement factor for energy content on flat top
    sycomore_data.f_Greenwald = data.Greenwaldfraction;       % Greenwald density fraction on flat top
    sycomore_data.ne_peak     = data.Centrene ./ data.Volaveragedne;
    if isfield(sycomore_data,'shot')
      sycomore_data.shot = sycomore_data.shot;                 % shot number, will be used to write data in UAL
      sycomore_data.run  = sycomore_data.run + 1;                 % run  number, will be used to write data in UAL
    else
      sycomore_data.shot = ceil(rand*1e4);                 % shot number, will be used to write data in UAL
      sycomore_data.run  = 1;                 % run  number, will be used to write data in UAL
    end
    if ~isfield(data,'Fractionnon_inductivecurrent');
	sycomore_data.f_ni = 0;                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
    elseif isnumeric(data.Fractionnon_inductivecurrent)
	sycomore_data.f_ni = data.Fractionnon_inductivecurrent;                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
    else
	sycomore_data.f_ni = 0;                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
    end
    sycomore_data.rimp = 0.1;               % ratio  between berillium and argon in core plasma
    sycomore_data.Recycling = 0.99;         % recycling  @ divertor 
    sycomore_data.zeff      = data.Zeff - 4 .* data.Z_2fraction ./ 100;          % line averaged Zeff without He ashes contribution
    sycomore_data.tau_He_o_tau_E_core = 5;  % ratio between core confinement time of He over energy confinement time.
    sycomore_data.rw = 0.7;                 % cyclotron radiation reflection coefficient
    sycomore_data.PNBI = data.NBIPower .* 1e6;              % NBI power used in Sycomore simulation (W, for ITER like default  case contains also ICRH power, if PNBI <  0, scale on ITER power).
    sycomore_data.fR_target = 1;            % position of outer target in unit of R0.
    sycomore_data.nbi_e_inj = 1e6;          % neutral beam energy injection.
    sycomore_data.signe = 1;                % sign of Bt . J_phi: must be provided by the machine description
    sycomore_data.flux_expansion = 3;       % divertor flux expansion 
    sycomore_data.target_orientation = 90;  % poloidal orientation of outer divertor target  (degrees)
    sycomore_data.S_factor = 3e-3;          % divertor flux spreading factor (m) 
    sycomore_data.cW               = 1e-5;  % tungsten concentration in core plasma   , must be provide by Sycomore 
    % There no enhancement of impurities concentration id divertor. have been checked  in soldiv module of Sycomore.
    sycomore_data.fzmax_div        = 0 .* data.Z_18fraction;     % Argon concentration enhancement in divertor (%)  , must be provide by Sycomore
    sycomore_data.f_DSOL             = 4;   % multiplicator applied to the SOL width when it is defined by a scaling law, for H mode only (Dsol_Hmode = factor_scale * Goldston scaling).
    sycomore_data.f_ne_LCFS          = 2;   % factor applied to edge scaling law for density:\nif > 0, ne_edge =  nea_factor * LCFS_denstity_scaling_law;\nif < 0,  ne_edge =  abs(nea_factor) * n_bar'electron density at LCFS (m^-3)  , must be provide by Sycomore
    sycomore_data.couplage         = 0.15   % coupling coefficient given the ratio between poloidal flux consumption in the plasam and in the central selenoid (CS) during breakdown to take into account dissipated flux in passive structure.
    sycomore_data.duration         = 7200;  % estimated duration of the shot in s
    % use constant P/R for allowed shine througth (and first orbit losses)
    % take < 5% of 33 MW for ITER @ 6.2m -> 2e5 W/m
    %sycomore_data.shine_through_limit =  2e5 .* data.Majorradius;  % limit of power lost in shine through (W)
    sycomore_data.shine_through_limit =  max(0.05 .* sycomore_data.PNBI , 2e5 .* data.Majorradius);  % limit of power lost in shine through (W)

    %  optionnal fields 
    sycomore_data.tokamak  = '';           % new UAL tokamak database name
    sycomore_data.user     = '';           % new UAL user database name
    sycomore_data.dataversion = '';        % selected a differente version of data (not the last one)
    sycomore_data.occurrence = '';         % input cpo scenario occurrence in Kepler (default = [])
    sycomore_data.path4metis_files = '';   % path for files save as a postprocessing of METIS (METIS data, figures, ...); if left empty, then no wrinting       
 
    %  structure details for sepa_option (with data example for ITER):
    sepa_option.rxup      = data.Uppertriangularity;     % upper triangularity (minor radius unit)
    sepa_option.zxup      = data.Upperelongation;     % upper altitude X point (minor radius unit)
    sepa_option.apup      = data.AngleLCMSupperX_outboard;         % upper separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amup      = data.AngleLCMSupperX_inboard;         % upper separatrix angle (-R,X) (HFS, degrees)
    sepa_option.ra        = data.Majorradius;       % major radius R0 (m) [6.2]
    sepa_option.za        = 0;         % altitude of the magnetic axis (m) [0.9]
    sepa_option.a         = data.Minorradius;         % minor radius (m) [2]
    sepa_option.rxdo      = data.Lowertriangularity;     % lower triangularity (minor radius unit)
    sepa_option.zxdo      = data.Lowerelongation;     % lower altitude X point (minor radius unit)
    sepa_option.apdo      = data.AngleLCMSlowerX_outboard;     % lower separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amdo      = data.AngleLCMSlowerX_inboard;     % lower separatrix angle (-R,X)  (HFS, degrees)
    sepa_option.b0        = sycomore_data.rb0 ./ data.RETF;      % magnetic field at R0
    sepa_option.delta     = data.RINPlasma - data.RETF;      % minimal distance between plasma and TF conductor 
 
    % metis call
    z0dinput = sycomore2metis(sepa_option,parameters_filename,sycomore_data,computation_mode);

end