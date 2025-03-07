% parameters
option = post.z0dinput.option;
option.neasser      = 1; 
option.natural      = 1;
option.fnbar_nat    = 1; 
option.fn0a             = 1;
option.fn0a_div         = 0.1;
%    option.dilution         = 1;          % effect not taken into account in Sycomore      
%    option.tau_limitation   =  'On';      % effect not taken into account in Sycomore but usefull for rampup and rampdown
%    option.ploss_exp        = 'max_power'; % best choice from DEMO1 and DEMO2 studies for Ploss computation
%    option.plhthr           = 'P_LCFS';    % used the real power crossing the LCFS
%    option.hysteresis       = 0;           % No hysteresis for the back transition to  L-mode 
%    option.fpped            = 1;          
%    option.hmore_pped       = 2;          
%    option.usepped_scl      = 2;         % use minmum pressure between scaling for predestal pressure and standard METIS rule Pped = K (W_H -W_L)
%    option.fintrinsic       = 0;         % machine dependent factor: intrinsic rotation must be tuned (used scaling in collisionality)         
option.angle_ece        = 180;        % optimisation of current drive by using HFS resonnance
option.sens              = 1;         % co current ECCD      
option.lhmode            = 5;          % LHCD channel used for 2nd ECCD system
option.etalh             = 1;          % 2nd eccd system in co-current
option.xlh               = 0;          % 2nd ECCD system position of maximum heat depostion
option.dlh               = 0.4000;     % 2nd ECCD system width of heat depostion profile
option.wlh               = 0;          % 2nd ECCD system position of maximum heat depostion
option.angle_ece2        = 180;        % optimisation of current drive by using HFS resonnance
% computation mode (2= fast without optimisation; 1 = fast with optimisation for testing; 0 = full)
%computation_mode = 2;
%computation_mode = 1;
computation_mode = 0;
% seach for flattop
pin =  post.z0dinput.cons.picrh + post.z0dinput.cons.plh + post.z0dinput.cons.pecrh + real(post.z0dinput.cons.pnbi)  + imag(post.z0dinput.cons.pnbi);
fact = post.z0dinput.cons.ip .* post.z0dinput.cons.nbar .* pin .*post.z0dinput.cons.hmore;
[n,x] = hist(fact,201);
ind_max = find(n == max(n));
ind_ok = find((fact >= min(x(ind_max) - mean(diff(x)))) & (fact <= max(x(ind_max) + mean(diff(x)))));
% mapping the data
sycomore_data.ip = mean(post.z0dinput.cons.ip(ind_ok));                % flat top plasma current (A).
sycomore_data.q95 = mean(post.zerod.q95(ind_ok));                % flat top safety factor
sycomore_data.rb0 = mean(post.z0dinput.geo.b0(ind_ok) .* post.z0dinput.geo.R(ind_ok));         % magnetic rigidity (flat top)
sycomore_data.available_flux = post.z0dinput.option.available_flux;
sycomore_data.device = sprintf('DEMO-R%d-a%d-RBt%d-Ip%d',ceil(mean(post.z0dinput.geo.R(ind_ok))*100), ...
                       ceil(mean(post.z0dinput.geo.a(ind_ok))*100), ...
                       ceil(sycomore_data.rb0),ceil(sycomore_data.ip/1e6)); % device name (used to generate file name and comments)
sycomore_data.scaling  = post.z0dinput.option.scaling;             % code for energy plasma content scaling law (same as METIS)

sycomore_data.H_H  = mean(post.z0dinput.cons.hmore(ind_ok));                 % enhancement factor for energy content on flat top
sycomore_data.f_Greenwald =  mean(post.z0dinput.cons.nbar(ind_ok) ./ 1e20 ./ ...
                             (post.z0dinput.cons.ip(ind_ok) ./ 1e6) .* pi .* post.z0dinput.geo.a(ind_ok) .^ 2);       % Greenwald density fraction on flat top
sycomore_data.ne_peak     = []; % kept mETIS parameters
if isfield(sycomore_data,'shot')
  sycomore_data.shot = sycomore_data.shot;                 % shot number, will be used to write data in UAL
  sycomore_data.run  = sycomore_data.run + 1;                 % run  number, will be used to write data in UAL
else
  sycomore_data.shot = ceil(rand*1e4);                 % shot number, will be used to write data in UAL
  sycomore_data.run  = 1;                 % run  number, will be used to write data in UAL
end
sycomore_data.f_ni = mean(post.zerod.ini(ind_ok) ./ post.zerod.ipar(ind_ok));                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
sycomore_data.rimp = post.z0dinput.option.rimp;               % ratio  between berillium and argon in core plasma
sycomore_data.Recycling = post.z0dinput.option.Recycling;         % recycling  @ divertor 
sycomore_data.zeff      = mean(post.z0dinput.cons.zeff(ind_ok));          % line averaged Zeff without He ashes contribution
sycomore_data.tau_He_o_tau_E_core = post.z0dinput.option.tauhemul;  % ratio between core confinement time of He over energy confinement time.
sycomore_data.rw = post.z0dinput.option.rw;                 % cyclotron radiation reflection coefficient
sycomore_data.PNBI = mean(real(post.z0dinput.cons.pnbi(ind_ok))  + imag(post.z0dinput.cons.pnbi(ind_ok)));              % NBI power used in Sycomore simulation (W, for ITER like default  case contains also ICRH power, if PNBI <  0, scale on ITER power).
sycomore_data.fR_target = post.z0dinput.option.fR_target;            % position of outer target in unit of R0.
sycomore_data.nbi_e_inj = max(post.z0dinput.option.einj,post.z0dinput.option.einj2);          % neutral beam energy injection.
sycomore_data.signe = post.z0dinput.option.signe;                % sign of Bt . J_phi: must be provided by the machine description
sycomore_data.flux_expansion = 3;       % divertor flux expansion 
sycomore_data.target_orientation = 90;  % poloidal orientation of outer divertor target  (degrees)
sycomore_data.S_factor = 3e-3;          % divertor flux spreading factor (m) 
sycomore_data.cW              =  post.z0dinput.option.cw_offset;  % tungsten concentration in core plasma   , must be provide by Sycomore 

% There no enhancement of impurities concentration id divertor. have been checked  in soldiv module of Sycomore.
sycomore_data.fzmax_div        = post.z0dinput.option.fzmax_div;     % Argon concentration enhancement in divertor (%)  , must be provide by Sycomore
sycomore_data.f_DSOL           = post.z0dinput.option.factor_scale;   % multiplicator applied to the SOL width when it is defined by a scaling law, for H mode only (Dsol_Hmode = factor_scale * Goldston scaling).
sycomore_data.f_ne_LCFS          = post.z0dinput.option.nea_factor;   % factor applied to edge scaling law for density:\nif > 0, ne_edge =  nea_factor * LCFS_denstity_scaling_law;\nif < 0,  ne_edge =  abs(nea_factor) * n_bar'electron density at LCFS (m^-3)  , must be provide by Sycomore
sycomore_data.couplage         = 0.15   % coupling coefficient given the ratio between poloidal flux consumption in the plasam and in the central selenoid (CS) during breakdown to take into account dissipated flux in passive structure.
sycomore_data.duration         = (2/3) .* max(post.z0dinput.cons.temps(ind_ok));  % estimated duration of the shot in s
% use constant P/R for allowed shine througth (and first orbit losses)
% take < 5% of 33 MW for ITER @ 6.2m -> 2e5 W/m
%sycomore_data.shine_through_limit =  2e5 .* data.Majorradius;  % limit of power lost in shine through (W)
sycomore_data.shine_through_limit =  max(0.05 .* sycomore_data.PNBI , 2e5 .* mean(post.z0dinput.geo.R(ind_ok)));  % limit of power lost in shine through (W)

%  optionnal fields 
sycomore_data.tokamak  = '';           % new UAL tokamak database name
sycomore_data.user     = '';           % new UAL user database name
sycomore_data.dataversion = '';        % selected a differente version of data (not the last one)
sycomore_data.occurrence = '';         % input cpo scenario occurrence in Kepler (default = [])
% for testing 
[p,fname] = fileparts(tempname);
mkdir(fname)
%
sycomore_data.path4metis_files = fname;   % path for files save as a postprocessing of METIS (METIS data, figures, ...); if left empty, then no wrinting       
 
%  structure details for sepa_option (with data example for ITER):
if isfield(post.z0dinput.exp0d,'Rsepa') && ~isempty(post.z0dinput.exp0d.Rsepa) 
    ind_lcfs = find(post.zerod.rm == max(post.zerod.rm),1);
    R = post.z0dinput.exp0d.Rsepa(ind_lcfs,:);
    Z =post.z0dinput.exp0d.Zsepa(ind_lcfs,:);
    lcfs_fname = fullfile(fname,'LCFS.mat');
    save(lcfs_fname,'R','Z');
    sepa_option.filename = lcfs_fname;
    sepa_option.rxup      = post.z0dinput.geo.d(ind_lcfs);     % upper triangularity (minor radius unit)
    sepa_option.zxup      = post.z0dinput.geo.K(ind_lcfs);     % upper altitude X point (minor radius unit)
    sepa_option.apup      = 0;         % upper separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amup      = 0;         % upper separatrix angle (-R,X) (HFS, degrees)
    sepa_option.ra        = post.z0dinput.geo.R(ind_lcfs);       % major radius R0 (m) [6.2]
    sepa_option.za        = 0;         % altitude of the magnetic axis (m) [0.9]
    sepa_option.a         = post.z0dinput.geo.a(ind_lcfs);         % minor radius (m) [2]
    sepa_option.rxdo      = sepa_option.rxup;     % lower triangularity (minor radius unit)
    sepa_option.zxdo      = sepa_option.zxup;     % lower altitude X point (minor radius unit)
    sepa_option.apdo      = 0;     % lower separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amdo      = 0;     % lower separatrix angle (-R,X)  (HFS, degrees)
    sepa_option.b0        = sycomore_data.rb0 ./ (sepa_option.ra - sepa_option.a - 1.4);      % magnetic field at R0
    sepa_option.delta     = 1.4;      % minimal distance between plasma and TF conductor
    disp('using LCFS given by points');
elseif isfield(post.z0dinput,'sepa_option')
    sepa_option = post.z0dinput.sepa_option;
    disp('using generated LCFS');
else
    sepa_option.rxup      = mean(post.z0dinput.geo.d(ind_ok));     % upper triangularity (minor radius unit)
    sepa_option.zxup      = mean(post.z0dinput.geo.K(ind_ok));     % upper altitude X point (minor radius unit)
    sepa_option.apup      = 0;         % upper separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amup      = 0;         % upper separatrix angle (-R,X) (HFS, degrees)
    sepa_option.ra        = mean(post.z0dinput.geo.R(ind_ok));       % major radius R0 (m) [6.2]
    sepa_option.za        = 0;         % altitude of the magnetic axis (m) [0.9]
    sepa_option.a         = mean(post.z0dinput.geo.a(ind_ok));         % minor radius (m) [2]
    sepa_option.rxdo      = sepa_option.rxup;     % lower triangularity (minor radius unit)
    sepa_option.zxdo      = sepa_option.zxup;     % lower altitude X point (minor radius unit)
    sepa_option.apdo      = 0;     % lower separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amdo      = 0;     % lower separatrix angle (-R,X)  (HFS, degrees)
    sepa_option.b0        = sycomore_data.rb0 ./ (sepa_option.ra - sepa_option.a - 1.4);      % magnetic field at R0
    sepa_option.delta     = 1.4;      % minimal distance between plasma and TF conductor 
    disp('using moment for LCFS');
end
% metis call
%z0dinput = sycomore2metis(sepa_option,option,sycomore_data,computation_mode);
%return
[z0dinput,post,option,sepa_option,sycomore_data,tsnapshot] = sycomore2metis(sepa_option,option,sycomore_data,computation_mode);

% comparison of METIS and Sycomore results 
ind_fusmax = find(post.zerod.temps >= tsnapshot,1); 
ind_fusmax_1d = find(post.profil0d.temps >= tsnapshot,1); 
disp('-------------------------------------------------------------------------------');
fprintf('Snap shot data @ t = %g s\n',tsnapshot);
ind_fus_max = find(post.zerod.pfus == max(post.zerod.pfus),1);
tend  = post.zerod.temps(ind_fus_max);
demo0d_2012_print;

% auxilliary data
rfan = (3.56e6 + 14.03e6) ./ 3.56e6 ;
padd = (post.z0dinput.cons.picrh + post.z0dinput.cons.pecrh + real(post.z0dinput.cons.pnbi) + imag(post.z0dinput.cons.pnbi) + post.zerod.pohm + post.z0dinput.cons.plh) / 1e6;
pin  = padd + post.zerod.pfus ./ 1e6;
ploss = pin - (post.zerod.pbrem + post.zerod.pcyclo + post.z0dinput.option.fprad .* post.zerod.prad + post.zerod.pioniz) ./ 1e6;
ip = post.zerod.ip ./ 1e6;
Bt = post.z0dinput.geo.b0;
ne = post.z0dinput.cons.nbar./ 1e19;
R  =  post.z0dinput.geo.R;
a  =  post.z0dinput.geo.a;
ep  = a ./ R;
K   = post.z0dinput.geo.K; 
Vp = post.zerod.vp;
Ka  = Vp ./ (2*pi^2.*R.*a.^2);
meff = post.zerod.meff;
% ITERH-98P(y,2)        
tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
	R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s    
	
% flush diary file 
diary off
diary on