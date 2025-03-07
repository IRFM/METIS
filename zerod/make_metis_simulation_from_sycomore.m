% name of Sycomore result file
[filename, pathname] = uigetfile('*.out', 'Select a Sycomore result file');
%sycomore_filename = '~/zineb/data/DEMO2014/Sycomore/DEMOp_fGW_0p8.txt';
%sycomore_filename = '~/zineb/zerod_data/Sycomore/sycomore_output_DEMO1_end.out';
if isempty(filename)
	return
else
	sycomore_filename = fullfile(pathname,filename);
end
% parameters
% name of file containig METIS parameters (leave empty to used default parameters)
parameters_filename = input('METIS parameters file (optionnal,leave empty to use default parameters)? ','s');
if isempty(parameters_filename)
		parameters_filename =[]; % to use default tuning
end		
% computation mode (2= fast without optimisation; 1 = fast with optimisation for testing; 0 = full)
%computation_mode = 2;
%computation_mode = 1;
%computation_mode = 0;
computation_mode = input('computing mode: 0= full; 1=  fast with optimisation (for testing); 2= fast without optimisation (just prepare METIS file) ? ');
switch computation_mode
case {0,1,2}
	%rien
otherwise
	return
end
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
  if isnumeric(data.Fractionnon_inductivecurrent)
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
if isnumeric(data.Fractionnon_inductivecurrent)
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
% for testing 
[p,fname] = fileparts(tempname);
mkdir(fname)
%
sycomore_data.path4metis_files = fname;   % path for files save as a postprocessing of METIS (METIS data, figures, ...); if left empty, then no wrinting       
 
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
%z0dinput = sycomore2metis(sepa_option,parameters_filename,sycomore_data,computation_mode);
%return
[z0dinput,post,option,sepa_option,sycomore_data,tsnapshot] = sycomore2metis(sepa_option,parameters_filename,sycomore_data,computation_mode);

% comparison of METIS and Sycomore results 
ind_fusmax = find(post.zerod.temps >= tsnapshot,1); 
ind_fusmax_1d = find(post.profil0d.temps >= tsnapshot,1); 
disp('-------------------------------------------------------------------------------');
fprintf('Snap shot data @ t = %g s\n',tsnapshot);
fprintf('Variable name\t\t\t\t\tSYCOMORE\t\t\t\t\tMETIS\n');     
fprintf('Plasma volume (m^3):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Plasmavolume,post.zerod.vp(ind_fusmax));
fprintf('Poloidal crosssection (m^2):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Poloidalcrosssection,post.zerod.sp(ind_fusmax));
fprintf('Separatrix area (m^2):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Separatrixarea,post.zerod.sext(ind_fusmax));
fprintf('Plasma current (MA):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Plasmacurrent,post.zerod.ip(ind_fusmax) ./ 1e6);
fprintf('Bt on axis (T):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Btonaxis,post.profil0d.fdia(ind_fusmax_1d,end) ./ post.profil0d.Raxe(ind_fusmax_1d,end));
fprintf('q95 (-):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.q95,post.zerod.q95(ind_fusmax));
fprintf('tauE (s):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.tauE,post.zerod.taue(ind_fusmax));

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
	
fprintf('H98 (a.u.):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.H98,post.zerod.taue(ind_fusmax)./tauh(ind_fusmax));
fprintf('Greenwald fraction (-):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Greenwaldfraction,post.zerod.nbar(ind_fusmax)./post.zerod.negr(ind_fusmax));
if ischar(data.Fractionnon_inductivecurrent)
    data.Fractionnon_inductivecurrent = NaN;
end
fprintf('Fraction non inductive current (%%):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Fractionnon_inductivecurrent,post.zerod.ini(ind_fusmax)./post.zerod.ipar(ind_fusmax));
fprintf('Fraction Boostrap current (%%):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.FractionBoostrapcurrent,100 .* post.zerod.iboot(ind_fusmax)./post.zerod.ipar(ind_fusmax));
fprintf('Normalized Beta (%%):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.NormalizedBeta,100 .* post.zerod.betan(ind_fusmax));
fprintf('Normalized thermal Beta (%%):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.NormalizedthermalBeta,100 .* post.zerod.betan(ind_fusmax) .* post.zerod.wth(ind_fusmax) ./ post.zerod.w(ind_fusmax));
fprintf('Poloidal thermal Beta (-):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.PoloidalthermalBeta,post.zerod.betap(ind_fusmax));
fprintf('Normalized pedestal radius (-):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Pedestalnorm_radius,post.profil0d.rmx(ind_fusmax_1d,end - 1) ./ post.profil0d.rmx(ind_fusmax_1d,end));
fprintf('Volume averaged ne (10^19 m^-3):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Volaveragedne/1e19, post.zerod.nem(ind_fusmax)/1e19);
fprintf('Centre ne (10^19 m^-3):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Centrene/1e19, post.zerod.ne0(ind_fusmax)/1e19);
fprintf('Pedestal ne (10^19 m^-3):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Pedestalne/1e19, post.zerod.neped(ind_fusmax)/1e19);
fprintf('LCFS ne (10^19 m^-3):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Edgene/1e19, post.zerod.nebord(ind_fusmax)/1e19);
fprintf('Volume averaged Te (keV):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.VolaveragedTe,post.zerod.tem(ind_fusmax)/1e3);
fprintf('Centre Te (keV):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.CentreTe,post.zerod.te0(ind_fusmax)/1e3);
fprintf('Pedestal Te (keV):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.PedestalTe,post.zerod.teped(ind_fusmax)/1e3);
fprintf('Edge Te (keV):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.EdgeTe,post.zerod.tebord(ind_fusmax)/1e3);
fprintf('<Ti/Te> (-):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.TioTe,post.zerod.tite(ind_fusmax));
fprintf('Z_2fraction (%%):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Z_2fraction,100 .* post.zerod.nhem(ind_fusmax) ./ post.zerod.nem(ind_fusmax));
if  post.z0dinput.option.zimp == 18
    fprintf('Z_18fraction (%%):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Z_18fraction,100 .* post.zerod.nimpm(ind_fusmax) ./ post.zerod.nem(ind_fusmax));
elseif post.z0dinput.option.zmax == 18
    fprintf('Z_18fraction (%%):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Z_18fraction,100 .* post.z0dinput.option.rimp .* post.zerod.nimpm(ind_fusmax) ./ post.zerod.nem(ind_fusmax));
end
fprintf('Zeff (%%):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Zeff,post.zerod.zeff(ind_fusmax));
fprintf('Q (-):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Q,rfan .* post.zerod.pfus(ind_fusmax) ./ 1e6 ./ padd(ind_fusmax));
fprintf('Fusion power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Fusionpower ,rfan .* post.zerod.pfus(ind_fusmax) ./ 1e6);
fprintf('NBI power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.NBIPower,real(post.zerod.pnbi(ind_fusmax) ./ 1e6) + imag(post.zerod.pnbi(ind_fusmax) ./ 1e6));
fprintf('ICRH power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.ICRHPower,post.zerod.picrh(ind_fusmax));
if post.z0dinput.option.lhmode == 5
  fprintf('ECRH power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.ECRHPower,post.zerod.plh(ind_fusmax) ./ 1e6 + post.zerod.pecrh(ind_fusmax));
  fprintf('LowerHybrid power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.LowerHybridPower,0);
else
  fprintf('ECRH power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',post.zerod.pecrh(ind_fusmax));
  fprintf('LowerHybrid power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.LowerHybridPower,data.ECRHPower,post.zerod.plh(ind_fusmax) ./ 1e6);
end
fprintf('Ohmic power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Ohmicpower,post.zerod.pohm(ind_fusmax) ./ 1e6);
fprintf('Bremsstrahlung power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Bremsstrahlungpower,post.zerod.pbrem(ind_fusmax) ./ 1e6);
fprintf('Synchrotron power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Synchrotronpower,post.zerod.pcyclo(ind_fusmax) ./ 1e6);
fprintf('Line radiation power (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Lineradiationpower,post.zerod.prad(ind_fusmax) ./ 1e6 + post.zerod.pioniz(ind_fusmax) ./ 1e6);
fprintf('Power through separatrix (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Powerthroughseparatrix,post.zerod.plhthr(ind_fusmax) ./ 1e6);
fprintf('L -> H transition power threshold (MW):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.L_Htransitionpowerthres_,post.zerod.plossl2h(ind_fusmax) ./ 1e6  + post.z0dinput.option.l2hmul);
fprintf('Thermal energy content (MJ):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Thermalenergycontent,post.zerod.wth(ind_fusmax) ./ 1e6);
fprintf('Qpeak divertor (MW/m^2):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Qpeak_divertor,post.zerod.peakdiv(ind_fusmax) ./ 1e6);
% flush diary file 
diary off
diary on
