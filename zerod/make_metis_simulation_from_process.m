% parameters
% name of file containig METIS parameters (leave empty to used default parameters)
%parameters_filename ='';
%parameters_filename =[]; % to use default tuning
% template for METIS parameters (to tune the composition)
z0dinput = sycomore2metis;
parameters_filename = z0dinput.option
% name of Sycomore result file
process_filename = '/donnees/JA132999/zineb/zerod_data/DEMO2015/Design_march_2015/DEMO2_Alternative_Design_-_March_2015_-__2MEMAF_v1_0.dat';
[process_filename,pf,idx]= uigetfile('*','Select a PROCESS .dat results file',process_filename);
if idx == 0
  return
end
process_filename = fullfile(pf,process_filename);
% computation mode (1 = fast for testing; 0 = full)
%computation_mode = 1;
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
      if data.plasma.(name).value ~= 0
	    label{end+1} = name_liste{ind}
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
sycomore_data.f_ni = data.plasma.faccd_.value + data.plasma.bootipf_.value;                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
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
sycomore_data.PNBI = data.hcd.pnbeam.value .* 1e6;              % NBI power used in Sycomore simulation (W, for ITER like default  case contains also ICRH power, if PNBI <  0, scale on ITER power).
sycomore_data.fR_target = 1;            % position of outer target in unit of R0.
sycomore_data.nbi_e_inj = data.hcd.enbeam.value .* 1e3;          % neutral beam energy injection.
sycomore_data.signe = 1;                % sign of Bt . J_phi: must be provided by the machine description
sycomore_data.flux_expansion = 3;       % divertor flux expansion 
sycomore_data.target_orientation = 90;  % poloidal orientation of outer divertor target  (degrees)
sycomore_data.S_factor = 3e-3;          % divertor flux spreading factor (m) 
% There no enhancement of impurities concentration id divertor. have been checked  in soldiv module of Sycomore.
sycomore_data.fzmax_div        = 100 .* max(0,(data.divertor.zeffso.value -  data.plasma.zeff.value) ./ parameters_filename.zmax .^ 2 - fzmax);     % Argon concentration enhancement in divertor (%)  , must be provide by Sycomore
sycomore_data.f_DSOL             = data.divertor.delw.value./3e-3;   % multiplicator applied to the SOL width when it is defined by a scaling law, for H mode only (Dsol_Hmode = factor_scale * Goldston scaling).
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
% for testing 
[p,fname] = fileparts(tempname);
mkdir(fname)
%
sycomore_data.path4metis_files = fname;   % path for files save as a postprocessing of METIS (METIS data, figures, ...); if left empty, then no wrinting       
 
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
sepa_option.delta     = data.plasma.rmajor.value - data.plasma.rminor.value - data.radial_build.TFcoilinboardleg.value(2,1);      % magnetic field at R0
 
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
fprintf('Poloidal thermal Beta (-):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Pedestalnorm_radius,post.profil0d.rmx(ind_fusmax_1d,end - 1) ./ post.profil0d.rmx(ind_fusmax_1d,end));
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
fprintf('Q (-):\t\t\t\t\t%g\t\t\t\t\t%g\n',data.Q,rfan .* post.zerod.pfus(ind_fusmax) ./ padd(ind_fusmax));
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
