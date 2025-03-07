% this function make auto run of METIS for Tore Supra shot with predefined parameters.
% the energy content is provided by Wdia mesurement.
% if there are available, TPROF electron density and temperature profiles are used.
%
% syntax:
%
% 	[post,tprofok] = autoTSmetis(shot,gas,time,filename,tprofonoff)
%
% input:
%	shot = shot number [5877,48361]
%
%       gas  = 1 = D & 2 = He;
%              If isempty or missing switch to automatic selection.
%
%       time = time slices for calculation (vector [nbt,1]).
%               If isempty or missing, use default time slices vector of METIS.
%               If is saclar, compute steady state solution for this selected time.
%
%	filename = optionnal file name in which METIS data are saved (without .mat).
%
%       tprofonoff = optional control of TPROF use.
%                    If = 0, METIS will run without TPROF profiles even if some TPROF run exsit.
%
% output:
%	post = standard METIS data structure
%
%       tprofok = flag set to 1 if TPROF data are used
%
function [post,tprofok] = autoTSmetis(shot,gas,time,filename,tprofonoff)

% default outputs
post  = [];
tprofok  = [];

% reset external data
rm_cs4m;

% inputs handling
if nargin < 1
    error('shot number must be provided !')
end
if nargin < 2
    gas = 1; % assume D
    % must be 2 for He
elseif isempty(gas)
    gas = 1; % assume D
    % must be 2 for He
end
workingpoint = NaN;
if nargin < 3
    time = []; % METIS select time slices
elseif length(time) == 1
    workingpoint = time;
    time = [];
end
if nargin < 4
    filename  = ''; % no nackup of METIS result
end
if nargin < 5
    tprofonoff = 1;
elseif isempty(tprofonoff)
    tprofonoff = 1;
end

% initialise METIS data
fprintf('Preparing METIS data set for Tore Supra shot %d\n',shot);
try
    z0dinput = zerod_init(1,shot,gas,time);
catch
    z0dinput = [];
end
if isempty(z0dinput)
    disp('No data for this shot');
    return
end

% set METIS parameters
z0dinput.option = map_option_ts(z0dinput.option);

% initialising  external data
if tprofonoff ~= 0
    % search for occurrence of TPROF
    tprofok = 0;
    for k=10:-1:1
        occ = k /10;
        occ(occ >= 1) = 0;
        fprintf('tying %g\n',shot + occ);
        [amin,tamin] = tsbase(shot + occ,'SPROFIP');
        if ~isempty(amin)
            tprofok = 1;
            break
        end
    end
    if tprofok == 1
        fprintf('TPROF data found at occurrence %d\n',10 .* occ);
        % reading density for TPROF
        [nefit,ttprof,rtprof,cert] = tsbase(shot + occ,'GPROFNEFIT');
        if ~isempty(nefit)
            % 1- Electron density
            NE_EXP.x     =  rtprof;
            NE_EXP.ne    =  nefit;
            NE_EXP.temps =  ttprof;
            setappdata(0,'NE_EXP',NE_EXP);
            % change configuration
            z0dinput.option.neasser = 0;
            z0dinput.option.natural = 0;
            % recompute nbar 
            nbar = trapz(rtprof,nefit,2);
            %figure;plot(ttprof,nbar,'r',z0dinput.cons.temps,z0dinput.cons.nbar,'b');
            nbar= interp1(ttprof,nbar,z0dinput.cons.temps,'nearest',NaN);
            nbar(~isfinite(nbar)) = z0dinput.cons.nbar(~isfinite(nbar));
            z0dinput.cons.nbar    = nbar + sqrt(-1) .* imag(z0dinput.cons.nbar);
            %keyboard
        end
        [tefit,ttprof,rtprof,cert] = tsbase(shot + occ,'GPROFTEFIT');
        if ~isempty(tefit)
            % 1- Electron temperature
            TE_EXP.x     =  rtprof;
            TE_EXP.te    =  tefit .* 1e3;
            TE_EXP.temps =  ttprof;
            setappdata(0,'TE_EXP',TE_EXP)
        end
        
    else
        % try to find a TPROF datafile on disk
        occ =  sum(abs('disk'));
        file = fullfile('/cgc_data/zineb/.upper/ts/data/',sprintf('%d',shot),'bile.mat.gz');
        if exist(file)
            tpf     = tempname;
            tpf_ext = sprintf('%s.mat.gz',tpf);
            copyfile(file,tpf_ext);
            unix(sprintf('gzip -d %s',tpf_ext));
            data = load(tpf);
            delete(sprintf('%s.mat',tpf));
            if ~isempty(data) && isfield(data,'rhofit') &&  ...
                    isfield(data,'nefit') && isfield(data,'tefit') && ...
                    ~isempty(data.rhofit) && ~isempty(data.tefit) && ~isempty(data.nefit)
                
                tprofok = 1;
                
                fprintf('TPROF data found in file %s\n',file);
                if ~isempty(data.nefit)
                    % 1- Electron density
                    NE_EXP.x     =  data.rhofit;
                    NE_EXP.ne    =  data.nefit;
                    NE_EXP.temps =  data.times;
                    setappdata(0,'NE_EXP',NE_EXP);
                    % change configuration
                    z0dinput.option.neasser = 0;
                end
                if ~isempty(data.tefit)
                    % 1- Electron temperature
                    TE_EXP.x     =  data.rhofit;
                    TE_EXP.te    =  data.tefit .* 1e3;
                    TE_EXP.temps =  data.times;
                    setappdata(0,'WORKING_POINT_TIME',workingpoint);
                    zassignin('base','z0dinput',z0dinput);
                    evalin('base','z0working_point;');
                    post = evalin('base','post');
                end
                if ~isempty(filename)
                    if tprofok == 1
                        metis_save(sprintf('%s_TS@%d_occ_%d',filename,shot,occ*10),post);
                    else
                        metis_save(sprintf('%s_TS@%d_notprof',filename,shot),post);
                    end
                end
                setappdata(0,'TE_EXP',TE_EXP)
            end
        end
    end
    if tprofok == 0
        disp('No TPROF data found !');
    end
end
% call METIS execution
if ~isfinite(workingpoint)
    zassignin('base','z0dinput',z0dinput);
    evalin('base','metis_run;');
    post = evalin('base','post');
else
    setappdata(0,'WORKING_POINT_TIME',workingpoint);
    zassignin('base','z0dinput',z0dinput);
    evalin('base','z0working_point;');
    post = evalin('base','post');
end
if ~isempty(filename)
    if tprofok == 1
        metis_save(sprintf('%s_TS@%d_occ_%d',filename,shot,occ*10),post);
    else
        metis_save(sprintf('%s_TS@%d_notprof',filename,shot),post);
    end
end

% default METIS paramaters for TS
function option = map_option_ts(option)

%option.gaz = 2;
option.neasser = 1;
option.Recycling = 0.7;
option.natural = 1;
option.ane = 11;
option.fn0a = 1;
option.fn0a_div = 0.1;
%
option.scaling = 4;
option.dilution = 1;
option.tau_limitation = 'On';
option.l2hscaling = 3;
option.pl2h_mass_charge = 1;
option.modeh = 0;
option.hysteresis = 0;
option.configuration = 3;
option.l2hslope =  0.5;
option.usepped_scl = 1;
option.taurotmul = 0;
option.fintrinsic = 0.2;
option.xiioxie = -4.5;
option.kishape = 0;
option.xieorkie = 0;
option.omega_shape = 0;
option.fstiff = 1;
option.ploss_exp = 'max_power';
option.xiioxie_ped = 0;
option.hmore_pped  = 2;
option.fpl2h_lim   = 2;
option.ki_expo     = 2;
option.plhthr      = 'P_LCFS';
option.grad_ped    = 3;
option.ode_pped    = 1;
option.adiabatic   = 1;

%
option.qdds = 1;
option.kidds = 3;
option.sitb = 2;
option.itb_sensitivity = 1;
option.itb_slope_max = 2;
option.smhd = 100;
option.tmhd = 0;
%
option.runaway = 0;
option.modeboot = 2;
%
option.li = 1;
option.breakdown = - 10;
option.berror = 0; % no breakdown simulation by default
option.L_eddy = 4.76e-4;
option.R_eddy = 4.1e-5;
option.C_eddy = 1;
option.B_eddy = 0.1;
option.I_eddy = 0;
option.p_prefill = 0.7e-03;
option.VV_volume = 25;
%
option.zeff = 0;
option.faccu = 0;
option.heat_acc = 0;
option.fne_acc = 0;
option.zmax = 8;
option.zimp = 6;
option.rimp = 0.3;
option.density_model ='minconv';


%
option.frad = 1;
option.matthews = -1;
option.fte_edge = 1;
option.gaunt = 1;
option.noncoronal = -1;
option.z_prad = 'Stangeby';
%
option.sol_lscale = 0;
option.eioniz     = 0;
option.fnesol     = 0;
option.sol_model  = 'scaling';
option.lcx = 1;
option.fcond = -1;
option.fmom = 0;
option.lambda_scale = 3;
option.sol_rad = 'decoupled';
option.fzmax_div = -1;
option.W_effect = 0;
option.cw_factor = 0;
option.cw_offset = 5e-5;
option.factor_scale = 1;
option.fpower = 0.6000;
option.fR_target = 2.7 / 2.94;
option.mach_corr = 1;

%
option.angle_ece = 90;
option.synergie  = 1;
option.sens      = 1;
option.eccdmul   = 1;
%
option.angle_nbi = 90;
option.rtang     = 2.85;
option.zext      = 0.3;
option.einj      = 500000;
option.nbicdmul  = 1;
option.nb_nbi    = 2;
option.e_shielding = 'Honda-Sauter';
option.drs1        = 0;
option.dzs1        = 0;
%
option.angle_nbi2 = 0;
option.rtang2 = 2.85;
option.zext2  = 0.1;
option.einj2  = 85000;
option.nbicdmul2 = 1;
option.drs2   = 0;
option.dzs2   = 0;
%
option.lhmode = 0;
option.etalh  = 0.8;
option.wlh = 32 .* 11e-3;
option.xlh = 0.3;
option.dlh = 0.2;
option.angle_ece2 = 90;
%option.npar0 = 2;
%option.npar_neg = -4;

% used as third injector
%  option.fwcd = 0;
%  option.mino = 'H';
%  option.cmin = 0.15;
%  option.nphi = 30;
%  option.freq = 55.5;
%  option.icrh_width = 0.7;
%  option.fact_mino  = 0;
option.orbit_width  = 1;
option.icrh_model = 'Dumont-Vu';
%
option.equi_ppar = 3;
option.signe = 1;
option.cronos_regul= 2;
option.available_flux =  9.8  - (-7.82);  %WB
%  option.machine = 'WEST';
option.evolution = 0;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UAL writing control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


option.init_output_cpo = 0;
option.restart= '';
option.coreprof= 1;
option.coretransp= 1;
option.coresource_lhcd= 1;
option.coresource_eccd= 1;
option.coresource_icrh= 1;
option.coresource_nbicd= 1;
option.coresource_fusion= 1;
option.coreneutrals= 1;
option.coresource_radiation= 1;
option.coresource_cyclotron= 1;
option.neoclassic= 1;
option.coresource_full= 1;
option.equilibrium= 1;
option.grid_equi= 0;
option.scenario_occurrence= '';
option.coreprof_occurrence= '';
option.coretransp_occurrence= '';
option.coreneutrals_occurrence= '';
option.neoclassic_occurrence= '';
option.equilibrium_occurrence= '';
option.coresources_occurrence= '';

% default value for IMAS
option.COCOS  = 11;
option.signe = 1;
option.orientation = -1;


