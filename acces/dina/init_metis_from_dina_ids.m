% Function: init_metis_from_dina_ids
% This function allows to read DINA data stored in IMAS format and to
% create set of input data for METIS. The function creates also cache data
% useful for METIS external data.
%
% syntax 1:
%   z0dinput = init_metis_from_dina_ids(metis_reference,uri,nbt)
%
% syntax 2: 
%   z0dinput = init_metis_from_dina_ids(metis_reference,shot,run,user,machine,version,backend,nbt)
%
% inputs:
%
%   metis_reference        = reference METIS file (with full path) providing METIS parameters (z0dinput.option)
%   uri                    = IMAS URI for the DINA IMAS dataset
%   nbt                    = number of time slices in METIS simulation
%   shot                   = shot number
%   run                    = run numner
%   user                   = user name
%   machine                = machine name
%   version                = main IMAS version (char)
%   backend                = backend number (MDS+ = 12, HDF5 = 13)
%
% outputs:
%   
%   z0dinput               = METIS input data structure
%
%
% for testing :
%
% z0dinput = init_metis_from_dina_ids('./certification/metis/reference_NTM_ITER.mat','imas:hdf5?path=/home/ITER/artaudj1/public/imasdb/iter/3/105092/1');
% z0dinput = init_metis_from_dina_ids('./certification/metis/reference_NTM_ITER.mat',105092,1);

function z0dinput = init_metis_from_dina_ids(metis_reference,shot,run,user,machine,version,backend,nbt)

if ischar(shot)
    mode_uri = true;
    uri = shot;    
    if (nargin < 3) || isempty(run)
        nbt = 301;
    else
        nbt = run;
    end
    elem = split(uri,'/');
    try
        run  = str2double(elem{end});
    catch
        run = 0;
    end
    if isempty(run) || ~isfinite(run)
        run = 0;
    end
    try
        shot = str2double(elem{end-1});
    catch
        shot = 999999;
    end
    if isempty(shot) || ~isfinite(shot)
        shot = 1;
    end
    try
        version = elem{end-2};
    catch
        version = '4';
    end
    if isempty(version) 
        version = '4';
    elseif length(version) > 1
        version = version(1);
    end
    try
        machine = elem{end-3};
    catch
        machine =  'iter';
    end
    if isempty(machine)
        machine = 'test';
    end
    if contains(elem{1},'hdf5')
        backend = 13;
    else
        backend = 12;
    end
    if isempty(backend) || ~isfinite(backend)
        backend = 13;
    end
    user = getenv('USER');
    for k=3:length(elem)
        if contains(elem{k},'public')
            user = elem{k-1};
        end
    end
    if isempty(user)
      user = getenv('USER');
    end 
else
    mode_uri = false;
    if nargin < 3
        error('syntax: z0dinput = init_metis_from_dina_ids(metis_reference,shot,run,user,machine,version,backend,nbt);');
    end
    if (nargin < 4) || isempty(user)
        user = getenv('USER');
    end
    if (nargin < 5) || isempty(machine)
        machine = 'iter';
    end
    if (nargin < 6) || isempty(version)
        version = '3';
    end
    if (nargin < 7) || isempty(backend)
        backend = 13;
    end
    if (nargin < 8) || isempty(nbt)
        nbt = 301;
    end
end



% load metis data
data = load_metis_imas(metis_reference);

% read IDSs
if mode_uri
    fprintf('reading IDS from %s\n',uri);
    idx = imas_open(uri);
else
    fprintf('reading IDS from %s/.../%s/%s/%d/%d\n',user,machine,version,shot,run);
    idx = imas_open_env_backend(shot,run,user,machine,version,backend);
end
disp('Reading summary')
sm = ids_get(idx,'summary');
disp('Reading equilibrium')
eq = ids_get(idx,'equilibrium');
disp('Reading pulse_schedule')
ps = ids_get(idx,'pulse_schedule');
disp('Reading *_profiles')
try
    cp = ids_get(idx,'core_profiles');
catch
    cp = ids_get(idx,'plasma_profiles');
end
disp('Reading *_sources')
try
    cs = ids_get(idx,'core_sources');
catch
    cs = ids_get(idx,'plasma_sources');
end
imas_close(idx);

% get insteresting time slices
selected_time = get_interresting_time_from_dina(sm,nbt);

% init METIS data
z0dinput = zerod_scalaire(unique(selected_time),5.3,6.2,2,1.8,0.45,15e6,1e20,0,20e6,40e6,33e6,1.62,0.2,1,0.5,0,180);
z0dinput.machine = machine;
z0dinput.shot    = shot;
noms = fieldnames(z0dinput.option);
for k=1:length(noms)
    if isfield(data.z0dinput.option,noms{k})
        z0dinput.option.(noms{k})  = data.z0dinput.option.(noms{k});
    end
end

% resample reference
setappdata(0,'NONAN_EX',true);
setappdata(0,'EXTRAPOLATION_NEAREST_EX',true);

% read references in pulse_schedule
if isfield(ps.flux_control,'ip')
    z0dinput.cons.ip = interp1_ex(ps.flux_control.time,abs(ps.flux_control.ip.reference),z0dinput.cons.temps,'linear');
else
    z0dinput.cons.ip = interp1_ex(ps.flux_control.time,abs(ps.flux_control.i_plasma.reference),z0dinput.cons.temps,'linear');
end
if isfield(sm.global_quantities,'ip') && ~isempty(sm.global_quantities.ip.value)
        z0dinput.cons.ip = interp1_ex(sm.time,abs(sm.global_quantities.ip.value),z0dinput.cons.temps,'linear');
end

% read ion density
t_ion = ps.density_control.time;
nem   = zeros(size(z0dinput.cons.temps));
zeff  = zeros(size(z0dinput.cons.temps));
meff  = zeros(size(z0dinput.cons.temps));
label_list = {};
for k=1:length(ps.density_control.ion)
    label = ps.density_control.ion{k}.name;
    if isempty(label)
        label = ps.density_control.ion{k}.label;
    end   
    ion.(label)  =  interp1_ex(t_ion,ps.density_control.ion{k}.n_i_volume_average.reference,z0dinput.cons.temps,'linear');
    zion.(label) = double(ps.density_control.ion{k}.element{1}.z_n);
    A.(label)    = double(ps.density_control.ion{k}.element{1}.a);
    nem          = nem + zion.(label) .* ion.(label);
    zeff         = zeff + zion.(label) .* zion.(label) .* ion.(label);
    meff         = meff + A.(label) .* ion.(label);
    label_list{k} = label;
end
if ~isfield(ion,'H')
    ion.H = zeros(size(z0dinput.cons.temps));
    zion.H  = 1;
    A.H = 1;
end
if ~isfield(ion,'D')
    ion.D = zeros(size(z0dinput.cons.temps));
    zion.D  = 1;
end
if ~isfield(ion,'T')
    ion.T = zeros(size(z0dinput.cons.temps));
    zion.T  = 1;
    A.T = 3;
end
if ~isfield(ion,'He3') 
    ion.He3 = zeros(size(z0dinput.cons.temps));
    zion.He3  = 2;
    A.He3 = 3;    
end
if ~isfield(ion,'He4') 
    ion.He4 = zeros(size(z0dinput.cons.temps));
    zion.He4  = 2;
    A.He4 = 4;    
end
if ~isfield(ion,'He') 
    ion.He = zeros(size(z0dinput.cons.temps));
    zion.He  = 2;
    A.He = 4;    
end
zeff  = zeff ./ nem;
meff  = meff ./ nem;

% compute line averaged density over volume averaged density
nbar_cp = NaN * ones(length(cp.profiles_1d),1);
nem_cp = NaN * ones(length(cp.profiles_1d),1);
if ~isempty(cp.time)
    time_cp =cp.time;
else
    time_cp = NaN * ones(length(cp.profiles_1d),1);
end
for k = 1:length(cp.profiles_1d)
    nbar_cp(k) = trapz(cp.profiles_1d{k}.grid.rho_tor_norm , ...
        cp.profiles_1d{k}.electrons.density);
    nem_cp(k) =  trapz(cp.profiles_1d{k}.grid.rho_tor_norm , ...
        cp.profiles_1d{k}.electrons.density .* cp.profiles_1d{k}.grid.volume) ./ ...
                trapz(cp.profiles_1d{k}.grid.rho_tor_norm , ...
                      cp.profiles_1d{k}.grid.volume);
   if ~ isempty(cp.profiles_1d{k}.time)
       time_cp(k) = cp.profiles_1d{k}.time;
   end
end
rap = min(10,max(0.1,interp1_ex(time_cp,max(1,nbar_cp) ./ max(1,nem_cp),z0dinput.cons.temps,'linear')));
z0dinput.cons.zeff = rap .* zeff;
z0dinput.cons.nbar = rap .* nem;

if ~isempty(sm.line_average.zeff.value)
            z0dinput.cons.zeff = interp1_ex(sm.time,abs(sm.line_average.zeff.value),z0dinput.cons.temps,'linear');
elseif ~isempty(sm.volume_average.zeff.value)
            z0dinput.cons.zeff = rap .* interp1_ex(sm.time,abs(sm.volume_average.zeff.value),z0dinput.cons.temps,'linear');
end
z0dinput.cons.zeff = min(10,max(1+sqrt(eps),z0dinput.cons.zeff));
if ~isempty(sm.line_average.n_e.value)
            z0dinput.cons.nbar = interp1_ex(sm.time,abs(sm.line_average.n_e.value),z0dinput.cons.temps,'linear');
elseif ~isempty(sm.volume_average.n_e.value)
            z0dinput.cons.nbar = rap .* interp1_ex(sm.time,abs(sm.volume_average.n_e.value),z0dinput.cons.temps,'linear');
end
z0dinput.cons.nbar = max(1e13,z0dinput.cons.nbar);

% power
if ~isempty(ps.ec.time) && ~isempty(ps.ec.power_launched.reference)
    z0dinput.cons.pecrh = interp1_ex(ps.ec.time,ps.ec.power_launched.reference,z0dinput.cons.temps,'linear');
else
    z0dinput.cons.pecrh = zeros(size(z0dinput.cons.temps));
end
if ~isempty(ps.lh.time) && ~isempty(ps.lh.power.reference)
    z0dinput.cons.plh = interp1_ex(ps.lh.time,ps.lh.power.reference,z0dinput.cons.temps,'linear');
else
    z0dinput.cons.plh = zeros(size(z0dinput.cons.temps));
end
if ~isempty(ps.ic.time) && ~isempty(ps.ic.power.reference)
    z0dinput.cons.picrh = interp1_ex(ps.ic.time,ps.ic.power.reference,z0dinput.cons.temps,'linear');
else
    z0dinput.cons.picrh = zeros(size(z0dinput.cons.temps));
end
if ~isempty(ps.nbi.time) && ~isempty(ps.nbi.power.reference)
    z0dinput.cons.pnbi = interp1_ex(ps.nbi.time,ps.nbi.power.reference,z0dinput.cons.temps,'linear');
else
    z0dinput.cons.pnbi = zeros(size(z0dinput.cons.temps));
end

% plasam composition
ion.D = max(ion.D,ion.H .* z0dinput.option.natural_nD_o_nH);
z0dinput.cons.iso = ion.T ./ max(1,ion.D);
% search for gaz
if any(z0dinput.cons.iso > 1e-3)
    z0dinput.option.gaz = 3;
    if (mean(ion.H ./ (ion.D + ion.T)) > 1e-2)
        z0dinput.option.mino = 'H';
        z0dinput.option.cmino = mean(ion.H ./ (ion.D + ion.T));
    else
        z0dinput.option.mino  = 'T';
        z0dinput.option.cmino = 1;
    end    
elseif (mean(ion.H ./ ion.D) > 10)
    z0dinput.option.gaz = 1;
    z0dinput.option.mino = 'D';
    z0dinput.option.cmino = mean(ion.D ./ ion.H);
else
    z0dinput.option.gaz = 2;
    if (mean(ion.H ./ ion.D) > 1e-2)
        z0dinput.option.mino = 'H';
        z0dinput.option.cmino = mean(ion.H ./ ion.D);
    else
        z0dinput.option.mino  = 'He3';
        z0dinput.option.cmino = 3e-2;
    end    
end

% search for He
nhe = ion.He3 + ion.He4 + ion.He;
if any(nhe>0)
    z0dinput.option.frhe0 = max(0,min(0.25,nhe ./ nem));
end

% search for main impurities
nimp_z1 = zeros(size(z0dinput.cons.temps));
z1 = 0;
nimp_z2 = zeros(size(z0dinput.cons.temps));
z2 = 0;
for k=1:length(label_list)
    if (zion.(label_list{k}) > 2 ) && (zion.(label_list{k}) < 50 ) 
        test1 =  mean(zion.(label_list{k}) .* ion.(label_list{k})) >  ...
                 mean(nimp_z1 .* z1);
        test2 =  mean(zion.(label_list{k}) .* ion.(label_list{k})) >  ...
                 mean(nimp_z2 .* z2);
        if test1 && test2
            if mean(nimp_z1 .* z1) > mean(nimp_z2 .* z2)
                nimp_z2 = ion.(label_list{k});
                z2      = zion.(label_list{k});
            else
                nimp_z1 = ion.(label_list{k});
                z1      = zion.(label_list{k});
            end
        elseif test1
            nimp_z1 = ion.(label_list{k});
            z1      = zion.(label_list{k});
        elseif test2
            nimp_z2 = ion.(label_list{k});
            z2      = zion.(label_list{k});
        end
    end
end

if mean(nimp_z1) > mean(nimp_z2)
    z0dinput.option.zimp = z1;
    z0dinput.option.zmax = z2;
    z0dinput.option.rimp = max(0,min(1,mean(nimp_z2 ./ nimp_z1)));
else
    z0dinput.option.zimp = z2;
    z0dinput.option.zmax = z1;
    z0dinput.option.rimp = max(0,min(1,mean(nimp_z1 ./ nimp_z2)));
end

% Tungsten
cW = ion.W ./ nem;
mat = cat(1,1e6 .*ones(size(z0dinput.cons.temps))', z0dinput.cons.pecrh', ...
            z0dinput.cons.plh',z0dinput.cons.picrh', ...
            real(z0dinput.cons.pnbi)' + imag(z0dinput.cons.pnbi)');
rep = lscov(mat',cW);
%
z0dinput.option.cw_offset = 1e6 .* rep(1);
z0dinput.option.cw_factor = 0;
z0dinput.option.cw_ecrh   = rep(2);
z0dinput.option.cw_lhcd   = rep(3);
z0dinput.option.cw_icrh   = rep(4);
z0dinput.option.cw_nbi1   = rep(5);
z0dinput.option.cw_nbi2   = rep(5);




% here ftnbi computation is missing: need NBI sources filled
index_nbi = [];
for k = 1:length(cs.source)
    switch cs.source{k}.identifier.index
        case 2
            index_nbi = k;
            break;
    end
end
% no yet suficient information
z0dinput.cons.ftnbi(:) = 0;

% xece
index_ecrh = [];
for k = 1:length(cs.source)
    switch cs.source{k}.identifier.index
        case 3
            index_ecrh = k;
            break;
    end
end
if ~isempty(index_ecrh)
    if ~isempty(cs.time)
        time_ec = cs.time;
    else
        time_ec = [];
    end
    xece = 0.2 * ones(length(cs.source{index_ecrh}.profiles_1d),1);
    for k=1:length(cs.source{index_ecrh}.profiles_1d)
        if ~isempty(cs.source{index_ecrh}.profiles_1d{k}.electrons.energy)
            [~,l] = max(cs.source{index_ecrh}.profiles_1d{k}.electrons.energy);
            if ~isempty(l)
                xece(k) = cs.source{index_ecrh}.profiles_1d{k}.grid.rho_tor_norm(l(1));
            end
        elseif ~isempty(cs.source{index_ecrh}.profiles_1d{k}.j_parallel)
            [~,l] = max(cs.source{index_ecrh}.profiles_1d{k}.j_parallel);
            if ~isempty(l)
                xece(k) = cs.source{index_ecrh}.profiles_1d{k}.grid.rho_tor_norm(l(1));
            end
        end
        if ~isempty(cs.source{index_ecrh}.profiles_1d{k}.time)
            time_ec(k) = cs.source{index_ecrh}.profiles_1d{k}.time;
        end
    end
    z0dinput.cons.xece = interp1_ex(time_ec,xece,z0dinput.cons.temps,'linear');
else
    z0dinput.cons.xece(:) = 0.2;
end

% geometry
dina_lcfs = {};
for l = 1:length(z0dinput.cons.temps)
    k = find(eq.time >= z0dinput.cons.temps(l),1);
    if isempty(k)
        k = length(eq.time);
    end
    while isempty(eq.time_slice{k}.boundary.outline.r(:))
        if k > (length(eq.time) / 2)
             k = k - 1;
        else
             k = k + 1;
        end
    end
    dina_lcfs{l}.time = z0dinput.cons.temps(l);
    dina_lcfs{l}.ip   = z0dinput.cons.ip(l);
    dina_lcfs{l}.R    = eq.time_slice{k}.boundary.outline.r(:);
    dina_lcfs{l}.Z    = eq.time_slice{k}.boundary.outline.z(:);

end
% generate  new LCFS for METIS
z0dinput = z0separatrix_fromfbetometis(z0dinput,dina_lcfs,'none',201,1);

% toroidal field
z0dinput.geo.b0 = eq.vacuum_toroidal_field.r0 ./ z0dinput.geo.R .* ...
                  interp1_ex(eq.time,eq.vacuum_toroidal_field.b0,z0dinput.cons.temps,'linear');
sign_b0 = sign(mean(z0dinput.geo.b0));
z0dinput.geo.b0 = sign_b0 .* z0dinput.geo.b0;
if isfield(ps.flux_control,'ip')
    z0dinput.option.signe = sign(mean(ps.flux_control.ip.reference)) .* sign(mean(eq.vacuum_toroidal_field.b0));
else
    z0dinput.option.signe = sign(mean(ps.flux_control.i_plasma.reference)) .* sign(mean(eq.vacuum_toroidal_field.b0));
end
z0dinput.option.orientation = sign(mean(eq.vacuum_toroidal_field.b0));

% try to get info from summary
% compute Ka
Kaeq = ones(1,length(eq.time_slice));
for k=1:length(eq.time_slice)
    Kaeq(k) = max(1,eq.time_slice{k}.global_quantities.area ./ pi ./ eq.time_slice{k}.boundary.minor_radius .^ 2);
end
Ka = interp1_ex(eq.time,Kaeq',z0dinput.cons.temps,'linear');
% ITERH-98P(y,2)
ploss = interp1_ex(sm.time,sgolayfilt(medfilt1(sm.global_quantities.power_loss.value,3),1,71),z0dinput.cons.temps,'linear') / 1e6;
pin   = z0dinput.cons.pecrh + z0dinput.cons.plh + z0dinput.cons.picrh + z0dinput.cons.pnbi;
pin   = real(pin) + imag(pin) + interp1_ex(sm.time,sgolayfilt(medfilt1(sm.global_quantities.power_ohm.value,3),1,71),z0dinput.cons.temps,'linear');
pin   = pin ./ 1e6;
ep       = z0dinput.geo.a ./ z0dinput.geo.R;
tau_98   = 56.2e-3  .* abs(z0dinput.cons.ip/1e6) .^ 0.93 .* ...
           z0dinput.geo.b0 .^ 0.15 .* ...
           (z0dinput.cons.nbar / 1e19) .^ 0.41 .*  ...
           ploss .^ -0.69 .* ...
           z0dinput.geo.R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s
% ITERH-96P(th)
tau_96  = 23e-3  .* abs(z0dinput.cons.ip/1e6) .^ 0.96 .* ...
          z0dinput.geo.b0 .^ 0.03 .* ...
          (z0dinput.cons.nbar / 1e19) .^ 0.4 .* pin .^ -0.73 .* ...
           z0dinput.geo.R .^ 1.83 .* z0dinput.geo.K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s
hmode = (interp1_ex(sm.time,double(sm.global_quantities.h_mode.value),z0dinput.cons.temps,'linear') > 0);
tau_E = double(hmode) .* tau_98 + double(~hmode) .* tau_96;
pth  = double(hmode) .* ploss + double(~hmode) .* pin;
tau_E_dina = sgolayfilt(medfilt1(sm.global_quantities.tau_energy.value,3),1,71);
z0dinput.cons.hmore = max(0.3,min(3,interp1_ex(sm.time,tau_E_dina,z0dinput.cons.temps,'linear') ./ max(eps,tau_E)));
% change ploss computation
z0dinput.option.ploss_exp = 'max_power';


% make data structure for external IMAS data in METIS
noms  = {'summary','equilibrium','pulse_schedule','core_profiles','plasma_profiles','core_sources','plasma_sources'};
alias = {'sm','eq','ps','cp','cp','cs','cs'};
ref.shot_num   = shot;
ref.run_num    = run;
ref.ual_occ    = '0';
ref.ual_tok    = machine;
ref.ual_user   = user;
ref.ual_ver    = version;
setappdata(0,'IMAS_REFERENCE4EXTERNAL_METIS_DATA',ref);
setappdata(0,'UAL_TOKAMAK',machine);
setappdata(0,'UAL_USER',user);
setappdata(0,'UAL_DATAVERSION',version);
setappdata(0,'UAL_BACKEND',13);
for k=1:length(noms)
    nom_cache = sprintf('IMAS_EXTERNAL_METIS_DATA_%s',noms{k});
    cache_data.ref = ref;
    cache_data.data = eval(alias{k});
    setappdata(0,nom_cache,cache_data);
end

% exp0d data
z0dinput.exp0d.RR      = interp1_ex(sm.time,sm.global_quantities.resistance.value,z0dinput.cons.temps,'linear');
% z0dinput.exp0d.ane: [300x1 double]
% z0dinput.exp0d.ane_actual: [300x1 double]
% z0dinput.exp0d.ape: [300x1 double]
% z0dinput.exp0d.ate: [300x1 double]
z0dinput.exp0d.betan = interp1_ex(sm.time,sm.global_quantities.beta_tor_norm.value,z0dinput.cons.temps,'linear');
z0dinput.exp0d.betaptot = interp1_ex(sm.time,sm.global_quantities.beta_pol.value,z0dinput.cons.temps,'linear');
%z0dinput.exp0d.betap  
%z0dinput.exp0d.d0: [300x1 double]
%               dsol: [300x1 double]
%     dwbpdt: [300x1 double]
%       dwdt: [300x1 double]
% dwmagtordt: [300x1 double]
%     dwthdt: [300x1 double]
%       ecrit_he: [300x1 double]
%  ecrit_he3_DDn: [300x1 double]
% ecrit_he4_DHe3: [300x1 double]
%   ecrit_he4_DT: [300x1 double]
%     ecrit_icrh: [300x1 double]
%      ecrit_nbi: [300x1 double]
%    ecrit_p_DDp: [300x1 double]
%   ecrit_p_DHe3: [300x1 double]
%    ecrit_t_DDp: [300x1 double]
% eddy_current: [300x1 double]
%          edgeflux: [300x1 double]
%    efficiency: [300x1 double]
%     einj_icrh: [300x1 double]
%       einj_lh: [300x1 double]
% einj_nbi_icrh: [300x1 double]
%      esup_fus: [300x1 double]
%     esup_icrh: [300x1 double]
%       esup_lh: [300x1 double]
%      esup_nbi: [300x1 double]
%         etalh: [300x1 double]
%        etalh0: [300x1 double]
%        etalh1: [300x1 double]
%  firstorb_nbi: [300x1 double]
% flux_edge_cor: [300x1 double]
% frac_pellet: [300x1 double]
%    fracmino: [300x1 double]
% frloss_icrh: [300x1 double]
%       frnbi: [300x1 double]
%      fwcorr: [300x1 double]
%        harm: [300x1 double]
%        hitb: [300x1 double]
%        hmhd: [300x1 double]
%      ialign: [300x1 double]
z0dinput.exp0d.ip = interp1_ex(sm.time,sm.global_quantities.ip.value,z0dinput.cons.temps,'linear');
sign_ip = sign(mean(z0dinput.exp0d.ip));
z0dinput.exp0d.ip = sign_ip .* z0dinput.exp0d.ip;
z0dinput.exp0d.iboot = sign_ip .* interp1_ex(sm.time,sm.global_quantities.current_bootstrap.value,z0dinput.cons.temps,'linear'); 
%      ieccd: [300x1 double]
%       ifus: [300x1 double]
%      ifwcd: [300x1 double]
%        ilh: [300x1 double]
%     inbicd: [300x1 double]
%     indice_inv: [300x1 double]
z0dinput.exp0d.ini  = sign_ip .* interp1_ex(sm.time,sm.global_quantities.current_non_inductive.value,z0dinput.cons.temps,'linear'); 
z0dinput.exp0d.icd  = z0dinput.exp0d.ini - z0dinput.exp0d.iboot;
z0dinput.exp0d.iohm = sign_ip .* interp1_ex(sm.time,sm.global_quantities.current_ohm.value,z0dinput.cons.temps,'linear'); 
%      ipar: [300x1 double]
%      irun: [300x1 double]
%      j0fus: [300x1 double]
%      jxfus: [300x1 double]
% kidds_evol: [300x1 double]
if isfield(sm.global_quantities,'li')
    z0dinput.exp0d.li = interp1_ex(sm.time,sm.global_quantities.li.value,z0dinput.cons.temps,'linear');
else
    z0dinput.exp0d.li = interp1_ex(sm.time,sm.global_quantities.li_3.value,z0dinput.cons.temps,'linear');
end
z0dinput.exp0d.meff = meff;
z0dinput.exp0d.modeh = double(hmode);
% mu0_nbi: [300x1 double]
%z0dinput.exp0d.n0a  
z0dinput.exp0d.n1m    = ion.H + ion.D + ion.T;
z0dinput.exp0d.nDm    = ion.D;
z0dinput.exp0d.nTm    = ion.T;
%                nb: [300x1 double]
z0dinput.exp0d.nbar = z0dinput.cons.nbar;
%          nbar_nat: [300x1 double]
%         ndd: [300x1 double]
% ndd_nbi_nbi: [300x1 double]
%  ndd_nbi_th: [300x1 double]
%      ndd_th: [300x1 double]
z0dinput.exp0d.ne0       = interp1_ex(sm.time,sm.local.magnetic_axis.n_e.value,z0dinput.cons.temps,'linear'); 
z0dinput.exp0d.nebord    = interp1_ex(sm.time,sm.local.separatrix.n_e.value,z0dinput.cons.temps,'linear');
z0dinput.exp0d.negr      = interp1_ex(sm.time,sm.global_quantities.greenwald_fraction.value,z0dinput.cons.temps,'linear') .* z0dinput.exp0d.nbar;
%             nelim: [300x1 double]
z0dinput.exp0d.nem = nem;
%z0dinput.exp0d.neped = interp1_ex(sm.time,sm.local.pedestal.n_e.value,z0dinput.cons.temps,'linear');
 % nhe4m: [300x1 double]
%z0dinput.exp0d.nhem   = ion.He;
%    ni0: [300x1 double]
% nibord: [300x1 double]
%    nim: [300x1 double]
z0dinput.exp0d.nimpm  = max(ion.O , ion.Be);
% niped: [300x1 double]
% nmino: [300x1 double]
% nsat: [300x1 double]
z0dinput.exp0d.nwm  = ion.W;
z0dinput.exp0d.pbrem = abs(interp1_ex(sm.time,sm.global_quantities.power_bremsstrahlung.value,z0dinput.cons.temps,'linear')); 
if ~isempty(sm.global_quantities.power_synchrotron.value)
    z0dinput.exp0d.pcyclo = abs(interp1_ex(sm.time,sm.global_quantities.power_synchrotron.value,z0dinput.cons.temps,'linear'));
end
%  pddfus: [300x1 double]
% peakdiv: [300x1 double]
z0dinput.exp0d.pecrh = z0dinput.cons.pecrh;
%      pei: [300x1 double]
%      pel: [300x1 double]
%  pel_fus: [300x1 double]
% pel_icrh: [300x1 double]
%  pel_nbi: [300x1 double]
%    pelec: [300x1 double]
%     peri: [300x1 double]
z0dinput.exp0d.pfus = abs(interp1_ex(sm.time,sm.fusion.power.value,z0dinput.cons.temps,'linear'));
% pfus_loss: [300x1 double]
%  pfus_nbi: [300x1 double]
%   pfus_th: [300x1 double]
% phiplasma: [300x1 double]
z0dinput.exp0d.picrh = z0dinput.cons.picrh;
% picrh_th: [300x1 double]
%z0dinput.exp0d.pin  = abs(interp1_ex(sm.time,sm.global_quantities.power_steady.value,z0dinput.cons.temps,'linear'));
%      pion: [300x1 double]
%  pion_fus: [300x1 double]
% pion_icrh: [300x1 double]
%  pion_nbi: [300x1 double]
%    pioniz: [300x1 double]
%  pioniz_i: [300x1 double]
%      piqj: [300x1 double]
%    piqnbi: [300x1 double]
z0dinput.exp0d.plh = z0dinput.cons.plh;
% plh_th: [300x1 double]
% plhrip: [300x1 double]
% plhthr: [300x1 double]
%   plim: [300x1 double]
z0dinput.exp0d.ploss = ploss .* 1e6;
% plossl2h: [300x1 double]
z0dinput.exp0d.pnbi = z0dinput.cons.plh;
%   pnbi_th: [300x1 double]
z0dinput.exp0d.pohm = abs(interp1_ex(sm.time,sm.global_quantities.power_ohm.value,z0dinput.cons.temps,'linear'));
% poynting: [300x1 double]
%     pped: [300x1 double]
%  ppedmax: [300x1 double]
%z0dinput.exp0d.prad = abs(interp1_ex(sm.time,sm.global_quantities.power_line.value,z0dinput.cons.temps,'linear'));
%   pradsol: [300x1 double]
% priptherm: [300x1 double]
z0dinput.exp0d.pth =  pth .* 1e6;
%    pttfus: [300x1 double]
%        pw: [300x1 double]
z0dinput.exp0d.q0  = abs(interp1_ex(sm.time,sm.local.magnetic_axis.q.value,z0dinput.cons.temps,'linear'));
z0dinput.exp0d.q95 = abs(interp1_ex(sm.time,sm.global_quantities.q_95.value,z0dinput.cons.temps,'linear'));
z0dinput.exp0d.qa = abs(interp1_ex(sm.time,sm.local.separatrix.q.value,z0dinput.cons.temps,'linear'));
% qeff: [300x1 double]
% qmin: [300x1 double]
% rm: [300x1 double]
%       rres: [300x1 double]
%     salpha: [300x1 double]
% salpha_he3: [300x1 double]
% salpha_he4: [300x1 double]
%   salpha_n: [300x1 double]
%   salpha_p: [300x1 double]
%   salpha_t: [300x1 double]
%       sext: [300x1 double]
%        slh: [300x1 double]
% sn0fr: [300x1 double]
%  snbi: [300x1 double]
%    sp: [300x1 double]
%   stf: [300x1 double]
z0dinput.exp0d.taue  = interp1_ex(sm.time,sm.global_quantities.tau_energy.value,z0dinput.cons.temps,'linear');
% taue_alt: [300x1 double]
%    tauee: [300x1 double]
%    tauei: [300x1 double]
%     tauh: [300x1 double]
%             z0dinput.exp0d.tauhe: [300x1 double]
% tauhe_h: [300x1 double]
% tauhe_l: [300x1 double]
%   tauii: [300x1 double]
z0dinput.exp0d.tauip = interp1_ex(sm.time,sm.global_quantities.tau_resistive.value,z0dinput.cons.temps,'linear');
% tauj: [300x1 double]
% taup: [300x1 double]
%   taus_he: [300x1 double]
% taus_icrh: [300x1 double]
%  taus_nbi: [300x1 double]
%    tauthl: [300x1 double]
z0dinput.exp0d.te0    = interp1_ex(sm.time,sm.local.magnetic_axis.t_e.value,z0dinput.cons.temps,'linear'); 
z0dinput.exp0d.tebord = interp1_ex(sm.time,sm.local.separatrix.t_e.value,z0dinput.cons.temps,'linear');
% telim: [300x1 double]
z0dinput.exp0d.tem = interp1_ex(sm.time,sm.volume_average.t_e.value,z0dinput.cons.temps,'linear');
% temps: [300x1 double]
%  teped: [300x1 double]
z0dinput.exp0d.tibord = interp1_ex(sm.time,sm.local.separatrix.t_i_average.value,z0dinput.cons.temps,'linear');
 % tiped: [300x1 double]
z0dinput.exp0d.tite = interp1_ex(sm.time,sm.volume_average.t_i_average.value,z0dinput.cons.temps,'linear') ./ max(eps,z0dinput.exp0d.tem);
% temps: [300x1 double]
z0dinput.exp0d.vloop = sign_ip .* interp1_ex(sm.time,sm.global_quantities.v_loop.value,z0dinput.cons.temps,'linear');
z0dinput.exp0d.edgeflux = - sign_ip .* interp1_ex(sm.time,cumtrapz(sm.time,sm.global_quantities.v_loop.value),z0dinput.cons.temps,'linear');
% vmes: [300x1 double]
z0dinput.exp0d.vp = interp1_ex(sm.time,sm.global_quantities.volume.value,z0dinput.cons.temps,'linear');
%       w: [300x1 double]
%     wbp: [300x1 double]
%    wdia: [300x1 double]
% wmagtor: [300x1 double]
%    wrad: [300x1 double]
%    wrlw: [300x1 double]
%    wrot: [300x1 double]
z0dinput.exp0d.wth = interp1_ex(sm.time,sm.global_quantities.energy_thermal.value,z0dinput.cons.temps,'linear');
z0dinput.exp0d.xeccd = z0dinput.cons.xece;
% xfus: [300x1 double]
% xitb: [300x1 double]
%  xlh: [300x1 double]
% xnbi: [300x1 double]
% xped: [300x1 double]
z0dinput.exp0d.xpoint = double(interp1_ex(sm.time,double(sm.boundary.type.value),z0dinput.cons.temps,'linear') > 0);
% xres: [300x1 double]
z0dinput.exp0d.zeff = z0dinput.cons.zeff;

% update orientation and sign
z0dinput.option.signe = sign_ip .* sign_b0;
z0dinput.option.orientation = sign_b0;

% new evaluation of hmore as the first one does not provide correct value
wth_dina = interp1_ex(sm.time,sgolayfilt(medfilt1(sm.global_quantities.energy_thermal.value,3),1,71),z0dinput.cons.temps,'linear');
z0dinput.cons.hmore = max(0.3,min(3, wth_dina ./ max(eps,interp1_ex(sm.time,tau_E_dina,z0dinput.cons.temps,'linear').* (pin * 1e6 - z0dxdt(wth_dina,z0dinput.cons.temps)))));

% compute poloidal flux offset
k = find(eq.time >= z0dinput.cons.temps(1),1);
if isempty(k)
    k= length(eq.time);
end
while(all(eq.time_slice{k}.profiles_1d.psi == 0))
    k = k + 1;
end
psi_offset = sign_ip .*  eq.time_slice{k}.profiles_1d.psi(end) / 2 /pi;
fprintf('Poloidal flux offset for METIS MUSCLE3 actor = %g\n',psi_offset);
fid = fopen('PSI_OFFSET','w');
fprintf(fid,'export PSI_OFFSET="%6.3g"\n',psi_offset);
fclose(fid);




disp('Data ready !')
