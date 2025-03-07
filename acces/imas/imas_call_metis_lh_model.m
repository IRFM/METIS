% function imas_call_metis_lh_model: psuedo IMAS actor to call METIS LHCD  model.
%
% syntax: 
%  1/ declaration of parameters:
%           [parameters,ref_model] = imas_call_metis_lh_model;
% 
%     Return parameters and a model for IMAS reference.
%     Also create xsd and xml models for IMAS (imas_call_metis_lh_model.xsd
%     & imas_call_metis_lh_model.xml).
%
% 2/ call of LHCD model
%
%  core_sources = imas_call_metis_lh_model(imas_ref_in,imas_ref_out,codeparam_filename, ... 
%                        equilibrium,core_profiles,lh_antennas,core_sources);
%
%    Call METIS LHCD model and return core_sources IDS. If imas_ref_out is
%    not empty, write core_sources IDS in IMAS database.
%
%    LH wave parameters (parallel index and frequency) and LH power are
%    read, first, in lh_antennas IDS. If this IDS is not available, 
%    parallel index and frequency are read from parameters and LH power
%    from core_sources.
%
% input:
%    imas_ref_in   = input reference structure following ref_model. Information for
%                    IMAS input data (tokamak, shot,run, ...)
%
%    imas_ref_out  = output reference structure following ref_model. Information for
%                    IMAS output data (tokamak, shot,run, ...)
%
%    codeparam_filename = code parameters. Either using parameters.valeur
%    format or imas_call_metis_lh_model.xml string or file. See METIS
%    documentation for parameters information.
%
% optional input:
%
% if IDS is available in Matlab caller workspace, the following IDS cans be
% provided in input : equilibrium,core_profiles,lh_antennas,core_sources;
% preventing the program to re-read it.
%
% output :
%
%   core_sources = core source IDS.
%
% example : see test_imas_call_metis_lh_model.m
%
% mailto: jean-francois.artaud@cea.fr
%



%test:
% [info,ref_model] = imas_call_metis_lh_model;
% option = info.valeur;
% imas_ref_in = ref_model;
% imas_ref_in.shot = 54719;
% imas_ref_in.run  = 20;
% imas_ref_in.tokamak  = 'west';
% imas_ref_in.occurrence.equilibrium  = 0;
% imas_ref_out = imas_ref_in;
% imas_ref_out.run = 21;
% option.gaz      = 2;
% option.etalh    = 0.8;
% option.npar0    = 1.8;
% option.freqlh    = 3.7;
% option.wlh      = 0.58;
% option.npar_neg = -4; 
% option.fupshift = 1;
% core_sources    = imas_call_metis_lh_model(imas_ref_in,imas_ref_out,option);

function [core_sources,ref_model] = imas_call_metis_lh_model(imas_ref_in,imas_ref_out,codeparam_filename,equilibrium,core_profiles,lh_antennas,core_sources)

% test input
if nargin <= 1
   % auto declaration mode + creation of template xlm and xsd 
   zs = zerod_param;
   % restriction
   list = {'gaz','etalh','wlh','npar0','npar_neg','freqlh','fupshift'};
   option = [];
   noms = fieldnames(zs);
   for k=1:length(noms)
       if isstruct(zs.(noms{k}))
           for l=1:length(list)
               if isfield(zs.(noms{k}),list{l})
                   option.(noms{k}).(list{l}) = zs.(noms{k}).(list{l});
               end
           end
       end
   end
   core_sources = option;
   if nargout > 1
       ref_model = complete_ref;
   end
   if nargin == 1
       ref_model = complete_ref;
       return;
   end
   % create xml and xsd template
   module_xsd_make('imas_call_metis_lh_model',1);
   % end of declaration
   return
end
if nargin < 3 
    error('syntax:  core_sources = imas_call_metis_lh_model(imas_ref_in,imas_ref_out,codeparam_filename,[equilibrium,core_profile,lh_antennas,core_sources])');
end
if nargin < 4
    equilibrium = [];
end
if nargin < 5
    core_profiles = [];
end
if nargin < 6
    lh_antennas = [];
end
if nargin < 7
    core_sources = [];
end

% complete imas_ref if needed
imas_ref_in = complete_ref(imas_ref_in);
if (nargout == 0) || ~isempty(imas_ref_out)
    imas_ref_out = complete_ref(imas_ref_out,imas_ref_in);
else
    imas_ref_out = [];
end
% read or decode code param
info   = imas_call_metis_lh_model(1);
option = info.valeur;
if  ~isempty(codeparam_filename) && isstruct(codeparam_filename)
    % codeparam_filename est une structure matlab
    noms = fieldnames(codeparam_filename);
    for k=1:length(noms)
        if isfield(option,noms{k})
            option.(noms{k}) = codeparam_filename.(noms{k});
        end
    end
    fprintf('imas_call_metis_lh_model using parameters from input :\n')
    disp(option)
    
elseif ~isempty(codeparam_filename)
    % codeparam_filename designe un fichier xml
    % lecture du fichier
    fid = fopen(codeparam_filename,'r');
    if fid > 0
        codeparam = char(fread(fid,Inf,'char')');
        fclose(fid);
        % add codeparam to options
        if ~isempty(codeparam)
            info  = xml_read(codeparam_filename);
            noms = fieldnames(info);
            for k=1:length(noms)
                if isfield(option,noms{k})
                    option.(noms{k}) = info.(noms{k});
                end
            end
        end
        fprintf('imas_call_metis_lh_model using parameters from %s :\n',codeparam_filename)
        disp(option)
    elseif ~isempty(codeparam_filename)
        
        % codeparam_filename est une chaine xml       
        tnp = tempname;
        fid = fopen(tnp,'w');
        fprintf(fid,'%s\n',codeparam_filename);
        fclose(fid);
        % add codeparam to options
        info  = xml_read(tnp);
        noms = fieldnames(info);
        for k=1:length(noms)
            if isfield(option,noms{k})
                option.(noms{k}) = info.(noms{k});
            end
        end
        fprintf('METIS4IMAS using parameters from input :\n')
        disp(option)
        delete(tnp);
    else
        error(sprintf('Unable to read codeparam file %s',codeparam_filename));
    end
    
end

% read input IDS if needed
if isempty(equilibrium) || isempty(core_profiles) ||  isempty(lh_antennas) ||  isempty(core_sources)
    disp('openning IMAS database:')
    idx = imas_open_env_backend(imas_ref_in.shot, imas_ref_in.run, imas_ref_in.user, imas_ref_in.tokamak, imas_ref_in.version, imas_ref_in.backend_id);
    if isempty(idx) || ~isfinite(idx) || (idx < 0)
        error('unable to open database');
    end
    if isempty(equilibrium)
        fprintf('-> reading equilibrium IDS occurrence %d\n',imas_ref_in.occurrence.equilibrium)
        switch imas_ref_in.occurrence.equilibrium
            case 0
                equilibrium = ids_get(idx,'equilibrium');
            otherwise
                equilibrium = ids_get(idx,sprintf('equilibrium/%d',imas_ref_in.occurrence.equilibrium));
        end
        
        if isempty(equilibrium) || isempty(equilibrium.time_slice)
            error('equilibrium IDS is empty');
        end
    end
    if isempty(core_profiles)
        fprintf('-> reading core_profiles IDS occurrence %d\n',imas_ref_in.occurrence.core_profiles)
        switch imas_ref_in.occurrence.core_profiles
            case 0
                core_profiles = ids_get(idx,'core_profiles');
            otherwise
                core_profiles = ids_get(idx,sprintf('core_profiles/%d',imas_ref_in.occurrence.core_profiles));
        end
        if isempty(core_profiles) || isempty(core_profiles.profiles_1d)
           error('core_profiles IDS is empty'); 
        end
    end
    if isempty(lh_antennas)
        fprintf('-> reading lh_antennas IDS occurrence %d\n',imas_ref_in.occurrence.lh_antennas)
        switch imas_ref_in.occurrence.lh_antennas
            case 0
                 lh_antennas = ids_get(idx,'lh_antennas');
            otherwise
                 lh_antennas = ids_get(idx,sprintf('lh_antennas/%d',imas_ref_in.occurrence.lh_antennas));
        end
        if isempty(lh_antennas) || isempty(lh_antennas.power.data)
           disp('lh_antennas IDS is empty'); 
           no_lh_antennas = true;
        else
           no_lh_antennas = false;
        end
    end
    % reading existing cores_sources
    if isempty(core_sources)
        fprintf('-> reading core_sources IDS occurrence %d\n',imas_ref_in.occurrence.core_sources)
        switch imas_ref_in.occurrence.core_sources
            case 0
                core_sources = ids_get(idx,'core_sources');
            otherwise
                core_sources = ids_get(idx,sprintf('core_sources/%d',imas_ref_in.occurrence.core_sources));
        end
        if isempty(core_sources) 
            disp('No existing core_sources IDS: initializing core_sources IDS');
            core_sources = ids_gen('core_sources');
        end
    end
    disp('closing IMAS database');
    imas_close(idx);
else
    if isempty(lh_antennas) || isempty(lh_antennas.power.data)
        disp('lh_antennas IDS is empty');
        no_lh_antennas = true;
    else
        no_lh_antennas = false;
    end
end    

% get time vector from core_profiles
if ~isempty(core_profiles.time)
    time_cp = core_profiles.time;
else
    time_cp = (1:length(core_profiles.profiles_1d))';
    for k=1:length(core_profiles.profiles_1d)
        time_cp(k) = core_profiles.profiles_1d{k}.time;
    end
end
% get time vector for equilibrium
if isempty(equilibrium.time)
    equilibrium.time = NaN * ones(length(equilibrium.time_slice),1);
end
if any(~isfinite(equilibrium.time))
    for kl=1:length(equilibrium.time_slice)
        if ~isempty(equilibrium.time_slice{kl}.time) && isfinite(equilibrium.time_slice{kl}.time)
            equilibrium.time(kl) = equilibrium.time_slice{kl}.time;
        end
    end
end
% count valid equilibrium
if length(equilibrium.code.output_flag) ~= length(equilibrium.time)
    if length(equilibrium.code.output_flag) > length(equilibrium.time)
        % this is a METIS case before bug correction
        equilibrium.code.output_flag = zeros(size(equilibrium.time));
    end 
end
nb_equi   = sum(double( (equilibrium.code.output_flag> -1) & isfinite(equilibrium.time)));
time_equi = NaN * ones(nb_equi,1);
% get equilibrium timed data
xli = [];
count = 1;
for k=1:length(equilibrium.time_slice)
    if (equilibrium.code.output_flag(k) > -1) && isfinite(equilibrium.time(k))
        % initialise data
        if isempty(xli)
            xli     = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            rin     = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            rout    = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            aminor  = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            Raxe    = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            epsi    = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            psi     = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            fdia    = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            q       = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            volume  = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            surface = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            area    = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            rho_tor = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            phi     = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            Rim     = ones(nb_equi,length(equilibrium.time_slice{k}.profiles_1d.psi));
            ip      = ones(nb_equi,1);
        end
        if ~isempty(equilibrium.time)
            time_equi(count) = equilibrium.time(k);
        else
            time_equi(count) = equilibrium.time_slice{k}.time;
        end
        rin(count,:)     = equilibrium.time_slice{k}.profiles_1d.r_inboard;
        rout(count,:)    = equilibrium.time_slice{k}.profiles_1d.r_outboard;
        aminor(count,:)  = (rout(count,:) - rin(count,:)) ./ 2;
        Raxe(count,:)    = (rout(count,:) + rin(count,:)) ./ 2;
        epsi(count,:)    = aminor(count,:) ./ Raxe(count,:);
        xli(count,:)     = aminor(count,:) ./ max(aminor(count,:));
        %
        % check psi orientation
        psi(count,:)     = - equilibrium.time_slice{k}.profiles_1d.psi / 2 / pi;      % 2*pi + sign ?
        %
        fdia(count,:)    = abs(equilibrium.time_slice{k}.profiles_1d.f);
        q(count,:)       = abs(equilibrium.time_slice{k}.profiles_1d.q);
        volume(count,:)  = equilibrium.time_slice{k}.profiles_1d.volume;
        surface(count,:) = equilibrium.time_slice{k}.profiles_1d.surface;
        area(count,:)    = equilibrium.time_slice{k}.profiles_1d.area;
        rho_tor(count,:) = equilibrium.time_slice{k}.profiles_1d.rho_tor;
        phi(count,:)     = equilibrium.time_slice{k}.profiles_1d.phi;
        Rim(count,:)     = equilibrium.time_slice{k}.profiles_1d.gm9;
        %
        ip(count)        = abs(equilibrium.time_slice{k}.global_quantities.ip);
        % next valid equilibrium
        count = count + 1;
    end
end

% compute epar from equi if not available in core_profiles
dpsidt_equi  = diff(psi,1,1) ./ (diff(time_equi,1,1) * ones(1,size(epsi,2)));
dpsidt_equi  = (interp1(time_equi(2:end) ,dpsidt_equi,time_equi,'linear',0) +  ...
               interp1(time_equi(1:end-1) ,dpsidt_equi,time_equi,'linear',0)) ./ 2;
%
ephi    = -Rim .*  dpsidt_equi;

% physical constants
phys = cphys;

% search for existing sources
index_sources = [];
for k=1:length(core_sources.source)
    if ~isempty(core_sources.source{k}.identifier.index)
        index_sources(end+1) = core_sources.source{k}.identifier.index;
    end
end
if isempty(index_sources)
    index_lh = 1;
    init_core_sources = true;
else
    index_lh = find(index_sources == 4);
    if isempty(index_lh)
        index_lh = length(core_sources.source) + 1;
        sources_lh_mem = [];
    else
        sources_lh_mem = core_sources.source{index_lh};
        sources_lh_mem.identifier.name  = 'lh_previous';
        sources_lh_mem.identifier.index = 99999;
        sources_lh_mem.identifier.description = sprintf('previous computation of LH: %s',sources_lh_mem.identifier.description);
        if  (core_sources.ids_properties.homogeneous_time == 1) && isempty(core_sources.time);
            sources_lh_mem_time = core_sources.time;
        else
            sources_lh_mem_time = NaN * ones(length(sources_lh_mem.global_quantities),1);
        end
        sources_lh_mem_plh  = NaN * ones(length(sources_lh_mem.global_quantities),1);
        for kl=1:length(sources_lh_mem.global_quantities)
            if ~isempty(sources_lh_mem.global_quantities{kl}.time)
                sources_lh_mem_time(kl) = sources_lh_mem.global_quantities{kl}.time;
            end
            sources_lh_mem_plh(kl) = sources_lh_mem.global_quantities{kl}.power;
        end
    end
    init_core_sources = false;
end

% initialisation of ouput IDS
if init_core_sources
    core_sources.ids_properties.homogeneous_time = 1;
    core_sources.ids_properties.comment ='LHCD source computed with METIS LHCD model';
    core_sources.ids_properties.source  = which(mfilename);
    core_sources.ids_properties.provider  = getenv('USER');
    core_sources.ids_properties.creation_date = sprintf('%s (julian date in second = %f)',datestr(now,'dd-mmm-yyyy HH:MM:SS'),clock2julday);
    % code info
    [repository,commit,version] = metis_info_imas;
    core_sources.code.name        = mfilename;
    core_sources.code.version     = version;
    core_sources.code.commit      = commit;
    core_sources.code.repository  = repository;
    %
    core_sources.time = time_equi;
    time_loop = time_equi;
else
    if  core_sources.ids_properties.homogeneous_time == 1
        time_loop = core_sources.time;
    else
        time_loop = time_equi;
    end
    core_sources.ids_properties.comment =sprintf('%s & added = LHCD source computed with METIS LHCD model',core_sources.ids_properties.comment);
    core_sources.ids_properties.source  = sprintf('%s & LHCD from %s',core_sources.ids_properties.source,which(mfilename));
    core_sources.ids_properties.provider  = sprintf('%s & LHCD from %s',core_sources.ids_properties.provider,getenv('USER'));
    core_sources.ids_properties.creation_date = sprintf('%s & LHCD at %s', ...
           core_sources.ids_properties.creation_date, ...
           sprintf('%s (julian date in second = %f)',datestr(now,'dd-mmm-yyyy HH:MM:SS'),clock2julday));
    % code info
    [repository,commit,version] = metis_info_imas;
    core_sources.code.name        = sprintf('%s & LHCD computed with %s',core_sources.code.name,mfilename);
    core_sources.code.version     = sprintf('%s & LHCD version =  %s',core_sources.code.version,version);
    core_sources.code.commit      = sprintf('%s & LHCD commit =  %s',core_sources.code.commit,commit);
    core_sources.code.repository  = sprintf('%s & LHCD repository =  %s',core_sources.code.repository,repository);
end
%
% declare LHCD source
%
core_sources.source{index_lh}.identifier.name  = 'lh';
core_sources.source{index_lh}.identifier.index = 4;
core_sources.source{index_lh}.identifier.description = 'LHCD source of power, current and momentum from METIS model';

% structure model
profiles_1d_model       = ids_allocate('core_sources','source/profiles_1d',1);
profiles_1d_model       = profiles_1d_model{1};
global_quantities_model = ids_allocate('core_sources','source/global_quantities',1);
global_quantities_model = global_quantities_model{1};


% call of the model
option_mem = option;
% loop on time
for k=1:length(time_loop)
    cons.temps    = time_loop(k);
    cons.ip       = ip(k);
    cons.plh      = 0;
    npar_peak     = 0;
    frequency     = 0;
    if no_lh_antennas
        cons.plh  = interp1(sources_lh_mem_time,sources_lh_mem_plh,cons.temps,'nearest',0); 
        npar_peak = option_mem.npar0;
        frequency = option_mem.freqlh * 1e9;
    else
        for l=1:length(lh_antennas.antenna)
            if ~isempty(lh_antennas.antenna{l}.power_forward.data)
                plh_1     = interp1(lh_antennas.antenna{l}.power_forward.time, ...
                    lh_antennas.antenna{l}.power_forward.data,time_loop(k),'nearest',0);
                cons.plh  = cons.plh + max(eps,plh_1);
                ind_peak = find(lh_antennas.antenna{l}.n_parallel_peak.data>=1);
                if length(ind_peak) < 3
                    npar_peak = npar_peak + option_mem.npar0 .* max(eps,plh_1);
                    disp('no data in n_parallel_peak')
                else
                    npar_peak = npar_peak + interp1(lh_antennas.antenna{l}.n_parallel_peak.time(ind_peak), ...
                        lh_antennas.antenna{l}.n_parallel_peak.data(ind_peak),time_loop(k),'nearest','extrap') .* max(eps,plh_1);
                end
                if isempty(lh_antennas.antenna{l}.frequency) || lh_antennas.antenna{l}.frequency < 1e6
                    disp('no data ifor wave frequency')
                    frequency = frequency + (option_mem.freqlh *1e9) .* max(eps,plh_1);
                else
                    frequency = frequency + lh_antennas.antenna{l}.frequency .* max(eps,plh_1);
                end
            end
        end
        npar_peak = npar_peak ./ cons.plh;
        frequency = frequency ./ cons.plh;
    end
    %
    % index in core profile
    ind_cp = find(time_cp >= time_loop(k),1);
    if isempty(ind_cp)
        ind_cp = length(time_cp);
    end
    % index equilibrium
    ind_eq = find(time_equi >= time_loop(k),1);
    if isempty(ind_eq)
        ind_eq = length(time_eq);
    end
    
    profil.xli    = xli(ind_eq,:);
    profil.Raxe   = Raxe(ind_eq,:);
    profil.epsi   = epsi(ind_eq,:);
    profil.fdia   = fdia(ind_eq,:);
    profil.qjli   = q(ind_eq,:);
    if ~isempty(core_profiles.profiles_1d{ind_cp}.electrons.density_thermal)
        profil.nep    = core_profiles.profiles_1d{ind_cp}.electrons.density_thermal(:)';
    else
        profil.nep    = core_profiles.profiles_1d{ind_cp}.electrons.density(:)';        
    end
    profil.tep    = core_profiles.profiles_1d{ind_cp}.electrons.temperature(:)';
    profil.rmx    = rho_tor(ind_eq,:);
    profil.spr    = pdederive(profil.xli,area(ind_eq,:),2,2,2,1);
    profil.vpr    = pdederive(profil.xli,volume(ind_eq,:),2,2,2,1);

    if ~isempty(core_profiles.profiles_1d{ind_cp}.zeff)
        profil.zeff   = core_profiles.profiles_1d{ind_cp}.zeff(:)';
    else
        profil.zeff   = core_profiles.global_quantities.z_eff_resistive(min(ind_cp,legnth(core_profiles.global_quantities.z_eff_resistive)));
    end
    % check sign !
    if ~isempty(core_profiles.profiles_1d{ind_cp}.e_field.parallel)       
            profil.epar   = core_profiles.profiles_1d{ind_cp}.e_field.parallel(:)';
    elseif  ~isempty(core_profiles.profiles_1d{ind_cp}.e_field_parallel)
            profil.epar   = core_profiles.profiles_1d{ind_cp}.e_field_parallel(:)';
    else
            profil.epar   = ephi(ind_eq,:);  % approximation
    end
    % problem ?
    profil.epar   = abs(profil.epar);
    %
    % complete option    
    option.freqlh      = frequency / 1e9;
    option.npar0       = npar_peak;
    
    % call model
    zz = zeros(size(xli(k,:)));
    if cons.plh > 1e4
        [time_out,plh_tot,ilh,x_out,plh,jlh,efficiency] = external_call_metis_lh_model(cons,profil,option);
        % momentum
        slh     = - sign(option.etalh) .* Raxe(ind_eq,end) .* option.npar0 ./ phys.c .* plh_tot .* (2 .* abs(option.etalh) - 1);
        rot_lh  = jlh  .* ((slh ./ ilh) * ones(1,size(jlh,2))) ./ 2 ./ pi ./ Raxe(ind_eq,end);
    else
        time_out      = time_loop(k);
        plh_tot       = 0;
        ilh           = 0;
        x_out         = xli(ind_eq,:);
        plh           = zz;
        jlh           = zz;
        efficiency    = 0;
        rot_lh        = zz;
    end
    
    
    % map IDS core source
    profiles_1d                              = profiles_1d_model;
	profiles_1d.grid.rho_tor_norm            = rho_tor(ind_eq,:) ./ max(rho_tor(ind_eq,:));
	profiles_1d.grid.rho_tor                 = rho_tor(ind_eq,:);
	profiles_1d.grid.psi                     = -2 .* pi .* psi(ind_eq,:);
	profiles_1d.grid.psi_magnetic_axis       = -2 .* pi .* psi(ind_eq,1);
	profiles_1d.grid.psi_boundary            = -2 .* pi .* psi(ind_eq,end);
    profiles_1d.grid.rho_pol_norm            = sqrt((psi(ind_eq,:)-psi(ind_eq,1)) / (psi(ind_eq,end)-psi(ind_eq,1)));
	profiles_1d.grid.volume                  = volume(ind_eq,:);
	profiles_1d.grid.area                    = area(ind_eq,:);
	profiles_1d.grid.surface                 = surface(ind_eq,:);

	profiles_1d.electrons.particles          = zz;
	profiles_1d.electrons.energy             = plh;
    profiles_1d.total_ion_energy             = zz;
    profiles_1d.momentum_tor                 = rot_lh;
    profiles_1d.j_parallel                   = -jlh;
    profiles_1d.conductivity_parallel        = zz;
    profiles_1d.time                         = time_out;
        
    % donnees complementaires
	profiles_1d.current_parallel_inside       = cumtrapz(xli(ind_eq,:),profil.spr .* jlh,2);
	profiles_1d.torque_tor_inside             = cumtrapz(xli(ind_eq,:),profil.vpr .* rot_lh,2);
	profiles_1d.total_ion_power_inside        = zz;
	profiles_1d.electrons.particles_inside    = zz;
	profiles_1d.electrons.power_inside        = cumtrapz(xli(ind_eq,:),profil.vpr .* plh,2);
    %            
	core_sources.source{index_lh}.profiles_1d{k}     = profiles_1d;

    % adding global data
    global_quantities                         = global_quantities_model;
    global_quantities.total_ion_power         = 0;
    global_quantities.torque_tor              = trapz(xli(ind_eq,:),profil.vpr .* rot_lh,2);
    global_quantities.current_parallel        = -ilh;
    global_quantities.time                    = time_out;
    global_quantities.electrons.particles     = 0;
    global_quantities.electrons.power         = plh_tot;
    global_quantities.power                   = plh_tot;
    %
    core_sources.source{index_lh}.global_quantities{k} = global_quantities;
    
end

% for debgug
if ~isempty(sources_lh_mem)
    core_sources.source{end+1} = sources_lh_mem;
end

% write IDS if requiered 
if ~isempty(imas_ref_out)
    fprintf('try to openning existing IMAS database:')
    try
        idx = imas_open_env_backend(imas_ref_out.shot, imas_ref_out.run, imas_ref_out.user, imas_ref_out.tokamak, imas_ref_out.version, imas_ref_out.backend_id);
    catch
        idx = [];
    end
    if ~isempty(idx) && (idx > 0)
         fprintf('reuse existing database\n')   
    else
         fprintf('no existing database entry -> creating a new database entry:')   
         try
             idx = imas_create_env_backend(imas_ref_out.shot, imas_ref_out.run, imas_ref_out.user, imas_ref_out.tokamak, imas_ref_out.version, imas_ref_out.backend_id);
         catch
             idx = [];
         end      
         if ~isempty(idx) && (idx > 0)
             fprintf('writing in new entry\n')
         else
            error('unable to create new entry in database'); 
         end
    end
    disp('-> wrinting core_sources IDS')
    switch imas_ref_out.occurrence.core_sources
        case 0
            ids_put(idx,'core_sources',core_sources);
        otherwise
            ids_put(idx,sprintf('core_sources/%d',imas_ref_out.occurrence.core_sources),core_sources);
    end
    disp('closing IMAS database');
    imas_close(idx);
end
% end of function





function ref = complete_ref(ref,ref_model)

% reference content
%   shot:    shot number.
%   run:     run number.
%   user:    owner of the database
%   tokamak: name of the machine
%   version: version of the datastructure
%   backend_id: backend identifier
if nargin <= 1
    ref_model.shot = NaN;
    ref_model.run  = 0;
    ref_model.user = getenv('USER');
    ref_model.tokamak = '';
    ref_model.version = '3';
    ref_model.backend_id = 13;
    ref_model.occurrence.equilibrium     = 0;
    ref_model.occurrence.core_profiles   = 0;
    ref_model.occurrence.lh_antennas     = 0;
    ref_model.occurrence.core_sources    = 0;
end

if (nargin == 0) || isempty(ref)
    ref = ref_model;
    return
end
    
noms = fieldnames(ref_model);
for k=1:length(noms)
   if ~isfield(ref,noms{k}) || isempty(ref.(noms{k}))
       ref.(noms{k}) = ref_model.(noms{k});
   end
end
ref.tokamak = strtrim(ref.tokamak);
ref.version = strtrim(ref.version);
ref.user    = strtrim(ref.user);

if isempty(ref.shot) || ~isfinite(ref.shot)
    error('shot number must be defined');
end
if isempty(ref.tokamak) 
    error('tokamak name must be defined');
end

function phys = cphys

% constante physique (phys)
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)



