function [error_flag,xsd] = metis4muscle3_make_idses(data_zerod,profil0d,z0dstruct,psi_offset,sawtooth,r0_norm)

% get IMAS dictionnary major version
IMAS_VERSION = get_imas_version;

% effect off ST crash
if ~isempty(sawtooth.next.delta_eta)
   delta_eta = inpterp1(sawtooth.next.x_ext,sawtooth.next.delta_eta,profil0d.xli,'linear','extrap');
   profil0d.eta = profil0d.eta + delta_eta;
   z0dstruct.profil.eta = z0dstruct.profil.eta + delta_eta;
   z0dstruct.profil0d.eta = z0dstruct.profil0d.eta + delta_eta;
end
if ~isempty(sawtooth.next.delta_jni)
   delta_jni = inpterp1(sawtooth.next.x_ext,sawtooth.next.delta_jni,profil0d.xli,'linear','extrap');
   profil0d.jni = profil0d.jni + delta_jni;
   z0dstruct.profil.jni = z0dstruct.profil.jni + delta_jni;
   z0dstruct.profil0d.jni = z0dstruct.profil0d.jni + delta_jni;    
end

% reset psi offset
data_zerod.edgeflux = data_zerod.edgeflux + 2 .* pi .* psi_offset; % in Vs !
z0dstruct.zerod.edgeflux = z0dstruct.zerod.edgeflux + 2 .* pi .* psi_offset;
z0dstruct.zs.edgeflux = z0dstruct.zs.edgeflux + 2 .* pi .* psi_offset;
profil0d.psi = profil0d.psi + psi_offset;
z0dstruct.profil.psi = z0dstruct.profil.psi + psi_offset;
z0dstruct.profil0d.psi = z0dstruct.profil0d.psi + psi_offset;


% start to manage COCOS in outputted IDSs
% backport for metis4imas
% default METIS in COCOS 7 with Ip and Bt positive
sigma_B0_in = 1;
sigma_Ip_in = 1;
factor_two_pi = 1;
sigma_bvac_r  = 1;
sigma_RphiZ_eff = 1;
% applied sign that is not completely defined in METIS
if isfield(z0dstruct.z0dinput.option,'orientation')
    sigma_B0_out = z0dstruct.z0dinput.option.orientation;
else
    sigma_B0_out = 1;
end
% disp('Signs are:')
% sigma_B0_out
sigma_Ip_out = sign(z0dstruct.z0dinput.option.signe .* sigma_B0_out);
% before COCOS transformation if any
sigma_B0_eff = sigma_B0_out;

COCOS_Sauter = false;
z0dstruct_loc = z0dstruct;
if ~isempty(which('ids_check_cocos'))
    switch  z0dstruct.z0dinput.option.COCOS_method
        case 'Sauter'
            COCOS_Sauter = true;
    end
end
if COCOS_Sauter
    % use Olivier Sauter tools for COCOS in this case
    %COCOS_Sauter = true;
        
    % change orientation and signe in z0dstruct if COCOS_Sauter to 1:
    % no modication of data before call of tool for COCOS transformation.
    z0dstruct_loc.z0dinput.option.signe = sigma_B0_in .* sigma_Ip_in;
    z0dstruct_loc.z0dinput.option.orientation = sigma_B0_in;
    
    % METIS data are in COCOS 7 with Ip and Bt positive !
    % neutral transformations if some field is not in COCOS 7 in METIS
    [data_zerod,profil0d,sigma_bvac_r,factor_two_pi,sigma_RphiZ_eff] = ...
      makecocos_imas_from_seven(data_zerod,profil0d,z0dstruct,7,sigma_Ip_in,sigma_B0_in);
  
    sigma_B0_eff = sigma_bvac_r;
    
             
elseif exist('profil0d','var') && exist('data_zerod','var')

  % METIS data are in COCOS 7 with Ip and Bt positive !
  [data_zerod,profil0d,sigma_bvac_r,factor_two_pi,sigma_RphiZ_eff] = ...
      makecocos_imas_from_seven(data_zerod,profil0d,z0dstruct,z0dstruct.z0dinput.option.COCOS,sigma_Ip_out,sigma_B0_out);
 
   sigma_B0_eff = sigma_bvac_r;
   
end 

% Structure for O. Sauter tools
cocos_struct.COCOS_Sauter = COCOS_Sauter;
% METIS native COCOS is 7:
% q > 0, Psi decreasing, dP/dPsi positive and phi is counted counter clock-wise.
% So (R,phi,Z) for computing BR  and BZ with sigma_Bp = -1 and  sigma_RZphi = 1.
% Theta from front is clockwise (sigma_rho_theta_phi = 1 and coordinate are
% (rho,theta,phi)
% By default with signe = 1 and orientation = 1, Bphi and the
% plasma current are in the same direction and Bphi is positive and counter
% clockwise the tokamak se from above.
% Psi is in Wb/rad so e_Bp = 0
cocos_struct.cocos_in     = 7; % METIS native
cocos_struct.cocos_out    = z0dstruct.z0dinput.option.COCOS;
cocos_struct.ipsign_out   = sign(z0dstruct.z0dinput.option.signe .* z0dstruct.z0dinput.option.orientation);
cocos_struct.b0sign_out   = z0dstruct.z0dinput.option.orientation;
cocos_struct.ipsign_in    = sigma_Ip_in; % METIS default
cocos_struct.b0sign_in    = sigma_B0_in; % METIS default
cocos_struct.error_bar    = 'none'; % no error bar in METIS output
cocos_struct.verbose      = z0dstruct.z0dinput.option.COCOS_verbose; % can be switch for debug
cocos_struct.check        = z0dstruct.z0dinput.option.COCOS_check; % can be switch for debug

disp('COCOS transform information:')
cocos_struct


%% CREATE MATLAB OUTPUT STRUCTURE
texte_diary = 'IDSes mapped from METIS4MUSCLE3';
error_flag = 0;
disp('Mapping IDS data structure pulse_schedule')
% this one is mandatory to get substructure code and prefill ids_properties
xsd.pulse_schedule = mappulse_schedule_imas(z0dstruct_loc,data_zerod);
disp('Mapping IDS data structure dataset_description')
if (z0dstruct.z0dinput.option.dataset_description == 1) || (z0dstruct.z0dinput.option.dataset_fair == 1)
    xsd.dataset_description = mapdataset_description_imas(z0dstruct_loc,data_zerod,texte_diary,error_flag,0,0,xsd.pulse_schedule.code);
    if z0dstruct.z0dinput.option.dataset_fair == 1
        xsd.dataset_fair        = map_dataset_fair(xsd.dataset_description,[]);
    end
end
if z0dstruct.z0dinput.option.summary == 1
    disp('Mapping IDS data structure summary')
    xsd.summary             = ids_gen('summary');
    xsd.summary.code        = xsd.pulse_schedule.code;
    xsd.summary = mapsummary_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.summary,0,0,sigma_B0_eff,sigma_bvac_r);
end
if z0dstruct.z0dinput.option.plasma_profiles == 1
    disp('Mapping IDS data structure plasma_profiles')
    %xsd.plasma_profiles       = ids_gen('plasma_profiles');
    xsd.plasma_profiles       =  merge_core_and_plasma('profiles');
    xsd.plasma_profiles.code  = xsd.pulse_schedule.code;
    xsd.plasma_profiles.ids_properties  =  xsd.pulse_schedule.ids_properties;
    xsd.plasma_profiles = mapcore_profiles_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.plasma_profiles,sigma_B0_eff);
    % alternative
    xsd.core_profiles = xsd.plasma_profiles;
end
if z0dstruct.z0dinput.option.plasma_transport == 1
    disp('Mapping IDS data structure plasma_transport')
    %xsd.plasma_transport      = ids_gen('plasma_transport');
    xsd.plasma_transport      = merge_core_and_plasma('transport');
    xsd.plasma_transport.code = xsd.pulse_schedule.code;
    xsd.plasma_transport.ids_properties  = xsd.pulse_schedule.ids_properties;
    xsd.plasma_transport = mapcore_transport_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.plasma_transport,sigma_B0_eff);
    % alternative
    xsd.core_transport = xsd.plasma_transport;
end
if z0dstruct.z0dinput.option.plasma_sources == 1
    %xsd.plasma_sources       = ids_gen('plasma_sources');
    xsd.plasma_sources       = merge_core_and_plasma('sources');
    xsd.plasma_sources.code  = xsd.pulse_schedule.code;
    xsd.plasma_sources.ids_properties  = xsd.pulse_schedule.ids_properties;
    model_source = xsd.plasma_sources.source{1};
    disp('Mapping IDS data structure plasma_sources (generic part)')
    xsd.plasma_sources = mapcore_sources_generic_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.plasma_sources,sigma_B0_eff);
    count_sources = 1;
    disp('Mapping IDS data structure plasma_sources for lhcd')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_lhcd_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    disp('Mapping IDS data structure plasma_sources for eccd')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_eccd_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    disp('Mapping IDS data structure plasma_sources for radiation')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_radiation_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    xsd.plasma_sources.source{count_sources} = mapcore_sources_brem_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    disp('Mapping IDS data structure plasma_sources for cyclotron')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_cyclotron_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    disp('Mapping IDS data structure plasma_sources for icrh')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_icrh_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    disp('Mapping IDS data structure plasma_sources for nbicd')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_nbicd_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    disp('Mapping IDS data structure plasma_sources for fusion')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_fusion_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    disp('Mapping IDS data structure plasma_sources ohm & boot')
    [xsd.plasma_sources.source{count_sources},xsd.plasma_sources.source{count_sources + 1}] = mapcore_sources_ohm_boot_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 2;
    disp('Mapping IDS data structure plasma_sources neutrals')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_neutral_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    disp('Mapping IDS data structure plasma_sources equipartition')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_equipartition_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    disp('Mapping IDS data structure plasma_sources for sum of sources')
    xsd.plasma_sources.source{count_sources} = mapcore_sources_full_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
    count_sources = count_sources + 1;
    % extra source for runaways
    if z0dstruct.z0dinput.option.runaway > 0
        disp('Mapping IDS data structure plasma_sources for runaway electron current')
        xsd.plasma_sources.source{count_sources} = mapcore_sources_runaways_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
    end
    % alternative
    xsd.core_sources = xsd.plasma_sources;
end
if z0dstruct.z0dinput.option.radiation == 1
    disp('Mapping IDS data structure radiation')
    xsd.radiation          = ids_gen('radiation');
    xsd.radiation.code     = xsd.pulse_schedule.code;
    xsd.radiation.ids_properties  = xsd.pulse_schedule.ids_properties;
    xsd.radiation = mapradiation_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.radiation,sigma_B0_eff);
end
switch z0dstruct.z0dinput.option.sol_model
    case '2_points'
        if z0dstruct.z0dinput.option.edge == 1
            disp('Mapping IDS data structure edge_profiles and edge_transport')
            [xsd.edge_profiles,xsd.edge_transport] = mapedge_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.pulse_schedule,sigma_B0_eff);
        end
        % the same data structure in METIS IDS is compatible with edge and
        % and plasma, just the serialization command parameter should
        % change.
        %xsd.plasma_profiles_1  = xsd.edge_profiles;
        %xsd.plasma_transport_1 = xsd.edge_transport;
end
if z0dstruct.z0dinput.option.numerics == 1
    disp('Mapping IDS data structure transport_solver_numerics')
    xsd.transport_solver_numerics                 = ids_gen('transport_solver_numerics');
    xsd.transport_solver_numerics.code            = xsd.pulse_schedule.code;
    xsd.transport_solver_numerics.ids_properties  =  xsd.pulse_schedule.ids_properties;
    xsd.transport_solver_numerics = maptransport_solver_numerics(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.transport_solver_numerics,0,0,sigma_B0_eff,sigma_bvac_r);
end
if z0dstruct.z0dinput.option.equilibrium == 1
    disp('Mapping IDS data structure equilibrium')
    xsd.equilibrium         = ids_gen('equilibrium');
    xsd.equilibrium.code    = xsd.pulse_schedule.code;
    xsd.equilibrium.ids_properties  = xsd.pulse_schedule.ids_properties;
    if isfield(xsd,'summary')
        xsd.equilibrium = mapequilibrium_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.equilibrium,xsd.summary, ...
            1,z0dstruct.z0dinput.option.COCOS,z0dstruct.z0dinput.option.equi_extrap,factor_two_pi,sigma_RphiZ_eff);
    else
        xsd.equilibrium = mapequilibrium_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.equilibrium,[], ...
            1,z0dstruct.z0dinput.option.COCOS,z0dstruct.z0dinput.option.equi_extrap,factor_two_pi,sigma_RphiZ_eff);
    end
end

% patch magnetic rigidity to use only one fixed value of r0 for the whole
% life of the actor
noms = fieldnames(xsd);
for k=1:length(noms)
    if isfield(xsd.(noms{k}),'vacuum_toroidal_field')
            xsd.(noms{k}).vacuum_toroidal_field.b0 =  xsd.(noms{k}).vacuum_toroidal_field.b0 .* ...
                 xsd.(noms{k}).vacuum_toroidal_field.r0 ./ r0_norm;
                 xsd.(noms{k}).vacuum_toroidal_field.r0 =  r0_norm;
    end       
end

% change all COCOS is xsd
noms = fieldnames(xsd);
for k=1:length(noms)
    xsd.(noms{k}) = change_cocos_OS(xsd.(noms{k}),noms{k},cocos_struct);
end






