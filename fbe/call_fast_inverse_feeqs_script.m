% call fast inverse FEEQS mode script for ITER or other device; to be updated
zineb_path;
root_feeqs = '/Applications/software/FEEQS.M';
set_reactor_feeqs_path(root_feeqs);
%
info = compute_reactor_fbe_inverse_fast;
option_feeqs_reactor = info.valeur;

% metis root
root_metis    = fileparts(which('metis'));   
%
option_feeqs_reactor.coils_description = fullfile(root_feeqs,'Projects/Reactors/iter_coils_fast_mode_alt.txt'); % some example for ITER
option_feeqs_reactor.metisfilename = fullfile(root_metis,'certification/metis/iter_2nbi_pfus.mat'); % METIS filename with path
option_feeqs_reactor.first_wall_name = fullfile(root_feeqs,'Projects/Reactors/iter_first_wall.txt'); % some example for ITER
option_feeqs_reactor.number_of_time_slices  = 11; % numeber of time slice computed, put a large number to compute all tim slices
% get som information from METIS file
data_metis = load_metis_imas(option_feeqs_reactor.metisfilename);
disp('------------------------------------------------------------')
if isfinite(data_metis.z0dinput.option.available_flux)
    fprintf('Declared available poloidal flux in METIS is %g Wb\n', 2 * pi * data_metis.z0dinput.option_feeqs_reactor.available_flux);
else
    disp('No available poloidal flux delcared in METIS:');
    disp('It will be preferable to set it in the parameter: z0dinput.option_feeqs_reactor.available_flux before running the simulation');
end
disp('-------------------------------------------------------------');
%
option_feeqs_reactor.first_time =  min(data_metis.z0dinput.cons.temps);
option_feeqs_reactor.last_time  =  max(data_metis.z0dinput.cons.temps); 
option_feeqs_reactor.reinforce  =  30; % weight to match extremal point in LCFS
option_feeqs_reactor.run_name   =  'test_1'; % to retreive the run when there are many
option_feeqs_reactor.plotonoff  =  0; % set graph display off/on (0 = no graph, 1 = some graphs and 2 = all graphs')
option_feeqs_reactor.movie_onoff = 0; % if = 1 generate a movie for presentation
% the best is to provide available poloidal flux !
%option_feeqs_reactor.flux_constant_mode = 'manual'; % if manual read z0dinput.option_feeqs_reactor.available_flux to compute initial magnetisation (i.e.  initial currents in central CS)
%option_feeqs_reactor.flux_constant_mode = 'auto'; % try to compute the intial magnetisation from simulation flux consumption; work only for complete scenario with ramp-up, flattop and ramp-down
option_feeqs_reactor.flux_constant_mode = 182 / 2/pi / 2; % in Wb/radian - set the value in Wb for premagnetisation
option_feeqs_reactor.flux_offset        = 0;   % to match some external computation a first time (for example from CREATE_NL)
option_feeqs_reactor.breakdown_flux     = 0.3000;  % flux consumption for the breakdown 
option_feeqs_reactor.weight             =  1; % regularisation term weight in cost function (sum(I^2) and dipole)
option_feeqs_reactor.mode_profile       = 'on'; % if = on, use P' and FF', otherwise use beta-p and li_3
option_feeqs_reactor.V_Force            =  'on'; % Compute (on) or not (off) volatages and forces esimation in postprocessing of FEEQS computation (need boboz to be compiled)

%
disp(' ');
% call FEESQ.M
dirmem = pwd;
fprintf('Initial directory is %s\n',dirmem);
switch option_feeqs_reactor.flux_constant_mode
    case 'manual'
        offset_flux = option_feeqs_reactor.flux_offset / 2 / pi;
    case 'auto'
        offset_flux = NaN;
    otherwise
        offset_flux = option_feeqs_reactor.flux_constant_mode;
end
% make a local copy of the project
projectpath = fileparts(which('compute_inverse4metis_fast'));
[~,projectname]  = fileparts(projectpath);
% temprary file
if isdir('/dev/shm')
    templocaldir = fullfile('/','dev','shm',sprintf('%s_%s',fileparts(tempname),getenv('USER')));
else
    templocaldir = tempname;
end

try
    output =  reactor_inverse4metis_fast(option_feeqs_reactor.coils_description,option_feeqs_reactor.metisfilename, ...
        option_feeqs_reactor.number_of_time_slices,option_feeqs_reactor.first_time,option_feeqs_reactor.last_time, ...
        option_feeqs_reactor.reinforce,option_feeqs_reactor.plotonoff, option_feeqs_reactor.movie_onoff,option_feeqs_reactor.run_name, ...
        option_feeqs_reactor.weight,offset_flux,dirmem,option_feeqs_reactor.breakdown_flux, ...
        option_feeqs_reactor.mode_profile,option_feeqs_reactor.V_Force,option_feeqs_reactor.first_wall_name);
    
    cd(dirmem);
catch
    keyboard
    cd(dirmem);
    error(lasterror);
end
option_mem = option_feeqs_reactor;
option_feeqs_reactor = output;
option_feeqs_reactor.flux_offset = 2 .* pi .* output.psi_off;
noms = fieldnames(option_mem);
for k=1:length(noms)
    if ~isfield(option_feeqs_reactor,noms{k})
        option_feeqs_reactor.(noms{k}) = option_mem.(noms{k});
    end
end

% final graphs
fprintf('Results has been stored in %s\n',output.output_name);
plotfigure_reactor_inverse4metis(output.output_name)

