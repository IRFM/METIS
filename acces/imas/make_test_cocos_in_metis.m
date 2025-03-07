% automatic test of COCOS coherence in METIS for all configurations
function make_test_cocos_in_metis(dataversion,moderefplus)

% clear workspace for test with working point
evalin('base','clear');

% get METIS root path
root = fileparts(which('metis'));

% report
if exist('report_test_cocos_in_metis.txt')
    delete('report_test_cocos_in_metis.txt');
end
diary off

% test availability of tools
if isempty(which('ids_check_cocos'))
    init_cocos_test_in_metis([],true);
end
% initialisation IMAS variables
if nargin < 1
    dataversion = '';
elseif isnumeric(dataversion)
    dataversion = sprintf('%d',dataversion);
end
%
[ual_tokamak,ual_user,ual_dataversion,ual_backend] = imasdb('test_cocos_imas','',dataversion,13);

if (nargin < 2) || isempty(moderefplus)
    moderefplus = false;
end
    

% remove existing reccord (to be sure to update non timed data)
home = fileparts(getenv('HOME'));
pathtodata = fullfile(home,ual_user,'public','imasdb',ual_tokamak,ual_dataversion(1),'99999','99');
if isdir(pathtodata) && ~isempty(dir(fullfile(pathtodata,'*.h5')))
   [s,t] = unix(sprintf('rm %s',fullfile(pathtodata,'*.h5')));
   if s~= 0
       error(t);
   end
end

% reference COCOS of METIS
% METIS native COCOS is 7:
% q > 0, Psi decreasing, dP/dPsi positive and phi is counted counter clock-wise.
% So (R,phi,Z) for computing BR  and BZ with sigma_Bp = -1 and  sigma_RZphi = 1.
% Theta from front is clockwise (sigma_rho_theta_phi = 1 and coordinate are
% (rho,theta,phi)
% By default with signe = 1 and orientation = 1, Bphi and the
% plasma current are in the same direction and Bphi is positive and counter
% clockwise the tokamak se from above.
% Psi is in Wb/rad so e_Bp = 0
COCOS_ref = 7;

% we need just one time slice
% random orientation and sign at initialisation
info = metis4imas(1);
options = info.valeur;
options.signe = 1;
options.orientation = 1;
options.COCOS = COCOS_ref;
options.COCOS_method = 'Sauter';


% make test imas 
% create IMAS data
[error_flag,output_data] = metis4imas(99999,99,'','test',[],options);

% initialisation
[error_flag,output_data] = metis4imas(99999,99,'','init',19,options);
options = pulse_schedule2option_imas(output_data.pulse_schedule);

% loop of configuration signs for Bt and Ip
for k=1:2
    switch k
        case 1
            options.signe = 1;
        otherwise
            options.signe = -1;
    end
    for l=1:2
        switch l
            case 1
                options.orientation = 1;
            otherwise
                options.orientation = -1;
        end
        for m=1:3
            switch m
                case 1
                    options.COCOS = 11;
                case 2
                    options.COCOS = 17;
                otherwise
                    options.COCOS = COCOS_ref;                    
            end
           
            % skip no reference mode
            if moderefplus && (options.COCOS ~= COCOS_ref) 
                continue;
            end
            
            % switch option in the loop
            %z0dstruct = getappdata(0,'IMAS_Z0DSTRUCT');
            %z0dstruct.z0dinput.option = options;
            %setappdata(0,'IMAS_Z0DSTRUCT',z0dstruct);
            % need to read and write pulse_schdule !
            
            % call METIS one time
            try
                [error_flag,output_data] = metis4imas(99999,99,'','one_time',98+k,options);
            catch
                keyboard
            end
            
            % get orientation from METIS parameters
            % Does not work as IMAS don't update non timed data
            % options_real= pulse_schedule2option_imas(output_data.pulse_schedule);
            options_real = options;
            
            % maybe there is a sign ?
            sign_bt = options_real.orientation;
            % sign of Bphi projected on plasma current direction
            sign_ip = options_real.signe .* sign_bt;
            

            
            % test results
            % loop on IDS
            idslist = fieldnames(output_data);
            diary on
            diary('report_test_cocos_in_metis.txt');
            disp(' ')
            disp('****************************************************************************');
            fprintf('COCOS = %d, case sign_bt = %d &  sign_ip= %d\n\n',options_real.COCOS,sign_bt,sign_ip)
            for k=1:length(idslist)
                ids_name = idslist{k};
                disp(' ')
                disp('----------------------------------------------------------------------------');
                fprintf('Testing IDS %s:\n',ids_name)
                switch ids_name
                    case 'core_profiles'
                        % remove new field not defined in ids_gen due to differences in model version
                        model = ids_gen(ids_name);
                        model_fields = fieldnames(model);
                        data  = output_data.(ids_name);
                        data_fields  = fieldnames(data);
                        for lk =1:length(data_fields)
                            if isempty(strmatch(data_fields{lk},model_fields,'exact'))
                                fprintf('Non existing field %s in IMAS model\n',data_fields{lk});
                                data = rmfield(data,data_fields{lk});
                            end
                        end
                        %
                        ids_generic_cocos_check(data,ids_name,options_real.COCOS,sign_ip,sign_bt);
                    case 'equilibrium'
                        % remove new field not defined in ids_gen due to differences in model version
                        model = ids_gen(ids_name);
                        model_fields = fieldnames(model);
                        data  = output_data.(ids_name);
                        data_fields  = fieldnames(data);
                        for lk =1:length(data_fields)
                            if isempty(strmatch(data_fields{lk},model_fields,'exact'))
                                fprintf('Non existing field %s in IMAS model\n',data_fields{lk});
                                data = rmfield(data,data_fields{lk});
                            end
                        end
                        %
                       ids_generic_cocos_check(data,ids_name,options_real.COCOS,sign_ip,sign_bt);
                    otherwise
                        % for the moment
                        disp('test not yet available');
                        continue;
                end
            end
            diary off
        end
    end
end

% End of part one.
% now use real METIS run
metisfile = fullfile(root,'certification','metis','ITER_rampup_ECCD.mat');
metis_load(metisfile)
time = evalin('base','z0dinput.cons.temps');
setappdata(0,'WORKING_POINT_TIME',time(end-1));
% complete parameters set
options = evalin('base','z0dinput.option');
info = metis4imas(1);
option = info.valeur;
noms = fieldnames(option);
for k=1:length(noms)
    if ~isfield(options,noms{k})
        options.(noms{k}) = option.(noms{k});
    end
end
options = secure_option_imas(options);
options.COCOS_method = 'Sauter';

% loop on signe, orientation  and COCOS
for k=1:2
    switch k
        case 1
            options.signe = 1;
        otherwise
            options.signe = -1;
    end
    for l=1:2
        switch l
            case 1
                options.orientation = 1;
            otherwise
                options.orientation = -1;
        end
        for m=1:3
            switch m
                case 1
                    options.COCOS = 11;
                case 2
                    options.COCOS = 17;
                otherwise
                    options.COCOS = COCOS_ref;                    
            end
            
            % skip no reference mode
            if moderefplus && (options.COCOS ~= COCOS_ref)
                continue;
            end
            
           % change option in the workspace
            zassignin('base','z0dinput.option',options);
            % run METIS working point
            evalin('base','z0working_point;');
            
            % test IDS COCOS
            post = evalin('base','post');
            
            diary('report_test_cocos_in_metis.txt');
            disp(' ')
            disp('****************************************************************************');
            fprintf('COCOS = %d, case sign_bt = %d &  sign_ip= %d\n\n',options.COCOS,options.signe,options.orientation)

            test_cocos_IMAS_METIS(post,ual_dataversion);
            diary off
            diary('void.txt');

        end
    end
end


% 
edit('report_test_cocos_in_metis.txt');
