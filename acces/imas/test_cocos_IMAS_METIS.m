% test COCOS coherence for one METIS dataset 
function test_cocos_IMAS_METIS(post,dataversion,targeted_cocos)

% test availability of tools
if isempty(which('ids_check_cocos'))
    init_cocos_test_in_metis([],true);
end
if nargin < 2
    dataversion = '';
elseif isnumeric(dataversion)
    dataversion = sprintf('%d',dataversion);
end
%
[ual_tokamak,ual_user,ual_dataversion,ual_backend] = imasdb('test_cocos_imas','',dataversion,13);
fprintf('Selected major data version is %s\n',ual_dataversion(1));
% 
% use test reccord
% occurrence.tokamak = ual_tokamak;
% occurrence.user    = ual_user;
% occurrence.dataversion = ual_dataversion;
% occurrence.occurrence = 0;
% occurrence.backend = ual_backend;
%
% add IMAS control parameters
info = metis4imas;
names = fieldnames(info.valeur);
for k=1:length(names)
    if ~isfield(post.z0dinput.option,names{k})
        post.z0dinput.option.(names{k}) = info.valeur.(names{k});
    end
end
if (nargin >= 3) && ~isempty(targeted_cocos)
    post.z0dinput.option.COCOS = targeted_cocos;
end
fprintf('Targeted COCOS is %d\n',post.z0dinput.option.COCOS);

% error_flag = metis4imas(99999,99,occurrence,post);
%
[error_flag,output_data] = metis4imas(99999,99,0,post,[]);
%
% get code param
options = pulse_schedule2option_imas(output_data.pulse_schedule);
fprintf('Obtained COCOS is %d\n',options.COCOS);

% get orientation from METIS parameters
% maybe there is a sign ?
sign_bt = options.orientation;
% sign of Bphi projected on plasma current direction
sign_ip = options.signe .* sign_bt;

% get list of created IDSs
%rep = dir(fullfile(getenv('HOME'),'public','imasdb',ual_tokamak,ual_dataversion(1),'99999','99','*.h5'));

% loop on IDS
diary on
idslist = fieldnames(output_data);
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
            for l =1:length(data_fields)
                if isempty(strmatch(data_fields{l},model_fields,'exact'))
                    fprintf('Non existing field %s in IMAS model\n',data_fields{l});
                    data = rmfield(data,data_fields{l});
                end
            end
            %
            ids_generic_cocos_check(data,ids_name,options.COCOS,sign_ip,sign_bt);
        case 'equilibrium'
            % remove new field not defined in ids_gen due to differences in model version
            model = ids_gen(ids_name);
            model_fields = fieldnames(model);
            data  = output_data.(ids_name);
            data_fields  = fieldnames(data);
            for l =1:length(data_fields)
                if isempty(strmatch(data_fields{l},model_fields,'exact'))
                    fprintf('Non existing field %s in IMAS model\n',data_fields{l});
                    data = rmfield(data,data_fields{l});
                end
            end
            %
            ids_generic_cocos_check(data,ids_name,options.COCOS,sign_ip,sign_bt);
        otherwise
            % for the moment
            disp('test not yet available');
            continue;
            
            fprintf('trying for a core_profiles like IDS\n')
            % remove/add new field not defined in ids_gen due to differences in model version
            model = ids_gen('core_profiles');
            model_fields = fieldnames(model);
            data  = output_data.(ids_name);
            data_fields  = fieldnames(data);
            for l =1:length(data_fields)
                if isempty(strmatch(data_fields{l},model_fields,'exact'))
                    fprintf('Non existing field %s in IMAS model\n',data_fields{l});
                    data = rmfield(data,data_fields{l});
                end
            end
            data_fields  = fieldnames(data);
            for l =1:length(model_fields)
                if isempty(strmatch(model_fields{l},data_fields,'exact'))
                    fprintf('Non existing field %s in IMAS data\n',model_fields{l});
                    data.(model_fields{l}) = model.(model_fields{l});
                end
            end
            %
            try
                ids_generic_cocos_check(data,ids_name,options.COCOS,sign_ip,sign_bt);
            catch
                disp('does not working with core_profiles like IDS');
            end
            fprintf('trying for a equilibrium like IDS\n')
            % remove/add new field not defined in ids_gen due to differences in model version
            model = ids_gen('core_profiles');
            model_fields = fieldnames(model);
            data  = output_data.(ids_name);
            data_fields  = fieldnames(data);
            for l =1:length(data_fields)
                if isempty(strmatch(data_fields{l},model_fields,'exact'))
                    fprintf('Non existing field %s in IMAS model\n',data_fields{l});
                    data = rmfield(data,data_fields{l});
                end
            end
            data_fields  = fieldnames(data);
            for l =1:length(model_fields)
                if isempty(strmatch(model_fields{l},data_fields,'exact'))
                    fprintf('Non existing field %s in IMAS data\n',model_fields{l});
                    data.(model_fields{l}) = model.(model_fields{l});
                end
            end
            try
                ids_generic_cocos_check(data,ids_name,options.COCOS,sign_ip,sign_bt);
            catch
                disp('does not working with equilibrium like IDS');
            end
    end
end

% that all !

