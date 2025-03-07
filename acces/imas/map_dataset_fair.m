function data = map_dataset_fair(dataset_description,data_fair_in,shot,run)



% model of IDS
try
    data = ids_gen('dataset_fair');
catch
    % not available in present version of IMAS
    data = [];
    return
end

data.ids_properties.comment          = dataset_description.ids_properties.comment;
data.ids_properties.homogeneous_time = 1;
data.ids_properties.source           = dataset_description.ids_properties.source;
data.ids_properties.provider         = dataset_description.ids_properties.provider;
data.ids_properties.creation_date    = dataset_description.ids_properties.creation_date;
data.time                            = dataset_description.time;
data.identifier                      = '';
if ~isempty(data_fair_in) && isfield(data_fair_in,'identifier');
    data.replaces                        = data_fair_in.identifier;
else
    data.replaces                        = '';
end
data.is_replaced_by                  = '';
data.valid                           = sprintf('%s/%s',datestr(now,29),datestr(now+99*365.25,29));   %YYYY-MM-DD/YYYY-MM-DD
data.rights_holder                   = '';
data.license                         = '';
data.is_referenced_by                = ''; 
data.ids_properties.version_put.data_dictionary       = dataset_description.dd_version;
data.ids_properties.version_put.access_layer          = dataset_description.dd_version;
data.ids_properties.version_put.access_layer_language = 'Matlab';
                
% try to get better definition of IMAS UAL version 
try
    % extraction of data
    rep = which('imas_open_env');
    [rep,l] = fileparts(rep);
    [rep,l] = fileparts(rep);
    [rep,info] = fileparts(rep);
    [dd_version,access_layer] = strtok(info,'-');
    data.ids_properties.version_put.data_dictionary       = dd_version;
    data.ids_properties.version_put.access_layer          = access_layer(2:end);
catch
    try
        data.ids_properties.version_put.data_dictionary  = load(fullfile(fileparts(which('litidss')),'noimas_installed','DATA4IMAS_WITHOUT_INFRASTRUCTURE'),'data_version');
        data.ids_properties.version_put.access_layer = 'noimas_installed';
    catch
         disp('unable to update IMAS UAL information');
    end
end

% some information identifier ?
try
    hn = z0hostname;
    if isappdata(0,'UAL_USER') && ~isempty(getappdata(0,'UAL_USER'))
        user = getappdata(0,'UAL_USER');
    else
        user = dataset_description.data_entry.user;
    end
    if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK'))
        machine = getappdata(0,'UAL_TOKAMAK');
    else
        machine = dataset_description.data_entry.machine;
    end
    if nargin < 3
        shot = dataset_description.data_entry.pulse;
    elseif isempty(shot)
        shot = dataset_description.data_entry.pulse;
    end   
    if nargin < 4
        run = dataset_description.data_entry.run;
    elseif isempty(run)
        run = dataset_description.data_entry.run;
    end   
    data.identifier  = sprintf('file://%s/~%s/public/imasdb/%s/%s/%d/ids_%d%4.4d.*',hn, ...
        user,machine, ...
        data.ids_properties.version_put.data_dictionary(1), ...
        fix(run/10000),shot,run);
    if ~isempty(data.replaces)
        data.is_replaced_by = data.identifier;
    end
    % update 
    while(sum(hn=='.') > 1)
        [rh,hn] = strtok(hn,'.'); 
    end
    data.rights_holder = upper(rh);
catch
    disp('unable to update identifier information');
end


%disp(data)

