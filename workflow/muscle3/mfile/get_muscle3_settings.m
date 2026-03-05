function settings_muscle3 = get_muscle3_settings(instance,list_of_muscle3_parameters,prefix_metis_option,prefix_external_data,list_exd)

if isempty(list_of_muscle3_parameters)
    settings_muscle3 = [];
else
    
    for k=1:length(list_of_muscle3_parameters)
        name    = list_of_muscle3_parameters{k}{1};
        type    = list_of_muscle3_parameters{k}{2};
        default = list_of_muscle3_parameters{k}{3};
        try
            data = instance.get_setting(name, type);
        catch
            data = [];
        end
        if isempty(data)
            fprintf('Warning: %s setting is not defined: set to default value\n');
            eval(sprintf('settings_muscle3.%s=default;',strrep(name,'_','.')));
        else
            %dbstop if all error
            switch type
                case 'str'
                    var  = char(data);
                otherwise
                    var  = double(data);
            end
            eval(sprintf('settings_muscle3.%s=var;',strrep(name,'_','.')));
            %dbclear all
        end
    end
end

% read optional parameters
if (nargin > 2) && ~isempty(prefix_metis_option)
    model = metis4imas(1);
    % get list of all settings
    listofsettings = cell(instance.list_settings);
    if ~isempty(listofsettings)
        for k=1:length(listofsettings)
            name = char(listofsettings{k});
            if startsWith(name,prefix_metis_option)
               metis_name = strrep(name,sprintf('%s_',prefix_metis_option),'');
               % for convinience and to be user friendly, metis parameter
               % can be design by is label
               if ~isfield(model.valeur,metis_name)
                   % search for label
                   nomsl = fieldnames(model.label);
                   for l=1:length(nomsl)
                      switch  model.label.(nomsl{l})
                          case metis_name
                              metis_name = nomsl{l};
                              break
                      end
                   end
               end
               if isfield(model.type,metis_name)
                   switch model.type.(metis_name)
                       case {'int', 'integer'}
                           data = double(instance.get_setting(name, 'int'));
                       case {'real','float'}
                           data = double(instance.get_setting(name, 'float'));                           
                       case {'string','char','str'}
                           data = char(instance.get_setting(name, 'str'));                          
                       otherwise
                           error('undefined type');
                   end
                   if ~isempty(data)
                        settings_muscle3.(prefix_metis_option).(strrep(name,sprintf('%s_',prefix_metis_option),'')) = data;
                   else
                        fprintf('Unvalid METIS parameter %s defined in settings\n', metis_name);   
                   end
               else
                   fprintf('Unvalid METIS parameter name %s defined in settings\n', metis_name);
               end
            end
        end
    end
end

% read external data settings
if (nargin > 4) && ~isempty(prefix_external_data)
    % get list of all settings
    listofsettings = cell(instance.list_settings);
    if ~isempty(listofsettings)
        for l=1:length(list_exd)
            ext_name = strrep(list_exd{l}{1},sprintf('%s_',prefix_external_data),'');
            settings_muscle3.(prefix_external_data).(ext_name).state    = list_exd{l}{4};
            settings_muscle3.(prefix_external_data).(ext_name).ids_name = list_exd{l}{2};
            settings_muscle3.(prefix_external_data).(ext_name).ext_name = list_exd{l}{3};
            for k=1:length(listofsettings)
                name = char(listofsettings{k});
                if startsWith(name,prefix_external_data)
                    if strcmp(name,list_exd{l}{1})
                        settings_muscle3.(prefix_external_data).(ext_name).state = double(instance.get_setting(name, 'int'));
                    end
                end
            end
        end
    end
end

% for later use in graphics display
instance_listofsettings = cell(instance.list_settings);
for k=1:length(instance_listofsettings)
    settings_muscle3.instance_listofsettings{k} = char(instance_listofsettings{k});
end


