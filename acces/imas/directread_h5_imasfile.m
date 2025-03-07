% helper function to access directly to h5 file for IMAS data
% syntax:
%    data = directread_h5_imasfile(path2file,ids_name,imas_data_name,[occurrence])
%
% inputs:
%   path2file       = full path to directory containing .h5 files or IMAS
%                     structure reference
%   ids_name        = name of the ids ('equilibrium')
%   imas_data_name  = name of the data as appears in Matlab without taking
%                     into account the cell or array of struct syntax
%                     ('time_slice.profiles_1d.dpressure_dpsi')
%  occurrence       = optional occurrence number
%
% output:
%   data = data retrieve from h5 file.
%          time is on the first coordinate.
%
function data = directread_h5_imasfile(path2file,ids_name,imas_data_name,occurrence)

% if IMAS reference
if isstruct(path2file)
    ref = path2file;
    if isfield(ref,'user')
        switch lower(ref.user)
            case 'imas_public'
                path2file = '/Imas_public';
            case 'imas_initiation'
                path2file = '/Imas_public/imas_initiation';
            case 'imas_public_continuous'
                path2file = '/Imas_public/imas_public_continuous';
            case 'imas_simulation'
                path2file = '/Imas_public/imas_simulation';
            case 'imas_static'
                path2file = '/Imas_public/imas_static';
            otherwise
                home_root = fileparts(getenv('HOME'));
                path2file = fullfile(home_root,strtrim(ref.user));
        end
    else
        path2file = getenv('HOME');
    end
    path2file = fullfile(path2file,'public','imasdb');
    if isfield(ref,'machine')
        path2file = fullfile(path2file,strtrim(ref.machine));
    else
        path2file = fullfile(path2file,strtrim(ref.tokamak));
    end
    if isfield(ref,'version')
        if isnumeric(ref.version)
            version = sprintf('%d',ref.version);
        else
            version = ref.version;
        end
    else
        version = '3';
    end
    path2file = fullfile(path2file,version(1));
    path2file = fullfile(path2file,sprintf('%d',ref.shot));
    if isfield(ref,'run')
        if isnumeric(ref.run)
            run = sprintf('%d',ref.run);
        else
            run = ref.run;
        end
    else
        run = '0';
    end
    path2file = fullfile(path2file,run);        
end

%
% generate H5 path
if (nargin > 3)  && ~isempty(occurrence) && (occurrence ~= 0)
    ids_name_occ = sprintf('%s_%d',strtrim(ids_name),occurrence);
    h5path = sprintf('/%s_%d/',strtrim(ids_name),occurrence);
    link_str = sprintf('%s_%d',strtrim(ids_name),occurrence);
else
    h5path   = sprintf('/%s/',strtrim(ids_name));
    ids_name_occ = ids_name;
    link_str = sprintf('%s',strtrim(ids_name));
end
sepa = '';
if any(imas_data_name == '/')
    sepa = '/';
elseif any(imas_data_name == '.')
    sepa = '.';
end
if isempty(sepa)
    %
    ids_model = ids_gen(ids_name);
    name1 = imas_data_name;
    if iscell(ids_model.(name1))
        h5path = sprintf('%s%s[]&',h5path,name1);
    else
        h5path = sprintf('%s%s&',h5path,name1);
    end
else
    %
    ids_model = ids_gen(ids_name);
    
    k= 999;
    while ~isempty(imas_data_name ) && (k > 0)
        k = k -1;
        [name1,imas_data_name] = strtok(imas_data_name,sepa);
        if iscell(ids_model.(name1))
            h5path = sprintf('%s%s[]&',h5path,name1);
            ids_model = ids_model.(name1){1};
        else
            h5path = sprintf('%s%s&',h5path,name1);
            ids_model = ids_model.(name1);
        end
        imas_data_name = imas_data_name(2:end);
    end
end
try
    data = h5read(fullfile(path2file,sprintf('%s.h5',ids_name_occ)),h5path(1:end-1));
catch
    % find is there is some link to follow
    info = h5info(fullfile(path2file,sprintf('%s.h5','master')));
    if isfield(info,'Links')
        link = info.Links;
        for k=1:length(link)
            if strcmp(strtrim(link(k).Name),link_str)
                try
                    new_run = split(link(k).Value{1},'/');
                catch
                    new_run = {};
                    tosplit = link(k).Value{1};
                    while ~isempty(tosplit)
                        [new_run{end+1},tosplit] = strtok(tosplit,'/');
                        tosplit = tosplit(2:end);
                    end
                end
                path2file = fullfile(fileparts(path2file),strtrim(new_run{2}));
                data = h5read(fullfile(path2file,sprintf('%s.h5',ids_name_occ)),h5path(1:end-1));
                break;
            end
        end
    end
end
if all(size(data) > 1) && isnumeric(data) && (length(size(data)) == 2)
    ss = size(data);
    ss = ss(:)';
    data = reshape(data.',cat(2,ss(end),ss(1:end-1)));
elseif all(size(data) > 1) && isnumeric(data)
    data = shiftdim(data,1);  
end
