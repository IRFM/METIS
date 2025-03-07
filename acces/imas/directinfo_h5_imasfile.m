% helper function to get information directly to h5 file for IMAS data
% syntax:
%    data = directinfo_h5_imasfile(path2file,keep_shape)
%
% inputs:
%   path2file       = full path to directory containing .h5 files or IMAS
%                     structure reference
%   keep_shape      = if true return also _AOS_SHAPE and _SHAPE in data name
%
% output:
%   data = information on IDSs (and master file)
%          Each field of data is a IDS (at the exeption of master).
%          Each data.(ids_name) substructure contains a field Size.
%          data.(ids_name).Size.(IMAS_data_name) is the size of the data.
%
% remark: 
%  1/ The occurrence is encoded in IDS name. Example:
%     equilibrium occurrence 0 = 'equilibrium'
%     equilibrium occurrence 1 = 'equilibrium_1'
%
%  2/ if a field do not exist, the corresponding data has not been produced.
%
%  3/ if the size is a matrix, the first column is the size of elements 
%     and the second is the  chunk size.
%
%
function data = directinfo_h5_imasfile(path2file,keep_shape)

% inputs
if nargin < 2
   keep_shape = false;
elseif isempty(keep_shape)
   keep_shape = false;
end   

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

% first: files in the directory
data = [];
lf = dir(fullfile(path2file,'*.h5'));
if isempty(lf)
    return;
end
for k= 1:length(lf)
    try
        info = h5info(fullfile(path2file,lf(k).name));
    catch
        info = [];
    end
    data.(strtok(lf(k).name,'.')) = info;
end

% next: if master.h5, follows links
noms = fieldnames(data);
for k=1:length(noms)
   if isfield(data.(noms{k}),'Links') && ~isempty(data.(noms{k}).Links)
       link = data.(noms{k}).Links;
       for l=1:length(link)
           try
               info = h5info(fullfile(path2file,link(l).Value{1}));
           catch
               info = [];
           end
           data.(link(l).Name) = info;
       end
   end
end

% get sizes
noms = fieldnames(data);
for k=1:length(noms)
   data.(noms{k}).Size = [];
   if isfield(data.(noms{k}),'Groups')  && ~isempty(data.(noms{k}).Groups)
        data.(noms{k}).Size = get_size_info(data.(noms{k}).Groups,keep_shape);
   end
end
    
    
function info_size = get_size_info(groups_in,keep_shape)


info_size = [];
if isempty(groups_in)
    return;
else
    for k=1:length(groups_in)
        loc_name = strrep(groups_in(k).Name,'/','');
        if isfield(groups_in(k),'Groups') && ~isempty(groups_in(k).Groups)
            if length(groups_in) == 1
                info_size.Groups = get_size_info(groups_in(k).Groups,keep_shape);
            else
                info_size.(loc_name).Groups = get_size_info(groups_in(k).Groups,keep_shape);
            end
        end
        if isfield(groups_in(k),'Datasets') && ~isempty(groups_in(k).Datasets)
            for l=1:length(groups_in(k).Datasets)
                data_name = strrep(groups_in(k).Datasets(l).Name,'&','__');
                data_name = strrep(data_name,'[','');
                data_name = strrep(data_name,']','');
                if ~isshape(data_name) || keep_shape
                    if isempty(groups_in(k).Datasets(l).Dataspace.Size)
                        if length(groups_in) == 1
                            info_size.(data_name) = numel(groups_in(k).Datasets(l).Dataspace);
                        else
                            info_size.(loc_name).(data_name) = numel(groups_in(k).Datasets(l).Dataspace);
                        end
                    else
                        if length(groups_in) == 1
                            info_size.(data_name) = groups_in(k).Datasets(l).Dataspace.Size;
                        else
                            info_size.(loc_name).(data_name) = groups_in(k).Datasets(l).Dataspace.Size;
                        end
                    end
                    if ~isempty(groups_in(k).Datasets(l).ChunkSize)
                        if length(groups_in) == 1
                            if (length(info_size.(data_name)) ~= length(groups_in(k).Datasets(l).ChunkSize)) || ...
                                    ~all(info_size.(data_name)(:) == groups_in(k).Datasets(l).ChunkSize(:))
                                info_size.(data_name) =  cat(2,info_size.(data_name)(:),groups_in(k).Datasets(l).ChunkSize(:));
                            end
                        else
                            if (length(info_size.(data_name)) ~= length(groups_in(k).Datasets(l).ChunkSize)) || ...
                                    ~all(info_size.(data_name)(:) == groups_in(k).Datasets(l).ChunkSize(:))
                                info_size.(loc_name).(data_name) = cat(2,info_size.(loc_name).(data_name)(:), groups_in(k).Datasets(l).ChunkSize(:));
                            end
                        end
                    end
                end
            end
        end
    end
end
    
function flag = isshape(name_in)

flag = false;
ind = strfind(name_in,'_SHAPE');
if ~isempty(ind)
    ind = ind(end);
    if length(name_in) == (ind + length('_SHAPE') - 1)
        flag = true;
    end
end


