function init_cocos_test_in_metis(local_imas_version,no_update)

% input
if (nargin < 1) || isempty(local_imas_version)
    local_imas_version = false;
end
if (nargin < 2) || isempty(no_update)
    no_update= false;
end

% metis root directory
root = fileparts(which('metis'));

% add path to the cocos test function
root_cocos = fullfile(root,'acces','imas','cocos_test');
if ~isdir(root_cocos)
    s = mkdir(root_cocos);
    if ~s
        error(sprintf('unable to make directory %s\n',root_cocos));
    end
    no_update= false;
end
addpath(root_cocos);

% test if METIS is owned by the user
[s,t] = unix(sprintf('ls -l %s',which('metis')));
if  s~= 0
    tobeinstalled = false;
    disp(t)
else
    t = strtok(t,'/');
    if  ~isempty(strfind(strrep(t,which('metis'),''),getenv('USER')))
        tobeinstalled = true;       
    else
        tobeinstalled = false;        
    end
end

% look if installation is needed in mode no_update
if no_update && tobeinstalled
    if isdir(fullfile(root_cocos,'cocos')) && isdir(fullfile(root_cocos,'cocos','matlab'))
        tobeinstalled = false; 
    end
end


% installation/update
if tobeinstalled
    
    % mv previous version of cocos tools
    if isdir(fullfile(root_cocos,'cocos'))
        if isdir(fullfile(root_cocos,'cocos.old'))
            [s,t] = unix(sprintf('rm -rf %s',fullfile(root_cocos,'cocos.old')));
            if s~= 0
                fprintf('Warning: unable to remove old version of cocos tools:\n%s\n',t);
            end
        end
        [s,t] = unix(sprintf('mv -f %s %s',fullfile(root_cocos,'cocos'),fullfile(root_cocos,'cocos.old')));
        if s == 0
            old = true;
        else
            error(t);
        end
    else
        old = false;
    end
    
    
    % clone last version of tools made by Olivier Sauter
    pwd_mem = pwd;
    cd(root_cocos);
    if ~isempty(which('gitclone'))
       try
            repo = gitclone('https://gitlab.epfl.ch/spc/cocos.git');
            if ~isempty(repo.GitFolder)
                s = 0;            
                t = repo.GitFolder;
            else
                s = 1;
                t = 'error cloning gitlab.epfl.ch/spc/cocos.git';
            end
        catch
            repo = [];
            t = lasterror;
            s = 1;
        end
        if isempty(repo)
            s = 1;
        end
    else
        [s,t] = unix('export path=''/bin:/usr/bin''; git -c http.sslVerify=false clone https://gitlab.epfl.ch/spc/cocos.git');
    end
    cd(pwd_mem);
    if s == 0
        fprintf('COCOS repository has been cloned:\n%s\n',t);
        addpath(fullfile(root_cocos,'cocos','matlab'));
        % remove .git directory
        [s,t] = unix(sprintf('rm -rf %s',fullfile(root_cocos,'.git')));
        if s ~= 0
            fprintf('unable to remove .git directory:\n%s\n',t);
        end
        [s,t] = unix(sprintf('rm -rf %s',fullfile(root_cocos,'cocos','.git')));
        if s ~= 0
            fprintf('unable to remove .git directory:\n%s\n',t);
        end
        % try to update .csv file with current IMAS version
        imas_root = getenv('IMAS_prefix');
        if isempty(imas_root)
            imas_root = fileparts(fileparts(which('imas_open_env_backend')));
        end
        % install file corresponding to in use IMAS version
        if ~isempty(imas_root) && local_imas_version
           disp('Try to install conversion rules (.csv files) corresponding to in use IMAS version');
           [s,t] = unix(sprintf('find %s  -name "grid_type_transformations_symbolic_table.csv" -print | grep "grid_type_transformations_symbolic_table.csv"',imas_root)); 
           if s == 0
               tt = tseparec(t);
               source = strtrim(tt(1,:));
               [smv,tmv] = unix(sprintf('mv -f %s %s',fullfile(root_cocos,'cocos','matlab','grid_type_transformations_symbolic_table.csv'), ...
                                                  fullfile(root_cocos,'cocos','matlab','grid_type_transformations_symbolic_table.csv.origin')));
               [s,t] = unix(sprintf('cp -pf %s %s',source,fullfile(root_cocos,'cocos','matlab','grid_type_transformations_symbolic_table.csv')));
               if s ~= 0
                   disp('unable to install grid_type_transformations_symbolic_table.csv coherent with present version of IMAS');
                   disp(t);
                   if smv == 0
                        unix(sprintf('mv -f %s %s',fullfile(root_cocos,'cocos','matlab','grid_type_transformations_symbolic_table.csv.origin'), ...
                                                  fullfile(root_cocos,'cocos','matlab','grid_type_transformations_symbolic_table.csv')));
                   else
                       disp(tmv);
                   end
               end               
           else
               disp('unable to install grid_type_transformations_symbolic_table.csv coherent with present version of IMAS');
           end
           %
           [s,t] = unix(sprintf('find %s  -name "ids_cocos_transformations_symbolic_table.csv" -print | grep "ids_cocos_transformations_symbolic_table.csv"',imas_root)); 
           if s == 0
               tt = tseparec(t);
               source = strtrim(tt(1,:));
               [smv,tmv] = unix(sprintf('mv -f %s %s',fullfile(root_cocos,'cocos','matlab','ids_cocos_transformations_symbolic_table.csv'), ...
                                                  fullfile(root_cocos,'cocos','matlab','ids_cocos_transformations_symbolic_table.csv.origin')));
               sc = copyfile(source,fullfile(root_cocos,'cocos','matlab','ids_cocos_transformations_symbolic_table.csv'),'f');
               if sc ~= 1
                   disp('unable to install ids_cocos_transformations_symbolic_table.csv coherent with present version of IMAS');
                   disp(t);
                   if smv == 0
                        unix(sprintf('mv -f %s %s',fullfile(root_cocos,'cocos','matlab','ids_cocos_transformations_symbolic_table.csv.origin'), ...
                                                  fullfile(root_cocos,'cocos','matlab','ids_cocos_transformations_symbolic_table.csv')));
                   else
                       disp(tmv);
                   end
               end               
           else
               disp('unable to install ids_cocos_transformations_symbolic_table.csv coherent with present version of IMAS');
           end
       end
        
    else
        fprintf('Unable to clone cocos repository:\n%s\n',t);
        if old
            [s,t] = unix(sprintf('cp -rpf %s %s',fullfile(root_cocos,'cocos.old'),fullfile(root_cocos,'cocos')));
            if s ~= 0
                error(t);
            else
                disp('Previous version restored')
                addpath(fullfile(root_cocos,'cocos','matlab'));
            end
        else
            error('unable to get COCOS tools');
        end
    end   
    
%    % open readme
%     try
%         edit(fullfile(root_cocos,'cocos','matlab','README.md'));
%         edit(fullfile(root_cocos,'cocos','matlab','README_tags'));
%     catch
%        % if no graphics
% in any case just display files to be less intrusif for GUI users.
        disp('==============================================================');
        disp('COCOS tools README:');
        type(fullfile(root_cocos,'cocos','matlab','README.md'));
        type(fullfile(root_cocos,'cocos','matlab','README_tags'));
        disp('==============================================================');
        disp(' ');
        
%    end


    
else
    addpath(fullfile(root_cocos,'cocos','matlab'));
    disp('Path to COCOS tools added');
end

% that all !



