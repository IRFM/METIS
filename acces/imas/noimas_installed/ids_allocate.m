  
%   Allocate the desired array of structure to the requested size.
%   
%   IDSname : the name of the IDS to modify.
%   AOSpath : path to the desired array of structure in the IDS ('/' separated following DD
%     convention)
%   size    : desired size.
%  
%   Example:
%      eq = ids_init('equilibrium');
%      eq.time_slice = ids_allocate('equilibrium','time_slice',5);
%      eq.time_slice{1}.profiles_2d = ids_allocate('equilibrium','time_slice/profiles_2d',1);
%  
%   See also : ids_init
% 
  
function aos_out =  ids_allocate(IDSname, AOSpath, sizeAos)

% init
persistent output
if isempty(output)
   load(fullfile(fileparts(which('litidss')),'noimas_installed','DATA4IMAS_WITHOUT_INFRASTRUCTURE'),'output');
end

model = output.(IDSname);
cellmode = false;
while(~isempty(AOSpath))
    [loc_path,AOSpath] = strtok(AOSpath,'/');
    if ~isempty(AOSpath)
        AOSpath = AOSpath(2:end);
    end
    if any(loc_path == '{')
        indc = find(loc_path == '{',1);
        loc_path(indc:end) = [];
    elseif any(loc_path == '(')
        indc = find(loc_path == '(',1);
        loc_path(indc:end) = [];
     end
    if iscell(model.(loc_path))
        model = model.(loc_path){1};
        cellmode = true;
    elseif isstruct(model.(loc_path))
        model = model.(loc_path)(1);
        cellmode = false;
    else
        model = model.(loc_path);
        cellmode = false;
    end
end   
aos_out  = [];
for k=1:sizeAos
    if cellmode
        aos_out{k} = model;
    else
        aos_out(k) = model;       
    end
end




