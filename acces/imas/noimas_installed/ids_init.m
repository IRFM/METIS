%   Generate an IDS of specified type with default field values.
%   
%   IDSname : name of the IDS to generate.
%  
%   NOTE: The array of structures in the resulting IDS will be
%   empty. This is not as in ids_gen where they are all of size 1. Use
%   ids_allocate to fill the array of structures.
%  
%   See also: ids_allocate
% 
function  ids_out = ids_init(IDSname)

% init
persistent output
if isempty(output)
   load(fullfile(fileparts(which('litidss')),'noimas_installed','DATA4IMAS_WITHOUT_INFRASTRUCTURE'),'output');
end

ids_out = output.(IDSname);
noms = fieldnames(ids_out);
for k=1:length(noms)
    if iscell(ids_out.(noms{k}))
        ids_out.(noms{k}) = {};
    end
end


  
