function out = IDS_list
% Return the list of existing IDSs
% This file was generated automatically
persistent idss_list
if isempty(idss_list)
   load(fullfile(fileparts(which('litidss')),'noimas_installed','DATA4IMAS_WITHOUT_INFRASTRUCTURE'),'idss_list');
end
out = idss_list;

