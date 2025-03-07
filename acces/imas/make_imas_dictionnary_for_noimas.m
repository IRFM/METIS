% generate data for ids_gen without infrastructure 
[output,idss_list,description,mixed] = litidss;
data_version = getenv('IMAS_VERSION');
save(fullfile(fileparts(which('litidss')),'noimas_installed','DATA4IMAS_WITHOUT_INFRASTRUCTURE'),'output','idss_list','description','mixed','data_version');

