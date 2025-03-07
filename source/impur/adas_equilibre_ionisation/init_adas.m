function init_adas

setappdata(0,'ADAS_ROOT','/Home/adas/adas');
setappdata(0,'ADAS_URL','open.adas.ac.uk');
setappdata(0,'ADAS_PREF_ORDER',[12,96,93,92,89,85,74,50]);

[p,f]=fileparts(which('init_adas'));
addpath(p);
addpath(fullfile(p,'import'));



