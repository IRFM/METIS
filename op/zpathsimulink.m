% script d'ajout des path pour simulink
function zpathsimulink
ch = which('date');
ind = max(findstr(ch,'matlab'))-1;
ch  = fullfile(ch(1:ind),'simulink');
addpath(fullfile(ch,'simulink'));
addpath(fullfile(ch,'simdemos'));
addpath(fullfile(ch,'fixedandfloat'));
addpath(fullfile(ch,'dee'));
addpath(fullfile(ch,'dastudio'));
addpath(fullfile(ch,'components'));
addpath(fullfile(ch,'blocks'));



