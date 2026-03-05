% This function generate documentation for outputed METIS IMAS data
% input filename are optionnal and can be used to create alernative IDS
% documenttaion.
function make_metis_imas_output_documentation(filename)

if (nargin == 0) || isempty(filename)
    filename = fullfile(fileparts(which('metis')),'certification','metis','ITER_rampup_ECCD.mat');
end

% load METIS data
evalin('base','clear post z0dinput');
data = load_metis_imas(filename);

% change parameters to get all IDSs  possible data
option = data.z0dinput.option;
option.nb_nbi   = 2;
option.zmax     = 18;
option.zmin     = 4;
option.W_effect = 1;
option.Sn_fraction = 0.5;
option.cw_offset   = 1e-6;
option.sol_rad  = 'decoupled';
option.sol_model = '2_points';
cons  = data.z0dinput.cons;
cons.picrh = cons.pecrh;
cons.plh   = cons.pecrh;
cons.nbi   = cons.pecrh + sqrt(-1) .* cons.pecrh;
z0dinput = data.z0dinput;
z0dinput.option = option;
z0dinput.cons   = cons;
zassignin('base','z0dinput',z0dinput);
evalin('base','metis_run');
post = evalin('base','post');

% generate IDS like data structure
[error_flag,output_data] = metis4imas(1,1,0,post,[]);

% read data description
[output,idss_list,description,mixed] = litidss;

% opent text file
filetxt = fullfile(fileparts(which('metis')),'doc','metis_imas_output_documentation.md');
fid = fopen(filetxt,'w');
% format it in Markdown
fprintf(fid,'List of IMAS data produced by METIS\n');
fprintf(fid,'===================================\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'Introduction\n');
fprintf(fid,'------------\n');
fprintf(fid,'\n');
fprintf(fid,'This document list all the possible IDS fields METIS can fill.\n');
fprintf(fid,'Some data are produced only when:\n');
fprintf(fid,'\t-two point model is activated (egde data).\n');
fprintf(fid,'\t-the model for tunsten is ativated.\n');
fprintf(fid,'\t-the divertor contains tin.\n');
fprintf(fid,'\t-plasma contains D & T.\n');
fprintf(fid,'\t-plasma contains D & He3.\n');
fprintf(fid,'\t-plasma contains H & B11.\n');
fprintf(fid,'\n');
fprintf(fid,'The IDS equilibrium provides two 2D grids: one in cartesian coordinates (R,Z) ala EFIT and one inverse grid in (rho,theta).\n');
fprintf(fid,'IDSs named *_1 are alternative IDS for edge data stored in occurrence 1.\n');
fprintf(fid,'All other IDSs are stored in occurrence 0.\n');

fprintf(fid,'\n');

% loop on produced IDSs
noms = fieldnames(output_data);
for k=1:length(noms)
    nomloc = strrep(noms{k},'_1','');
    nomloc = strrep(nomloc,'_2','');
    nomloc = strrep(nomloc,'_3','');    
    len = fprintf(fid,'IDS %s\n',noms{k});
    fprintf(fid,char(abs('-') * ones(len,1)));
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    scan_list_variable_imas(output_data.(nomloc),description.(nomloc),noms{k},fid);
    fprintf(fid,'\n');
    fprintf(fid,'\n');
end

fclose(fid);





