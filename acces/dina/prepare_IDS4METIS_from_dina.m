% This function prepare automatically a METIS IMAS dataset from a DINA IMAS dataset.
%
% syntax:
%  [post,error_flag] = prepare_IDS4METIS_from_dina(dina_uri,output_uri,metis_reference,metis_mode,nbt)
%
% input:
%     dina_uri          = URI for the inputed DINA IMAS dataset
%     output_uri        = URI of the METIS dataset created in IMAS format
%     metis_reference   = METIS reference file (including path) for the
%                         configuration of the models in the METIS simulation
%     metis_mode        = either 'interpretative' or 'predictive'.
%                         In interpretative mode METIS use kinetic profiles
%                         from DINA dataset (Te,Ti,ne). In interpretative,
%                         all kinetic fields are computed using METIS
%                         model.
%     nbt               = number of time slices in the simulation (default=301) 
%
% output (for debug):
%    post       = METIS results
%    error_flag = status of IMAS IDSs creation
%
%for testing
% prepare_IDS4METIS_from_dina('imas:hdf5?path=/home/ITER/artaudj1/public/imasdb/iter/3/105092/1', ...
%                'imas:hdf5?path=/home/ITER/artaudj1/public/imasdb/iter/3/105092/1001',...
%               './certification/metis/reference_NTM_ITER.mat','interpretative',301);
%
% prepare_IDS4METIS_from_dina('imas:hdf5?path=/home/ITER/artaudj1/public/imasdb/iter/3/105092/1', ...
%                'imas:hdf5?path=/home/ITER/artaudj1/public/imasdb/iter/3/105092/1002',...
%               './certification/metis/reference_NTM_ITER.mat','predictive',301);
%
function [post,error_flag] = prepare_IDS4METIS_from_dina(dina_uri,output_uri,metis_reference,metis_mode,nbt)

% display some context information
path

%  input handeling
if nargin < 3 
    error('post = prepare_IDS4METIS_from_dina(dina_uri,output_uri,metis_reference,metis_mode,nbt)')
elseif (nargin < 4) || isempty(nbt)
    nbt = 301;
end

fprintf('reading DINA data from %s\n',dina_uri);
fprintf('Genrating METIS data in %s\n',output_uri);
fprintf('METIS reference dataset is %s\n',metis_reference);
fprintf('Running mode = %s\n',metis_mode);
fprintf('Number of time slices = %d\n',nbt);


% remove external data
rm_cs4m;

% read DINA data and make a METIS data set
z0dinput = init_metis_from_dina_ids(metis_reference,dina_uri,nbt);
zassignin("base","z0dinput",z0dinput);

% generate external data if needed
switch metis_mode
    case 'interpretative'
     list = {'Ne','Te','Ti','Momentum','Zeff','LHCD','ECRH','ICRH','NBI','NBI2','PLINE','RUNAWAY'};
     for k=1:length(list)
         try
             cs4m_imas(list{k},'no_inter');
             fprintf('Adding external data %s for METIS\n',list{k});
         catch
             fprintf('No data for %s for METIS\n',list{k});
         end
     end
end

% run metis
evalin('base','metis_run;');

% save metis run for debug
evalin('base',sprintf('metis_save(''METIS_from_DINA_%s_%d_%s'');', ...
               z0dinput.option.machine,z0dinput.option.shot,metis_mode));

% generate IDSs
post = evalin('base','post;');
error_flag = metis4imas(output_uri,0,0,post);
close all
