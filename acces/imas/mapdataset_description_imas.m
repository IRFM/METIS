%% ------------------------------------------------------------
%% CONVERT DATA FROM METIS TO DATASET_DESCRIPTION IDS (MAPPING)
%% ------------------------------------------------------------
function dataset_description = mapdataset_description_imas(z0dstruct,data_zerod,texte_diary,error_flag,run,shot,code)

%% generate empty ids
dataset_description = ids_gen('dataset_description');

%% DEFINE  CODE AND IDS_PROPERTIES STRUCTURES
dataset_description.code = code; 
dataset_description.output_flag = error_flag;
dataset_description.ids_properties.homogeneous_time = 1;
dataset_description.ids_properties.comment = texte_diary;
dataset_description.ids_properties.source    = 'METIS';
if ispc
  dataset_description.ids_properties.provider  = getenv('USERNAME');
else
  dataset_description.ids_properties.provider  = getenv('USER');
end
dataset_description.ids_properties.creation_date = sprintf('%s (julian date in second = %f)',datestr(now,'dd-mmm-yyyy HH:MM:SS'),clock2julday);

%% PARENT 
if isfield(z0dstruct.z0dinput,'dataset_description')
	dataset_description.parent_entry = z0dstruct.z0dinput.dataset_description.data_entry;
end

%% MACHINE
if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK'))
    dataset_description.data_entry.machine = getappdata(0,'UAL_TOKAMAK');
else
	itminfo = getenv('MDSPLUS_TREE_BASE_0');
	if ~isempty(itminfo)
        	fst ={};
        	item =' ';
        	while(~isempty(itminfo) && ~isempty(item))
			[itminfo,item] = fileparts(itminfo);
			fst{end+1} = item;
        	end
        	dataset_description.data_entry.machine = fst{strmatch('imasdb',fst,'exact') -1};
	else
		dataset_description.data_entry.machine = z0dstruct.z0dinput.machine;
	end
end

%% USER
dataset_description.data_entry.user = getenv('USER');	

%% SHOT NUMBER 
if isempty(shot)
   dataset_description.data_entry.pulse = real(shot)
else
  dataset_description.data_entry.pulse = real(z0dstruct.z0dinput.shot(1));
end

%% RUN NUMBER
if ~isempty(run)
  dataset_description.data_entry.run = real(run);
end

% machine 
dataset_description.data_entry.machine = z0dstruct.z0dinput.machine;

% pulse type 
dataset_description.data_entry.pulse_type = 'simulation';

%% TIME stamp for the simulation
dataset_description.time = clock2julday;

%% imas_version
if ~isempty(getenv('IMAS_VERSION'))
	dataset_description.imas_version = getenv('IMAS_VERSION');
elseif isappdata(0,'UAL_DATAVERSION')
	dataset_description.imas_version = getappdata(0,'UAL_DATAVERSION');
else
    p = which('imas_open_env');
    while ~isempty(strfind(p,'ual'))
      [p,v] = fileparts(p);
    end
    [p,v] = fileparts(p);
    rep = dir(p);
    if ~isempty(rep)
      list = {};
      for k=1:length(rep)
	if isempty(strmatch(rep(k).name,{'.','..'},'exact'))
	  list{end+1} = rep(k).name;
	end
      end
      dataset_description.imas_version = sprintf('>= %s', list{1});
      ual_dataversion = getenv('IMAS_VERSION');
      setappdata(0,'UAL_DATAVERSION',ual_dataversion);
    end
end

%% dd_version -> where is the information ?
if isfield(dataset_description,'imas_version')
  dataset_description.dd_version = dataset_description.imas_version;
end
 
%% simulation substructure
if isappdata(0,'root')
	dataset_description.simulation.comment_before = sprintf('simulation using METIS located at %s',getappdata(0,'root'));
end
if isappdata(0,'METIS_INTERFACE_TITLE')
	dataset_description.simulation.comment_after = getappdata(0,'METIS_INTERFACE_TITLE');
end
dataset_description.simulation.time_begin   = min(z0dstruct.z0dinput.cons.temps);
dataset_description.simulation.time_end     = max(z0dstruct.z0dinput.cons.temps);
dataset_description.simulation.time_restart = min(data_zerod.temps);
dataset_description.simulation.time_current = max(data_zerod.temps);
if length(data_zerod.temps) > 1
	dataset_description.simulation.time_step = mean(diff(data_zerod.temps));
else
	indicet = find(z0dstruct.z0dinput.cons.temps < data_zerod.temps(1),1);
	if ~isempty(indicet)
		dataset_description.simulation.time_step = abs(z0dstruct.z0dinput.cons.temps(indicet) - data_zerod.temps(1));
	end
end
 
%% not avilable from METIS  : time_begun and ime_ended

% workflow
if isfield(z0dstruct,'simout') && ~isempty(z0dstruct.simout)
	dataset_description.simulation.workflow = 'METIS embedded in Simulink workflow';
else
	dataset_description.simulation.workflow = 'METIS alone';
end

