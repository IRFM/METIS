function set_feeqs_path_jt60sa(folder_name)

if nargin == 0
    folder_name ='';
end  
if ~isempty(folder_name)  && ~exist(fullfile(folder_name,'start_up_Unix.m'))
	  warning('invalid path to FEEQS.M');
end

% test if FEEQS is in the path
if isempty(which('start_up_Unix.m'))

  if isempty(folder_name)
      list_path{1} = '.';	
      switch getenv('USER')
       case 'JA132999'
	          list_path{end+1} = '/ZONE_TRAVAIL/JA132999/feeqs.m/FEEQS.M';
	          list_path{end+1} = '/ZONE_TRAVAIL/JA132999/feeqs.m/FEEQS.M.Git/feeqs/FEEQS.M/';
      end
      if isappdata(0,'PATH2FEEQS')
	list_path{end+1} =  getappdata(0,'PATH2FEEQS');
      end
      list_path{end+1} = '/ZONE_TRAVAIL/JA132999/public/FEEQS.M';
      list_path{end+1} = '/donnees/JA132999/public/FEEQS.M';
      list_path{end+1} = '/afs/eufus.eu/g2itmdev/user/g2jfa/public/FEEQS.M';

      %
      for k=1:length(list_path)
	if exist(fullfile(list_path{k},'start_up_Unix.m'))
	    folder_name = list_path{k};
        break
	end
      end
  end
  if isempty(folder_name)
	folder_name = uigetdir;
  end
  if isempty(folder_name) || ~exist(fullfile(folder_name,'start_up_Unix.m'))
	  error('this is not a path to FEEQS.M');
  end

  % set FEEQS environement for JT-60SA
  addpath(folder_name);
  dirmem = pwd;
  cd(folder_name)
  evalin('base','start_up_Unix');
  cd(dirmem);
  addpath(fullfile(folder_name,'Projects','JT60_SA'));
  addpath(fullfile(folder_name,'Projects','JT60_SA','Lib'));
  addpath(fullfile(folder_name,'Projects','JT60_SA','Data'));
  %% addpath(fullfile(folder_name,'Projects','JT60_SA','YAML'));
  addpath(fullfile(folder_name,'Projects','Reactors'));

  % add path to fast code
  addpath(fullfile(folder_name,'Projects','FastCoilCurrentIdentification'));
  addpath(fullfile(folder_name,'Projects','FastCoilCurrentIdentification','Lib'));
  


end  
