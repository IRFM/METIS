function [machine,user,ver,backend] = which_MDSdatabase_metis_imas




% ---------------------------------------------------------------------
% From LUKE - Function that check which UAL MDS+ database is connected
% by Y. Peysson CEA-IRFM <yves.peysson@cea.fr> and Joan Decker
% Modified by M. Schneider to work with IMAS environment (imasdb)
% ---------------------------------------------------------------------

%hdf5base    = getenv('HDF5_BASE');
mdstreebase = getenv('MDSPLUS_TREE_BASE_0');

%if isempty(hdf5base) | isempty(mdstreebase),
if isempty(mdstreebase)
    ver = getenv('IMAS_VERSION');
    machine = '';
    user = '';
    if isempty(ver)
	warning('Environment variables for IMAS MDS+/HDF5 database are not set up.');
   end
   return
end

is = 1;
str = {};
%remain = hdf5base;
remain = mdstreebase;
while true
  [strtemp,remain] = strtok(remain,'/');	
  if isempty(strtemp)
    break
  end
  str{is} = strtemp;
  is = is + 1;
end

%% Private MDS+ database (may be own by any IMAS user)
if strcmp(str{1},'home'),
  ver = str{length(str)-1};
  machine = str{length(str)-2};
  user = str{2};
%% Public MDS+ database
else
  ver = str{length(str)-1};
  machine = str{length(str)-2};
  user = str{1};
end	

backend = getenv('IMAS_AL_BACKEND');

