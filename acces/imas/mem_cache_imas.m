%% ----------------------------------
%% FUNCTION FOR UAL CACHE MANAGEMENT
%% ----------------------------------
function [idx,state] = mem_cache_imas(idx)

% not available anymore in current Matlab API
state = 0;
return



if isempty(import) && ~isappdata(0,'JAVAINTERFACESET') && isempty(getappdata(0,'JAVAINTERFACESET'))
  import ualmemory.javainterface.*;
end
setappdata(0,'JAVAINTERFACESET','done')

if nargin  == 0
  idx = [];
  state = 0;
  setappdata(0,'imas_enable_mem_cache_imas',0);
  return
elseif isempty(idx)
  state = 0;
  setappdata(0,'imas_enable_mem_cache_imas',0);
  return
end

%% MECHANISM NOT USED UNTIL A SYNCRHONIZATION WITH KEPLER IS SET UP
setappdata(0,'imas_enable_mem_cache_imas',0);
state = 0;
return

try
  UALAccess.enableMemCache(idx)
  setappdata(0,'imas_enable_mem_cache_imas',1);
  state = 1;
catch
  fprintf('Unable to turn on UAL memory cache:\n%s\n',lasterr);
  setappdata(0,'imas_enable_mem_cache_imas',0);
  state = 0;
end
