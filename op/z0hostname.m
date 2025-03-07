function name = z0hostname
[ret, name] = system('hostname');
if ret ~= 0,
  if ispc
    name = getenv('COMPUTERNAME');
    else
    name = getenv('HOSTNAME');
  end
end
name = lower(name);
name(name < 32) = [];
name = strtrim(name);