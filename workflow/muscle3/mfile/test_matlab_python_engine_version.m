% test python version from module
function test_matlab_python_engine_version
py_info = pyenv;
[version_major,version_minor] = strtok(char(py_info.Version),'.');
version_minor = strtrim(version_minor(2:end));
if str2double(version_major) < 3
    python_version_ok = false;
elseif str2double(version_major) > 3
    python_version_ok = true;
elseif str2double(version_minor) < 6
    python_version_ok = false;
else
    python_version_ok = true;    
end
if ~python_version_ok
    disp('Wrong Python version !')
    [s,t] = unix('which python');
    disp('Launch Matlab separately and run the command:')
    fprintf('pyenv(''Version'',''%s'')\n',strtrim(t));
    disp('and exit from Matlab')
    error('Wrong Python version used in Matlab');
end
