function Out = ids_gen(IDSname)
% Out = ids_gen(IDSname)
% Create full Matlab IDS structures
% 
% IDSname : name of the IDS to create

% YB, YF

persistent output

% check number of arguments
if (nargin ~= 1)
    error('Bad number of input arguments (1: IDS name) .');
end

% check validity of the IDS name
if ~ischar(IDSname);
    error('First input argument must be a string (IDS name)')
end
%IDSname = lower(IDSname);
if isempty(strmatch(IDSname, IDS_list,'exact'))
    msgerror = strcat(IDSname,' is not a valid IDS name!!');
    error(msgerror);
end

if isempty(output)
    load(fullfile(fileparts(which('litidss')),'noimas_installed','DATA4IMAS_WITHOUT_INFRASTRUCTURE'),'output');
end
% create a Java object and convert it to a Matlab structure
Out = output.(char(IDSname));
