function finish
dbclear all
disp('Exit from METIS actor has been requested');
disp('Saving Matlab variables for caller context and base workspace');
caller_fun_name = sprintf('meti4muscle3_caller_%s_context.mat',evalin('caller','mfilename'));
if exist(caller_fun_name,'file')
    evalin('caller',sprintf('save(''%s'',''-append'');',caller_fun_name));
else
    evalin('caller',sprintf('save(''%s'');',caller_fun_name));
end
evalin('base','save(''meti4muscle3_workspace_context.mat'');');
disp('Matlab appdata are:')
getappdata(0);

disp('Matlab variables in METIS MUSCLE3 actor are:');
evalin('caller','whos');
disp('Matlab variables in workspace are:');
evalin('base','whos');


% manage profiler results
profiler_status = profile('status');
switch profiler_status.ProfilerStatus
    case 'on'
         disp('Preparing profiler report');
         profsave(profile('info'),'metis4muscle3_profile_results')
         profile report
end
