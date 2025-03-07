% Matlab script to call METIS through the Matlab actor 
% open logfile for debugging :
diary(sprintf('metis4itm_%d_%d_%d.txt',shot,run,occurrence_in));
pwd
disp('start metis4itm');
fprintf('shot = %d; run = %d; occurrence_in = %d; time = %d; method = %s; interpolation = %d\n', ...
shot,run,occurrence_in,time,method,interpolation_method);
disp('code param input = ');
disp(code_param_xml);
%save(sprintf('metis4itm_input_%d_%d_%d',shot,run,occurrence_in));
%disp('save input data of metis4itm');
% initialize UAL environnement :
if ~exist('metis4itm.m')
        try
        	import_ual;
        catch
	disp('UAL script configuration uavaible,try alternative method')
                      javaaddpath( [getenv('UAL') '/jar/ualjava.jar'] )
                      addpath( [getenv('UAL') '/matlabinterface'] )
                      addpath( [getenv('MDSPLUS_DIR') '/matlab'] )
        end
        addpath /pfs/home/jfa/public/metis4itm/trunk
        zineb_path
        disp('define Matlab path for metis4itm');
else
   exist('metis4itm.m')
   which metis4itm
end
% discard time if value < -999999
% that allows to control from which the value of the time is read.
% if < -999999, time vetcor is read in CPO scenario
% otherwise it is an input of Matlab actor.
if time <= -999999
        time = [];
end

% initialise output
flag_out = -1;
occurrence_out = 0;
output =-1 ;

% call to METIS with the help of metis4itm interface
disp('call metis4itm');
try
    flag_out  = metis4itm(shot,run,occurrence_in,method,time,code_param_xml,interpolation_method)
catch
	struct_lasterror = lasterror
	save contexte_error_matlab_actor
	struct_lasterror.message
end

% traitement of error
if  flag_out ~= 0
	occurrence_out = 0;
	output =-1 ;
        disp('error executing metis4itm');
	disp('end of metis4itm');
	diary off
else 
	% handling of actor output
	disp('retreive data of metis4itm');
	z0dstruct = getappdata(0,'ITM_Z0DSTRUCT');
	if ~isempty(z0dstruct.z0dinput.option.scenario_occurrence)
		occurrence_out = z0dstruct.z0dinput.option.scenario_occurrence;
	else
		occurrence_out = 0;
	end
	%disp('save output data of metis4itm');
	%save(sprintf('metis4itm_output_%d_%d_%d',shot,run,occurrence_in));
        if isempty(time)
		output =z0dstruct.zerod.temps(end)
        else
		output =time;
        end
	disp('end of metis4itm');
	diary off
end

