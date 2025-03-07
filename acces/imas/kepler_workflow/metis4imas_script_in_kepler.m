% Matlab script to call METIS through the Matlab actor 
% open logfile for debugging :
diary(sprintf('metis4imas_%d_%d_%d.txt',shot,run,occurrence_in));
pwd
disp('start metis4imas');
fprintf('shot = %d; run = %d; occurrence_in = %d; time = %d; method = %s; interpolation = %d\n', ...
shot,run,occurrence_in,time,method,interpolation_method);
disp('code param input = ');
disp(code_param_xml);
% detection of initialisation
if ~ exist('init_state','var')
   init_state = 1;
   if strcmp(method,'one_time')
         method = 'init';
   end 
elseif isempty(init_state)
   init_state = 1;
   if strcmp(method,'one_time')
         method = 'init';
   end 
else
   init_state = 0;
end
switch method
case 'init'
  if isappdata(0,'IMAS_Z0DSTRUCT')
  	rmappdata(0,'IMAS_Z0DSTRUCT');
  end
end
save(sprintf('metis4imas_input_%d_%d_%d',shot,run,occurrence_in));
disp('save input data of metis4itm');
% initialize UAL environnement :
if not(exist('metis4imas.m'))
%       try
%       	import_ual;
%      catch
%	disp('UAL script configuration uavaible,try alternative method')
%                      javaaddpath( [getenv('UAL') '/jar/ualjava.jar'] )
%                      addpath( [getenv('UAL') '/matlabinterface'] )
%                      addpath( [getenv('MDSPLUS_DIR') '/matlab'] )
%       end
 %       addpath /home/ITER/artaudj1/public/METIS/metis
        addpath(fileparts(fileparts(fileparts(fileparts(which('metis4imas_script_in_kepler'))))))
        zineb_path
        disp('defining Matlab path for metis4imas');
        path
else
   exist('metis4imas.m')
   which metis4imas
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
%iteration = 0;

% call to METIS with the help of metis4imas interface
% flush
diary off
!sync
diary on
disp('call metis4imas');
try
        flag_out  = metis4imas(shot,run,occurrence_in,method,time,code_param_xml,interpolation_method)
catch
	struct_lasterror = lasterror
	save contexte_error_matlab_actor
	struct_lasterror.message
end
disp('end metis4imas');
% flush
diary off
!sync
diary on
% traitement of error
if  flag_out ~= 0
	occurrence_out = 0;
	output =-1 ;
        disp('error executing metis4imas');
	disp('end of metis4imas');
	diary off
else 
	% handling of actor output
	disp('retreive data of metis4imas');
	z0dstruct = getappdata(0,'IMAS_Z0DSTRUCT');
        if isempty(z0dstruct)			
        	appdata_zero = getappdata(0);
 		struct_lasterror = lasterror
		struct_lasterror.message
       		save(sprintf('metis4itm_context_%d_%d_%d',shot,run,occurrence_in));
        	z0dstruct
        	z0dstruct.z0dinput
        	z0dstruct.z0dinput.option
                error(lasterror);
	elseif ~isempty(z0dstruct.z0dinput.option.summary_occurrence)
		occurrence_out = z0dstruct.z0dinput.option.summary_occurrence;
	else
		occurrence_out = 0;
	end
	disp('save output data of metis4imas');
	save(sprintf('metis4itm_output_%d_%d_%d',shot,run,occurrence_in));
	% flush
	diary off
        !sync
	diary on
	if isempty(time)
		output = z0dstruct.zerod.temps(end)
	else
		output =time;
	end
	disp('end of metis4imas');
	diary off
end
