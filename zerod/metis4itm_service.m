% METIS4ITM_SERVICE: this is a service (unix like) that pilote metis4itm through file exchange
% the compile version must be used
% environnements variables are used to configure the service:
%
%	KEPLER_CMD4METIS      = path + root on file name used as input file (default = cmd4metis)
%       METIS_OUTPUT4KEPLER   = path + root on file name used as output file (default = metis4kepler)
%       METIS_SERVICE_TIMEOUT = time out in s (maximum sleep time before recieve a command, default = 300)
%
% each file have a rank (command number), so the name of input/ouput file are <string>@<command_number>
% the command number is increase of one unity after each call
% the command 0 reset the service (at any time)
% 
% the command file must contains any of these lines :
%
%  		exit_flag = 0
%  		shot = 1
%  		run  = 1
%  		occurrence =
%  		method =
%  		time = 10
%  		codeparam_xml =
%  		interpolation_method =
%  		tokamak =
%  		user =
%  		dataversion =
%  		shot_ref =
%  		run_ref  =
%
%  the output file contains:
%
%
%  		output = 10
%  		occurrence_out = 0
%  		flag_out = 0
%
function metis4itm_service

% service metis4itm 
if isdeployed 
    dir
    if exist('Metis_splash_screen3.png')
        try
             splash('Metis_splash_screen3','png');
             drawnow
        catch
          f = errordlg(lasterr, 'splash error');
          waitfor(f);
        end
    else
         f = errordlg('unable to find image', 'splash error');
         waitfor(f);
    end
    maxNumCompThreads('automatic')
else
	% root du programme et initilisation du path
	root = getappdata(0,'root');
	if isempty(root)
		addpath /pfs/work/jfa/Project/metis4itm/trunk
		zineb_path
		disp('define Matlab path for metis4itm');
		zineb_path;
		root = getappdata(0,'root');
	end
end

% version
[zver,zdate]        = zinebversion;

clear functions

setappdata(0,'langue_cronos','anglais');
setappdata(0,'uicrossref',[]) ; % securite
rmappdata(0,'uicrossref') ;

if exist('import_ual')
	try 
		import_ual;
        catch
		disp('UAL script configuration uavaible,try alternative method')
                javaaddpath( [getenv('UAL') '/jar/ualjava.jar'] )
                addpath( [getenv('UAL') '/matlabinterface'] )
                addpath( [getenv('MDSPLUS_DIR') '/matlab'] )
 	end
end
%link with UAL
try
	rep = javaclasspath('-all');
	for l=1:length(rep)
		sl = rep{l}; 
		ind = findstr(sl,'ualjava.jar');
		if ~isempty(ind)
			[f,s,e] = fileparts(fileparts(fileparts(sl)));
			%warning off
			%addpath('/afs/efda-itm.eu/project/switm/');
			%addpath(sprintf('/afs/efda-itm.eu/project/switm/ual/%s%s/matlabinterface',s,e));
			%warning on            
            		setappdata(0,'UALVERSION',sprintf('%s%s',s,e));
            		% use :  any(getappdata(0,'UALVERSION') >'4.08b') to test if
            		% ual version is greater          
		end
	end
end
%  addpath('/afs/efda-itm.eu/project/switm/');
%  addpath('/afs/efda-itm.eu/project/switm/ual/4.07b/matlabinterface');
if isempty(import)
	import ualmemory.javainterface.*;
end

% récupere le nom du fichier d'entree
fnin = getenv('KEPLER_CMD4METIS');
if isempty(fnin)
	fnin='cmd4metis';
end
fnout = getenv('METIS_OUTPUT4KEPLER');
if isempty(fnout)
	fnout='metis4kepler';
end
time_out = getenv('METIS_SERVICE_TIMEOUT');
if isempty(time_out)
	time_out= 300; %s
else
	try
		time_out = str2num(time_out);
	catch
		time_out= 300; %s
	end
end
pwd
disp('start metis4itm');

% boucle avec attente d'appel a METIS4ITM
exit_flag  = 0;
file_number = 0;
shot = [];
run_shot  = [];
occurrence_in = [];
method ='';
time = NaN;
codeparam_filename= '';
interpolation_method='';
tokamak='';
user ='';
dataversion = '';
shot_ref = [];
run_ref  = [];
next_file_number = file_number;

while (exit_flag == 0) && (file_number < 1e9)

	[exit_flag,shot,run_shot ,occurrence_in,method,time,codeparam_xml,interpolation_method,tokamak, ...
        user,dataversion,shot_ref,run_ref,next_file_number] =  ...
        read_metis_service_input(file_number,fnin,time_out);

	if exit_flag
		break
	end

        % gestion des entrees
        if ~isempty(tokamak)
		occurrence_in.occurrence  = occurrence_in;
		occurrence_in.tokamak     = tokamak;
		occurrence_in.user	  = user;
		occurrence_in.dataversion = dataversion;
	end
	if ~isempty(shot_ref)
		if isempty(run_ref)
			shot      = shot + sqrt(-1) .* shot_ref;
			run_shot  = run_shot + sqrt(-1)  .* run_shot;
		else
			shot      = shot + sqrt(-1) .* shot_ref;
			run_shot  = run_shot + sqrt(-1)  .* run_ref;
		end
	end
	fprintf('shot = %d; run = %d; occurrence_in = %d; time = %d; method = %s; interpolation = %d\n', ...
	shot,run_shot,occurrence_in,time,method,interpolation_method);
	disp('code param input = ');
	disp(codeparam_xml);
	disp('call metis4itm');
	try
		flag_out  = metis4itm(shot,run_shot,occurrence_in,method,time,codeparam_xml,interpolation_method)
	catch
		struct_lasterror = lasterror
		save contexte_error_matlab_metis4itm_actor
		struct_lasterror.message
	end
	if  flag_out ~= 0
		occurrence_out = 0;
		output =-1 ;
		disp('error executing metis4itm');
		disp('end of metis4itm');
	else 
		% lecture de l'occurrence de sortie
		disp('retreive data of metis4itm');
		z0dstruct = getappdata(0,'ITM_Z0DSTRUCT');
                if isempty(z0dstruct)
			occurrence_out = 0;
			output =-1 ;
		elseif ~isempty(z0dstruct.z0dinput.option.scenario_occurrence)
			occurrence_out = z0dstruct.z0dinput.option.scenario_occurrence;
			output =time;
		else
			occurrence_out = 0;
			output =time;
		end
		% ecriture du fichier de sortie
        	write_metis_service_output(file_number,fnout,output,occurrence_out,flag_out);
	end
        % mise ajour du numero de fichier
	file_number = next_file_number;
end
disp('end of metis4itm');
diary off
if isdeployed 
	exit;
end


function [exit_flag,shot,run,occurrence,method,time,codeparam_xml,interpolation_method, ...
          tokamak,user,dataversion,shot_ref,run_ref,next_file_number] =  ...
                       read_metis_service_input(file_number,fnin,time_out)

% prototype de sortie
exit_flag = 0;
shot = 1;
run  = 1;
occurrence = '';
method = '';
time= 0;
codeparam_xml ='';
interpolation_method = [];
tokamak='';
user='';
dataversion='';
shot_ref = [];
run_ref  = [];
next_file_number =NaN;

% gestion des entrees
fid        = [];
reinit     = 0;
dt         = 0.1;
exit_flag  = 0; 
fnin_num0   = sprintf('%s@%d',fnin,0);
if exist(fnin_num0,'file')
	reinit =1;
        next_file_number = 1;
	fid = fopen(fnin_num0,'r');
else
	fnin_num = sprintf('%s@%d',fnin,file_number);
        fprintf('METIS4ITM_service : waiting for new command (#%d in file %s)\n',file_number,fnin_num);
	while(time_out > 0)
		if exist(fnin_num,'file')
			fid = fopen(fnin_num,'r');
        		next_file_number =  file_number + 1;
			break
		elseif exist(fnin_num0,'file')
			reinit =1;
        		next_file_number = 1;
			fid = fopen(fnin_num0,'r');
			fnin_num = fnin_num0;
			break			
		end
		pause(dt);
		time_out = time_out - dt;
	end
end
if isempty(fid)
	disp('METIS4ITM_SERVICE : time out');
        exit_flag  = 1; 
        return
elseif fid < 0
	disp('METIS4ITM_SERVICE : error opening input file');	
        exit_flag  = 1; 
        return
end

% lecture
rawdata = textscan(fid,'%s%s','delimiter', '=');
fclose(fid);
delete(fnin_num);


% boucle sur le cellule
if isempty(rawdata)
	return
end
for k=1:length(rawdata{1})
	nom = rawdata{1}{k};
        nom(nom <= ' ') = [];
	switch lower(nom)
	case 'exit_flag'
		try 
			r = str2num(rawdata{2}{k});
			if ~isempty(r)
				exit_flag = fix(r);
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'shot'
		try 
			r = str2num(rawdata{2}{k});
			if ~isempty(r)
				shot = fix(r);
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'run'
		try 
			r = str2num(rawdata{2}{k});
			if ~isempty(r)
				run = fix(r);
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'occurrence'
		try 
			r = deblank(rawdata{2}{k});
			if ~isempty(r)
				occurrence = r;
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'method'
		try 
			r = deblank(rawdata{2}{k});
			if ~isempty(r)
				method = r;
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'time'
		try 
			r = str2num(rawdata{2}{k});
			if ~isempty(r)
				time = r;
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'codeparam_xml'
		try 
			r = deblank(rawdata{2}{k});
			if ~isempty(r)
				codeparam_xml = r;
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'interpolation_method'
		try 
			r = str2num(rawdata{2}{k});
			if ~isempty(r)
				interpolation_method = fix(r);
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'tokamak'
		try 
			r = deblank(rawdata{2}{k});
			if ~isempty(r)
				tokamak = r;
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'user'
		try 
			r = deblank(rawdata{2}{k});
			if ~isempty(r)
				user = r;
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'dataversion'
		try 
			r = deblank(rawdata{2}{k});
			if ~isempty(r)
				dataversion = r;
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'shot_ref'
		try 
			r = str2num(rawdata{2}{k});
			if ~isempty(r)
				shot_ref = fix(r);
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	case 'run_ref'
		try 
			r = str2num(rawdata{2}{k});
			if ~isempty(r)
				run_ref = fix(r);
			end
		catch
			fprintf('METIS4ITM_SERVICE : unable to interpret field name %s in input file %s (value = %s)\n', ...
                                 nom,fnin_num,rawdata{2}{k});			
		end
	otherwise
		fprintf('METIS4ITM_SERVICE : unknown field name %s in input file %s\n',nom,fnin_num);
	end
end


function write_metis_service_output(file_number,fnout,output,occurrence_out,flag_out)

% ouverture du fichier temporaire
ftout = tempname;
fnout_num = sprintf('%s@%d',fnout,file_number);
fid = fopen(ftout,'w');
if fid < 0
	disp('METIS4ITM_SERVICE : error opening output file');	
        return
end
% ecriture
fprintf(fid,'output = %g\n',output);
if ischar(occurrence_out)
	fprintf(fid,'occurrence_out = %s\n',occurrence_out);
else
	fprintf(fid,'occurrence_out = %d\n',occurrence_out);
end
fprintf(fid,'flag_out = %d\n',flag_out);
fclose(fid);
% deplacement (empeche de lire un fichier partiellement ecrit)
movefile(ftout,fnout_num);

