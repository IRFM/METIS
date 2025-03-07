% SIMMETIS is an interface betwen Simulink and METIS (s-function)
%-------------------------------------------------------------------------
% fichier simmetis.m ->  simcronos, mdlInitializeSizes,mdlDerivatives,
%                         mdlOutputs ,mdlTerminate
%
% use : zerodevolution
% contain global variables, this is reserved names :
%     metis parameters:       Z0DINPUT_MEM
%     metis output memory :   Z0DSTRUCT_MEM
%     time offset between simulink and metis : timeoffsetS
%     fisrt computation flag :  premiercalcul
%     present indice  :  indice_courant
%     last computed scalars output : ZEROD_MEM
%     last cimputed profiles : PROFILE0D_MEM
%     list of selected output :ACTIVE_OUTPUT

% function Matlab 7 :
%
% This fnction is a s-function for simulink and
% do not be use directely in matlab. Simmetis is
% an interface betwen Simulink and METIS.
% The function use the METIS data set prescribe in parameter.
%
% the last computed data are store in matlab workspace in :
%       zerod_simmetis        = scalar datas
%	profil0d_simmetis     = profiles 
%
% the list of inputs tag are store in "list_input_metis_sfunction"
%
% syntaxe  :
%
%     [sys,x0,str,ts] = simmetis(t,x,u,flag,metis_file,list_output,save_file,verbose,input_option,z0dstruct_init);
%
% entrees :
%
%     t,s,u,flag    =  simulink s-function standard data set
%
%     metis_file    =  complete name of METIS file load a initialisation. 
%                      The metis configuration parameters are read in this file. 
%                      All reference that are not prescribed in simulink are interpolate from reference value 
%                      contains in this file. if empty used z0dinput form the workspace 
%                      (that allow to tune the simulation from the Metis interface)
%
%     list_output   = cell array that defined the output of this s-functions. if empty, the output is wth.
%                     example for scalar: list_output = {'w','vloop,'te0'}. 
%                     The output can be any scalar output or points of profiles of metis.
%                     the possible tag are store in "list_output_metis_sfunction_max" for scalar and
%                     in "list_profiles_metis_sfunction_max" for profiles. 
%                     For profiles, the syntax is <name>@<matlab subscript>. 
%                     example of profile tag : 'qjli@1:5:21'. The 'end' method is not available.  
%                     The selected list of output are store in "list_output_metis_sfunction".
%
%     save_file     = result filename where all data are store at the end of simulation, optionnel
%
%     verbose       = if = 1, display information during run
%
%     input_option  = optionnal cell list of METIS option parameter used as input of the S-function (default : {'npar0'})
%		      warning : if you change this parameter, you must also update the simulink mutiplexer 
%		      that are before S-function call in the Simulink models.
%                     To handle reference for second NBI, keywords 'pnbi2', for NBI power, and 'ftnbi2', for beam composition, 
%                     can be used as input parameter name. 
%
%     z0dstruct_init = optionnal restart internal state. Same structure than z0dstruct variable backuped inside the result file (see save_file)
% 
% sorties :
%
%     standart s-function data set, see simulink documentation
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 31/05/2005.
%
%
% liste des modifications :
%
%
%-------------------------------------------------------------------------
%
function [sys,x0,str,ts] = simmetis(t,x,u,flag,metis_file,list_output,save_file,verbose,input_option,z0dstruct_init)


if nargin < 5
	error('input file must be defined')
end
if nargin < 6
	list_output ={};
end
if nargin < 7
	save_file='';
end
if nargin < 8
	verbose = 0;
end
if nargin < 9
	input_option = {};
end
if nargin < 10
	z0dstruct_init = [];
end



% simmetis input order
input_order = {'iso','xece','ip','nbar','picrh','plh','pnbi','pecrh','hmore', ...
               'zeff','ftnbi','flux','a','R','K','d','b0','z0'};

% This example S-function illustrates how to create a variable
% step block in Simulink.  This block implements a variable step
% delay in which the first input is delayed by an amount of time
% determined by the second input.
%
%     dt      = u(2)
%     y(t+dt) = u(t)
%

% variable persistantes matlab/Cronos
% parametre de metis
global Z0DINPUT_MEM
% donnees de metis (accumulation)
global Z0DSTRUCT_MEM
% offset de temps entre simulink et metis
%global timeoffsetS
% flag de premier calcul pour metis
global premiercalcul
% indice courant dans les donnees d'entrees
%global indice_courant
% donnees scalaires courantes
global ZEROD_MEM
% profils  courants
global PROFILE0D_MEM
% liste des outputs
global ACTIVE_OUTPUT

% securite simulink
sys = [];
x0  = [];
str = [];
ts  = [];


% debug
if verbose ==1
  	fprintf('SIMMETIS : flag = %d,temps = %g\n',flag,t);
	t 
	x
	u
end

switch flag

  case 0
    %initialisation de cronos
    [sys,x0,str,ts,z0dinput,list_output] = mdlInitializeSizes(metis_file,list_output,verbose,input_order,input_option); % Initialization
    Z0DINPUT_MEM = z0dinput;
    ACTIVE_OUTPUT = list_output;
    premiercalcul = 1;

  case 3
    if premiercalcul == 1
       [sys,z0dstruct,zs,profli] = mdlOutputs(t,x,u,z0dstruct_init,Z0DINPUT_MEM,ACTIVE_OUTPUT,verbose);
       Z0DSTRUCT_MEM = z0dstruct;
       ZEROD_MEM = zs;
       PROFILE0D_MEM = profli;
       premiercalcul  = 0;
       %assignin('base','zerod_simmetis',zs);
       %assignin('base','profil0d_simmetis',profli);
      
    else
       % mise a jour
       [sys,z0dstruct,zs,profli] = mdlOutputs(t,x,u,Z0DSTRUCT_MEM,Z0DINPUT_MEM,ACTIVE_OUTPUT,verbose);      
       Z0DSTRUCT_MEM = z0dstruct;
       ZEROD_MEM = zs;
       PROFILE0D_MEM = profli;
       %assignin('base','zerod_simmetis',zs);
       %assignin('base','profil0d_simmetis',profli);
   end
   
  case 9
    % fin de la simulation
    if premiercalcul == 0
    	% il est possible de demander ici la sauvegarde des donnees de cronos
    	sys = mdlTerminate(t,x,u,Z0DSTRUCT_MEM,Z0DINPUT_MEM,metis_file,ACTIVE_OUTPUT,save_file,verbose);
    end
  
  case {1,2,4}
  	sys = [];
  otherwise
    error(['Unhandled flag = ',num2str(flag)]); % Error handling
end


if verbose == 1
	disp('end of simmetis :')
end

% End of simcronos.
%==============================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the
% S-function.
%==============================================================
%
function [sys,x0,str,ts,z0dinput,list_output] = mdlInitializeSizes(metis_file,list_output,verbose,input_order,input_option)

% sauve le model 
save_system(bdroot);


%chargement du fichier de reference
if ~isempty(metis_file)
	data = load(metis_file);
	if isfield(data,'post')
		z0dinput = data.post.z0dinput;
	elseif isfield(data,'z0dinput')
		z0dinput = data.z0dinput;
	else
		z0dinput = data;
	end
else
	z0dinput = evalin('base','z0dinput');
end
% suppression de la separatrice et Xdur
if isfield(z0dinput.exp0d,'Rsepa')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
end
if isfield(z0dinput.exp0d,'Zsepa')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
end
		
if isfield(z0dinput.exp0d,'XDURt')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURt');
end
if isfield(z0dinput.exp0d,'XDURx')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURx');
end
if isfield(z0dinput.exp0d,'XDURv')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURv');
end

% creation de la liste des consignes d'entrees
list_input = cat(1,fieldnames(z0dinput.cons),fieldnames(z0dinput.geo));
list_input(strmatch('temps',list_input,'exact')) = [];
list_input(strmatch('sp',list_input,'exact')) = [];
list_input(strmatch('vp',list_input,'exact')) = [];
list_input(strmatch('sext',list_input,'exact')) = [];

% test entrees de metis et de la sfunction
if length(input_order) ~= length(list_input)
	error('input simmetis and Metis missmatch');
end
for k=1:length(input_order)
	if isempty(strmatch(input_order{k},list_input,'exact'))
		error(sprintf('undefined input %s in Metis input',input_order{k}));	
	end
end

% parametre des options de Metis
if ~isempty(input_option)
	list_option = fieldnames(z0dinput.option);
        % extention au deuxieme injecteur
        list_option{end+1} = 'pnbi2';
        list_option{end+1} = 'ftnbi2';
        % boucel sur les options
	for k=1:length(input_option)
		if isempty(strmatch(input_option{k},list_option,'exact'))
			error(sprintf('undefined input %s in Metis input',input_option{k}));	
		end
	end
	input_order = cat(2,input_order,input_option);
end

% exportation dans le workspace
assignin('base','list_input_metis_sfunction',input_order);
z0dinput.list_input = input_order;

% creation de la liste des sorties scalaires
list_output_max = fieldnames(zero1t);
assignin('base','list_output_metis_sfunction_max',list_output_max);
z0dinput.list_output_max = list_output_max;

% creation de la liste des sorties profils
list_output_prof = fieldnames(z0dprofinfo);
assignin('base','list_profil_metis_sfunction_max',list_output_prof);
z0dinput.list_output_prof = list_output_prof;

% test de compatibilite + comptage
nboutput = 0;
if isempty(list_output)
	list_output{1} = 'wth';
end 
for k=1:length(list_output)
     outc = list_output{k}; 
     if any(outc == '@')
     	[tag,sub] = strtok(outc,'@');
	if isempty(strmatch(tag,list_output_prof,'exact'))
     		error(sprintf('unknown output profile %s',outc))
	else
		nbp = eval(sub(2:end));
		nboutput = nboutput + length(nbp);
	end
     else
     	if isempty(strmatch(outc,list_output_max,'exact'))
     		error(sprintf('unknown output signal %s',outc))
	else
		nboutput = nboutput + 1;
     	end
     end
end
% creation de la liste des sorties
assignin('base','list_output_metis_sfunction',list_output);
z0dinput.list_output = list_output;  
z0dinput.nboutput    = nboutput;    
z0dinput.nbinput     = length(input_order);    

% affichage de la configuration de Metis
disp('Metis configuration :')
disp(z0dinput.option)


if verbose == 1
	disp('Input :')
	list_input
	disp(' ')
	disp('Possible scalar output :')
	zero1t
	disp('Possible profile output :')
	z0dprofinfo
	disp(' ')
	disp('Sleected output :')
	list_output
end

% ceci est une initialisation pour simulink uniquement
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumInputs      = length(input_order);
sizes.NumOutputs     = nboutput;
sizes.DirFeedthrough = 1;  % flag=4 requires direct feedthrough
                           % if input u is involved in
                           % calculating the next sample time
                           % hit.
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
%
% Initialize the initial conditions.
%
x0 = []; % pas d'etat initial a exporter
%
% Set str to an empty matrix.
%
str = [];
%
% Initialize the array of sample times.
%
ts = [-1 0];      % variable sample time
% End of mdlInitializeSizes.

%
%==============================================================
% mdlOutputs
% Return the block outputs.
%==============================================================
%
function   [sys,z0dstruct,zs,profli] = mdlOutputs(time,x,u,z0dstruct,z0dinput,list_output,verbose,input_order)     

% ici on recupere les donnees en provenance de simulink
% this is an exemple of input from simulink
% ici on recupere les donnees en provenance de simulink
% this is an exemple of input from simulink
% separation des champs
geo_nom = fieldnames(z0dinput.geo);
cons_nom = fieldnames(z0dinput.cons);
for k=1:length(z0dinput.list_input)
	nomc = z0dinput.list_input{k};
	if isempty(z0dstruct)
		% dans ce cas les option ne sont pas a  modifier
		if ~isempty(strmatch(nomc,geo_nom,'exact'))
			geo1t.(nomc) = interp1(z0dinput.cons.temps,z0dinput.geo.(nomc),time,'linear',z0dinput.geo.(nomc)(1));
		elseif ~isempty(strmatch(nomc,cons_nom,'exact'))
			cons1t.(nomc) = interp1(z0dinput.cons.temps,z0dinput.cons.(nomc),time,'linear',z0dinput.cons.(nomc)(1));	
		end	
		
	elseif isfinite(u(k))
		if ~isempty(strmatch(nomc,geo_nom,'exact'))
			geo1t.(nomc) = u(k);
		elseif ~isempty(strmatch(nomc,cons_nom,'exact'))
			cons1t.(nomc) = u(k);	
                elseif ~isempty(strmatch(nomc(1:end-1),cons_nom,'exact')) && (nomc(end) == '2')
                       % extention au deuxieme injecteur
		       cons1t.(nomc(1:end-1)) = cons1t.(nomc(1:end-1)) + sqrt(-1) .* u(k);	
		else
		       z0dinput.option.(nomc) = u(k);
		end
%  	elseif isempty(z0dstruct)
%  		% dans ce cas les option ne sont pas a  modifier
%  		if ~isempty(strmatch(nomc,geo_nom,'exact'))
%  			geo1t.(nomc) = interp1(z0dinput.cons.temps,z0dinput.geo.(nomc),time,'linear',z0dinput.geo.(nomc)(1));
%  		elseif ~isempty(strmatch(nomc,cons_nom,'exact'))
%  			cons1t.(nomc) = interp1(z0dinput.cons.temps,z0dinput.cons.(nomc),time,'linear',z0dinput.cons.(nomc)(1));	
%  		end	
	else
		% dans ce cas les option ne sont pas a modifier
		if ~isempty(strmatch(nomc,geo_nom,'exact'))
			geo1t.(nomc) = interp1(z0dinput.cons.temps,z0dinput.geo.(nomc),time,'linear',z0dinput.geo.(nomc)(end));
		elseif ~isempty(strmatch(nomc,cons_nom,'exact'))
			cons1t.(nomc) = interp1(z0dinput.cons.temps,z0dinput.cons.(nomc),time,'linear',z0dinput.cons.(nomc)(end));	
		end
	end
end
cons1t.temps = time;

if verbose == 1
 	cons1t
 	geo1t	
end

% recherche de la connexion avec un equilibre
try
	sepa1t = evalin('base','SEPARATRIX4METIS');		
	%appel de zerod evolution (sans donnees experimental et avec separatrice pour le moment)
	[zs,profli,z0dstruct] = zerodevolution(z0dstruct,z0dinput.option,time,cons1t,geo1t,[],sepa1t);
        assignin('base','zerod_simmetis',zs);
        assignin('base','profil0d_simmetis',profli);
        assignin('base','cons_simmetis',cons1t);
        assignin('base','geo_simmetis',geo1t);
        assignin('base','sepa_simmetis',sepa1t);
catch
	%appel de zerod evolution (sans donnees experimental et sans separtarice pour le moment)
	[zs,profli,z0dstruct] = zerodevolution(z0dstruct,z0dinput.option,time,cons1t,geo1t);
        assignin('base','zerod_simmetis',zs);
        assignin('base','profil0d_simmetis',profli);
        assignin('base','cons_simmetis',cons1t);
        assignin('base','geo_simmetis',geo1t);
end

if verbose == 1
 	zs
end

% ici les donnees en provenance de cronos doivent etre retournees a simulink
sys = NaN .* ones(1,z0dinput.nboutput);
nbc = 1;
for k=1:length(list_output)
     outc = list_output{k};
     if any(outc == '@')
     	[tag,sub] = strtok(outc,'@');
	pc        = profli.(tag);
	indc      = eval(sub(2:end));
	sys(1,nbc + (1:length(indc)) - 1) = pc(indc);
	nbc = nbc + length(indc);
     else
	sys(1,nbc) = zs.(list_output{k});
	nbc = nbc + 1;
     end
end
% end mdlOutputs


% si necessaire
function sys = mdlDerivatives(t,x,u,z0dstruct,z0dinput,list_output,verbose,input_order)

% reservation memoire
sys = NaN .* ones(1,length(list_output));
% boucle sur les sorties
nbc = 1;
for k=1:length(list_output)
     outc = list_output{k};
     if any(outc == '@')
     	[tag,sub] = strtok(outc,'@');
	pc        = profli.(tag);
	indc      = eval(sub(2:end));
	sys(1,nbc+(1:length(indc)) - 1) =  diff(pc((end - 1):end,indc),1,1) ./  ...
	                               (diff(z0dstruct.profil.temps((end - 1):end)) * ones(1,length(indc)));
	nbc = nbc + length(indc);
     else
	sys(1,nbc) = diff(z0dstruct.zs.(list_output{k})((end - 1):end)) ./ diff(z0dstruct.zs.temps((end - 1):end));
	nbc = nbc + 1;
     end
end

if verbose == 1
	disp('dUdt =')
	sys
end

% End of mdlGetTimeOfNextVarHit.


function  sys = mdlTerminate(t,x,u,z0dstruct,z0dinput,metis_file,list_output,save_file,verbose)

% sauvegarde dans le workspace pour utiliser les plot de metis
try
	post = evalin('base','post');
catch
	post = [];
end
% 
% ajout du nom du system et du system
%
z0dstruct.z0dinput.system.name     = strtok(gcs,filesep);
z0dstruct.z0dinput.system.fullname = which(z0dstruct.z0dinput.system.name);
% sauvvegarde temporaire du systeme
%[s,z0dstruct.z0dinput.system.mdl]  = unix(sprintf('cat %s',z0dstruct.z0dinput.system.fullname));
fid = fopen(z0dstruct.z0dinput.system.fullname,'r');
if fid > 0
    if ~isempty(strfind(z0dstruct.z0dinput.system.fullname,'.slx'))
        z0dstruct.z0dinput.system.mdl = fread(fid,Inf,'uint8')';        
    else
        z0dstruct.z0dinput.system.mdl = char(fread(fid,Inf,'char')');
    end
	fclose(fid);
else
	mdl_loc ='unable to read system text file';
end			

% donnees d'entree du modele
z0dstruct.z0dinput.system.z0dinput = z0dinput;
%
post.zerod            = z0dstruct.zs;
post.profil0d         = z0dstruct.profil;
post.z0dinput         = z0dstruct.z0dinput;

assignin('base','post',post);

% sauvegarde du resultat 
if ~isempty(save_file)
	save(save_file,'z0dstruct','z0dinput','metis_file','list_output');
	if verbose  == 1
		fprintf('save Metis simulation results in %s \n',save_file);
	end
end

% on ne retourne rien vers simulink
sys = [];

% End of mdlTerminate

