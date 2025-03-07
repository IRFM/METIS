function option_feeqs_jt60sa = compute_jt60sa_magnetisation(option_feeqs_jt60sa)

% test if FEEQS is in the path
set_feeqs_path_jt60sa;


% GUI generation
if nargin < 1

    valeur.ics2ka     = 19.07;   % CS2 initila current 
    type.ics2ka       = 'float';
    borne.ics2ka      = [0,20];  
    defaut.ics2ka     = 19.07;
    info.ics2ka       = 'Magnetisation: initial current in CS2 coil (kA)';
    
    valeur.mesh_type = 'standard';     
    type.mesh_type   = 'string';
    borne.mesh_type  = {'test','standard','fine'};  
    defaut.mesh_type = 'standard';
    info.mesh_type   = 'Number of elements in the 2D mesh (test is reserved for testing and must be not used for production)';    
    
    valeur.output_file_name 	    = '';
    type.output_file_name           = 'string';                    
    borne.output_file_name          = '';  
    defaut.output_file_name         = '';                
    info.output_file_name           = 'matfile name where results are stored (if empty, no file will be created)';
    
    valeur.plotonoff    = 2;      
    type.plotonoffs     = 'integer';
    borne.plotonoff     = {0,1,2};  
    defaut.plotonoff    = 2;
    info.plotonoff      = 'Set level of displayed graphs:  0 = no graphs; 1 = results only; 2 = all graphs (debug)';
    
    valeur.regularisation    = 1e-14;   
    type.regularisation      = 'float';
    borne.regularisation     = [1e-8,1e-20];  
    defaut.regularisation    = 1e-14;
    info.regularisation      = 'weight of regularisation term (sum(I^2))';

    valeur.R_ini    = 2.8;   
    type.R_ini      = 'float';
    borne.R_ini     = [1.44,4.2];  
    defaut.R_ini    = 2.8;
    info.R_ini      = 'Central postition of null zone (m)';

    valeur.a_ini    = 1.1;   
    type.a_ini      = 'float';
    borne.a_ini     = [0.2,2.4];  
    defaut.a_ini    = 1.1;
    info.a_ini      = 'radius of null zone (m)';


    interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
    interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

    option_feeqs_jt60sa.valeur     = valeur;
    option_feeqs_jt60sa.type       = type;
    option_feeqs_jt60sa.borne      = borne;
    option_feeqs_jt60sa.defaut     = defaut;
    option_feeqs_jt60sa.info       = info;
    option_feeqs_jt60sa.interface  = interface;
    
    option_feeqs_jt60sa.description = 'Computation of pre magnetisation before breakdown for JT-60SA';   % description (une ligne) de la fonction

    option_feeqs_jt60sa.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
    option_feeqs_jt60sa.gui      ='';                             % nom de l'interface graphique specifique si elle existe
    option_feeqs_jt60sa.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

    % end of GUI form declaration
    return

end

% call FEESQ.M
dirmem = pwd;
% make a local copy of the project
projectpath = fileparts(which('compute_breakdown_current'));
[~,projectname]  = fileparts(projectpath);
% temprary file 
if isdir('/dev/shm')
  templocaldir = fullfile('/','dev','shm',sprintf('%s_%s',fileparts(tempname),getenv('USER')));
else
  templocaldir = tempname;
end
if ~isdir(templocaldir)
  mkdir(templocaldir);
end
unix(sprintf('cp -rpf %s %s',projectpath,templocaldir));
try
  cd(fullfile(templocaldir,projectname));
  [coil_currents,CurrentsInt,psi_ini_out,B_null,cost,regul] = compute_breakdown_current(option_feeqs_jt60sa.ics2ka,option_feeqs_jt60sa.mesh_type,option_feeqs_jt60sa.plotonoff, ...
                           option_feeqs_jt60sa.regularisation,option_feeqs_jt60sa.R_ini,option_feeqs_jt60sa.a_ini);
  cd(dirmem);
catch
  cd(dirmem);
  error(lasterror);
end
% remove local temporary project directory
unix(sprintf('rm  -rf %s',fullfile(templocaldir,projectname)));
% save
if ~isempty(option_feeqs_jt60sa.output_file_name)    
    save(option_feeqs_jt60sa.output_file_name,'option_feeqs_jt60sa','coil_currents','CurrentsInt','psi_ini_out','B_null','cost','regul');
end
% output
premagnetisation_jt60sa.option_feeqs_jt60sa = option_feeqs_jt60sa;
premagnetisation_jt60sa.coil_currents = coil_currents;
premagnetisation_jt60sa.CurrentsInt   = CurrentsInt;
premagnetisation_jt60sa.psi_ini_out   = psi_ini_out;
premagnetisation_jt60sa.B_null        = B_null;
premagnetisation_jt60sa.cost          = cost;
premagnetisation_jt60sa.regul         = regul;
option_feeqs_jt60sa = premagnetisation_jt60sa;
