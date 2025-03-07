% fonction assistant d'appel du 0d cronos
function zcall0d(action)

%#function zcall0d_fct
% initilisation in expert mode if not defined
if isappdata(0,'GUI_GLOBAL_BASIC_MODE')
    if getappdata(0,'GUI_GLOBAL_BASIC_MODE') == 0
        z0dbasic_mode_switcher(2);
        defmode_init = 2;
    else
        z0dbasic_mode_switcher(1);
        defmode_init = 1;
    end
elseif isdeployed
    defmode_init = 1;
elseif strcmp(upper(getenv('METIS_MODE')),'EXPERT')
    defmode_init = 2;
elseif strcmp(upper(getenv('METIS_MODE')),'STANDARD')
    defmode_init = 1;
else
    user = getenv('USER');
    if length(user) > 2
        if isempty(str2num(user(1:2))) && ~isempty(str2num(user(3:end)))
            defmode_init = 2;
        else
            defmode_init = 1;
        end
    else
        defmode_init = 1;
    end
end

if nargin == 0
	action = 'init';
elseif isempty(action);
	action = 'init';
end	
% choix de la langue
%langue =  lower(getappdata(0,'langue_cronos'));
langue =  'anglais';
setappdata(0,'langue_cronos',langue);

% initilisation selon presence ou nom de donnees cronos
if strcmp(action,'init')
    
    if ~isappdata(0,'GUI_GLOBAL_BASIC_MODE')
        switch defmode_init
            case 1
                z0dbasic_mode_switcher(1);
            otherwise
                z0dbasic_mode_switcher(2);
        end
    end
    
    if  ~isempty(dir('reference_scenario_jt60-sa'))
        mode_jt60_sa = 2;
        txt_jt60_sa = 'JT-60SA';
    elseif isdir(fullfile(fileparts(which('jt60sawall.mat')),'reference_scenario_jt60-sa'))
        mode_jt60_sa = 1;
        txt_jt60_sa = 'JT-60SA';
    else
        mode_jt60_sa = 0;
        txt_jt60_sa = 'JT-60SA scenario generator';
    end
    
    % detection de Cronos
    itemliste = {'Tore Supra','WEST','JET','DIIID','ASDEX-U','COMPASS','EAST',txt_jt60_sa,'TCV','Reactors (ITER & DEMO)','Other','load'};
    %yc = evalin('base','exist(''param'',''var'') & exist(''data'',''var'')');
    %if yc == 1
    itemliste{end+1} = 'Cronos';
    cronos_num = length(itemliste);
    %else
    %    cronos_num = 999991;
    %end
    
    
    % detection UAL
    if ~isempty(getenv('HDF5_BASE')) | ~isempty(getenv('MDSPLUS_TREE_BASE_0')) | ~isempty(getenv('UAL'))
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
                    % use :  any(getappdata(0,'UALVERSION') < '4.08b')
                    
                end
            end
        end
        %  addpath('/afs/efda-itm.eu/project/switm/');
        %  addpath('/afs/efda-itm.eu/project/switm/ual/4.07b/matlabinterface');
        if isempty(import) && ~isappdata(0,'JAVAINTERFACESET') && isempty(getappdata(0,'JAVAINTERFACESET'))
            import ualmemory.javainterface.*;
        end
        setappdata(0,'JAVAINTERFACESET','done')
    end
    % Detection of IMAS presence
    if exist('imas_open_env')
        setappdata(0,'IMAS_EXIST','YES');
        itemliste{end+1} = 'IMAS';
        imas_num = length(itemliste);
    else
        if isappdata(0,'IMAS_EXIST')
            rmappdata(0,'IMAS_EXIST');
        end
        imas_num = 999993;
    end
    
    if isappdata(0,'UALVERSION') && ~isempty(getappdata(0,'UALVERSION'))
        itemliste{end+1} = 'UAL';
        ual_num = length(itemliste);
    else
        ual_num = 999992;
    end
    
    itemliste{end+1} = 'System code';
    system_code_num = length(itemliste);
    
    itemliste{end+1} = 'ST40';
    st40_code_num = length(itemliste);
    
    
    itemliste{end+1} = 'Cancel';
    
    switch langue
        case 'francais'
            header = 'Quel source de donnees ?';
        otherwise
            header = 'Which data source ?';
    end
    
    Button = menu(header,itemliste);
    
    if Button == 0
        return
    end
    
    switch Button
        case length(itemliste)
            disp('initialisation cancelled')
        case imas_num
            disp('Reading IMAS UAL data')
            evalin('base','z0dinput = zerod_init(118);');
        case ual_num
            disp('reading UAL data')
            evalin('base','z0dinput = zerod_init(117);');
        case cronos_num
            yc = evalin('base','exist(''param'',''var'') & exist(''data'',''var'')');
            if yc ~= 1
                evalin('base','zuiload');
            end
            disp('reading CRONOS data')
            evalin('base','z0dinput = zerod_init(0);');
        case system_code_num
            disp('Reading System code output')
            evalin('base','z0dinput = zerod_init(119);');
        case st40_code_num
            disp('Reading data for ST40 Tokamak')
            evalin('base','z0dinput = zerod_init(201);');       
        case 1
            disp('access to Tore Supra data')
            evalin('base','z0dinput = zerod_init(1);');
        case 2
            disp('Preparation of WEST simulation')
            %evalin('base','z0dinput = zerod_init(11);');
            rep = menu('WEST simulation: what kind of simulation ?','Scenario generator','Reference simulation','Before shot','After shot','Interpretative');
            switch rep
                case 1
                    if preset_west_scenario_choice ~= 0
                        evalin('base','zwaitfor(zuicreefunform(''westmetissimulation'',''west_param'',1,0,''z0dinput=westmetissimulation(west_param);''),''visible'');');
                    end
                case 2
                    pwd_mem =pwd;
                    if isdir(fullfile(fileparts(which('west_wall.mat')),'reference_scenario_WEST'))
                        cd(fullfile(fileparts(which('west_wall.mat')),'reference_scenario_WEST'));
                    elseif ~isempty(dir('reference_scenario_WEST'))
                        try
                            cd('reference_scenario_WEST');
                            pwd_mem	= '..';
                        catch
                            pwd;
                        end
                    end
                    evalin('base','z0dinput = zerod_init(-2);');
                    drawnow;
                    zcall0d_fct('radio_load');
                    cd(pwd_mem);
                case 3
                    evalin('base','z0dinput = zerod_init(13);');
                    
                case 4
                    evalin('base','z0dinput = zerod_init(12);');
                    
                otherwise
                    %evalin('base','[z0dinput,post,cp,equi] = autoWESTmetis;');
                    evalin('base','zwaitfor(zuicreefunform(''autoWESTmetis2'',''autoWESTmetis_option'',1,0,''[z0dinput,post,cp,equi] = autoWESTmetis2(autoWESTmetis_option);''),''visible'');');

            end
        case 3
            disp('access to JET data')
            evalin('base','z0dinput = zerod_init(2);');
            %
            % Read TRANSP data if needed
            if isappdata(0,'TRANSP_RUN') && isfinite(getappdata(0,'TRANSP_RUN'))
                try
                    evalin('base','transp4metis(z0dinput.shot,getappdata(0,''TRANSP_RUN''),getappdata(0,''TRANSP_USER''));');
                catch
                   disp('Error accessing TRANSP data');
                   disp(lasterr);
                end
                
            end
        case 4
            disp('access to DIII-D data')
            evalin('base','z0dinput = zerod_init(3);');
        case 5
            disp('access to ASDEX-U data')
            evalin('base','z0dinput = zerod_init(4);');
        case 6
            disp('access to COMPASS data')
            evalin('base','z0dinput = zerod_init(5);');
        case 7
            disp('access to EAST data')
            evalin('base','z0dinput = zerod_init(6);');
        case 8
            switch mode_jt60_sa
                case {1,2}
                    rep = menu('Kind of simulation ?','JT-60SA scenario generator','Reference simulation files', ...
                        'Pre magnetisation computation before breakdown','Fast estimation of current in coils', ...
                        'Free boundary equilibrium in inverse static mode','Free boundary equilibrium in inverse evolutive mode', ...
                        'Plot results of free boundary equilibrium in inverse mode');
                    switch rep
                        case 1
                            if preset_jt60sa_scenario_choice ~= 0
                                evalin('base','zwaitfor(zuicreefunform(''jt60sametissimulation'',''jt60sa_param'',1,0,''z0dinput=jt60sametissimulation(jt60sa_param);''),''visible'');');
                            end
                        case 2
                            if mode_jt60_sa == 1
                                pwd_mem =pwd;
                                cd(fullfile(fileparts(which('jt60sawall.mat')),'reference_scenario_jt60-sa'));
                            elseif ~isempty(dir('reference_scenario_jt60-sa'))
                                try
                                    cd('reference_scenario_jt60-sa');
                                catch
                                    pwd;
                                end
                            end
                            evalin('base','z0dinput = zerod_init(-2);');
                            drawnow;
                            zcall0d_fct('radio_load');
                            if mode_jt60_sa == 1
                                cd(pwd_mem);
                            else
                                try
                                    cd('..');
                                catch
                                    pwd
                                end
                            end
                        case 3
                            try
                                evalin('base','zwaitfor(zuicreefunform(''compute_jt60sa_magnetisation'',''option_feeqs_jt60sa'',1,0,''premagnetisation_jt60sa = compute_jt60sa_magnetisation(option_feeqs_jt60sa);''),''visible'');')
                            catch
                                if isempty(which('start_up_Unix.m'))
                                    warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
                                else
                                    errordlg(lasterr, 'Error calling FEEQS.M');
                                end
                            end
                        case 4
                            try
                                evalin('base','zwaitfor(zuicreefunform(''compute_jt60sa_fbe_inverse_fast'',''option_feeqs_jt60sa_fast'',1,0,''option_feeqs_jt60sa_fast = compute_jt60sa_fbe_inverse_fast(option_feeqs_jt60sa_fast);''),''visible'');')
                            catch
                                if isempty(which('start_up_Unix.m'))
                                    warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
                                else
                                    errordlg(lasterr, 'Error calling FEEQS.M');
                                end
                            end
                        case 5
                            try
                                evalin('base','zwaitfor(zuicreefunform(''compute_jt60sa_fbe_inverse'',''option_feeqs_jt60sa'',1,0,''[option_feeqs_jt60sa,fbe_lcfs] = compute_jt60sa_fbe_inverse(option_feeqs_jt60sa);''),''visible'');')
                            catch
                                if isempty(which('start_up_Unix.m'))
                                    warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
                                else
                                    errordlg(lasterr, 'Error calling FEEQS.M');
                                end
                            end
                            
                         case 6
                            try
                                evalin('base','zwaitfor(zuicreefunform(''compute_jt60sa_fbe_inverse_evol'',''option_feeqs_jt60sa'',1,0,''[option_feeqs_jt60sa,fbe_lcfs] = compute_jt60sa_fbe_inverse_evol(option_feeqs_jt60sa);''),''visible'');')
                            catch
                                if isempty(which('start_up_Unix.m'))
                                    warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
                                else
                                    errordlg(lasterr, 'Error calling FEEQS.M');
                                end
                            end
                            
                       case 7
                            try
                                evalin('base','set_feeqs_path_jt60sa;plotfigure_inverse4metis;')
                            catch
                                if isempty(which('start_up_Unix.m'))
                                    warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
                                else
                                    errordlg(lasterr, 'Error calling FEEQS.M');
                                end
                            end
                            
                        otherwise
                            evalin('base','z0dinput = zerod_init(-2);');
                            
                    end
                otherwise
                    evalin('base','zwaitfor(zuicreefunform(''jt60sametissimulation'',''jt60sa_param'',1,0,''z0dinput=jt60sametissimulation(jt60sa_param);''),''visible'');');
            end
        case 9
            disp('access to TCV data')
            evalin('base','z0dinput = zerod_init(7);');
        case 10
            rep = menu('Kind of simulation ?','Reactor scenario generator','Fast esitimation of current in coils', ...
                         'Plot results of free boundary equilibrium in inverse mode','Cancel');

            switch rep
               case 1
                    evalin('base','zwaitfor(zuicreefunform(''reactormetissimulation'',''reactor_option'',1,0,''z0dinput=reactormetissimulation(reactor_option);''),''visible'');');
               case 2
                          try
                                evalin('base','set_reactor_feeqs_path;zwaitfor(zuicreefunform(''compute_reactor_fbe_inverse_fast'',''option_feeqs_reactor'',1,0,''option_feeqs_reactor = compute_reactor_fbe_inverse_fast(option_feeqs_reactor);''),''visible'');')
                            catch
                                if isempty(which('start_up_Unix.m'))
                                    warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
                                else
                                    errordlg(lasterr, 'Error calling FEEQS.M');
                                end
                            end
               
               case 3
                            try
                                evalin('base','set_reactor_feeqs_path;plotfigure_reactor_inverse4metis;')
                            catch
                                if isempty(which('start_up_Unix.m'))
                                    warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
                                else
                                    errordlg(lasterr, 'Error calling FEEQS.M');
                                end
                            end
               
               otherwise
                   disp('action cancelled');
               end
        case 11
            evalin('base','z0dinput = zerod_init(-1);');
        otherwise
            evalin('base','z0dinput = zerod_init(-2);');
            drawnow;
            zcall0d_fct('radio_load');
    end
    if evalin('base','exist(''z0dinput'',''var'')') && evalin('base','isempty(z0dinput)')
        [hform,hui] = zuiformhandle('zeroda');
        if ishandle(hform)
            zuiformvisible(hform);
        else
            evalin('base','metis');
        end
        return
    end
    
elseif ~strcmp(action,'deja')
    evalin('base','z0dinput = zerod_init;');
end

% creation de l'interface    
% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('zeroda');
if ishandle(hform)
        zuiformvisible(hform);
        z0dbasic_mode_switcher(get(findobj(zuiformhandle('zeroda'),'tag','popup_advance'),'value'));
	return
end


% 1eres lignes : Nom du fichier de travail
% --------------------------------------
form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'titre','text@full',' Metis ',50,'',[]};
form{length(form)+1} = {col1} ;

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% alignement jump
% ---------------
colj = {'void','jump','void',[],''} ;

% 1ere ligne
% ----------
col1 = {'radio_param'   ,'radio@left' ,'Parameters' ,0,'edit Metis parameters'};
col2 = {'radio_load'   ,'radio' ,'Load' ,0,'load Metis data and parameters'};
col3 = {'radio_save'   ,'radio' ,'Save' ,0,'save Metis data and parameters'};
col3bis = {'radio_export'   ,'radio' ,'Export' ,0,'Export METIS data to IMAS data model stored in .mat file'};
col4 = {'radio_ref'   ,'radio' ,'Create reference' ,0,'Create a reference data set from Metis results'};
col5 = {'radio_audit'   ,'radio' ,'Compare' ,0,'Compare results with reference'};
col6 = {'radio_pdf'   ,'radio' ,'PDF ouput' ,0,'make a PDF of all figures of Metis simulator'};
col7 = {'radio_closeall'   ,'radio' ,'Close' ,0,'Close all figures'};
col8 = {'popup_advance'   ,'popup@right' ,'Standard mode|Expert mode' ,defmode_init,'Select the mode of the interface: if = Standard mode, only main features are available; if = Expert mode, all features are available'};
form{length(form)+1} = {col1,colj,col2,col3,col3bis,colj,col4,col5,colj,col6,colj,colj,col7,colj,col8};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
col1 = {'titre2','text@full',' Waveforms & data edition ',50,'',[]};
form{length(form)+1} = {col1} ;
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% 2eme ligne
% ----------
col1 = {'radio_ip'   ,'radio' ,'Ip   ' ,0,'edit plasma current waveform Ip (A)'};
col2 = {'radio_nbar'     ,'radio' ,'Nbar     ' ,0,'edit Nbar waveform (m^-3)'};
col3 = {'radio_zeff'   ,'radio' ,'Zeff   ' ,0,'edit average Zeff waveform'};
col4 = {'radio_xece'     ,'radio' ,'Xecrh  ' ,0,'edit ECRH radial deposition position waveform'};
col5 = {'radio_b0'     ,'radio' ,'B0    ' ,0,'edit toroidal field waveform B0 (@ R0, T)'};
col6 = {'radio_iso'     ,'radio' ,'nT/nD    ' ,0,sprintf('edit waveform for ratio for:\n    - n_T / n_D, if gas = 3;\n    - n_He3 / n_D, if gas = 5;\n    - n_B11 / n_H, if gas = 11')};
col7 = {'radio_flux'     ,'radio' ,'Flux    ' ,0,'edit waveform for poloidal edge flux'};
col8 = {'radio_make_flux'     ,'radio' ,'Create flux' ,0,'Create poloidal flux waveform from output of previous run'};
if exist('metis2luke')
    col9 = {'radio_luke'     ,'radio','External Luke' ,0, ...
	    'open menu to select external data take from LUKE data set to be used in METIS and that will not be computer in METIS' };
    form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8,col9};
else
    col9 = {'radio_void_8'     ,'radio' ,'' ,0,''};
    form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8};
end

% 3eme ligne
% ----------
col1 = {'radio_ecrh'   ,'radio' ,'ECRH   ' ,0,'edit ECRH power waveform (W)'};
col2 = {'radio_icrh'    ,'radio' ,'ICRH    ' ,0,'edit ICRH power waveform (W)'};
col3 = {'radio_lh'   ,'radio' ,'LH   ' ,0,'edit LH  power waveform (W)'};
col4 = {'radio_nbi' ,'radio' ,'NBI ' ,0,'edit NBI power waveform (W)'};
col5 = {'radio_hmore' ,'radio' ,'H factor' ,0,'edit time confinement multiplication factor (for any confiment mode Ohmic,L,H ...)'};
col6 = {'radio_ftnbi' ,'radio' ,'FT_NBI ' ,0,sprintf('%s\n%s\n%s\n%s','edit NBI tritium power fraction waveform: P_NBI_T = FT_NBI * P_NBI and P_NBI_D = (1 - FT_NBI) * P_NBI (for DT plasma) and:', ...
        'if gas  = 1 or 4,  edit NBI hydrogen power fraction waveform: P_NBI_H = FT_NBI * P_NBI and P_NBI_D = (1 - FT_NBI) * P_NBI', ...
        'if gas  = 5, edit NBI He3 power fraction waveform: P_NBI_He3 = FT_NBI * P_NBI and P_NBI_D = (1 - FT_NBI) * P_NBI', ...
        'if gas  = 11, edit NBI boron power fraction waveform: P_NBI_B11 = FT_NBI * P_NBI and P_NBI_H = (1 - FT_NBI) * P_NBI')}; 
col7 = {'radio_temps'     ,'radio' ,'Time edition' ,0,'Edition of the time vector of METIS'};
%if isdeployed
%    col8 = {'radio_free2'     ,'radio' ,'' ,0,''};
%else
    col8 = {'radio_remove'     ,'radio' ,'Clear external' ,0,'Removes external data used in METIS'};
%end
col9 = {'radio_void_9'     ,'radio' ,'' ,0,''};
if exist('metis2luke')
  form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8,col9};
else
  form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8};
end
% 4eme ligne
% ----------
col1 = {'radio_r0'       ,'radio' ,'R0   ' ,0,'edit plasma major radius waveform (m)'};
col2 = {'radio_z0'       ,'radio' ,'z0   ' ,0,'edit plasma vertical position (m)'};
col3 = {'radio_a'     ,'radio' ,'a   ' ,0,'edit plasma minor radius waveform  (m)'};
col4 = {'radio_k'     ,'radio' ,'K    ' ,0,'edit plasma elongation waveform  (b/a @ LCFS)'};
col5 = {'radio_d'      ,'radio' ,'d  ' ,0,'edit plasma up/down averaged triangularity waveform (@ LCFS)'};
col6 = {'radio_sepa'   ,'radio' ,'Separatrix' ,0,'compute complete LCMS'};
col7 = {'radio_noise'     ,'radio' ,'de-noising' ,0,'de-noising METIS waveforms'};
%if isdeployed
    col8 = {'radio_free3'     ,'radio' ,'' ,0,''};
%else
  col8 = {'radio_external'     ,'radio' ,'External data' ,0, ...
	  'open menu to select external data take from CRONOS data set to be used in METIS and that will not be computer in METIS'};
%end
col9 = {'radio_void_10'     ,'radio' ,'' ,0,''};
if exist('metis2luke')
  form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8,col9};
else
  form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8};
end
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
col1 = {'titre4','text@full',' Command ',50,'',[]};
form{length(form)+1} = {col1} ;
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'radio_go'   ,'radio' ,'Run METIS' ,0,'run of the Metis simulator'};
col2 = {'radio_fast'   ,'radio' ,'Run METIS in fast mode' ,0,'run of Metis simulator, fast mode for scenario preparation'};
col3 = {'radio_workingpoint'   ,'radio' ,'Operation point' ,0,' Compute operation point for one time slice (steady state, all d/dt = 0, only dPsi/dt #0, design for reactor studies)'};
col4 = {'radio_hyb'   ,'radio' ,'Fit of LH efficiency & Wdia' ,0,'search of the best value of LH efficiency by fitting the poloidal flux consummation and used measured plasma energy contents'};
col5 = {'radio_evolution'   ,'radio' ,'Evolution' ,0,'Compute Metis simulation using evolution mode'};
col6 = {'radio_evolution_restart'   ,'radio' ,'Restart' ,0,'Restart Metis simulation using evolution mode for interrupted simulation'};
form{length(form)+1} = {col1,colj,col2,colj,col3,colj,col4,colj,col5,col6};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
col1 = {'titre3','text@full','Visualisation ',50,'',[]};
form{length(form)+1} = {col1} ;
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% 5 ieme lign = boutton important
col1 = {'radio_scenar'        ,'radio' ,'Overview' ,0,'Overview of the shot'};
col2 = {'radio_digest'   ,'radio' ,'simulation summary' ,0,'Print summary parameter table and plot some graphs'};
col3 = {'radio_prof'        ,'radio' ,'Profiles' ,0,'plot main profiles'};
col4 = {'radio_2deq'        ,'radio' ,'2D equi.' ,0,'2D Metis equilibrium'};
col5 = {'radio_dataplot0'        ,'radio' ,'Data browser' ,0,'open complete METIS data browser'};
col6 = {'fig2pub'        ,'radio@right' ,'Fig2Pub' ,0,'tranform figure in set of figures ready to be used un document (select figure in menu)'};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6};


% ligne vide pour separer
sepa ={'separation_comm','frame','',-10,''};
form{length(form)+1} = {sepa};

% 6eme ligne
% z0plotconv.m        
% ----------
% col1 = {'radio_sc'       ,'radio' ,'over view' ,0,'plot 0D main results'};
col1 = {'radio_p'     ,'radio' ,'power' ,0,'Plot of 0D powers versus time'};
col2 = {'radio_e'        ,'radio' ,'energy' ,0,'Plot of 0D energy quantities versus time'};
col3 = {'radio_c'        ,'radio' ,'confinement' ,0,'Plot of 0D confinment quantities versus time'};
col4 = {'radio_t'   ,'radio' ,'temperature' ,0,'Plot of 0D temperatures versus time'};
col5 = {'radio_n'    ,'radio' ,'density' ,0,'Plot of 0D densities versus time'};
col6 = {'radio_lhacc'          ,'radio' ,'LH wave' ,0,'Lower Hybrid propagation and absorption diagrams'};
col7 = {'radio_ddsts'          ,'radio' ,'Sawtooth' ,0,'Sawtooth period, NTM triggering and control'};
%col8 = {'radio_qlkANN_k'     ,'radio' ,'QlkANNk' ,0,'Plot flux from Qualikiz neural network with kinetic electrons'};
%col9 = {'radio_qlkz_std_wp'   ,'radio' ,'OP + QLKZ' ,0,'Compute operating point for one time slice using kinetic profiles predicted with Qualikiz \n(steady state, all d/dt = 0, only dPsi/dt #0,\ncreates external data for Te, Ti and Ne that can be reused;\n external data must be cleared before next call)'};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8,col9};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,colj};

% 6eme ligne
% ----------
col1 = {'radio_j'     ,'radio' ,'current' ,0,'Plot of 0D currents versus time'};
col2 = {'radio_eq'         ,'radio' ,'equilibrium' ,0,'Plot of 0D equilibrium quantities versus time'};
col3 = {'radio_lhe'          ,'radio' ,'LH efficiency ' ,0,'Plot of 0D LH efficiency'};
col4 = {'radio_geo'        ,'radio' ,'geometry' ,0,'Plot of geometrical quantities versus time'};
col5 = {'radio_star'        ,'radio' ,'nu* & rho*' ,0,'plot nu star, rho star and other quantities useful for transport, stability analysis or scenario tuning'};
%col6 = {'radio_nbijet'        ,'radio' ,'NBI JET' ,0,'comparison to TRANSP and PENCIl data  only for JET shot'};
col6 = {'radio_nbijet'        ,'radio' ,'JET graphs' ,0,'Graphs dedicated to JET tokamak'};
col7 = {'radio_conv'        ,'radio' ,'convergence' ,0,'Plot convergence indicator of the 0D'};
%col8 = {'radio_qlkANN_k_wp'   ,'radio' ,'OP + QlkANNk' ,0,'Compute operating point for one time slice using kinetic profiles predicted with Qualikiz neural network\n(steady state, all d/dt = 0, only dPsi/dt #0,\ncreates external data for Te, Ti and Ne that can be reused;\n external data must be cleared before next call)'};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,colj};

% 7eme ligne
% ----------
col1 = {'radio_ts'     ,'radio' ,'Neutrons' ,0,'plot neutron rates versus time for DD and DT reactions'};
col2 = {'radio_er'         ,'radio' ,'Er' ,0,'coarse estimation of radial electric field'};
col3 = {'radio_plotrad'     ,'radio' ,'Radiation' ,0,'radiative power summary'};
col4 = {'radio_density'    ,'radio' ,'Ne & Te exp.' ,0,'Plot density and temperature profiles computed by METIS and the experimental data for JET and Tore Supra'};
col5 = {'radio_gaz'        ,'radio' ,'Gas balance ' ,0,'plot gas balance'};
col6 = {'radio_breakdown'   ,'radio' ,'Breakdown' ,0,'Graph dedicated to breakdown and burn through visualisation'};
%col7 = {'radio_export_cpos'    ,'radio' ,'CPOS' ,0,'re-computes CPOS and 2D equilibrium data then exports it in a mat file'};
%col8 = {'radio_coher0d1d'   ,'radio' ,'Coherence' ,0,'coherence between 0d data and 1D data; mismatch sign an problem of convergence'};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8};
col7 = {'radio_coher0d1d'   ,'radio' ,'Coherence' ,0,'coherence between 0d data and 1D data; mismatch sign an problem of convergence'};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,colj};

% 8eme ligne
% ----------
col1 = {'radio_l2h'     ,'radio' ,'L->H' ,0,'plot graphs about L-> H transition'};
col2 = {'radio_fluxpol'         ,'radio' ,'Flux Consumption' ,0,'Consumed poloidal flux (infinite solenoid hypotesis)'};
col3 = {'radio_shine'    ,'radio' ,'Shine through' ,0,'plot NBI shine through estimation'};
col4 = {'radio_scaling'        ,'radio' ,'HH & HL' ,0,'confinement factors'};
col5 = {'radio_fusion_power'    ,'radio' ,'Fusion' ,0,'Total fusion power (DT,DD and TT with neutrons contribution and not only fraction heating the plasma)'};
col6 = {'radio_divertor'      ,'radio' ,'Divertor' ,0,'sketch of power depostion on divertor target '};
col7 = {'radio_2pts'        ,'radio' ,'2 points ' ,0,'plot two points model results'};
%col7 = {'radio_ramp'      ,'radio' ,'ramp 2pts ' ,0,'plot two points model results with ramp of density'};
%col8 = {'radio_cost'    ,'radio' ,'cost' ,0,'device construction cost (estimation)'};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8};
col8 = {'radio_cost'    ,'radio' ,'cost' ,0,'device construction cost (estimation)'};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8};

%  % 91eme ligne
%  % ----------
%  col1 = {'radio_void91'      ,'radio' ,'' ,0,''};
%  col2 = {'radio_void92'      ,'radio' ,'' ,0,''};
%  col3 = {'radio_void93'      ,'radio' ,'' ,0,''};
%  col4 = {'radio_void94'      ,'radio' ,'' ,0,''};
%  col5 = {'radio_void95'      ,'radio' ,'' ,0,''};
%  col6 = {'radio_void96'      ,'radio' ,'' ,0,''};
%  col7 = {'radio_void97'      ,'radio' ,'' ,0,''};
%  col8 = {'radio_void98'      ,'radio' ,'' ,0,''};
%  form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8};

% post processing new section
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
col1 = {'titre2','text@full',' Post processing ',50,'',[]};
form{length(form)+1} = {col1} ;
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% first line for post processing commands
col1 = {'radio_ramp'      ,'radio' ,'ramp 2pts ' ,0,'plot two points model results with ramp of density'};
col2 = {'radio_export_cpos'    ,'radio' ,'CPOS' ,0,'re-computes CPOS and 2D equilibrium data then exports it in a mat file'};
col3 = {'radio_qlkANN_k'     ,'radio' ,'QlkANNk' ,0,'Plot flux from Qualikiz neural network with kinetic electrons'};
col4 = {'radio_qlkANN_k_wp'   ,'radio' ,'OP + QlkANNk' ,0,sprintf('Compute operating point for one time slice using kinetic profiles predicted with Qualikiz neural network\n(steady state, all d/dt = 0, only dPsi/dt #0,\ncreates external data for Te, Ti and Ne that can be reused;\n external data must be cleared before next call)')};
col5 = {'radio_qlkz_std_wp'   ,'radio' ,'OP + QLKZ' ,0,sprintf('Compute operating point for one time slice using kinetic profiles predicted with Qualikiz \n(steady state, all d/dt = 0, only dPsi/dt #0,\ncreates external data for Te, Ti and Ne that can be reused;\n external data must be cleared before next call)')};
col6 = {'radio_coils_currents'   ,'radio' ,'Coil currents' ,0,sprintf('Compute esitimation of current in poloidal coils using fast inverse module of FEEQS.M\n(FEEQS.M is not part of METIS and should be installed separately; not suitable for tokamak with iron core).')};
col7 = {'radio_plot_coils_currents'   ,'radio' ,'Plot coils currents' ,0,sprintf('Plot results of current in poloidal coils computed with fast inverse module of FEEQS.M\n(FEEQS.M is not part of METIS and should be installed separately).')};
col8 = {'radio_eqdsk'   ,'radio' ,'EQDSK file' ,0,sprintf('Computation of EQDSK files using fast fixe boundary equilibrium solver of FEEQS.M;\nextrapolation outside the lCFS will not included real x-point and is valid only far from external current sources\n(FEEQS.M is not part of METIS and should be installed separtely)')};
col9 = {'radio_equi_gs'   ,'radio' ,'G-S' ,0,sprintf('Creation of an external equilibrium structure for METIS by recomputing Grad-Shafranov equation\nbased on fast fixed boundary mode of FEEQS.M code and using METIS data as input of FEEQS.M\n(FEEQS.M is not part of METIS and should be installed separtely)')};
if isdeployed
    form{length(form)+1} = {col1,col2};
else
    form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7,col8,col9};
end
% separation
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% 11eme ligne : Bouton Quit
% ------------------------
comm{1}={'btn_quit','radio@left','Quit',0,'Close this window'};
comm{2}={'btn_init','radio@center','Initialisation',0,'Initialisation of a new METIS simulation'};
comm{3}={'aide','radio@right','User''s guide',0,'User''s guide (requires PDF reader installed)'};

hout=zuicreeform('Metis : Fast tokamak simulator','zeroda','zcall0d_fct','',form,comm);


if isappdata(0,'METIS_INTERFACE_TITLE')
      [hfig,h] = zuiformhandle('zeroda');
      set(hfig,'name',getappdata(0,'METIS_INTERFACE_TITLE'));
end        
z0dbasic_mode_switcher(get(findobj(zuiformhandle('zeroda'),'tag','popup_advance'),'value'));

