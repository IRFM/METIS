function zcall0d_fct(action)

% variable global importante pour les test
% utilisation du mexfile si disponible	 
if isappdata(0,'MEXSOLVER_IN_METIS')
	mexsolver = getappdata(0,'MEXSOLVER_IN_METIS');
else	
	mexsolver =  [];
end
if isempty(mexsolver)
	repmex=which(strcat('mexpde1dsolver.',mexext));
	if ~isempty(repmex)
		mexsolver = 1;
	else
		mexsolver = 0;	
	end
	setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);
end  


if nargin ==0
	action = ' ';
end
% disp('callback : ')
%disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('zeroda');
if ~isempty(hfig)
    ho = gco;
else
    ho = [];
end
% information pour l'assistant
zuicr(hfig,action);

% recherche du tooltip
info = get(ho,'tooltip'); 

% variables d'entre de l'editeur de consignes zuieditcons
if ~strcmp(action,'btn_init') & ~strcmp(action,'radio_save') & ~strcmp(action,'radio_load')
   try
       x           = evalin('base','z0dinput.cons.temps') ;
   catch
       x           = [];
   end
   texte_x     = 'time (s)' ;
   var_x       = 'void';
   canal       = 1 ;
   code_retour = '' ;
   liste_ref   = {'Ip     ','Flux','Nbar','Gaspuff','Zeff','B0','ECRH','ICRH', ...
        	  'LH','NBI1','NBI2','Fact H','FT_NBI1','FT_NBI2','R0','z0','a','K','d','Iso','Xece','Empty'} ;
   var_ref     = {  ...
                    {'z0dinput.cons.temps','z0dinput.cons.ip','-'},
                    {'z0dinput.cons.temps','z0dinput.cons.flux','-'},
                    {'z0dinput.cons.temps','z0dinput.cons.nbar','-','real'},
                    {'z0dinput.cons.temps','z0dinput.cons.nbar','-','imag'},
                    {'z0dinput.cons.temps','z0dinput.cons.zeff','-'},
                    {'z0dinput.cons.temps','z0dinput.geo.b0','-'},
                    {'z0dinput.cons.temps','z0dinput.cons.pecrh','-'},
                    {'z0dinput.cons.temps','z0dinput.cons.picrh','-'},
                    {'z0dinput.cons.temps','z0dinput.cons.plh','-'},
                    {'z0dinput.cons.temps','z0dinput.cons.pnbi','-','real'},
                    {'z0dinput.cons.temps','z0dinput.cons.pnbi','-','imag'},
                    {'z0dinput.cons.temps','z0dinput.cons.hmore','-'},
                    {'z0dinput.cons.temps','z0dinput.cons.ftnbi','-','real'},
                    {'z0dinput.cons.temps','z0dinput.cons.ftnbi','-','imag'},
                    {'z0dinput.cons.temps','z0dinput.geo.R','-'},
                    {'z0dinput.cons.temps','z0dinput.geo.z0','-'},
                    {'z0dinput.cons.temps','z0dinput.geo.a','-'},
                    {'z0dinput.cons.temps','z0dinput.geo.K','-'},
                    {'z0dinput.cons.temps','z0dinput.geo.d','-'},
                    {'z0dinput.cons.temps','z0dinput.cons.iso','-'},
                    {'z0dinput.cons.temps','z0dinput.cons.xece','-'},
		    {'[]','[]',''}  } ;
end

% reroutage si IMAS  et UAL
if isappdata(0,'IMAS_EXIST') && ~isempty(getappdata(0,'IMAS_EXIST')) && isappdata(0,'UALVERSION') && ~isempty(getappdata(0,'UALVERSION'))
    switch lower(action)
    case 'radio_save'
        rep = menu('Save to (Data destination)? ', ...
                         'File', 'IMAS','ITM','Cancel');
	
        switch rep,
        case 4
            if ishandle(ho)
                set(ho,'value',0);
                return
            end	
        case 2
              action = 'imas_save';
        case 3
              action = 'ual_save';
        end

    end  
elseif isappdata(0,'IMAS_EXIST') && ~isempty(getappdata(0,'IMAS_EXIST'))
    % reroutage si IMAS 
    switch lower(action)
    case 'radio_save'
        ButtonName = questdlg('Save to ?', ...
                         'Data destination', ...
                         'File', 'IMAS', 'Cancel', 'File');
        switch ButtonName,
        case 'Cancel'
            if ishandle(ho)
                set(ho,'value',0);
                return
            end	
        case 'IMAS'
              action = 'imas_save';
        end

    end
    % reroutage si UAL
elseif isappdata(0,'UALVERSION') && ~isempty(getappdata(0,'UALVERSION'))
    switch lower(action)
    case 'radio_save'
        ButtonName = questdlg('Save to ?', ...
                         'Data destination', ...
                         'File', 'UAL', 'Cancel', 'File');
        switch ButtonName,
        case 'Cancel'
            if ishandle(ho)
                set(ho,'value',0);
                return
            end	
        case 'UAL'
              action = 'ual_save';
        end

    end
end    
% selon ation
switch lower(action)
case 'btn_init'
        [hfig,hui] = zuiformhandle('zerod');
	zuifaitfunaction('close',hfig);
	zcall0d('init');
	if ishandle(ho)
        	set(ho,'value',0);
	end	   
        if isappdata(0,'METIS_INTERFACE_TITLE')
		[hfig,h] = zuiformhandle('zeroda');
		set(hfig,'name',getappdata(0,'METIS_INTERFACE_TITLE'));
        end
case {'btn_quit','close'}
	delete(hfig);
        if isdeployed
             exit
        end
case 'popup_advance'
	z0dbasic_mode_switcher(get(ho,'value'))
case 'radio_param'
        % test if data exist
        if evalin('base','~exist(''z0dinput'')')
            warning('No METIS data available in the workspace !');
            return
        end
        %
        if isappdata(0,'IMAS_EXIST') && ~isempty(strmatch(getappdata(0,'IMAS_EXIST'),'YES','exact'))
                    zuicreefunform('metis4imas','z0dinput.option',1,[],'',1);
        elseif isappdata(0,'UALVERSION') && ~isempty(getappdata(0,'UALVERSION'))
                    zuicreefunform('metis4itm','z0dinput.option',1,[],'',1);
        else
                    zuicreefunform('zerod','z0dinput.option',1,[],'',1);
        end
        
        set(ho,'value',0);
case 'radio_load'
    
    
    langue      =  lower(getappdata(0,'langue_cronos'));
    switch langue
        case 'francais'
            [file,path]=uigetfile('*.mat','Nom du fichier a charger ?');
        otherwise
            [file,path]=uigetfile('*.mat','load : choose a file ?');
    end
    drawnow
    if ~ischar(file)
        set(ho,'value',0);
        return
    end
    filename = strcat(path,file);
    data = load(filename,'post');
    if isempty(data) ||  (length(fieldnames(data)) == 0)
        data = load(filename,'z0dinput');
        fullresult = 0;
    else
        fullresult = 1;
    end
    if isempty(data) ||  (length(fieldnames(data)) == 0)
        switch langue
            case 'francais'
                warning('Ce fichier ne contient pas de donnees pour Metis')
            otherwise
                warning('This is not a Metis file');
        end
        return
    elseif fullresult == 1
        data =data.post;
    end
    if ~isfield(data,'zerod') && (fullresult == 1)
        switch langue
            case 'francais'
                warning('Ce fichier ne contient pas de donnees pour Metis')
            otherwise
                warning('This is not a Metis file');
        end
        return
    end
    if ~isfield(data,'z0dinput')
        switch langue
            case 'francais'
                warning('Ce fichier ne contient pas de donnees pour Metis')
            otherwise
                warning('This is not a Metis file');
        end
        return
    end
    
    % compatibilite entre version
    if isfield( data.z0dinput,'exp') & ~isfield( data.z0dinput,'exp0d')
        data.z0dinput.exp0d = data.z0dinput.exp;
    end
    
    % mise a niveau des anciennes versions
    model = zerod_init(-2,data.z0dinput.shot,data.z0dinput.option.gaz,data.z0dinput.cons.temps);
    
    % option
    noms = fieldnames(model.option);
    for k = 1:length(noms)
        if ~isfield(data.z0dinput.option,noms{k})
            data.z0dinput.option.(noms{k}) = model.option.(noms{k});
        end
    end
    % info
    noms = fieldnames(model.info);
    for k = 1:length(noms)
        if ~isfield(data.z0dinput.info,noms{k})
            data.z0dinput.info.(noms{k}) = model.info.(noms{k});
        end
    end
    % zsinfo
    noms = fieldnames(model.zsinfo);
    for k = 1:length(noms)
        if ~isfield(data.z0dinput.zsinfo,noms{k})
            data.z0dinput.zsinfo.(noms{k}) = model.zsinfo.(noms{k});
        end
    end
    % profinfo
    if isfield(data.z0dinput,'profinfo')
        noms = fieldnames(model.profinfo);
        for k = 1:length(noms)
            if ~isfield(data.z0dinput.profinfo,noms{k})
                data.z0dinput.profinfo.(noms{k}) = model.profinfo.(noms{k});
            end
        end
    else
        noms = fieldnames(model.profinfo);
        for k = 1:length(noms)
            data.z0dinput.profinfo.(noms{k}) = model.profinfo.(noms{k});
        end
    end
    % consignes
    noms = fieldnames(model.cons);
    for k = 1:length(noms)
        if ~isfield(data.z0dinput.cons,noms{k})
            data.z0dinput.cons.(noms{k}) = model.cons.(noms{k});
        end
    end
    % geo
    noms = fieldnames(model.geo);
    for k = 1:length(noms)
        if ~isfield(data.z0dinput.geo,noms{k})
            data.z0dinput.geo.(noms{k}) = model.geo.(noms{k});
        end
    end
    
    %zerod data
    if fullresult == 1
        noms = fieldnames(model.zsinfo);
        for k = 1:length(noms)
            if ~isfield(data.zerod,noms{k})
                data.zerod.(noms{k}) = NaN .* data.z0dinput.cons.temps;
            end
        end
    end
    
    % profil 0d
    if isfield(data,'profil0d')
        if isfield(data.profil0d,'temps')
            prnan = NaN .* ones(length(data.profil0d.temps),21);
        else
            prnan = NaN .* ones(length(data.z0dinput.cons.temps),21);
        end
        if ~isfield(data.profil0d,'xli')
            data.profil0d.xli = linspace(0,1,21);
        end
        
        noms = fieldnames(model.profinfo);
        for k = 1:length(noms)
            if ~isfield(data.profil0d,noms{k})
                data.profil0d.(noms{k}) = prnan;
            end
        end
    else
        prnan = NaN .* ones(length(data.z0dinput.cons.temps),21);
        noms = fieldnames(model.profinfo);
        for k = 1:length(noms)
            data.profil0d.(noms{k}) = prnan;
        end
        data.profil0d.xli   = linspace(0,1,21);
        data.profil0d.temps = data.z0dinput.cons.temps;
    end
    
    % mise a jour de la structure experimentale vide
    noms = fieldnames(data.z0dinput.zsinfo);
    if ~isfield(data.z0dinput,'exp0d')
        data.z0dinput.exp0d=[];
    end
    exp0d  = data.z0dinput.exp0d;
    if isfield(exp0d,'temps')
        texp = exp0d.temps;
    else
        texp = data.z0dinput.cons.temps;
        exp0d.temps = texp;
    end
    nbt  = length(texp);
    %
    vtnan = NaN .* ones(nbt,1);
    
    for k = 1:length(noms)
        nomc = noms{k};
        if isfield(exp0d,nomc)
            var = getfield(exp0d,nomc);
            if length(var) ~= nbt
                disp('dimension mismatch')
                var = mean(var(isfinite(var))) .* ones(nbt,1);
                exp0d = setfield(exp0d,nomc,var);
            else
                % si donnnees non valides
                fnan = imag(var);
                var  = real(var);
                ind  = find(fnan~=0 & var == 0);
                if  ~isempty(ind)
                    var(ind) = NaN;
                end
                exp0d  = setfield(exp0d,nomc,var);
            end
        else
            exp0d = setfield(exp0d,nomc,vtnan);
        end
    end
    % backward compatibilty (bug shortcut)
    if isfield(exp0d,'flux')
        if ~isfield(exp0d,'edgeflux')
            exp0d.edgeflux = exp0d.flux;
        elseif all(~isfinite(exp0d.edgeflux))
            exp0d.edgeflux = exp0d.flux;
        end
        exp0d = rmfield(exp0d,'flux');
    end
    data.z0dinput.exp0d = exp0d;
    
    % securite sur le donnees cronos
    if  data.z0dinput.mode_exp  == 0
        try
            rep = evalin('base','isstruct(data)');
        catch
            data.z0dinput.mode_exp  = -100;
        end
    end
    
    zassignin('base','z0dinput',data.z0dinput);
    
    
    if fullresult == 0
        set(ho,'value',0);
        z0dinterfacetitle(filename);
        warndlg('Warning: this data set contain only METIS input','METIS load');
        return
    end
    zassignin('base','post.zerod',data.zerod);
    zassignin('base','post.z0dinput',data.z0dinput);
    zassignin('base','post.profil0d',data.profil0d);
    
    if isfield(data,'lukeinmetis')
        zassignin('base','post.lukeinmetis',data.lukeinmetis);
    end
    
    if isfield(data,'simout')
        zassignin('base','post.simout',data.simout);
        zassignin('base','simout',data.simout);
    end
    
    % mode evolution depuis l'interface
    if ~isfield(data,'z0dstruct_inter')
        % rien
    elseif ~isempty(data.z0dstruct_inter) && isfield(data.z0dstruct_inter,'profil')
        prnan = NaN .* ones(length(data.z0dstruct_inter.profil.temps),21);
        noms = fieldnames(model.profinfo);
        for k = 1:length(noms)
            if ~isfield(data.z0dstruct_inter.profil,noms{k})
                data.z0dstruct_inter.profil.(noms{k}) = prnan;
            end
        end
        zassignin('base','post.z0dstruct_inter',data.z0dstruct_inter);
        zassignin('base','z0dstruct',data.z0dstruct_inter);
    end
    
    
    % cas du couplage avec simulink
    if isfield(data.z0dinput,'system')
        vv = ver('simulink');
        if length(vv) > 0
            % restauration du model simulink
            chem = which(data.z0dinput.system.name);
            if isempty(chem)
                fprintf('Simulink system %s does not exist in Matlab path\n',data.z0dinput.system.name);
                if isdir(fileparts(data.z0dinput.system.fullname))
                    fprintf('Simulink system %s will be created in %s \n',data.z0dinput.system.name, ...
                        fileparts(data.z0dinput.system.fullname));
                    newname = data.z0dinput.system.fullname;
                    try
                        addpath(fileparts(data.z0dinput.system.fullname));
                    end
                else
                    fprintf('Simulink system %s will be created in %s \n',data.z0dinput.system.name,pwd);
                    [void,newname,ext] = fileparts(data.z0dinput.system.fullname);
                    newname = strcat(newname,ext);
                    try
                        addpath(pwd)
                    end
                end
                [fid,mess] = fopen(newname,'w');
                if fid >= 3
                    fprintf(fid,'%s',data.z0dinput.system.mdl);
                    fclose(fid);
                    drawnow
                    void=which(data.z0dinput.system.name);
                    try
                        evalin('base',data.z0dinput.system.name);
                    end
                else
                    error(mess)
                end
            else
                %[s,mdl_loc]  = unix(sprintf('cat %s',chem));
                fid = fopen(chem,'r');
                if fid > 0
                    mdl_loc = char(fread(fid,Inf,'char')');
                    fclose(fid);
                else
                    mdl_loc ='';
                end
                if strcmp(mdl_loc,data.z0dinput.system.mdl)
                    % model identique ouverture
                    try
                        evalin('base',data.z0dinput.system.name);
                    end
                else
                    fprintf('Simulink system %s as some differences with system saved in Metis data\n',data.z0dinput.system.name);
                    tage = num2str(fix(datenum(clock)*1e5));
                    reste = 62 - length(tage);
                    reste = min(length(data.z0dinput.system.name),reste);
                    newname =strcat(data.z0dinput.system.name(1:reste),tage);
                    fprintf('Simulink system %s will be created in same directory\n',newname);
                    if ~isempty(strfind(chem,'.slx'))
                        [fid,mess] = fopen(fullfile(fileparts(chem),strcat(newname,'.slx')),'w');
                    else
                        [fid,mess] = fopen(fullfile(fileparts(chem),strcat(newname,'.mdl')),'w');
                    end
                    if fid >= 3
                        if ~isempty(strfind(chem,'.slx'))
                            fwrite(fid,data.z0dinput.system.mdl,'uint8');
                        else
                            fprintf(fid,'%s',data.z0dinput.system.mdl);
                        end
                        fclose(fid);
                        drawnow
                        void=which(newname);
                        try
                            evalin('base',newname);
                        catch
                            fprintf('Simulink system %s can''t be opened\n',newname);
                        end
                    else
                        error(mess)
                    end
                end
            end
            zassignin('base','z0dinput',data.z0dinput.system.z0dinput);
            
        end
    end
    % test s'il y a des donnees externe
    if isfield(data,'appdata')
        noms = fieldnames(data.appdata);
        externaldata = 0;
        for k = 1:length(noms)
            if findstr(noms{k},'_EXP')
                fprintf('Loading external data in METIS for %s\n',strtok(noms{k},'_'));
                setappdata(0,noms{k},data.appdata.(noms{k}));
                externaldata = 1;
            end
        end
        if externaldata
            fprintf('METIS have some external data\n You can remove the external data with the instruction "rmappdata(0,''DATA_NAME'')"\n')
            fprintf('to obtain the list of data : getappdat(0); the METIS external data are name "XXXX_EXP"\n')
            help('external_data_rule_for_METIS')
        end
    else
        % information sur les donnees externes
        noms = fieldnames(getappdata(0));
        text_warn = '';
        for k = 1:length(noms)
            if findstr(noms{k},'_EXP')
                fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
                text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
            end
        end
        if  ~isempty(text_warn)
            warndlg(sprintf('Pay attention: external data from previous simulation in use\n%s',text_warn), ...
                    'Pay attention: external data from previous simulation use');
        end
    end
    
    % restrore sepa option if available
    if isfield(data.z0dinput,'sepa_option')
        zassignin('base','sepa_option',data.z0dinput.sepa_option);
        if isdeployed
            evalin('base','sepa_option');
        end
    end
    
    
    
    set(ho,'value',0);
    z0dinterfacetitle(filename);
case 'radio_save'
	try
		post = evalin('base','post');
	catch
		post = [];
	end
        langue      =  lower(getappdata(0,'langue_cronos'));
	if isempty(post)
        	switch langue
        	case 'francais'
			warning('Pas de donnee a sauver')
        	otherwise
			warning('No data to be saved');
        	end
        	set(ho,'value',0);
		return
	end
	if ~isfield(post,'zerod')
        	switch langue
        	case 'francais'
			warning('Pas de donnee a sauver')
        	otherwise
			warning('No data to be saved');
        	end
        	set(ho,'value',0);
		return
	end
	if ~isfield(post,'z0dinput')
        	switch langue
        	case 'francais'
			warning('Pas de donnee a sauver')
        	otherwise
			warning('No data to be saved');
        	end
        	set(ho,'value',0);
		return
	end
	data.zerod    = post.zerod;
	data.z0dinput = post.z0dinput;
	data.profil0d = post.profil0d;
	% sauvegarde des donnees externes
	noms = fieldnames(getappdata(0));
	for k = 1:length(noms)
		if findstr(noms{k},'_EXP')
			data.appdata.(noms{k}) = getappdata(0,noms{k});
		end
	end

	% sauvegarde de la sortie standard du model simulink
	try
		simout = evalin('base','simout');
	catch
		simout = [];
	end
   	zassignin('base','post.simout',simout);
	data.simout = simout;
	
	
	% sauvegarde  des donnees de la simulation en mode evolution
	try
		z0dstruct = evalin('base','z0dstruct');
	catch
		z0dstruct = [];
	end
	zassignin('base','post.z0dstruct_inter',z0dstruct);
	data.z0dstruct_inter = z0dstruct;
    
        % donnees LUKE  dans METIS
        if isfield(post,'lukeinmetis')
   	     data.lukeinmetis = post.lukeinmetis;
        end   


	% ajout des donnees de certication
	root = getappdata(0,'root');
	if isempty(root)
		zineb_path;
		root = getappdata(0,'root');
	end
	% structure d'information 
	data.info_test_metis.date        = clock;
    if ~ispc
        [s,t] = unix('uname -a');
        if s == 0
            data.info_test_metis.machine = t;
        else
            error(sprintf('error executing ''uname -a'' (%s)',t));
        end
        [s,t] = unix('env');
        if s == 0
            data.info_test_metis.env = t;
        else
            error(sprintf('error reading environnement variables  (%s)',t));
        end
    else
        data.info_test_metis.machine =getenv('COMPUTERNAME');
        [s,t] = system('set')
        if s == 0
            data.info_test_metis.env = t;
        else
            error(sprintf('error reading environnement variables  (%s)',t));
        end
    end
	[data.info_test_metis.version,data.info_test_metis.date_version]  = zinebversion;
	data.info_test_metis.matlab_version  = version;
	data.info_test_metis.toolbox_version = ver;
	data.info_test_metis.path = matlabpath;
	if ispc
		data.info_test_metis.user = getenv('USERNAME');
	else
		data.info_test_metis.user = getenv('USER');
	end
	data.info_test_metis.root = root;
	data.info_test_metis.mexsolver = mexsolver;
    
    % save cooling ratse to be able to know which one has been used for the
    % simulation.
    try
        load('Lz_zave.mat','tabmat')
        data.tabmat = tabmat;
    catch
        disp('metis_save: unable to read cooling rate');
    end

	post = data;
	clear data
        switch langue
        case 'francais'
	    	[file,path]=uiputfile('*.mat','Nom du fichier de sauvegarde ?');
        otherwise
		[file,path]=uiputfile('*.mat','Save : choose a file ?');
        end
	drawnow
	if ~ischar(file)        	
		set(ho,'value',0);
		return
	end
	filename = fullfile(path,file);
	save(filename,'post');
        set(ho,'value',0);
        z0dinterfacetitle(filename);
        
case 'ual_save'
    
	try
		post = evalin('base','post');
		z0dinput = evalin('base','z0dinput');
		post = copy_ual_parameters_2_output(post,z0dinput);
	catch
		post = [];
	end
        langue      =  lower(getappdata(0,'langue_cronos'));
	if isempty(post)
        	switch langue
        	case 'francais'
			warning('Pas de donnee a sauver')
        	otherwise
			warning('No data to be saved');
        	end
		set(ho,'value',0);
		return
	end
	if ~isfield(post,'zerod')
        	switch langue
        	case 'francais'
			warning('Pas de donnee a sauver')
        	otherwise
			warning('No data to be saved');
        	end
		set(ho,'value',0);
		return
	end
	if ~isfield(post,'z0dinput')
        	switch langue
        	case 'francais'
			warning('Pas de donnee a sauver')
        	otherwise
			warning('No data to be saved');
        	end
		set(ho,'value',0);
		return
	end
	data.zerod    = post.zerod;
	data.z0dinput = post.z0dinput;
	data.profil0d = post.profil0d;
	% sauvegarde des donnees externes
	noms = fieldnames(getappdata(0));
	for k = 1:length(noms)
		if findstr(noms{k},'_EXP')
			data.appdata.(noms{k}) = getappdata(0,noms{k});
		end
	end

         if isfield(post,'lukeinmetis')
   	     data.lukeinmetis = post.lukeinmetis;
         end   

	% sauvegarde de la sortie standard du model simulink
	try
		simout = evalin('base','simout');
	catch
		simout = [];
	end
   	zassignin('base','post.simout',simout);
	data.simout = simout;

	% ajout des donnees de certication
	root = getappdata(0,'root');
	if isempty(root)
		zineb_path;
		root = getappdata(0,'root');
	end
	% structure d'information 
	data.info_test_metis.date        = clock;
    if ~ispc
        [s,t] = unix('uname -a');
        if s == 0
            data.info_test_metis.machine = t;
        else
            error(sprintf('error executing ''uname -a'' (%s)',t));
        end
        [s,t] = unix('env');
        if s == 0
            data.info_test_metis.env = t;
        else
            error(sprintf('error reading environnement variables  (%s)',t));
        end
    else
        data.info_test_metis.machine =getenv('COMPUTERNAME');
        [s,t] = system('set')
        if s == 0
            data.info_test_metis.env = t;
        else
            error(sprintf('error reading environnement variables  (%s)',t));
        end
    end
	[data.info_test_metis.version,data.info_test_metis.date_version]  = zinebversion;
	data.info_test_metis.matlab_version  = version;
	data.info_test_metis.toolbox_version = ver;
	data.info_test_metis.path = matlabpath;
    if ispc
            data.info_test_metis.user = getenv('USERNAME');
    else
            data.info_test_metis.user = getenv('USER');
    end 
    
    data.info_test_metis.root = root;
    data.info_test_metis.mexsolver = mexsolver;
    post = data;
    clear data

    if isfield(post.z0dinput,'run')
        run_num = num2str(post.z0dinput.run + 1);
    else
        run_num = '1';
    end
    if isfield(post.z0dinput,'occ')
        occ = post.z0dinput.occ;
    else
        occ = '';
    end
    if isfield(post.z0dinput,'tokamak')
        ual_tok = post.z0dinput.tokamak;
    else
        ual_tok = '';
    end
    if isfield(post.z0dinput,'user')
        ual_user = post.z0dinput.user;
    else
        ual_user = '';
    end
    if isfield(post.z0dinput,'dataversion')
        ual_ver = post.z0dinput.dataversion;
    else
        ual_ver = '';
    end
    prompt={'shot number :                         ','Run :' ,'Occurrence:','Tokamak:','User:','Version:'};
    def={num2str(post.z0dinput.shot),run_num,occ,ual_tok,ual_user,ual_ver};
    dlgTitle='Saving data to UAL';
    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    drawnow
    if isempty(answer)
       set(ho,'value',0);
       return
    end
    shot  = str2num(answer{1});
    run   = str2num(answer{2});
    occ   = answer{3};
    ual_tok   = answer{4};
    ual_user   = answer{5};
    ual_ver    = answer{6};
    if ~isempty(ual_tok) && ~isempty(ual_tok(ual_tok > ' '))
	  ual_occ = occ;
          clear occ
          occ.tokamak = ual_tok(ual_tok>' ');
          occ.user    = ual_user(ual_user> ' ');
          occ.dataversion = ual_ver(ual_ver> ' ');
          occ.occurrence = ual_occ(ual_occ> ' ');
    end
    error_flag = metis4itm(shot,run,occ,post);
    post.z0dinput.run = run;
    if isstruct(occ)
    	post.z0dinput.tokamak = occ.tokamak;
    	post.z0dinput.user = occ.user;
    	post.z0dinput.dataversion = occ.dataversion;
    	post.z0dinput.occ = occ.occurrence;
    else
    	post.z0dinput.occ = occ;
    end
    post.z0dinput.shot = shot;
    post.z0dinput.option.shot = shot;
    zassignin('base','post',post);
    set(ho,'value',0);
    txt = 'Metis : Fast tokamak simulator';
    if isstruct(occ)
     	txt = sprintf('%s (%s@%d for run = %d and occ = %s)',txt,post.z0dinput.machine,shot,run,occ.occurrence);
   else
    	txt = sprintf('%s (%s@%d for run = %d and occ = %s)',txt,post.z0dinput.machine,shot,run,occ);
    end
    setappdata(0,'METIS_INTERFACE_TITLE',txt);

    [hfig,h] = zuiformhandle('zeroda');
    if isempty(hfig)
      return
    end

    set(hfig,'name',txt);


case 'imas_save'
    
	try
	  post = evalin('base','post');
	  z0dinput = evalin('base','z0dinput');
	  post = copy_ual_parameters_2_output(post,z0dinput);
	catch
	  post = [];
	end
	if isempty(post)
	  warning('No data to be saved');
	  set(ho,'value',0);
	  return
	end
	if ~isfield(post,'zerod')
	  warning('No data to be saved');
	  set(ho,'value',0);
	  return
	end
	if ~isfield(post,'z0dinput')
	  warning('No data to be saved');
	  set(ho,'value',0);
	  return
	end
	data.zerod    = post.zerod;
	data.z0dinput = post.z0dinput;
	data.profil0d = post.profil0d;
	% sauvegarde des donnees externes
	noms = fieldnames(getappdata(0));
	for k = 1:length(noms)
	  if findstr(noms{k},'_EXP')
	    data.appdata.(noms{k}) = getappdata(0,noms{k});
	  end
	end

	if isfield(post,'lukeinmetis')
	  data.lukeinmetis = post.lukeinmetis;
	end   

	% sauvegarde de la sortie standard du model simulink
	try
	  simout = evalin('base','simout');
	catch
	  simout = [];
	end
   	zassignin('base','post.simout',simout);
	data.simout = simout;

	% ajout des donnees de certication
	root = getappdata(0,'root');
	if isempty(root)
	  zineb_path;
	  root = getappdata(0,'root');
	end
	% structure d'information 
	data.info_test_metis.date        = clock;
	if ~ispc
	  [s,t] = unix('uname -a');
	  if s == 0
            data.info_test_metis.machine = t;
	  else
            error(sprintf('error executing ''uname -a'' (%s)',t));
	  end
	  [s,t] = unix('env');
	  if s == 0
            data.info_test_metis.env = t;
	  else
            error(sprintf('error reading environnement variables  (%s)',t));
	  end
	else
	  data.info_test_metis.machine =getenv('COMPUTERNAME');
	  [s,t] = system('set')
	  if s == 0
            data.info_test_metis.env = t;
	  else
            error(sprintf('error reading environnement variables  (%s)',t));
	  end
	end
	[data.info_test_metis.version,data.info_test_metis.date_version]  = zinebversion;
	data.info_test_metis.matlab_version  = version;
	data.info_test_metis.toolbox_version = ver;
	data.info_test_metis.path = matlabpath;
	if ispc
	  data.info_test_metis.user = getenv('USERNAME');
	else
	  data.info_test_metis.user = getenv('USER');
	end 
    
	data.info_test_metis.root = root;
	data.info_test_metis.mexsolver = mexsolver;
	post = data;
	clear data

	if isfield(post.z0dinput,'run')
	  run_num = num2str(post.z0dinput.run + 1);
	else
	  run_num = '1';
	end
	if isfield(post.z0dinput,'occ')
	  occ = post.z0dinput.occ;
	else
	  occ = '';
	end
	if isfield(post.z0dinput,'tokamak')
	  ual_tok = lower(post.z0dinput.tokamak);
	elseif isfield(post.z0dinput,'machine')
	  ual_tok = lower(post.z0dinput.machine);
	elseif isappdata(0,'UAL_TOKAMAK')
	  ual_tok = getappdata(0,'UAL_TOKAMAK');
	else
	  ual_tok = '';
	end
	if isfield(post.z0dinput,'user')
	  ual_user = post.z0dinput.user;
	elseif isappdata(0,'UAL_USER')
	  ual_user = getappdata(0,'UAL_USER');
	elseif ispc
	  ual_user = getenv('USERNAME');
	else
	  ual_user = getenv('USER');
	end 
	if isfield(post.z0dinput,'dataversion')
	  ual_ver = post.z0dinput.dataversion;
	elseif isappdata(0,'UAL_DATAVERSION')
	  ual_ver = getappdata(0,'UAL_DATAVERSION');
	else
	  ual_ver = getenv('IMAS_VERSION');
	end
	if isfield(post.z0dinput,'backend')
	    ual_backend =post.z0dinput.backend;
	elseif isappdata(0,'UAL_BACKEND')
	    ual_backend = getappdata(0,'UAL_BACKEND');
	elseif ~isempty(getenv('IMAS_AL_BACKEND'))
	    ual_backend = getenv('IMAS_AL_BACKEND');
    else
        ual_backend = '';
	end
	if isnumeric(ual_backend)
	    ual_backend = num2str(ual_backend);
	end
	if ~isempty(strfind(strtrim(upper(post.z0dinput.machine)),'WEST'))
	    time_shift = tsbase(post.z0dinput.shot,'rignitron');
	    if ~isempty(time_shift)
		mode_tf = true;
		prompt={'Comment :','shot number :','Run :' ,'Occurrence:','Tokamak:','User:','Version:','Time shift','Backend'};
	        def={'                                                                           ',num2str(post.z0dinput.shot),run_num,occ,ual_tok,ual_user,ual_ver,time_shift,ual_backend};
	    else
	        mode_tf = false;
		prompt={'Comment :','shot number :','Run :' ,'Occurrence:','Tokamak:','User:','Version:','Backend'};
		def={'                                                                           ',num2str(post.z0dinput.shot),run_num,occ,ual_tok,ual_user,ual_ver,ual_backend};	    
	    end
	else
	    mode_tf = false;
	    prompt={'Comment :','shot number :','Run :' ,'Occurrence:','Tokamak:','User:','Version:','Backend'};
	    def={'                                                                           ',num2str(post.z0dinput.shot),run_num,occ,ual_tok,ual_user,ual_ver,ual_backend};
	end
	dlgTitle='Saving data to UAL';
	lineNo=1;
	answer=zinputdlg(prompt,dlgTitle,lineNo,def);
	drawnow
	if isempty(answer)
	  set(ho,'value',0);
	  return
	end
	comment  = strtrim(answer{1});
	shot  = str2num(answer{2});
	run   = str2num(answer{3});
	occ   = strtrim(answer{4});
	ual_tok   = strtrim(answer{5});
	ual_user   = strtrim(answer{6});
	ual_ver    = strtrim(answer{7});
	ual_backend  = strtrim(answer{end});
	if ~isempty(ual_tok) && ~isempty(ual_tok(ual_tok > ' '))
	  ual_occ = occ;
          clear occ
          occ.tokamak = ual_tok(ual_tok>' ');
          occ.user    = ual_user(ual_user> ' ');
          occ.dataversion = ual_ver(ual_ver> ' ');
          occ.occurrence = ual_occ(ual_occ> ' ');
          occ.backend = ual_backend(ual_backend> ' ');
	end
	% if available, time shifting
	post_imas = post;
	if length(answer) > 8
	    time_shift = str2num(answer{8});
	    if ~isempty(time_shift)
	      post_imas.zerod.temps = post_imas.zerod.temps + time_shift;
              post_imas.z0dinput.cons.temps = post_imas.z0dinput.cons.temps + time_shift;
              post_imas.z0dinput.exp0d.temps = post_imas.z0dinput.exp0d.temps + time_shift;
              post_imas.profil0d.temps = post_imas.profil0d.temps + time_shift;
 	    end	
	end
	
	% set comment 
	if ~isempty(comment)
	    if isappdata(0,'METIS_INTERFACE_TITLE')
		      setappdata(0,'METIS_INTERFACE_TITLE',sprintf('%s (%s)',comment,getappdata(0,'METIS_INTERFACE_TITLE')));
	    else
		      setappdata(0,'METIS_INTERFACE_TITLE',comment);
	    end
        end
	% write imas data
	error_flag = metis4imas(shot,run,occ,post_imas);
	post.z0dinput.run = run;
	if isstruct(occ)
	  post.z0dinput.tokamak = occ.tokamak;
	  post.z0dinput.user = occ.user;
	  post.z0dinput.dataversion = occ.dataversion;
	  post.z0dinput.occ = occ.occurrence;
	  post.z0dinput.backend = occ.backend;
	else
	  post.z0dinput.occ = occ;
	end
	post.z0dinput.shot = shot;
	post.z0dinput.option.shot = shot;
	zassignin('base','post',post);
	set(ho,'value',0);
	txt = 'Metis : Fast tokamak simulator';
	if isstruct(occ)
	  txt = sprintf('%s (%s@%d for run = %d and occ = %s)',txt,occ.tokamak,shot,run,occ.occurrence);
	else
	  txt = sprintf('%s (%s@%d for run = %d and occ = %s)',txt,post.z0dinput.machine,shot,run,occ);
	end
	setappdata(0,'METIS_INTERFACE_TITLE',txt);
	
	% if tokamak is west, ask for post processing: computation of synthetic daignostic
	post_func = sprintf('metis_post_process_%s',strtrim(lower(post.z0dinput.machine)));
	if exist(post_func);
	    if isstruct(occ)
		evalin('base',sprintf('%s(%d,%d,''%s'',''%s'',''%s'',''%s'',''%s'');',post_func,shot,run,occ.occurrence,occ.user,occ.tokamak,occ.dataversion,occ.backend));
	    end
	end
	
	[hfig,h] = zuiformhandle('zeroda');
	if isempty(hfig)
	  return
	end

	set(hfig,'name',txt);

case 'radio_export'
	% export data in mat file using IMAS data model
	try
		post = evalin('base','post');	  
		z0dinput = evalin('base','z0dinput');
		post = copy_ual_parameters_2_output(post,z0dinput);
	catch
		post = [];
	end
	if isempty(post)
		warning('No data to be exported');
		return
	end
	if ~isfield(post,'zerod')
		warning('No data to be exported');
		return
	end
	if ~isfield(post,'z0dinput')
		warning('No data to be exported');
 		return
	end		
        [file,path]=uiputfile('*.mat','Eport METIS data in IMAS data model stored in .mat file : choose a file ?');
	drawnow
	if ~ischar(file)
		set(ho,'value',0);
		return
	end
	imas_version = NaN;
	filename = fullfile(path,file);
	if isfield(post.z0dinput,'run')
	    run = post.z0dinput.run;
	else
            run = 0;
	end
    if isfield(post.z0dinput,'tokamak')
        tokamak = post.z0dinput.tokamak;
    else
        tokamak = post.z0dinput.machine;
    end
    shot = post.z0dinput.shot;
    [error_flag,exported_data] = metis4imas(shot,run,'',post,[]);
    if error_flag ~=  0
        error('unable to create Matlab IMAS like datastructure');
    else
        % update fair data
        if isfield(exported_data,'dataset_fair') && ~isempty(exported_data.dataset_fair)
            exported_data.dataset_fair.replaces       = exported_data.dataset_fair.identifier;
            hn = z0hostname;
            if filename(1) == filesep
                exported_data.dataset_fair.identifier     = sprintf('file://%s%s',hn,filename);
            else
                exported_data.dataset_fair.identifier     = sprintf('file://%s/%s',hn,filename);
            end           
            exported_data.dataset_fair.is_replaced_by = exported_data.dataset_fair.identifier;
        end
        lastwarn('','');
        save(filename,'exported_data','imas_version','tokamak','shot','run','filename','-v7');
        [~, msgid] = lastwarn;
        if strcmp(msgid,'MATLAB:save:sizeTooBigForMATFile')
            disp('change save format to 7.3');
            save(filename,'exported_data','imas_version','tokamak','shot','run','filename','-v7.3');
        end
    end
    set(ho,'value',0);

case 'radio_ref'
        % cas du mode standalone
        if isdeployed 
	    try
	      evalin('base','jeux1.post = post;');
	      set(ho,'value',0);
 	    end
            return
        end
	try
		post = evalin('base','post');
	catch
		post = [];
	end
        langue      =  lower(getappdata(0,'langue_cronos'));
	if isempty(post)
        	switch langue
        	case 'francais'
			warning('Pas de donnee a sauver')
        	otherwise
			warning('No data to be saved');
        	end
		return
	end
	if ~isfield(post,'zerod')
        	switch langue
        	case 'francais'
			warning('Pas de donnee a sauver')
        	otherwise
			warning('No data to be saved');
        	end
		return
	end
	if ~isfield(post,'z0dinput')
        	switch langue
        	case 'francais'
			warning('Pas de donnee a sauver')
        	otherwise
			warning('No data to be saved');
        	end
		return
	end
	evalin('base','jeux1.post = post;z0dlistparam(post);');

        set(ho,'value',0);

%  case 'radio_workspace'
%              workspace;
%              set(ho,'value',0);

case 'radio_audit'
	evalin('base','audit_metis_runs');
        set(ho,'value',0);    
        
case 'radio_pdf'
	rapfile = z0rapport;
        langue      =  lower(getappdata(0,'langue_cronos'));
        switch langue
        case 'francais'
	    	helpdlg(sprintf('document %s.pdf cree',rapfile),'Le document est disponible');
        otherwise
	    	helpdlg(sprintf('the repport name is %s.pdf',rapfile),'PDF ready !');
        end
       set(ho,'value',0);
       
case 'radio_closeall'
    list = findobj(0,'type','figure');
    main = findobj(0,'type','figure','tag','zeroda');
    for k=1:length(list)
      if list(k) ~= main  
    	close(list(k));
      end
    end
    set(ho,'value',0);
    
case 'radio_ip'
	nom        = 'z0dinput.cons.ip' ;
	y          = evalin('base',nom) ;
	texte_y    = 'Ip' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);

case 'radio_flux'
	nom        = 'z0dinput.cons.flux' ;
	y          = evalin('base',nom) ;
	texte_y    = 'Flux' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);

case 'radio_nbar'	  
	nom        = 'z0dinput.cons.nbar' ;
	y          = evalin('base',nom) ;
        %if any(imag(y))

	    texte_y    = 'Nbar' ;
	    var_y      = nom;
	    texte_prop = '' ;
	    var_prop   = '' ;
	    zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
			'nbar',liste_ref,var_ref,texte_prop,var_prop) ;

	    texte_y    = 'Gas puff' ;
	    var_y      = nom;
	    texte_prop = '' ;
	    var_prop   = '' ;
	    zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
			'gaspuff',liste_ref,var_ref,texte_prop,var_prop) ;
        
	%else    
	%    texte_y    = 'Nbar' ;
	%    var_y      = nom;
	%    texte_prop = '' ;
	%    var_prop   = '' ;
	%    zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	%		code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	%end
        set(ho,'value',0);	   
		    
case 'radio_iso'	  
    gas = evalin('base','z0dinput.option.gaz');	
    if gas == 5
        nom        = 'z0dinput.cons.iso' ;
        y          = evalin('base',nom) ;
        texte_y    = 'nhe3/nD' ;
        var_y      = nom;
        texte_prop = '' ;
        var_prop   = '' ;
        zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    'nhe3nd',liste_ref,var_ref,texte_prop,var_prop) ;
        
        texte_y    = 'nT/nD' ;
        var_y      = nom;
        texte_prop = '' ;
        var_prop   = '' ;
        zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    'ntnd',liste_ref,var_ref,texte_prop,var_prop) ;
   else
        nom        = 'z0dinput.cons.iso' ;
        y          = real(evalin('base',nom));
        texte_y    = 'Iso' ;
        var_y      = nom;
        texte_prop = '' ;
        var_prop   = '' ;
        zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
            code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);
    end
		    
case 'radio_zeff'	  
	nom        = 'z0dinput.cons.zeff' ;
	y          = evalin('base',nom) ;
	texte_y    = 'Zeff' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);	   
		    
case 'radio_xece'
	nom        = 'z0dinput.cons.xece' ;
	y          = evalin('base',nom) ;
	texte_y    = 'xece' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);
		    
case 'radio_hmore'	  
	nom        = 'z0dinput.cons.hmore' ;
	y          = evalin('base',nom) ;
	texte_y    = 'Fact H' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);

    case 'radio_ftnbi'
        nb_nbi = evalin('base','z0dinput.option.nb_nbi');
        nom        = 'z0dinput.cons.ftnbi' ;
        y          = evalin('base',nom) ;
        if nb_nbi == 2
            
            texte_y    = 'Ftnbi injector 1';
            var_y      = nom;
            texte_prop = '' ;
            var_prop   = '' ;
            zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                'idn',liste_ref,var_ref,texte_prop,var_prop) ;
            
            texte_y    = 'Ftnbi injector 2';
            var_y      = nom;
            texte_prop = '' ;
            var_prop   = '' ;
            zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                'idn2',liste_ref,var_ref,texte_prop,var_prop) ;
        else
            texte_y    = 'Ftnbi' ;
            var_y      = nom;
            texte_prop = '' ;
            var_prop   = '' ;
            zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        end
        set(ho,'value',0);

case 'radio_ecrh'	  
	nom        = 'z0dinput.cons.pecrh' ;
	y          = evalin('base',nom) ;
	texte_y    = 'ECRH' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);	   
		    
case 'radio_icrh'	  
	nom        = 'z0dinput.cons.picrh' ;
	y          = evalin('base',nom) ;
	texte_y    = 'ICRH' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);	   
		    
case 'radio_lh'	  
	nom        = 'z0dinput.cons.plh' ;
	y          = evalin('base',nom) ;
	texte_y    = 'LH' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);	   
		    
case 'radio_nbi'
        nb_nbi = evalin('base','z0dinput.option.nb_nbi');	  
	nom        = 'z0dinput.cons.pnbi' ;
	y          = evalin('base',nom) ;
        if nb_nbi == 2

	  texte_y    = 'NBI injector 1';
	  var_y      = nom;
	  texte_prop = '' ;
	  var_prop   = '' ;
	  zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
		      'idn',liste_ref,var_ref,texte_prop,var_prop) ;

	  texte_y    = 'NBI injector 2';
	  var_y      = nom;
	  texte_prop = '' ;
	  var_prop   = '' ;
	  zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
		      'idn2',liste_ref,var_ref,texte_prop,var_prop) ;
	else
	  texte_y    = 'NBI' ;
	  var_y      = nom;
	  texte_prop = '' ;
	  var_prop   = '' ;
	  zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
		      code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	end	    
        set(ho,'value',0);	   
case 'radio_r0'	  
	nom        = 'z0dinput.geo.R' ;
	y          = evalin('base',nom) ;
	texte_y    = 'R0' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    'LCFS',liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);
	%zclear0dsepa;

case 'radio_z0'	  
	nom        = 'z0dinput.geo.z0' ;
	y          = evalin('base',nom) ;
	texte_y    = 'z0' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    'LCFS',liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);
	%zclear0dsepa;

case 'radio_a'	  
	nom        = 'z0dinput.geo.a' ;
	y          = evalin('base',nom) ;
	texte_y    = 'a' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    'LCFS',liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);	   
	%zclear0dsepa;

case 'radio_k'	  
	nom        = 'z0dinput.geo.K' ;
	y          = evalin('base',nom) ;
	texte_y    = 'K' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    'LCFS',liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);	   
	%zclear0dsepa;

case 'radio_d'	  
	nom        = 'z0dinput.geo.d' ;
	y          = evalin('base',nom) ;
	texte_y    = 'd' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    'LCFS',liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);	   
	%zclear0dsepa;

case 'radio_b0'	  
	nom        = 'z0dinput.geo.b0' ;
	y          = evalin('base',nom) ;
	texte_y    = 'B0' ;
	var_y      = nom;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
        set(ho,'value',0);	   
		    
case 'radio_noise'

       evalin('base','z0dinput = z0dnoise(z0dinput);z0dinput.denoising = ''yes'';');
       set(ho,'value',0);	   
	  
case 'radio_temps'

       evalin('base','z0djivaro;');
       set(ho,'value',0);	   
	  
case 'radio_make_flux'

       evalin('base','z0dinput.cons.flux = post.zerod.edgeflux ./ 2 ./ pi;z0dinput.exp0d.edgeflux = post.zerod.edgeflux;');
       set(ho,'value',0);
	  
case 'radio_external'

       evalin('base','cs4m;');
       set(ho,'value',0);	   
	  
case 'radio_luke'

       evalin('base','cs4m_luke;');
       set(ho,'value',0);	   
	  
case 'radio_remove'

       evalin('base','rm_cs4m(1);');
       set(ho,'value',0);	   
	  
case 'radio_go'	
	hdlg = z0dpatience('full');
        close_0plot;
	drawnow

	% information sur les donnees externes
	noms = fieldnames(getappdata(0));
    text_warn = '';
	for k = 1:length(noms)
		if findstr(noms{k},'_EXP')
			fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
			text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
		end
    end
    if isdeployed && ~isempty(text_warn)
        warndlg(text_warn,'Pay attention: external data used');
    end

    % pour eviter les incoherence sur le mode restart
    evalin('base','clear  z0dstruct');
        
	% overwrite parameters whith reference parameters if defined
	evalin('base','z0dinput.option = z0doverwriteparam(z0dinput.option);');

	% sort du mode evolution
	evalin('base','z0dinput.option.evolution = 0;');

	% recopie machine dans option
	evalin('base','zerod_machine_name;');
	        
	evalin('base','[post.zerod,void,post.profil0d] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);post.z0dinput = z0dinput;z0plotsc;', ...
		'errordlg(lasterr,''Error during the run of Metis simulator !'');');
	delete(hdlg)
        set(ho,'value',0);	   

case 'radio_fast'
	hdlg = z0dpatience('fast');
        close_0plot;
	drawnow
	
	% information sur les donnees externes
	% information sur les donnees externes
	noms = fieldnames(getappdata(0));
    text_warn = '';
	for k = 1:length(noms)
		if findstr(noms{k},'_EXP')
			fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
			text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
		end
    end
    if isdeployed && ~isempty(text_warn)
        warndlg(text_warn,'Pay attention: external data used');
    end


        % pour eviter les incoherence sur le mode restart
        evalin('base','clear  z0dstruct');

        % overwrite parameters whith reference parameters if defined
	evalin('base','z0dinput.option = z0doverwriteparam(z0dinput.option);');

	% sort du mode evolution
	evalin('base','z0dinput.option.evolution = 0;');

	% recopie machine dans option
	evalin('base','zerod_machine_name;');
	        
        evalin('base','[post.zerod,void,post.profil0d] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);post.z0dinput = z0dinput;z0plotsc;', ...
		'errordlg(lasterr,''Error during the run of Metis simulator !'');');
	delete(hdlg)
        set(ho,'value',0);

case 'radio_hyb'
	hdlg = z0dpatience('fit');
        close_0plot;
	drawnow
	
	% information sur les donnees externes
	% information sur les donnees externes
	noms = fieldnames(getappdata(0));
    text_warn = '';
	for k = 1:length(noms)
		if findstr(noms{k},'_EXP')
			fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
			text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
		end
    end
    if isdeployed && ~isempty(text_warn)
        warndlg(text_warn,'Pay attention: external data used');
    end


        % pour eviter les incoherence sur le mode restart
        evalin('base','clear  z0dstruct');

        % overwrite parameters whith reference parameters if defined
	evalin('base','z0dinput.option = z0doverwriteparam(z0dinput.option);');

	% sort du mode evolution
	evalin('base','z0dinput.option.evolution = 0;');

	% recopie machine dans option
	evalin('base','zerod_machine_name;');
	        
	evalin('base','z0dfitlheta;', ...
		'errordlg(lasterr,''Error during the run of Metis simulator !'');');
	delete(hdlg)
        set(ho,'value',0);
        
case 'radio_evolution'	
 	hdlg = z0dpatience('full');
        close_0plot;
	drawnow

	% information sur les donnees externes
	% information sur les donnees externes
	noms = fieldnames(getappdata(0));
    text_warn = '';
	for k = 1:length(noms)
		if findstr(noms{k},'_EXP')
			fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
			text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
		end
    end
    if isdeployed && ~isempty(text_warn)
        warndlg(text_warn,'Pay attention: external data used');
    end


	% overwrite parameters whith reference parameters if defined
	evalin('base','z0dinput.option = z0doverwriteparam(z0dinput.option);');

	% sort du mode evolution
	evalin('base','z0dinput.option.evolution = 0;');

	% recopie machine dans option
	evalin('base','zerod_machine_name;');
	err_flag = 0;        
	evalin('base','zerodrunevolution;z0plotsc;','err_flag = 1');
	if err_flag == 1
		err_flag2 = 0;        
		evalin('base','post.z0dinput = z0dstruct.z0dinput;post.zerod = z0dstruct.zs;post.profil0d = z0dstruct.profil;','err_flag2 = 1;');
		if err_flag2 == 1
                       errordlg(lasterr,'Error during the run of Metis simulator with evolution mode!');
                       return
		else	
			warndlg('Simulation in evolution modeinterrupted by the user','Interrupted');
		end
	end
        set(ho,'value',0);	   

case 'radio_evolution_restart'	
 	hdlg = z0dpatience('full');
        close_0plot;
	drawnow

%  	% detecte si c'est possible
%  	err_flag = 0;        
%  	evalin('base','length(z0dstruct.zs.temps);','err_flag = 1;');
%          if err_flag == 1	
%  		set(ho,'value',0);
%  		errordlg('Restart mode works only with evolution run','Restart');
%  		return;
%          end
	% information sur les donnees externes
	% information sur les donnees externes
	noms = fieldnames(getappdata(0));
    text_warn = '';
	for k = 1:length(noms)
		if findstr(noms{k},'_EXP')
			fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
			text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
		end
    end
    if isdeployed && ~isempty(text_warn)
        warndlg(text_warn,'Pay attention: external data used');
    end


	% overwrite parameters whith reference parameters if defined
	evalin('base','z0dinput.option = z0doverwriteparam(z0dinput.option);');

	% sort du mode evolution
	evalin('base','z0dinput.option.evolution = 0;');

	% recopie machine dans option
	evalin('base','zerod_machine_name;');
	err_flag = 0;        
	evalin('base','zerodrunevolution_restart;z0plotsc;','err_flag = 1;');
	if err_flag == 1
		err_flag2 = 0;        
		evalin('base','post.z0dinput = z0dstruct.z0dinput;post.zerod = z0dstruct.zs;post.profil0d = z0dstruct.profil;','err_flag2 = 1;');
		if err_flag2 == 1
                       errordlg(lasterr,'Error during the run of Metis simulator with evolution mode!');
                        return
		else	
			warndlg('Simulation in evolution modeinterrupted by the user','Interrupted');
		end
	end
        set(ho,'value',0);	
        
case 'radio_workingpoint'	

 	hdlg = z0dpatience('full');
        close_0plot;
	drawnow

	evalin('base','z0working_point', ...
		'errordlg(lasterr,''Error during the run of Metis simulator !'');');
	delete(hdlg)
        set(ho,'value',0);	   

case 'radio_sc'
      evalin('base','z0plotsc;');
      set(ho,'value',0);

case 'radio_p' 
      evalin('base','z0plotp;');
      set(ho,'value',0);	   

case 'radio_e' 
      evalin('base','z0plote;');
      evalin('base','z0plot_bilan_energy;');
      set(ho,'value',0);	   

case 'radio_t' 
      evalin('base','z0plott;');
      set(ho,'value',0);	   

case 'radio_n' 
      evalin('base','z0plotn;');
      set(ho,'value',0);	
      
case 'radio_fusion_power' 
      evalin('base','compute_total_fusion_power(post);');
      set(ho,'value',0);	   

case 'radio_c' 
      evalin('base','z0plotc;');
      set(ho,'value',0);	   

case 'radio_eq' 
      evalin('base','z0ploteq;');
      evalin('base','z0plotstationnary');
      set(ho,'value',0);	   

case 'radio_lhe' 
      evalin('base','z0plotlh;');
      set(ho,'value',0);	   

case 'radio_conv' 
      evalin('base','z0plotconv;');
      liste_of_convergence_symbols
      set(ho,'value',0);	   

case 'radio_geo' 
      evalin('base','z0plotgeo;');
      set(ho,'value',0);	   

case 'radio_j' 
      evalin('base','z0plotj;');
      set(ho,'value',0);	   

case 'radio_ts' 
      evalin('base','z0dtsneutron;');
      set(ho,'value',0);

case 'radio_ddsts' 
      evalin('base','z0plot_sawtooth_info;');
      try
	  evalin('base','z0plotdds;');
      end
      set(ho,'value',0);	   

case 'radio_scenar'
      evalin('base','z0plotsc;');
      set(ho,'value',0);

case 'radio_prof'
      evalin('base','z0profview;');
      set(ho,'value',0);

case 'radio_2deq'

      movout = menu('What do you want ?','Create a movie of 2D plasma equilibrium and state', ...
                    'Display time evolution of 2D plasma equilibrium and state', ...
                    'Display simple 2D flux surfaces','Cancel');
      drawnow
      switch movout
      case 1
	  evalin('base','filename = z0dmovie(post);');
      case 2
	  evalin('base','z0dmovie(post);');
      case 3
          hwait = msgbox('processing the data','please wait');
	  evalin('base','z0plot2dequi(post)'); 
	  delete(hwait)
      end 
      set(ho,'value',0);

case 'radio_density'
      evalin('base','z0densview;');
      set(ho,'value',0);

case 'radio_cost'
      if ~isdeployed
	    evalin('base','disp(''---------------------------- '');[cost,info,scale,ref] = ztokcost(post);');
      else
	    warndlg('Cost function is not implemented in this version','Windows standalone version');
      end
      set(ho,'value',0);


case 'radio_er'
      evalin('base','z0ploter(post.zerod,post.z0dinput.geo,post.z0dinput.cons,post.z0dinput.option,post.profil0d);');
      set(ho,'value',0);

case 'radio_shine'
      evalin('base','z0plotshine;');
      set(ho,'value',0);

case 'radio_star'
      % for backward compatibilty
      evalin('base','z0plotnustar_fun(post);');
      set(ho,'value',0);

case 'radio_l2h'
      evalin('base','z0plotl2h;');
      set(ho,'value',0);

case 'radio_fluxpol'
      evalin('base','z0plotflux;');
      set(ho,'value',0);

case 'radio_lhacc'
      evalin('base','z0plotlhacc;');
      set(ho,'value',0);

case 'radio_nbijet'
      movout = menu('Select a graph','JET NBI (comparison to TRANSP and PENCIl data  only for JET shot)', ...
                    'Bilan ramp-up','JET sawteeth','JET Bolometry','Li beam','Cancel');
      drawnow
      switch movout
      case 1
	  evalin('base','z0plotnbi;');
      case 2
	  evalin('base','plot_bilan_ramp_up_jet;');
      case 3 
	  evalin('base','plot_st_jet;');
      case 4
      	  evalin('base','z0plot_jet_bolo;');
      case 5
       	  evalin('base','z0plot_jet_libeam;')    	  
      end
      set(ho,'value',0);

case 'radio_scaling'
      evalin('base','z0plot_test_scaling;');
      set(ho,'value',0);

case 'radio_gaz'
      evalin('base','metis_gaz;');
      set(ho,'value',0);

case 'radio_2pts'
      evalin('base','z0plot2points;');
      set(ho,'value',0);

case 'radio_plotrad'
      evalin('base','z0plotrad;');
      set(ho,'value',0);

case 'radio_digest'
      evalin('base','z0plot_shot_digest(post);');
      set(ho,'value',0);

case 'radio_coher0d1d'
      evalin('base','z0plotcoherence_0D_1D;');
      evalin('base','test_psi_lcfs_evol(post);');
      set(ho,'value',0);

case 'radio_divertor'
      evalin('base','z0plotdivertor;');
      set(ho,'value',0);

case 'radio_ramp'
      evalin('base','test_SOL_2points;');
      set(ho,'value',0);

case lower('radio_qlkANN_k')
     %evalin('base','z0qlkANN_kin_e_2018(post.z0dinput.option,post.zerod,post.profil0d);');
     evalin('base','zuicreefunform(''gui_qlkANN10D'',''qlkparams'',1,0,''[chie_loop,chii_loop,chii_neo,D_loop,V_loop] =  gui_qlkANN10D(qlkparams,post);'',''visible'');');
     set(ho,'value',0);

case lower('radio_qlkANN_k_wp')
      evalin('base','z0working_point_qlk_nn_10D;');
      set(ho,'value',0);
     
case lower('radio_qlkz_std_wp')
      evalin('base','z0working_point_qlkz_std;');
      set(ho,'value',0);
      
case  'radio_coils_currents'     
        try
            evalin('base','set_reactor_feeqs_path;zwaitfor(zuicreefunform(''compute_reactor_fbe_inverse_fast'',''option_feeqs_reactor'',1,0,''option_feeqs_reactor = compute_reactor_fbe_inverse_fast(option_feeqs_reactor);''),''visible'');')
        catch
            if isempty(which('start_up_Unix.m'))
                warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
            else
                errordlg(lasterr, 'Error calling FEEQS.M');
            end
        end
        set(ho,'value',0);

case 'radio_plot_coils_currents'
        try
            evalin('base','set_reactor_feeqs_path;plotfigure_reactor_inverse4metis;')
        catch
            if isempty(which('start_up_Unix.m'))
                warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
            else
                errordlg(lasterr, 'Error calling FEEQS.M');
            end
        end
        set(ho,'value',0);
        
case  'radio_eqdsk'     
        try
            evalin('base','set_reactor_feeqs_path;zwaitfor(zuicreefunform(''compute_eqdsk_feeqs_based'',''option_feeqs_eqdsk'',1,0,''option_feeqs_eqdsk = compute_eqdsk_feeqs_based(option_feeqs_eqdsk);''),''visible'');')
        catch
            if isempty(which('start_up_Unix.m'))
                warndlg('FEEQS.M is not available on this computer.', 'FEEQS.M');
            else
                errordlg(lasterr, 'Error calling FEEQS.M');
            end
        end
        set(ho,'value',0);
        
case 'radio_equi_gs'
      evalin('base','make_external_equi_with_feeqs(post);');
      set(ho,'value',0);
            
case 'radio_breakdown'
      evalin('base','z0taue_burnthrough(post.z0dinput.option,post.z0dinput.geo,post.z0dinput.cons,post.zerod,post.zerod.tauthl);');
      set(ho,'value',0);

case 'radio_dataplot0'
	% disp('appel de zdataplot');
	hdp = findobj(0,'type','figure','tag','zdataplot_metis');
	if ishandle(hdp)
		figure(hdp);
	else
		zdataplot_metis;
	end
	drawnow
        set(ho,'value',0);
 
case 'radio_sepa'
      rep = menu('LCFS source ? ','LCFS generator','Import from CREATE','Import from FEEQS','Cancel');	
      switch rep,
      case 1
            h=zuicreefunform('z0dsepanew2','sepa_option',1,0,'z0dinput=z0separatrix(z0dinput,sepa_option,1,1);');
      case 2
	    evalin('base','import_create_lcfs_data');
      case 3
	    evalin('base','import_feeqs_lcfs_data');
      end
      set(ho,'value',0);
case 'radio_export_cpos'
      evalin('base','metis_equi2d2matfile;');
      set(ho,'value',0); 	
case 'aide'

     disp(pwd)
     s = 0;
     t = '';
     if ispc
          % mettre ici la commande d'ouverture de l'aide pour windows
          winopen ('Howto_METIS_final.pdf');
     elseif ismac
          [s,t] = system(sprintf('open /Applications/Preview.app %s',fullfile(fileparts(which('zero1t')),'Howto_METIS_final.pdf')));
     elseif exist(fullfile(fileparts(which('zero1t')),'Howto_METIS_final.pdf'))
	if system('which acroread') == 0
	      [s,t] = system(sprintf('acroread %s &',fullfile(fileparts(which('zero1t')),'Howto_METIS_final.pdf')));
	elseif system('which okular') == 0
	      [s,t] = system(sprintf('okular %s &',fullfile(fileparts(which('zero1t')),'Howto_METIS_final.pdf')));
	elseif system('which evince') == 0
	      [s,t] = system(sprintf('evince %s &',fullfile(fileparts(which('zero1t')),'Howto_METIS_final.pdf')));
	end
     else
	if system('which acroread') == 0
	      [s,t] = system('acroread METIS_HELP.pdf &');
	elseif system('which okular') == 0
	      [s,t] = system('okular METIS_HELP.pdf &');
	elseif system('which evince') == 0
	      [s,t] = system('evince METIS_HELP.pdf &');
	end
      end
      if s~= 0
          disp(t);
      end
      set(ho,'value',0);

case 'fig2pub'
      hlf = findobj(0,'type','figure','HandleVisibility','on');
      item = {};
      hlfv = {};
      for lk = 1:length(hlf)
          if strmatch(get(hlf(lk),'tag'),{'zeroda','section_form_zerod'},'exact')
	      % fenetre METIS pas dans le menu
	  elseif ~isempty(get(hlf(lk),'name')) &&  ~isempty(get(hlf(lk),'tag'))
	      item{end+1} = sprintf('#%d title = %s (%s)',fix(double(hlf(lk))),get(hlf(lk),'name'),get(hlf(lk),'tag'));
	      hlfv{end+1} = hlf(lk);
	  elseif ~isempty(get(hlf(lk),'name')) 
	      item{end+1} = sprintf('#%d title = %s',fix(double(hlf(lk))),get(hlf(lk),'name'));
	      hlfv{end+1} = hlf(lk);
	  elseif ~isempty(get(hlf(lk),'tag'))
	      item{end+1} = sprintf('#%d  (%s)',fix(double(hlf(lk))),get(hlf(lk),'tag'));
	      hlfv{end+1} = hlf(lk);
	  else
	      item{end+1} = sprintf('#%d',fix(double(hlf(lk))));	  
	      hlfv{end+1} = hlf(lk);
	  end
      end
      if length(item) > 0
	K = menu('Select a figure number',item);
	if K > 0
	  layout_nice(hlfv{K});
	  drawnow
	end
      end
      set(ho,'value',0);

otherwise
      warning('action not taken into account')
      if strcmp(get(ho,'type'),'uicontrol') && strcmp(get(ho,'style'),'radio')
            set(ho,'value',0);
      end
end

try 
  if strcmp(get(gco,'type'),'uicontrol') && strcmp(get(gco,'style'),'radio')
	set(gco,'value',0');
  end
end
