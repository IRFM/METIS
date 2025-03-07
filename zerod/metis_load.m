function metis_load(filename)

data = load(filename,'post');
fullresult = 1;
if isempty(data) ||  (length(fieldnames(data)) == 0)
    data = load(filename,'z0dinput');
    fullresult = 0;
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

zassignin('base','z0dinput',data.z0dinput);


if fullresult == 0
   z0dinterfacetitle(filename);
   disp('Warning: this data set contain only METIS input');
   return
end
zassignin('base','post.zerod',data.zerod);
zassignin('base','post.z0dinput',data.z0dinput);
zassignin('base','post.profil0d',data.profil0d);

% donnees LUKE dans METIS
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
fprintf('METIS data file : %s loaded\n',filename);
z0dinterfacetitle(filename);

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
            %fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
            text_warn =sprintf('%susing external data in METIS for %s\n',text_warn,strtok(noms{k},'_'));
        end
    end   
    if  ~isempty(text_warn)
        warning(sprintf('Pay attention: external data from previous simulation in use\n%s',text_warn));
    end
end

% restrore sepa option if available
if isfield(data.z0dinput,'sepa_option')
  zassignin('base','sepa_option',data.z0dinput.sepa_option);
  if isdeployed
	evalin('base','sepa_option');
  end
end

