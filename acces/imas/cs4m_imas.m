% cette fonction prend les donnees de cronos pour en faire des donnees externes de METIS
% elle peut etre utilisee pour reprendre des donnees experimentales
function  cs4m_imas(choice)

if nargin == 0
	choice ='';
end

% preselectind database
% recover tokamak name and IMAS version
if ~isappdata(0,'UAL_TOKAMAK') || ~isappdata(0,'UAL_DATAVERSION') || ~isappdata(0,'UAL_USER')
    imasdb;
end

% choice for external data
if isempty(choice)
  % menu to select field
  itemliste = {'Ne','Te','Ti','LHCD','ECRH','ICRH','NBI','NBI2','PLINE','RUNAWAY','Momentum','Zeff','Clear cache','Cancel'};
  header    = sprintf('Copy of IMAS data in METIS external data structure\nWhich data or source must be copied ?');

  while(isempty(strmatch(choice,itemliste)) | isempty(choice))
	Button = menu(header,itemliste);
	
	if (Button == 0) || (Button == length(itemliste))
		return

        elseif Button == (length(itemliste) -  1)
            % clear cache
            ap0  = getappdata(0);
            noms = fieldnames(ap0);
            for k = 1:length(noms)
	      if ~isempty(strfind(noms{k},'IMAS_EXTERNAL_METIS_DATA_'))
		  rmappdata(0,noms{k});
	      end
            end
            if isappdata(0,'IMAS_REFERENCE4EXTERNAL_METIS_DATA_EQUILIBRIUM')
	       rmappdata(0,'IMAS_REFERENCE4EXTERNAL_METIS_DATA_EQUILIBRIUM');
            end
            if isappdata(0,'IMAS_REFERENCE4EXTERNAL_METIS_DATA')
	       rmappdata(0,'IMAS_REFERENCE4EXTERNAL_METIS_DATA');
            end
            return
	end
	choice = itemliste{Button};
  end
end


% open and read IDSs
if isappdata(0,'IMAS_REFERENCE4EXTERNAL_METIS_DATA')
	ref = getappdata(0,'IMAS_REFERENCE4EXTERNAL_METIS_DATA');
    run_num   = num2str(ref.run_num);
	ual_occ   = ref.ual_occ;
	ual_tok   = ref.ual_tok;
	ual_user  = ref.ual_user;
	ual_ver   = ref.ual_ver;
    shot_num  = num2str(ref.shot_num);

elseif evalin('base','exist(''post'',''var'')')
	post = evalin('base','post');

	if isfield(post.z0dinput,'run')
	  run_num = num2str(post.z0dinput.run + 1);
	else
	  run_num = '0';
	end
	if isfield(post.z0dinput,'occ')
	  ual_occ = post.z0dinput.occ;
	else
	  ual_occ = '';
	end
	if isfield(post.z0dinput,'machine')
	  ual_tok = post.z0dinput.machine;
    elseif  isappdata(0,'UAL_TOKAMAK') 
      ual_tok = getappdata(0,'UAL_TOKAMAK');
	else
	  ual_tok = '';
	end
	if isfield(post.z0dinput,'user')
	  ual_user = post.z0dinput.user;
    elseif isappdata(0,'UAL_USER')
      ual_user = getappdata(0,'UAL_USER');
	else
	  ual_user = getenv('USER');
	end
	if isfield(post.z0dinput,'dataversion')
	  ual_ver = post.z0dinput.dataversion;
    elseif isappdata(0,'UAL_DATAVERSION')
      ual_ver = getappdata(0,'UAL_DATAVERSION');
	else
	  ual_ver = '';
	end
    shot_num = num2str(post.z0dinput.shot);
else
	run_num = '0';
	ual_occ = '';
	if  isappdata(0,'UAL_TOKAMAK') 
      ual_tok = getappdata(0,'UAL_TOKAMAK');
	else
	  ual_tok = '';
    end
    if isappdata(0,'UAL_USER')
      ual_user = getappdata(0,'UAL_USER');
	else
	  ual_user = getenv('USER');
    end
    if isappdata(0,'UAL_DATAVERSION')
      ual_ver = getappdata(0,'UAL_DATAVERSION');
	else
	  ual_ver = '';
    end
    if evalin('base','exist(''z0dinput'',''var'')')
		z0dinput = evalin('base','z0dinput');
        shot_num = num2str(z0dinput.shot);
    else
        shot_num = '?';
    end
end
prompt={'shot number :','Run :' ,'Occurrence:','Tokamak:','User:','Version:'};
def={shot_num,run_num,ual_occ,ual_tok,ual_user,ual_ver};
dlgTitle='Reading data from UAL';
lineNo=1;
answer=zinputdlg(prompt,dlgTitle,lineNo,def);
drawnow
if isempty(answer)
	return
end
ref.shot_num   = str2num(answer{1});
ref.run_num    = str2num(answer{2});
ref.ual_occ    = strtrim(answer{3});
ref.ual_tok    = strtrim(answer{4});
ref.ual_user   = strtrim(answer{5});
ref.ual_ver    = strtrim(answer{6});

% check field filling
if isempty(ref.shot_num)
    disp('No shot number !');
    return
end
if  isempty(ref.run_num)
   ref.run_num = 0;
end
if isempty(ref.ual_tok) || isempty(ref.ual_user) || isempty(ref.ual_ver)
   % recover tokamak name and IMAS version 
   if ~isappdata(0,'UAL_TOKAMAK') || ~isappdata(0,'UAL_DATAVERSION') || ~isappdata(0,'UAL_USER')
       imasdb;
   end
   if isempty(ref.ual_tok)
       ref.ual_tok = getappdata(0,'UAL_TOKAMAK');
   end
   if isempty(ref.ual_user)
       ref.ual_user = getappdata(0,'UAL_USER');
   end
   if isempty(ref.ual_ver)
       ref.ual_ver = getappdata(0,'UAL_DATAVERSION');
   end
end

setappdata(0,'IMAS_REFERENCE4EXTERNAL_METIS_DATA',ref);

% time offset for WEST public database
if ~isempty(strmatch(ref.ual_user,{'imas_public'},'exact')) &&  ...
       ~isempty(strmatch(ref.ual_tok,{'west'},'exact'))
       timeoffset =  tsbase(ref.shot_num,'rignitron');
       if isempty(timeoffset)
	  timeoffset = 32;
       end
else
       timeoffset = 0;
end
% access to IDS
switch choice
case {'Ne','Te','Ti','Momentum','Zeff'}
    idss_list = {'core_profiles','equilibrium'};
case {'LHCD','ECRH','ICRH','NBI','NBI2','PLINE'}
    idss_list = {'core_sources','equilibrium'};
case 'RUNAWAY'
    idss_list = {'core_profiles','core_sources','equilibrium'};
otherwise 
  return
end
% search in cache
output = [];
idss_list_cache = idss_list;
idss_list = {};
for k = 1:length(idss_list_cache)
    nom_cache = sprintf('IMAS_EXTERNAL_METIS_DATA_%s',idss_list_cache{k});
    if isappdata(0,nom_cache);
	cache_data = getappdata(0,nom_cache);
	if same(cache_data.ref,ref)
	    fprintf('retreiving IDS %s from memory cache\n',idss_list_cache{k});
	    output.(idss_list_cache{k}) = cache_data.data;	    
	else
	   idss_list{end+1} = idss_list_cache{k};
	end
    else
	idss_list{end+1} = idss_list_cache{k};
    end
end

if ~isempty(idss_list)
  if ~isempty(ref.ual_occ) || ~isempty(ref.ual_tok) || ~isempty(ref.ual_user) ||  ~isempty(ref.ual_ver)
      output_db = litidss(idss_list,ref.shot_num,ref.run_num,ref.ual_user,ref.ual_tok,ref.ual_ver,ref.ual_occ);
  else
      output_db = litidss(idss_list,ref.shot_num,ref.run_num);
  end
  noms = fieldnames(output_db);
  for k=1:length(noms)
    nom_cache = sprintf('IMAS_EXTERNAL_METIS_DATA_%s',noms{k});
    cache_data.ref = ref;
    cache_data.data = output_db.(noms{k});
    setappdata(0,nom_cache,cache_data);
    output.(noms{k}) = output_db.(noms{k});
  end
end
% test if there is some data
switch choice
case {'Ne','Te','Ti','Momentum','Zeff'}
    if isempty(output.core_profiles.time) && isempty(output.core_profiles.profiles_1d{1}.time)
        disp('No core_profiles data available')
	return
    end
case {'LHCD','ECRH','ICRH','NBI','NBI2','PLINE','RUNAWAY'}
    if isempty(output.core_sources.time) && isempty(output.core_sources.profiles_1d{1}.time)
        disp('No core_sources data available')
	return
    end
otherwise 
  return
end

% test if data is available for equilibrium ; otherwise retry with new shot reference
if isempty(output.equilibrium.time) && (isempty(output.equilibrium.time_slice{1}.time) || (output.equilibrium.time_slice{1}.time == -9.0000e+40))
    prompt={'shot number :','Run :' ,'Occurrence:','Tokamak:','User:','Version:'};
    def={shot_num,run_num,ual_occ,ual_tok,ual_user,ual_ver};
    dlgTitle='Retry to read data from UAL for Equilibrium IDS';
    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    drawnow
    if isempty(answer)
	return
    end
    ref_equi.shot_num   = str2num(answer{1});
    ref_equi.run_num    = str2num(answer{2});
    ref_equi.ual_occ    = deblank(answer{3});
    ref_equi.ual_tok    = deblank(answer{4});
    ref_equi.ual_user   = deblank(answer{5});
    ref_equi.ual_ver    = deblank(answer{6});
    setappdata(0,'IMAS_REFERENCE4EXTERNAL_METIS_DATA_EQUILIBRIUM',ref_equi);
    nom_cache = sprintf('IMAS_EXTERNAL_METIS_DATA_%s','equilibrium');
    if isappdata(0,nom_cache);
    	cache_data = getappdata(0,nom_cache);
	if same(cache_data.ref,ref)
	    fprintf('retreiving IDS equilibrium from memory cache\n');
	    output.equilibrium = cache_data.data;	    
	elseif ~isempty(ref_equi.ual_occ) || ~isempty(ref_equi.ual_tok) || ~isempty(ref_equi.ual_user) ||  ~isempty(ref_equi.ual_ver)
	    output_db = litidss('equilibrium',ref.shot_num,ref.run_num,ref.ual_user,ref.ual_tok,ref.ual_ver,ref.ual_occ);
	    output.equilibrium = output_db.equilibrium;
	else
	    output_db = litidss('equilibrium',ref.shot_num,ref.run_num);
	    output.equilibrium = output_db.equilibrium;
	end
    end
end
clear output_db cache_data

% test if data is available for equilibrium
if isempty(output.equilibrium.time) && (isempty(output.equilibrium.time_slice{1}.time) || (output.equilibrium.time_slice{1}.time == -9.0000e+40))
    % ask fo METIS equilibrium in replacement
    if evalin('base','exist(''post'',''var'')')
        ButtonName = questdlg('Do you want use METIS equilibrium to map external data ?','Metis equilibrium','Yes', 'No','No');
        switch ButtonName
            case 'Yes'
                post              = evalin('base','post');
                time_equi         = post.profil0d.temps;
                amin              = interp1(post.z0dinput.cons.temps,post.z0dinput.geo.a,time_equi,'linear','extrap');
                rho_tor_norm_equi = cell(length(time_equi),1);
                amin_equi         = cell(length(time_equi),1);
                x_equi            = cell(length(time_equi),1);
                for k=1:length(time_equi)
                    rho_tor_norm_equi{k} = post.profil0d.rmx(k,:) ./ max(post.profil0d.rmx(k,:),[],2);
                    amin_equi{k}         = amin(k) * post.profil0d.xli;
                    x_equi{k}            = post.profil0d.xli;
                end
            otherwise
                disp('No equilibrium data available')
                return
        end
    else
        disp('No equilibrium data available')
        return
    end
else
    
    % coordinate translation from rho_tor to x METIS
    if output.equilibrium.ids_properties.homogeneous_time == 1
        time_equi	= output.equilibrium.time;
    else
        time_equi	= [];
    end
    if isfield(output.equilibrium,'time_slice')
        for k = 1:length(output.equilibrium.time_slice)
            if output.equilibrium.ids_properties.homogeneous_time ~= 1
                time_equi(k) = output.equilibrium.time_slice{k}.time;
            end
            rho_tor_norm_equi{k} = output.equilibrium.time_slice{k}.profiles_1d.rho_tor_norm;
            amin_equi{k}         = (output.equilibrium.time_slice{k}.profiles_1d.r_outboard - output.equilibrium.time_slice{k}.profiles_1d.r_inboard) ./ 2;
            x_equi{k}               =  amin_equi{k} ./ max(amin_equi{k});
        end
    else
        rho_tor_norm_equi = output.equilibrium.profiles_1d.rho_tor;
        rho_tor_norm_equi = rho_tor_norm_equi ./ (max(rho_tor_norm_equi,[],2) * ones(1,size(rho_tor_norm_equi,2)));
        amin_equi         = (output.equilibrium.profiles_1d.r_outboard - output.equilibrium.profiles_1d.r_inboard) ./ 2;
        x_equi            =  amin_equi ./ (max(amin_equi,[],2) *ones(1,size(amin_equi,2)));
    end
    % filtre sur le temps mauvais
    if iscell(rho_tor_norm_equi)
        indbad = [];
        for k=1:length(rho_tor_norm_equi)
            if all((rho_tor_norm_equi{k} == 0) | ~isfinite(rho_tor_norm_equi{k})) || (sum(double(rho_tor_norm_equi{k} == 0)) > 1)
                indbad(end+1) = k;
            end
        end
        if ~isempty(indbad)
            time_equi(indbad) = [];
            rho_tor_norm_equi{indbad} = [];
            amin_equi{indbad} = [];
            x_equi{indbad} = [];
        end
        
    else
        indbad = find(all((rho_tor_norm_equi == 0) | ~isfinite(rho_tor_norm_equi),2)| (sum(double(rho_tor_norm_equi == 0),2) > 1));
        if ~isempty(indbad)
            time_equi(indbad) = [];
            rho_tor_norm_equi(indbad,:) = [];
            amin_equi(indbad,:) = [];
            x_equi(indbad,:) = [];
        end
    end
end

switch choice 
case 'Ne'
	% 1- Electron density 
	[validity,data_out,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_profiles,{'electrons','density'});
	if validity == 1
	  NE_EXP.x     =  x_out;
	  NE_EXP.ne    =  data_out;
	  NE_EXP.temps =  time_out - timeoffset;
	  setappdata(0,'NE_EXP',NE_EXP);
	else
		error('CS4M @ Ne : no valid data');	
	end
	
case 'Te'
	% 1- Electron temperature 
	[validity,data_out,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_profiles,{'electrons','temperature'});
	if validity == 1
	  TE_EXP.x     =  x_out;
	  TE_EXP.te    =  data_out;
	  TE_EXP.temps =  time_out - timeoffset;
	  setappdata(0,'TE_EXP',TE_EXP);
	else
		error('CS4M @ Te : no valid data');
	end

case 'Ti'
	% 1- Ion temperature
	[validity,data_out,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_profiles,'t_i_average');
	if validity == 1
	  TI_EXP.x     =  x_out;
	  TI_EXP.ti    =  data_out;
	  TI_EXP.temps =  time_out - timeoffset;
	  setappdata(0,'TI_EXP',TI_EXP);
	else
		error('CS4M @ Ti : no valid data');
	end
	
case 'Momentum'
	[validity,data_out,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_profiles,'momentum_tor');
	if validity == 1
	  [validity,r2i] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.equilibrium,{'profiles_1d','gm1'});
	  % Mtor is proportionnal to density with a factor that is quite constant in radial direction
	  [validity,density] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_profiles,{'electrons','density'});
	  %rtor = Mtor .* profli.omega ./ profli.r2i;
	  omega = data_out ./ (density ./ max(density(:))) .* r2i;
	  omega = omega ./ max(eps,max(abs(omega(:))));
	  VTOR_SHAPE_EXP.x     =  x_out;
	  VTOR_SHAPE_EXP.omega =  omega;
	  VTOR_SHAPE_EXP.temps =  time_out - timeoffset;
	  setappdata(0,'VTOR_SHAPE_EXP',VTOR_SHAPE_EXP);
	else
		error('CS4M @ Momentum : no valid data');
	end

case 'Zeff'	
	% 1- Ion temperature
	[validity,data_out,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_profiles,'zeff');
	if validity == 1
	  ZEFF_SHAPE_EXP.x     =  x_out;
	  ZEFF_SHAPE_EXP.zeff  =  data_out;
	  ZEFF_SHAPE_EXP.temps =  time_out - timeoffset;
	  setappdata(0,'ZEFF_SHAPE_EXP',ZEFF_SHAPE_EXP);
	else
		error('CS4M @ Zeff : no valid data');
	end	
	
case 'LHCD'
	% 4- LHCD shape (scale on 0D data)
	[validity,plh,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,{'electrons','energy'},{'LH','LHCD'});
        if validity == 1
	  [validity,jlhcd,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,'j_parallel',{'LH','LHCD'});
	  LH_SHAPE_EXP.x     =  x_out;
	  LH_SHAPE_EXP.jlh   =  jlhcd;
	  LH_SHAPE_EXP.plh   =  plh;
	  LH_SHAPE_EXP.temps =  time_out - timeoffset;
	  setappdata(0,'LH_SHAPE_EXP',LH_SHAPE_EXP);
	else
		error('CS4M @ LHCD : no valid data');
	end
	
case 'ECRH'
	[validity,peccd,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,{'electrons','energy'},{'EC','ECRH','ECCD'});
        if validity == 1
	  [validity,jeccd,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,'j_parallel',{'EC','ECRH','ECCD'});
	% 4- LHCD shape (scale on 0D data)
	ECCD_SHAPE_EXP.x       =  x_out;
	ECCD_SHAPE_EXP.jeccd   =  jeccd;
	ECCD_SHAPE_EXP.peccd   =  peccd;
	ECCD_SHAPE_EXP.temps   =  time_out - timeoffset;
	setappdata(0,'ECCD_SHAPE_EXP',ECCD_SHAPE_EXP);
	else
		error('CS4M @ ECRH : no valid data');
	end
	
case 'ICRH'
	% 5- ICRH shape (scale on 0D data)
	[validity,pel,x_out,time_out]    = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,{'electrons','energy'},{'IC','ICRH','FW','FWCD'});
        if validity == 1
	  [validity_js,js,x_out,time_out]   = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,'j_parallel',{'IC','ICRH','FW','FWCD'});
	  [validity_ion,pion,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,'total_ion_energy',{'IC','ICRH','FW','FWCD'});
	  ICRH_SHAPE_EXP.x     =  x_out;
	  if validity_js  == 1
	      ICRH_SHAPE_EXP.jfwcd   = js;
	  else
	      ICRH_SHAPE_EXP.jfwcd   =  zeros(size(pel));
	  end
	  if validity_ion == 1
	    ICRH_SHAPE_EXP.pfw   =  zeros(size(pel));
	    ICRH_SHAPE_EXP.pel   =  pel;
	    ICRH_SHAPE_EXP.pion   =  pion;
	  else
	    ICRH_SHAPE_EXP.pel   =  zeros(size(pel));
	    ICRH_SHAPE_EXP.pion   =  zeros(size(pel));
	    ICRH_SHAPE_EXP.pfw   =  pel;
	  end
	  ICRH_SHAPE_EXP.temps =  time_out - timeoffset;
	  setappdata(0,'ICRH_SHAPE_EXP',ICRH_SHAPE_EXP);
	else
		error('CS4M @ ICRH : no valid data')
	end
	
case 'NBI'
	% 6- NBI shape (scale on 0D data)
	[validity,pel,x_out,time_out]    = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,{'electrons','energy'},{'NBI','NBICD'});
        if validity == 1
	  [validity_js,js,x_out,time_out]   = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,'j_parallel',{'NBI','NBICD'});
	  [validity_ion,pion,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,'total_ion_energy',{'NBI','NBICD'});
	  NBICD_SHAPE_EXP.x     =  x_out;
	  if validity_js  == 1
	      NBICD_SHAPE_EXP.jnbicd   = js;
	  else
	      NBICD_SHAPE_EXP.jnbicd   =  zeros(size(pel));
	  end
	  NBICD_SHAPE_EXP.pel     =  pel;
	  NBICD_SHAPE_EXP.pion      =  pion;
	  NBICD_SHAPE_EXP.temps     =  time_out - timeoffset;
	  setappdata(0,'NBICD_SHAPE_EXP',NBICD_SHAPE_EXP);
	else
		error('CS4M @ NBICD : no valid data')
	end
	
	
case 'NBI2'
      error('CS4M @ NBI2 : not yet implemented')

case 'PLINE'
	% 7- line radiation
	[validity,prad,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,{'electrons','energy'},{'LINE','LINERADIATION','PLINE','PRAD'});
        if validity == 1
	  PLINE_EXP.x         =  x_out;
	  PLINE_EXP.prad      =  abs(prad);
	  PLINE_EXP.temps     =  time_out - timeoffset;
	  setappdata(0,'PLINE_EXP',PLINE_EXP)
	else
		error('CS4M @ PLINE : no valid data')
	end
	
case 'RUNAWAY'
        % Runaway electrons current
	[validity,jrun,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_sources,'j_parallel',{'RUNAWAY','RUNAWAYS'});
        if validity == 1
	  [validity,rho_tor] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.equilibrium,{'profiles_1d','rho_tor'});
	  [validity,spr] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.equilibrium,{'profiles_1d','darea_drho_tor'});
	  [validity,johm] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,output.core_profiles,'j_ohmic');
	  for k=1:length(time_out)
	      irun(k) = trapz(rho_tor(k,:),spr(k,:) .* jrun(k,:),2);
	      iohm(k) = trapz(rho_tor(k,:),spr(k,:) .* johm(k,:),2);
	  end
	  RUNAWAY_EXP.temps  = time_out - timeoffset;
	  RUNAWAY_EXP.x      = x_out;
	  RUNAWAY_EXP.jrun   = jrun;
	  RUNAWAY_EXP.irun   = irun;
	  RUNAWAY_EXP.iohm   = iohm; 
	  RUNAWAY_EXP.ip     = interp1(output.core_profiles.time,output.core_profiles.global_quantities.ip,time_out,'linear','extrap');
	  setappdata(0,'RUNAWAY_EXP',RUNAWAY_EXP)
	else
		error('CS4M @ Runaway electrons : no valid data')
	end
	
otherwise
	% cas cancel
	return
end

function [validity,data_out,x_out,time_in] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,ids_in,nom_in,source_names)

% recover time slices vector
if ids_in.ids_properties.homogeneous_time == 1
    time_in = ids_in.time;
    if  ~isempty(time_in)
	homogene_time = true;
    else
	homogene_time = false;
    end
else
    homogene_time = false;
end
% comment for graph
comment = ids_in.ids_properties.comment;

% for core_sources selection of the source
index_source = NaN;
if nargin > 5
  for k=1:length(ids_in.source)
    name = deblank(upper(ids_in.source{k}.identifier.name));
    if ~isempty(name)
      if ~isempty(strmatch(name,source_names,'exact'))
	  index_source = k;
	  fprintf('Selected source index = %d\n',index_source);
	  comment = sprintf('%s (source = %s @ %d)',comment,name,index_source);
	  break
      end
    end
  end

  if isfinite(index_source)
    ids_in = ids_in.source{index_source};
  else
    error(sprintf('source %s not found in core_sources IDS',source_names));
  end
  
end

% output status
validity = 1;

% recover time slices vector
if ~homogene_time
  if isfield(ids_in,'profiles_1d')
    for k=1:length(ids_in.profiles_1d)
      time_in(k) = ids_in.profiles_1d{k}.time;
    end
  else
    for k=1:length(ids_in.time_slice)
      time_in(k) = ids_in.time_slice{k}.time;
    end  
  end
end

% check if all equilibrium 1d profiles have the same length
timechange = 0;
if iscell(x_equi)
  ll = length(x_equi{1});
  x_homogene = 1;
  for k = 2:length(time_equi)
    if (length(x_equi{k}) ~= ll) && all(x_equi{k} == x_equi{1})
	x_homogene	= 0;
    end
  end
else
  x_homogene = 1;
  timechange = 1;
end
if x_homogene == 1
  if iscell(x_equi)
    x_out = x_equi{1};
  else
    x_out = x_equi(1,:);
  end
else
  x_out = linspace(0,1,21);
end
data_out = NaN .* ones(length(time_in),length(x_out));
% change of coordinate
if timechange == 1
    for k = 1:length(time_in)
        if iscell(time_in)
            indice_equi = find(time_equi >= time_in{k},1);
        else
            indice_equi = find(time_equi >= time_in(k),1);
        end
        if isempty(indice_equi)
            indice_equi = length(time_equi);
        end
        if isfield(ids_in,'profiles_1d')
            loc_data =  ids_in.profiles_1d;
        elseif isfield(ids_in,'time_slice')
            loc_data =  ids_in.time_slice;
        else
            error('this data structure is not yet described in this function');
        end
        if ~iscell(loc_data)
            if isfield(loc_data,'grid')
                loc_rho = loc_data.grid.rho_tor_norm(k,:);
                loc_rho = loc_rho(loc_rho >= 0);
                x_local  = interp1(rho_tor_norm_equi(indice_equi,:),x_equi(indice_equi,:),loc_rho,'pchip','extrap');
            else
                loc_rho = loc_data.profiles_1d.rho_tor_norm(k,:);
                loc_rho = loc_rho(loc_rho >= 0);
                x_local  = interp1(rho_tor_norm_equi(indice_equi,:),x_equi(indice_equi,:),loc_rho,'pchip','extrap');
            end
        else
            if isfield(loc_data{1},'grid')
                x_local  = interp1(rho_tor_norm_equi(indice_equi,:),x_equi(indice_equi,:),loc_data{k}.grid.rho_tor_norm(loc_data{k}.grid.rho_tor_norm >= 0),'pchip','extrap');
            else
                x_local  = interp1(rho_tor_norm_equi(indice_equi,:),x_equi(indice_equi,:),loc_data{k}.profiles_1d.rho_tor_norm(loc_data{k}.profiles_1d.rho_tor_norm >= 0),'pchip','extrap');
            end
        end
        if ischar(nom_in)
            data_local = loc_data{k}.(nom_in);
        else
            data_local = loc_data{k}.(nom_in{1});
            validity = 0;
            for l=2:length(nom_in)
                nom_validity = sprintf('%s_validity',nom_in{l});
                if isfield(data_local,nom_validity)
                    validity = data_local.(nom_validity);
                    if validity == -999999999
                        validity = 0;
                    end
                end
                data_local = data_local.(nom_in{l});
            end
            if any(validity < 0) && isnumeric(data_local)
                data_local(:) = NaN;
            end
        end
        try
            if all(isfinite(data_local(:)))
                data_out(k,:) = interp1(x_local,data_local,x_out,'pchip','extrap');
            else
                if iscell(time_in)
                    fprintf('invalid data at t = %g s\n',time_in{k});
                else
                    fprintf('invalid data at t = %g s\n',time_in(k));
                end
            end
        catch
            if iscell(time_in)
                fprintf('invalid data at t = %g s\n',time_in{k});
            else
                fprintf('invalid data at t = %g s\n',time_in(k));
            end
        end
    end
else
    for k = 1:length(time_in)
        if iscell(time_in)
            indice_equi = find(time_equi >= time_in{k},1);
        else
            indice_equi = find(time_equi >= time_in(k),1);
        end
        if isfield(ids_in,'profiles_1d')
            loc_data =  ids_in.profiles_1d{k};
        elseif isfield(ids_in,'time_slice')
            loc_data =  ids_in.time_slice{k};
        else
            error('this data structure is not yet described in this function');
        end
        if isfield(loc_data,'grid')
            x_local  = interp1(rho_tor_norm_equi{indice_equi},x_equi{indice_equi},loc_data.grid.rho_tor_norm(loc_data.grid.rho_tor_norm >= 0),'pchip','extrap');
        else
            x_local  = interp1(rho_tor_norm_equi{indice_equi},x_equi{indice_equi},loc_data.profiles_1d.rho_tor_norm(loc_data.profiles_1d.rho_tor_norm >= 0),'pchip','extrap');
        end
        if ischar(nom_in)
            data_local = loc_data.(nom_in);
        else
            validity = 0;
            data_local = loc_data.(nom_in{1});
            for l=2:length(nom_in)
                nom_validity = sprintf('%s_validity',nom_in{l});
                if isfield(data_local,nom_validity)
                    validity = data_local.(nom_validity);
                    if validity == -999999999
                        validity = 0;
                    end
                end
                data_local = data_local.(nom_in{l});
            end
            if any(validity < 0) && isnumeric(data_local)
                data_local(:) = NaN;
            end
        end
        try
            if all(isfinite(data_local(:)))
                data_out(k,:) = interp1(x_local,data_local,x_out,'pchip','extrap');
            else
                if iscell(time_in)
                    fprintf('invalid data at t = %g s\n',time_in{k});
                else
                    fprintf('invalid data at t = %g s\n',time_in(k));
                end
            end
        catch
            if iscell(time_in)
                fprintf('invalid data at t = %g s\n',time_in{k});
            else
                fprintf('invalid data at t = %g s\n',time_in(k));
            end
        end
    end
end
if all(~isfinite(data_out(:)))
  validity = 0;
else
  indbad = find(all(~isfinite(data_out),2));
  if ~isempty(indbad) 
    data_out(indbad,:) = [];
    time_in(indbad)    = [];
  end
  validity = 1;
end
figure;
if length(time_in) == 1
  plot(x_out,data_out);
else
  zplotprof(gca,time_in,x_out,data_out);
end
xlabel('rho_{tor_norm}')
if iscell(nom_in)
  ylabel(sprintf('%s/%s',nom_in{end-1},nom_in{end}))
else
  ylabel(nom_in)
end  
title(comment);
