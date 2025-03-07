% fillCPOs4jet : reads JET data and coreprof CPO to make a METIS simulation using
%-------------------------------------------------------------
% fonction Matlab 12
%
%
% This function reads JET data and coreprof CPO to make a METIS simulation using experimental
% enegiecontents and, for time slices where data are available, fit of experimental Te, Ti and Ne.
%
% syntax : 
%
% 	fillCPOs4jet(shot,run_in,run_out,user,data_version,parameters)
%
%
% input :
%
%
%     shot                    = shot number (integer >0)
%
%     run_in                  = run number corresponding to the wanted coreprof data set(integer >0)
%				
%     run_out                 = run number in which METIS results is stored (integer >0)
%				
%     user                    = optional,UAL user database name for input coreprof data set
%
%     dataversion             = optional,select a different version of data than default
%
%     parameters              = optional, name and path to the file containing METIS user defined set of parameters.
%                               
%
% output :
%
%     CPOs reccords  and two files: logfile & METIS data file
%
% test:
% fillCPOs4jet(73344,13,117,'thia','4.10b')


% fonction ecrite par J-F Artaud
% version SVN (created the 07/12/2015)
%-----------------------------------------------------------------------
%
function fillCPOs4jet(shot,run_in,run_out,user,data_version,parameters)

% filename for diary and METIS output file 
filename_out = sprintf('fillCPOs4jet@%dfrom%dto%d',shot,run_in,run_out);
fprintf('results will be stored in %s.log and %s.mat files\n',filename_out,filename_out);
diary(strcat(filename_out,'.log'));
diary on

% add mds+ for Matlab
if isempty(which('mdsconnect'))
  % add access to JET data
  MDSPLUS_DIR = getenv('MDSPLUS_DIR');
  MDSPLUS_DIR ='/afs/ipp-garching.mpg.de/itm/switm/mdsplus/4.0_helsinki/mdsplus';
  rep = dir(fullfile(MDSPLUS_DIR,'java/classes/*.jar'));
  for k=1:length(rep)
      javaaddpath(fullfile(MDSPLUS_DIR,'java','classes',rep(k).name),'-end');
  end
  addpath(fullfile(MDSPLUS_DIR,'matlab'));
end

% prepare METIS data for JET
evalin('base',sprintf('z0dinput = zerod_init(2,%d,1);',shot));
% set parameters for experimental fit
evalin('base','z0dinput.option.scaling = 4;');
evalin('base','z0dinput.option.xiioxie = 3;');

% set optionnal file for paramaters
if nargin > 5
    evalin('base',sprintf('z0dinput.option.reference_parameters=''%s'';',parameters));
end

% run metis 1 for corrdinates change, is simpler than using equilibrium  datan
% information sur les donnees externes
noms = fieldnames(getappdata(0));
for k = 1:length(noms)
	if findstr(noms{k},'_EXP')
		fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
	end
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

% memorisation of the run data
evalin('base','jeux1.post = post;')



% acces au donnees UAL
%[output,cpos_list] = litcpos({'equilibrium','coreprof'},shot,run_in,user,'jet',data_version);
[output,cpos_list] = litcpos({'coreprof'},shot,run_in,user,'jet',data_version);

% create external data
post     = evalin('base','post');
t0d      = post.profil0d.temps;
x0d      = post.profil0d.xli;
rhon_0d  = post.profil0d.rmx ./ (max(post.profil0d.rmx,[],2) * ones(size(x0d)));
% remind to change raidial coordinate to Lao one.
% 1- Electron density 
[out,tout] = prep_ext(t0d,x0d,rhon_0d,post.profil0d.nep,output.coreprof.time,output.coreprof.rho_tor_norm, ...
               output.coreprof.ne.value);
ylabel('N_e');

if ~isempty(out)	       
	NE_EXP.temps =  t0d;
 	NE_EXP.x     =  x0d;
 	NE_EXP.ne    =  out;
	setappdata(0,'NE_EXP',NE_EXP)
	fprintf('Using experimental data for N_e from t= %g to t= %g \n',min(tout),max(tout));
else
	disp('no valid measurement for Ne')
end

% 2- Electron temperature 
[out,tout] = prep_ext(t0d,x0d,rhon_0d,post.profil0d.tep,output.coreprof.time,output.coreprof.rho_tor_norm, ...
               output.coreprof.te.value);
ylabel('T_e');
if ~isempty(out)	       
 	TE_EXP.temps =  t0d;
 	TE_EXP.x     =  x0d;
 	TE_EXP.te    =  out;
	setappdata(0,'TE_EXP',TE_EXP)
	fprintf('Using experimental data for T_e from t= %g to t= %g \n',min(tout),max(tout));
else
	disp('no valid measurement for Te')
end

%
% 3- ion temperature
[out,tout] = prep_ext(t0d,x0d,rhon_0d,post.profil0d.tip,output.coreprof.time,output.coreprof.rho_tor_norm, ...
               output.coreprof.ti.value);
ylabel('T_i');
if ~isempty(out)	       
 	TI_EXP.temps =  t0d;
 	TI_EXP.x     =  x0d;
 	TI_EXP.ti    =  out;
	setappdata(0,'TI_EXP',TI_EXP)
	fprintf('Using experimental data for T_i from t= %g to t= %g \n',min(tout),max(tout));
else
	disp('no valid measurement for Ti')
end

% 10- Line radiation
% 	PLINE_EXP.temps =  time vector [N,1]
% 	PLINE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	PLINE_EXP.prad  =  density of line radiative power [N,M] in W/m^3
%	setappdata(0,'PLINE_EXP',PLINE_EXP)
%
% 11- Toroidal rotation 
% 	VTOR_SHAPE_EXP.temps =  time vector [N,1]
% 	VTOR_SHAPE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	VTOR_SHAPE_EXP.omega =  profile shape of solid like plasma toroidal rotation [N,M] if possible in turn/s
%	setappdata(0,'VTOR_SHAPE_EXP',VTOR_SHAPE_EXP)
%
% 12- Effective charge profile shape (the line averaged value is given by the reference) 
% 	ZEFF_SHAPE_EXP.temps =  time vector [N,1]
% 	ZEFF_SHAPE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	ZEFF_SHAPE_EXP.zeff  =  profile shape of effective charge without He ashes coming from fusion reaction and without W contribution
%	setappdata(0,'ZEFF_SHAPE_EXP',ZEFF_SHAPE_EXP)
%

% run metis 2 with eternal data coming from measurments
% information sur les donnees externes
noms = fieldnames(getappdata(0));
for k = 1:length(noms)
	if findstr(noms{k},'_EXP')
		fprintf('using external data in METIS for %s\n',strtok(noms{k},'_'));
	end
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


% write outpu cpo from metis
evalin('base',sprintf('error_flag = metis4itm(%d,%d,'''',post);',shot,run_out));


% report
evalin('base',sprintf('metis_save %s',filename_out));
diary off

function [out,tout] = prep_ext(t0d,x0d,rhon_0d,in_0d,time,rho_norm,in_cpo)

% valide time slice
ind_ok = find(all((in_cpo >= -9.0000e+38) & (rho_norm >= -9.0000e+38),2));
if isempty(ind_ok)
	out = [];
	tout = [];
	return
else
	out = in_0d;
end
% time slices replaced by measurement
ind_rep = find((t0d>= min(time(ind_ok))) & (t0d <= max(time(ind_ok))));
% time resampling
inter     = interp1(time(ind_ok),in_cpo(ind_ok,:),t0d(ind_rep),'nearest','extrap');
rho_inter = interp1(time(ind_ok),rho_norm(ind_ok,:),t0d(ind_rep),'nearest','extrap');
% resampling on space
for k=1:length(ind_rep)
  out(ind_rep(k),:) = pchip(rho_inter(k,:),inter(k,:),rhon_0d(ind_rep(k),:)); 
end
%out(ind_rep,:) = griddata(time(ind_ok)*ones(1,size(rho_norm,2)),rho_norm(ind_ok,:),in_cpo(ind_ok,:), ...
%                          t0d(ind_rep)*ones(1,size(rhon_0d,2)),rhon_0d(ind_rep,:),'linear');

%figure;plot(rhon_0d(ind_rep,:)',in_0d(ind_rep,:)','b',rhon_0d(ind_rep,:)',out(ind_rep,:)','or',rho_norm(ind_ok,:)',in_cpo(ind_ok,:)','k');drawnow
figure
zplotprof(gca,t0d(ind_rep),rhon_0d(ind_rep,:),in_0d(ind_rep,:),'color','b');
zplotprof(gca,t0d(ind_rep),rhon_0d(ind_rep,:),out(ind_rep,:),'color','r');
zplotprof(gca,time(ind_ok),rho_norm(ind_ok,:),in_cpo(ind_ok,:),'color','k','linestyle',':');
xlabel('rho_tor_norm')
legend('METIS first run','METIS external data','Coreprof');
drawnow
 
tout = t0d(ind_rep);


