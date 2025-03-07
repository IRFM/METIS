% this function make auto run of METIS for  WEST shot with predefined
% parameters using external data for Te, Ti and Ne.
% the energy content is provided by Wdia mesurement and current drive
% efficiency of LHCD is adapted to match poloidal flux consumption.
%
% sybtax:
%
% in command line:
%
% 	[z0dinput,post,cp,equi] =
% 	autoWESTmetis(shot,gas,nhonhpnd,frhe0,density_shape,fcupper,time,filename,cp);
%
% with GUI:
%
%  [z0dinput,post,cp,equi] = autoWESTmetis;
%
% input:
%	shot = WEST shot number (> 50000)
%
%       gas  = 1 = D & 2 = He;
%              If empty or missing switch to automatic selection.
%
%       nhonhpnd   = content in hydrogen of the  plasma (n_H / (n_H  + n_D).
%                    if empty, set to 0.05.
%
%       frhe0      = fraction of He in plasam (if empty, set to 0)
%
%       density_shape = source of density profile ('tanh' or 'NICE')
%
%       fcupper       = fraction of maximum cupper density identified from spectroscopy
%
%       time = time slices for calculation.
%               If isempty or missing, use default time slices vector of METIS.
%               If is saclar, compute steady state solution for this selected time.
%
%	    filename = optionnal file name in which METIS data are saved (without .mat).
%
%       cp, equi  = optionnal precomputed core_profiles data and equilibrium data
%
% output:
%   z0dinput = METIS input data structure
%
%	post     = standard METIS data structure
%   
%   cp       = completed core_profiles IDS
%
%   equi     = time swapperd equilibrium IDS
%
%
function [z0dinput,post,cp,equi] = autoWESTmetis2(autoWESTmetis_option,time,filename,cp,equi)

persistent tabmat
% nouvelle donnees ADAS
if isempty(tabmat)
    try
        load('Lz_zave.mat')
    end
end

% default outputs
z0dinput    = [];
post        = [];
if nargin < 2
    time = []; % METIS select time slices
elseif length(time) == 1
    workingpoint = time;
    time = [];
end
if nargin < 3
    filename  = ''; % no nackup of METIS result
end

if nargin < 4
  cp          = [];
end
if nargin < 5
  equi        = [];
end
% reset external data
rm_cs4m;

% inputs handling    
workingpoint = NaN;
if nargin < 1
    % test access to tools_dc
    if exist('run_gui_dc') ~= 2
        warndlg('WEST data access requires IMAS framework and tools_dc functions. You must load the module tools_dc before launching METIS. This feature is not yet available under Windows or Mac OS/X','Access to data not avialable');
        z0dinput = [];
        return
    end
    try
        numshot = fix(evalin('base','param.from.shot.num'));
    catch
        try
            numshot = tsdernier_choc;
            rignitron = tsbase(numshot,'rignitron');
            if isempty(rignitron)
                numshot = numshot - 1;
            end
        catch
            numshot = 54178;
        end
    end
    if isempty(numshot)
      numshot = 54178;
    end
    %
    z0dinput = gui_compute_energy_west(numshot);
    % set specific option for METIS
    z0dinput.valeur.metis_mode = 'off';
    z0dinput.valeur.Wshape = 'density';
    z0dinput.valeur.rescale_te_with_wmhd = 'mixed';

    % end of GUI creation phase
    return;
    
else
    % backward compatibility
    shot      = abs(autoWESTmetis_option.shot);
    gas       = autoWESTmetis_option.gas_type;
    nhonhpnd  = autoWESTmetis_option.nHonZ1;
    frhe0     = autoWESTmetis_option.fr_he;
    density_shape     = autoWESTmetis_option.density_source;
    fcupper   = autoWESTmetis_option.fcupper;
    
    if ~isfield(autoWESTmetis_option,'no_cp_computation')
        autoWESTmetis_option.no_cp_computation = false;
    end

    % compatibility with tools_dc
    if ~autoWESTmetis_option.no_cp_computation
        zassignin('base','compute_energy_west_option',autoWESTmetis_option);
    end
end

% try to use equilibrium time slices by default
equi = imas_west_get(shot,'equilibrium',0,1);
if  ~isempty(equi) &&  length(equi.time) >= 11
    t_igni = tsbase(shot,'rignitron');
    if ~isempty(t_igni)
        time   = equi.time(:) - t_igni;
        disp('Using time slices from equilibrium IDS');
    end
end

% initialise METIS data
fprintf('Preparing METIS data set for WEST shot %d\n',shot);
try
    switch autoWESTmetis_option.metis_full
        case 'on'
             z0dinput = zerod_init(12,-abs(real(autoWESTmetis_option.shot)),gas,time);
       otherwise
            z0dinput = zerod_init(12,autoWESTmetis_option.shot,gas,time);
    end
catch
    z0dinput = [];
end
if isempty(z0dinput)
    disp('No data for this shot');
    return
end

% set METIS parameters
z0dinput.option = change_option_west(z0dinput.option);

% use Hard-X ray inversion if is available
if ~isfield(autoWESTmetis_option,'force_metis_lh_model')
    autoWESTmetis_option.force_metis_lh_model = 'off';
end
switch autoWESTmetis_option.force_metis_lh_model
    case 'on'
        disp('Internal METIS model for LHCD source shape is forced');
    otherwise
        if isfield(z0dinput.exp0d,'XDURt') && isfield(z0dinput.exp0d,'XDURx') && isfield(z0dinput.exp0d,'XDURv')
            indok            = find(all(z0dinput.exp0d.XDURv>=0,2));
            if length(indok) >3
                fprintf('dlh = %g & wlh = %g\n',z0dinput.option.dlh,z0dinput.option.wlh);
                z0dinput.option.dlh = 0;
                z0dinput.option.wlh = 0;
                disp('Switching to Hard-X ray inversion profile for LHCD source shape');
            end
        end
end

% initialising  external data
if nargin < 5 || isempty(cp)|| isempty(equi)
     [cp,equi,zeff_res_data] = compute_energy_west(shot,autoWESTmetis_option,autoWESTmetis_option.density_source);
elseif ~autoWESTmetis_option.no_cp_computation
     [cp,equi,zeff_res_data] = compute_energy_west(shot,autoWESTmetis_option,autoWESTmetis_option.density_source,cp);
else
     try
         zeff_res_data = evalin('base','zeff_res_data');
     catch
         zeff_res_data = [];
     end
     if isempty(zeff_res_data)
         keyboard
         [~,~,zeff_res_data] = compute_energy_west(shot,autoWESTmetis_option,autoWESTmetis_option.density_source,cp);
     end
end
if isempty(cp)
   error('No IMAS data available')
else
   zassignin('base','core_profiles',cp); 
end

% reading rignitron for time offset
timeoffset = tsbase(shot,'rignitron')

% make data from equilibrium 
rho_tor_norm_equi = equi.profiles_1d.rho_tor;
rho_tor_norm_equi = rho_tor_norm_equi ./ (max(rho_tor_norm_equi,[],2) * ones(1,size(rho_tor_norm_equi,2)));
amin_equi         = (equi.profiles_1d.r_outboard - equi.profiles_1d.r_inboard) ./ 2;
x_equi            =  amin_equi ./ (max(amin_equi,[],2) *ones(1,size(amin_equi,2)));
time_equi         = equi.time';

% sort good and evil
indgood  = [];
for k=1:length(time_equi)
    if all(diff(x_equi(k,:)) > 0) && all(diff(rho_tor_norm_equi(k,:)) > 0)
        indgood(end+1) = k;
    end
end
if isempty(indgood)
     error('No valid equilibrium time slice');
else
  rho_tor_norm_equi  = rho_tor_norm_equi(indgood,:);
  amin_equi          = amin_equi(indgood,:);
  x_equi             = x_equi(indgood,:);
  time_equi          = time_equi(indgood);
end

% 1- Electron density
[validity,data_out,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,cp,{'electrons','density'});
if validity == 1
    NE_EXP.x     =  x_out;
    NE_EXP.ne    =  max(1e13,data_out);
    NE_EXP.temps =  time_out - timeoffset;
    setappdata(0,'NE_EXP',NE_EXP);
else
    error('CS4M @ Ne : no valid data');
end

% updating line avaraged density reference from electron density profiles
if isappdata(0,'NE_EXP')
    nbar_new  = trapz(NE_EXP.x,NE_EXP.ne,2);
    nbar_ref  = interp1(z0dinput.cons.temps,real(z0dinput.cons.nbar),NE_EXP.temps,'linear',NaN);
    nbar_new(~isfinite(nbar_new)) = nbar_ref(~isfinite(nbar_new));
    %nbar_new(nbar_new < 1e13) = nbar_ref(nbar_new < 1e13);
    nbar_mem = real(z0dinput.cons.nbar);
    gas_puff = imag(z0dinput.cons.nbar);
    z0dinput.cons.nbar = interp1(NE_EXP.temps,nbar_new,z0dinput.cons.temps,'linear',NaN);
    z0dinput.cons.nbar(~isfinite(z0dinput.cons.nbar)) = nbar_mem(~isfinite(z0dinput.cons.nbar));
    z0dinput.cons.nbar(z0dinput.cons.temps < 0.5) = nbar_mem(z0dinput.cons.temps < 0.5);
    z0dinput.cons.nbar = max(1e13,z0dinput.cons.nbar) + sqrt(-1) .* gas_puff;
    figure;
    plot(NE_EXP.temps,nbar_ref,NE_EXP.temps,nbar_new,z0dinput.cons.temps,real(z0dinput.cons.nbar),'r');
    
    % updating also experimental value in exp0d
    ne0_new  = interp1(NE_EXP.temps,NE_EXP.ne(:,1),z0dinput.cons.temps,'linear',NaN);
    z0dinput.exp0d.ne0(isfinite(ne0_new)) = ne0_new(isfinite(ne0_new));
    
    nea_new  = interp1(NE_EXP.temps,NE_EXP.ne(:,end),z0dinput.cons.temps,'linear',NaN);
    z0dinput.exp0d.nebord(isfinite(nea_new)) = nea_new(isfinite(nea_new));
    
end

% 2- Electron temperature
[validity,data_out,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,cp,{'electrons','temperature'});
if validity == 1
    TE_EXP.x     =  x_out;
    TE_EXP.te    =  max(sqrt(eps),data_out);
    TE_EXP.temps =  time_out - timeoffset;
    setappdata(0,'TE_EXP',TE_EXP);
    
    % updating also experimental value in exp0d
    te0_new  = interp1(TE_EXP.temps,TE_EXP.te(:,1),z0dinput.cons.temps,'linear',NaN);
    z0dinput.exp0d.te0(isfinite(te0_new)) = te0_new(isfinite(te0_new));
    
    tea_new  = interp1(TE_EXP.temps,TE_EXP.te(:,end),z0dinput.cons.temps,'linear',NaN);
    z0dinput.exp0d.tebord(isfinite(tea_new)) = tea_new(isfinite(tea_new));
    
else
    error('CS4M @ Te : no valid data');
end

% 3- Ion temperature
[validity,data_out,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,cp,'t_i_average');
if validity == 1
    TI_EXP.x     =  x_out;
    TI_EXP.ti    =  max(sqrt(eps),data_out);
    TI_EXP.temps =  time_out - timeoffset;
    setappdata(0,'TI_EXP',TI_EXP);
    
else
    error('CS4M @ Ti : no valid data');
end

%4 - Wth
% filtrer Inf et NaN
indok = find(isfinite(cp.global_quantities.energy_diamagnetic) & (cp.global_quantities.energy_diamagnetic > -9e40));
z0dinput.exp0d.wth  = max(eps,interp1(cp.time(indok) - timeoffset,cp.global_quantities.energy_diamagnetic(indok), ...
    z0dinput.cons.temps,'nearest','extrap'));
% ne remplacer que si vide
indok = find(isfinite(equi.global_quantities.w_mhd) & (equi.global_quantities.w_mhd > -9e40));
z0dinput.exp0d.w = max(eps,interp1(equi.time(indok) - timeoffset,equi.global_quantities.w_mhd(indok), ...
    z0dinput.cons.temps,'linear','extrap'));  
%z0dinput.exp0d.wth  = z0dinput.exp0d.w;
indok = find(isfinite(equi.global_quantities.li_3) & (equi.global_quantities.li_3 > -9e40));
z0dinput.exp0d.li = max(eps,interp1(equi.time(indok) - timeoffset,equi.global_quantities.li_3(indok), ...
    z0dinput.cons.temps,'linear','extrap'));  


%5 -Zeff
zeff     = NaN * ones(length(cp.time),length(cp.profiles_1d{1}.zeff));
nWprof   = NaN * ones(length(cp.time),length(cp.profiles_1d{1}.zeff));
nNprof   = NaN * ones(length(cp.time),length(cp.profiles_1d{1}.zeff));
nCuprof   = NaN * ones(length(cp.time),length(cp.profiles_1d{1}.zeff));
neprof   = NaN * ones(length(cp.time),length(cp.profiles_1d{1}.zeff));
teprof   = NaN * ones(length(cp.time),length(cp.profiles_1d{1}.zeff));
rho      = NaN * ones(length(cp.time),length(cp.profiles_1d{1}.zeff));
vol_prof      = NaN * ones(length(cp.time),length(cp.profiles_1d{1}.zeff));
zeff_lav = NaN * ones(length(cp.time),1);
nW_lav   = NaN * ones(length(cp.time),1);
ne_lav   = NaN * ones(length(cp.time),1);
nN_lav   = NaN * ones(length(cp.time),1);
nCu_lav  = NaN * ones(length(cp.time),1);
delta_zeff = zeros(length(cp.time),1);
for k=1:length(cp.time);
    zeff(k,:) = cp.profiles_1d{k}.zeff;
    neprof(k,:) = cp.profiles_1d{k}.electrons.density;
    teprof(k,:) = cp.profiles_1d{k}.electrons.temperature;
    rho(k,:)  = cp.profiles_1d{k}.grid.rho_tor_norm;
    vol_prof(k,:)      = cp.profiles_1d{k}.grid.volume;
    zeff_lav(k)            = trapz(rho(k,:) ,zeff(k,:),2) ./ trapz(rho(k,:),ones(size(zeff(k,:))),2);
    ne_lav(k)              = max(1e13,trapz(rho(k,:) ,neprof(k,:),2) ./ trapz(rho(k,:),ones(size(zeff(k,:))),2));
    for l=1:length(cp.profiles_1d{k}.ion)
       switch  cp.profiles_1d{k}.ion{l}.label
           case 'W'
               nWprof(k,:) = cp.profiles_1d{k}.ion{l}.density_thermal;
               nW_lav(k)   = max(0,trapz(rho(k,:) ,nWprof(k,:),2) ./ trapz(rho(k,:),ones(size(zeff(k,:))),2));
               % computation of delta Zeff for METIS
               zave          = z0wavez(teprof(k,:));
               dzeffp        = nWprof(k,:) ./ max(1,max(nWprof(k,:) .* zave,neprof(k,:))) .* zave .^ 2;
               dzeffp(~isfinite(dzeffp)) = 0;
               dzeffp        = max(0,min(74,dzeffp));
               delta_zeff(k) = max(0,min(3,trapz(rho(k,:) ,dzeffp,2) ./ trapz(rho(k,:),ones(size(dzeffp)),2)));
               index_w = l;
               cp.profiles_1d{k}.nwp = nWprof(k,:);
           case 'Cu'
               nCuprof(k,:) = cp.profiles_1d{k}.ion{l}.density_thermal_error_upper;
               nCu_lav(k)  = max(0,trapz(rho(k,:) ,nCuprof(k,:),2) ./ trapz(rho(k,:),ones(size(zeff(k,:))),2));
           case 'N'
               nNprof(k,:) = cp.profiles_1d{k}.ion{l}.density_thermal;
               nN_lav(k)  = max(0,trapz(rho(k,:) ,nNprof(k,:),2) ./ trapz(rho(k,:),ones(size(zeff(k,:))),2));
       end
    end

end
% propagate option slection
switch autoWESTmetis_option.zeff_mode
    case {'resistive','resitive'}
        indok = find(isfinite(zeff_res_data.zeff_res_av));
        if isempty(indok)
            disp('No valide resistive Zeff available');
            z0dinput.cons.zeff  = 3 * ones(size(z0dinput.cons.temps));
        elseif length(indok) == 1
            z0dinput.cons.zeff  = zeff_res_data.zeff_res_av * ones(size(z0dinput.cons.temps));
       else
            z0dinput.cons.zeff  = interp1(cp.time(indok) - timeoffset,zeff_res_data.zeff_res_av(indok),z0dinput.cons.temps,'nearest','extrap');
        end
        disp('using resistive Zeff in METIS');
    case 'bremsstrahlung'
        indok = find(isfinite(zeff_lav));
        z0dinput.cons.zeff  = interp1(cp.time(indok) - timeoffset,zeff_lav(indok),z0dinput.cons.temps,'nearest','extrap');
        if all(abs(mean(zeff_lav(indok)) - zeff_lav(indok)) < 1e-6)
            disp('using resistive Zeff in METIS');
       else
            disp('using visible bremsstrahlung Zeff in METIS');
        end
    case 'Matthews'
       indok = find(isfinite(zeff_res_data.zeff_matthews_av));
       z0dinput.cons.zeff  = interp1(cp.time(indok) - timeoffset,zeff_res_data.zeff_matthews_av(indok),z0dinput.cons.temps,'nearest','extrap');        
       disp('using Matthews Zeff in METIS');
    otherwise
        error('Option not yet implemented');
end
%if length(indok)>= 3
%    delta_zeff = interp1(cp.time(indok) - timeoffset,delta_zeff(indok),z0dinput.cons.temps,'nearest','extrap');
%    z0dinput.cons.zeff  = max(z0dinput.cons.zeff, 1 + delta_zeff) - delta_zeff;
%end
z0dinput.exp0d.zeff = z0dinput.cons.zeff;

% cupper
indok = find(isfinite(nCu_lav) & isfinite(nN_lav) & (imag(nCu_lav) == 0) & (imag(nN_lav) == 0) & isfinite(nW_lav) & (imag(nW_lav) == 0));
if (length(indok) >=3) & (fcupper > 0)

    % from calibration using departure to Matthews scaling law
    % before : 20190403
    %rcu_mean = 5.1180e+34;
    %rcu_std  = 7.3182e+34;
    % since 20190403 :
    %rcu_mean = 3.8954e+34;
    %rcu_std  = 4.2036e+34;
    % with OH part of the shot and resitive Zeff (20190426):
    %rcu_mean = 2.0341e+34;
    %rcu_std =  8.2880e+33;
    % correction � 1 sigma 
    %fcupper = fcupper .* (rcu_mean - rcu_std) ./ rcu_mean;

    % W cooling rate and unscaled radiative power
    tzw  = tabmat.W.data(:,1)';
    lzw  = tabmat.W.data(:,2)';
    indok  = find((tzw > 1e-4) & (tzw <= 1e2));
    tzw = tzw(indok);
    lzw = lzw(indok);
    lzw  = cat(2,1e-38,lzw,lzw(end) + 5.355e3 .* 74  .*(sqrt(105)-sqrt(tzw(end))) .* 1e-28);
    tzw  = cat(2,1e-4,tzw,105);
    rzw   = reshape(10 .^ pchip(log10(tzw),log10(lzw),log10(teprof(:)/1e3)),size(teprof));
    pradwp_ref   = nWprof .* neprof .* rzw;
    pradwp       = sum(diff(vol_prof,1,2) .* (pradwp_ref(:,1:end-1) + pradwp_ref(:,2:end)) ./ 2,2) / 1e6;

    %cupper radiative power
    tzcu  = tabmat.Cu.data(:,1)';
    lzcu  = tabmat.Cu.data(:,2)';
    indok  = find((tzcu > 1e-4) & (tzcu <= 1e2));
    tzcu = tzcu(indok);
    lzcu = lzcu(indok);
    lzcu  = cat(2,1e-38,lzcu,lzcu(end) + 5.355e3 .* 74  .*(sqrt(105)-sqrt(tzcu(end))) .* 1e-28);
    tzcu  = cat(2,1e-4,tzcu,105);
    rzcu   = reshape(10 .^ pchip(log10(tzcu),log10(lzcu),log10(teprof(:)/1e3)),size(teprof));
    pradcup_ref   = fcupper .* nCuprof .* neprof .* rzcu;
    pradcup       = sum(diff(vol_prof,1,2) .* (pradcup_ref(:,1:end-1) + pradcup_ref(:,2:end)) ./ 2,2) / 1e6;
    
    % automatic correction of fcupper to prevent to high value
    %pradwnave = trapz(cp.time(indok),pradwp(indok) .* ne_lav(indok),1) ./ trapz(cp.time(indok),ne_lav(indok),1);
    %pradcunave = trapz(cp.time(indok),pradcup(indok) .* ne_lav(indok),1) ./ trapz(cp.time(indok),ne_lav(indok),1);
    %fcupper    = max(0,min(1,pradwnave ./ max(eps,pradcunave)));
    fprintf('fcupper effective value is %g\n',fcupper);
    %pradcup_ref   = fcupper .* nCuprof .* neprof .* rzcu;
    %pradcup       = sum(diff(vol_prof,1,2) .* (pradcup_ref(:,1:end-1) + pradcup_ref(:,2:end)) ./ 2,2) / 1e6;

    % correction W density
    rapwcor_raw = min(1,max(0,max(0,pradwp - pradcup) ./ max(eps,pradwp)));
    rapwcor = ones(size(cp.time));
    rapwcor(indok) = rapwcor_raw(indok);
    figure;plot(cp.time,rapwcor);drawnow
    rapwcor1 = trapz(cp.time,rapwcor.* ne_lav,1) ./ trapz(cp.time,ne_lav,1);
    fprintf('Tungsten density change by a factor %g\n',rapwcor1);
    nW_lav  = rapwcor .* nW_lav;
    nCuc = min(1,interp1(cp.time(indok) - timeoffset,fcupper .* nCu_lav(indok) ./ max(1,nN_lav(indok)),z0dinput.cons.temps,'nearest','extrap'));
    rimp = trapz(z0dinput.cons.temps,nCuc .* real(z0dinput.cons.nbar)) ./ trapz(z0dinput.cons.temps,real(z0dinput.cons.nbar));
    if isfinite(rimp) & (rimp > 0)
        z0dinput.option.rimp = rimp;
        z0dinput.option.zimp = 7;
        z0dinput.option.zmax = 29;
        disp('Cupper indicator available: switching to Nitrogen and Cupper impurities')
        fprintf('rimp set to %g\n',z0dinput.option.rimp);
        % specific treatment of Zeff in METIS
        delta_zeff = (1 - interp1(cp.time(indok) - timeoffset,rapwcor(indok) .* delta_zeff(indok),z0dinput.cons.temps,'nearest','extrap')) ;
        %
        z0dinput.cons.zeff  = z0dinput.cons.zeff + delta_zeff;
        z0dinput.exp0d.zeff = z0dinput.cons.zeff;
    end
end


%6- tune W content
switch autoWESTmetis_option.metis_mode
    case 'off'
        fprintf('Using identified tungsten density profile as external data for METIS;\nnW_shape as been identitied using the method: %s\n',autoWESTmetis_option.Wshape);
        [validity,data_out,x_out,time_out] = rhotornorm2xli(time_equi,rho_tor_norm_equi,x_equi,cp,'nwp');
        if validity == 1
            NW_SHAPE_EXP.x     =  x_out;
            NW_SHAPE_EXP.nwp   =  max(0,data_out);
            NW_SHAPE_EXP.temps =  time_out - timeoffset;
            setappdata(0,'NW_SHAPE_EXP',NW_SHAPE_EXP);
        else
            error('CS4M @ NW_SHAPE : no valid data');
        end
        
    otherwise
        disp('using fixed peaking factor for tungsten density profile (n_W = C(t) * n_e ^ gamma');
end


z0dinput.option.cw_ecrh   = 0;
z0dinput.option.cw_factor = 0;
z0dinput.option.cw_icrh   = 0;  
z0dinput.option.cw_lhcd   = 0;
z0dinput.option.cw_nbi1   = 0;
z0dinput.option.cw_nbi2   = 0;
z0dinput.option.cw_offset = 0; 
  
if ~isfield(cp,'Wshape_factor_metis') || (cp.Wshape_factor_metis  == 1)
   indok = find(isfinite(nW_lav) & (imag(nW_lav) == 0) & (nW_lav < (0.1 .* ne_lav))); 
   nwc = interp1(cp.time(indok) - timeoffset,nW_lav(indok) ./ ne_lav(indok),z0dinput.cons.temps,'nearest','extrap');
else
   indok = find(isfinite(nWprof(:,end)) & (imag(nWprof(:,end)) == 0) & (nWprof(:,end) < (0.1 .* neprof(:,end)))); 
   nwc = interp1(cp.time(indok) - timeoffset,nWprof(indok,end) ./ neprof(indok,end),z0dinput.cons.temps,'nearest','extrap');
end
indfit = find((z0dinput.cons.temps > (min(z0dinput.cons.temps) +1)) & (z0dinput.cons.temps < (max(z0dinput.cons.temps) -0.5)));
if isempty(indfit)
    indfit = 1:length(z0dinput.cons.temps);
end
%compute linear reponse
plh   = z0dinput.cons.plh ./ 1e6;
picrh = z0dinput.cons.picrh ./ 1e6;
v1    = ones(size(z0dinput.cons.temps));
if any(plh > 0.1) && any(picrh > 0.1)
    rep = cat(2,v1(indfit),plh(indfit),picrh(indfit)) \ nwc(indfit);
    z0dinput.option.cw_offset = rep(1);
    z0dinput.option.cw_lhcd   = rep(2);
    z0dinput.option.cw_icrh   = rep(3);
    
elseif any(plh > 0.1)
    rep = cat(2,v1(indfit),plh(indfit)) \ nwc(indfit);
    z0dinput.option.cw_offset = rep(1);
    z0dinput.option.cw_lhcd   = rep(2);
elseif any(picrh > 0.1)
    rep = cat(2,v1(indfit),picrh(indfit)) \ nwc(indfit);
    z0dinput.option.cw_offset = rep(1);
    z0dinput.option.cw_icrh   = rep(2);
else
    if length(indfit) >= 7
        z0dinput.option.cw_offset = min(sgolayfilt(nwc(indfit),1,5));
    else
        z0dinput.option.cw_offset = mean(nwc(indfit));
    end
end

if ~isfield(cp,'Wshape_factor_metis') || (cp.Wshape_factor_metis  == 1)
  fprintf('n_W_bar / n_e_bar = %6.3e  + %6.3e * P_{LH,MW} + %6.3e * P_{ICRH,MW}\n',z0dinput.option.cw_offset,z0dinput.option.cw_lhcd,z0dinput.option.cw_icrh);
else
  fprintf('n_W_LCFS / n_e_LCFS = %6.3e  + %6.3e * P_{LH,MW} + %6.3e * P_{ICRH,MW}\n',z0dinput.option.cw_offset,z0dinput.option.cw_lhcd,z0dinput.option.cw_icrh);
end
figure;
plot(z0dinput.cons.temps,nwc,'b',z0dinput.cons.temps(indfit),nwc(indfit),'.b' ,z0dinput.cons.temps, ...
     z0dinput.option.cw_offset + z0dinput.option.cw_lhcd .* plh + z0dinput.option.cw_icrh .* picrh,'r-.');
title(sprintf('Adjustement tungsten concentration for WEST shot #%d',abs(autoWESTmetis_option.shot)));
set(gca,'ylim',[0,max(nwc(indfit))]);
legend({'data','selected','fit'});
xlabel('time (s)');
ylabel('C_W');
drawnow

% ICRH
if autoWESTmetis_option.force_nHonZ1 == 1
    z0dinput.option.cmin  = 1 ./ ( 1./ nhonhpnd - 1);
end
% Helium
%z0dinput.option.frhe0 = frhe0;  ALREADY SET IN ZEROD_INIT_WEST
switch autoWESTmetis_option.gas_type
    case 'automatic'
        switch z0dinput.option.gaz
            case 1
                disp('initial METIS data have been set to hydrogen main ion')
            case 2
                disp('initial METIS data have been set to deuterium main ion')
            case 4
                disp('initial METIS data have been set to helium main ion')
        end
        fprintf('Helium fraction has been set to %g\n',z0dinput.option.frhe0);
    case 'D'
        disp('Plasma forced to be in deuterium');
        z0dinput.option.frhe0 =  cp.frhe0 ./ (1 + 2 .* cp.frhe0);
        z0dinput.option.gaz   = 2;
    case'He'
        disp('Plasma forced to be in helium');
        z0dinput.option.frhe0 =  0;
        z0dinput.option.gaz   = 4;
end
switch autoWESTmetis_option.fr_he_mode
    case 'automatic'
        % rien
    case 'prescribed'
        z0dinput.option.frhe0 = autoWESTmetis_option.fr_he ./  (1 + 2 .* autoWESTmetis_option.fr_he);
        disp('Helium fraction has been forced');
end
if z0dinput.option.gaz   ~= 4
    fprintf('Helium fraction has been set to %g\n',z0dinput.option.frhe0);
else
    z0dinput.option.frhe0 =  0;
    fprintf('Helium fraction has been reset to %g for full helium shot\n',z0dinput.option.frhe0);    
end

% nW shape if available
if isfield(cp,'Wshape_factor_metis')
    z0dinput.option.fne_acc = cp.Wshape_factor_metis;
end

% cocos vrification
switch autoWESTmetis_option.cocos_onoff
    case 'off'
        z0dinput.option.COCOS_check  = 'off';
        z0dinput.option.COCOS_method = 'Native';
end

% call METIS execution
if ~isfinite(workingpoint)
    zassignin('base','z0dinput',z0dinput);
    plh = z0dinput.cons.plh;
    indplh  = find(plh >0.1e5);
    if isempty(indplh)
        disp('Now launching METIS in full computation mode');
	evalin('base','metis_run;');
    else
        disp('Now launching METIS in computation mode optimized for LHCD efficiency fitting');
	evalin('base','metis_fitlh;');   
    end
    post = evalin('base','post');
else
    setappdata(0,'WORKING_POINT_TIME',workingpoint);
    zassignin('base','z0dinput',z0dinput);
    evalin('base','z0working_point;');
    post = evalin('base','post');
end
if ~isempty(filename)
    if tprofok == 1
        metis_save(sprintf('%s_TS@%d_occ_%d',filename,shot,occ*10),post);
    else
        metis_save(sprintf('%s_TS@%d_notprof',filename,shot),post);
    end
end

% default METIS paramaters for TS
function option = change_option_west(option)

% fir Wdia
option.scaling = 4;
% No runaway
option.runaway = 0;
% nW isomorphe to ne
option.fne_acc = 1;
% fit LHCD efficientcy
option.lhmode = 1;
% no modification of the density
option.neasser = 0;
% ionisation equilibrium
option.noncoronal = 0;
% ratio of impurities
option.rimp = 0.1;
% model tail like ?
% la source est trop centralee et ne rend pas compte des exp�riences.
%option.upshiftmode = 'newmodel + tail';
% No modification of density
option.natural = 0;


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
            if any(validity < -1) && isnumeric(data_local)
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
            if any(validity < -1) && isnumeric(data_local)
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


