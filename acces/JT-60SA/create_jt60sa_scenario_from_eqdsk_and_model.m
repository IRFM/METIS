% this function creating a new starting point for METIS simulation for
% JT-60SA using at set of EQDSK files names by time in ms and a METIS file
% simulation use as a model for parametring the new simulation:
%
% z0dinput = create_jt60sa_scenario_from_eqdsk_and_model(gfiles_dirname,metis_file_name)
%
% gfiles_dirname = path to directory containing the EQDSK files
%
% metis_file_name = full path and name of METIS file simulation
%
% z0dinput = new METIS input data structure
%
function z0dinput = create_jt60sa_scenario_from_eqdsk_and_model(dirname,metis_file)

if nargin < 2 
    error('z0dinput = create_jt60sa_scenario_from_eqdsk_and_model(gfiles_dirname,metis_file_name);')
end

if isempty(which('get_eqdsk_data_for_scenario'))
   addpath(fullfile(fileparts(which('create_jt60sa_scenario_from_eqdsk_and_model')),'gfile'));
end

% read metis file
 metis_load(metis_file)
% get z0dinput model
z0dinput_model = evalin('base','z0dinput');

% read data from eqdsk files
[time,ip,R_LCFS,Z_LCFS,RB0,Rsepa,Zsepa,R0,Z0,aminor,K,d,psi,psin,q,ptot,fdia,pprim,ffprim] = get_eqdsk_data_for_scenario(dirname);

% make new z0dinput
z0dinput = z0dinput_model; % to get all the fields before modification
%
z0dinput.machine = sprintf('JT-60SA-%d',fix(max(ip/1e5)));
z0dinput.option.machine = z0dinput.machine;
z0dinput.shot    = fix(max(ip/1e5));
z0dinput.option.shot = z0dinput.shot;

% synchronized time coordinate
[tip,ft]    = ftip(ip);
[tipm,ftm]  = ftip(z0dinput.cons.ip);
% make new time vector for model
tip_ru  = tip(1:(find(ft,1)-1));
tip_ru  = tip_ru - min(tip_ru);
tip_ru  = tip_ru ./ max(tip_ru);
tip_ft  = tip(find(ft));
tip_ft  = tip_ft - min(tip_ft);
tip_ft  = tip_ft ./ max(tip_ft);
tip_rd  = tip((max(find(ft))+1):end);
tip_rd  = tip_rd - min(tip_rd);
tip_rd  = tip_rd ./ max(tip_rd);
tipm_ru  = tipm(1:(find(ftm,1)-1));
tipm_ru  = tipm_ru - min(tipm_ru);
tipm_ru  = tipm_ru ./ max(tipm_ru);
tipm_ft  = tipm(find(ftm));
tipm_ft  = tipm_ft - min(tipm_ft);
tipm_ft  = tipm_ft ./ max(tipm_ft);
tipm_rd  = tipm((max(find(ftm))+1):end);
tipm_rd  = tipm_rd - min(tipm_rd);
tipm_rd  = tipm_rd ./ max(tipm_rd);
ind_ru   = (1:(find(ftm,1)-1));
ind_ft   = find(ftm);
ind_rd   = ((max(find(ftm))+1):length(ftm));
% new index ru
ind_new_ru = interp1(tipm_ru,ind_ru,tip_ru,'nearest','extrap');
ind_new_ru(1) = ind_ru(1);
ind_new_ru(end) = ind_ru(end);
ind_new_ft = interp1(tipm_ft,ind_ft,tip_ft,'nearest','extrap');
ind_new_ft(1) = ind_ft(1);
ind_new_ft(end) = ind_ft(end);
ind_new_rd = interp1(tipm_rd,ind_rd,tip_rd,'nearest','extrap');
ind_new_rd(1) = ind_rd(1);
ind_new_rd(end) = ind_rd(end);
ind_new  =cat(1,ind_new_ru(:),ind_new_ft(:),ind_new_rd(:));
%
% update structure data
%
noms = fieldnames(z0dinput.cons);
for k=1:length(noms)
    z0dinput.cons.(noms{k}) = z0dinput_model.cons.(noms{k})(ind_new);
end
z0dinput.cons.temps = time;
z0dinput.cons.ip    = ip;
z0dinput.cons.flux  = -(psi(:,end) - psi(1,end))/2/pi .* sign(psi(end,end) - psi(1,end));
z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
noms = fieldnames(z0dinput.exp0d);
for k=1:length(noms)
    z0dinput.exp0d.(noms{k}) = z0dinput_model.exp0d.(noms{k})(ind_new);
end
z0dinput.exp0d.temps = time;
z0dinput.exp0d.ip    = ip;
z0dinput.exp0d.Rsepa = Rsepa;
z0dinput.exp0d.Zsepa = Zsepa - Z0(:) * ones(1,size(Zsepa,2));
z0dinput.exp0d.q0    = q(:,1);
for k=1:length(time)
    z0dinput.exp0d.q95(k) = interp1(psin(k,:),q(k,:),0.95,'linear','extrap');
end
z0dinput.exp0d.qa    = q(:,end); 
z0dinput.exp0d.qeff  = q(:,end-1);
z0dinput.exp0d.qmin  = min(q,[],2);
z0dinput.exp0d.edgeflux = 2 .* pi .* z0dinput.cons.flux;
%
% fill geo substructure
%
z0dinput.geo.a    = aminor(:);
z0dinput.geo.R    = R0(:);
z0dinput.geo.K    = K(:);
z0dinput.geo.d    = d(:);
z0dinput.geo.b0   = RB0(:) ./ z0dinput.geo.R;
z0dinput.geo.z0   = Z0;
z0dinput.geo.sp   = [];
z0dinput.geo.vp   = [];
z0dinput.geo.sext = [];




function [tip,ft] = ftip(ip)

ft   = double(abs((ip - max(ip)) ./ (max(ip) - min(ip))) < 5e-2);
tip  = cumsum(1+ abs(diff(cat(1,ip -min(ip),0))));
tip  = tip - min(tip);
tip  = tip ./ max(tip);

