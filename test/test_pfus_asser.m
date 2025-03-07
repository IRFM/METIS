% script pour le calcul de metis avec zerodevolution 
% test for pfus asservissement
% conversion for gas
pam3 = 2.6558e+20;  % particles, atoms or molecules for 1 Pa.m^3
% load reference simulation
datafile = fullfile(fileparts(which('metis')),'certification','metis','reference_NTM_ITER.mat');
data           = load_metis_imas(datafile);
z0dinput       = data.z0dinput;
% securite
if isfield(z0dinput.exp0d,'Rsepa') && isempty(z0dinput.exp0d.Rsepa)
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
end
if isfield(z0dinput.exp0d,'XDURx') && isempty(z0dinput.exp0d.XDURx)
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURx');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURt');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURv');
end
% initialisation
cons1t  = zerod_get1t(z0dinput.cons,1);
geo1t   = zerod_get1t(z0dinput.geo,1);    % the density nbar is provider for the first time
exp0d1t = zerod_get1t(z0dinput.exp0d,1);
[zs,profil,z0dstruct] = zerodevolution([],z0dinput.option,z0dinput.cons.temps(1),cons1t,geo1t,z0dinput.exp0d,exp0d1t,exp0d1t);
%
% here is controller initialisation
%
tmin = 120;
tmax = z0dinput.cons.temps(end);
dt = 1;
tlist = tmin:dt:tmax;
pfus_ref = 80e6;
tau_e    = 3; 
pnbi_max = 3 .* 16.5e6;
inte_diff = 0;
gain = 3;
% start of loop onf time slices
for k_index = 1:length(tlist)
    time  = tlist(k_index);
	cons1t  = zerod_get_one_time(z0dinput.cons,z0dinput.cons.temps,time,'linear','extrap');
	geo1t   = zerod_get_one_time(z0dinput.geo,z0dinput.cons.temps,time,'linear','extrap');
	exp0d1t = zerod_get_one_time(z0dinput.exp0d,z0dinput.cons.temps,time,'linear','extrap');
    %
    % here is the density controller
    %
    prop      = pfus_ref - zs.pfus;
    inte_diff = inte_diff + prop * dt ; % 2 from molecule to electron
    feed      = real(cons1t.pnbi) + imag(cons1t.pnbi);
    pnbi      = pnbi_max .* max(0,tanh(gain .* (prop + inte_diff/ tau_e + feed) ./ pnbi_max));
    disp([prop,inte_diff,feed,pnbi])
    cons1t.pnbi  = (1 + sqrt(-1)) ./ 2  .* pnbi;
    %
    % call of METIS in evolution mode
	[zs,profil,z0dstruct] = zerodevolution(z0dstruct,z0dinput.option,time,cons1t,geo1t,z0dinput.exp0d,exp0d1t,exp0d1t);
    disp([pfus_ref,zs.pfus]);
end

post.z0dinput = z0dstruct.z0dinput;
post.zerod    = z0dstruct.zs;
post.profil0d = z0dstruct.profil;
