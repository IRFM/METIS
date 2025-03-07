% just load any METIS simulation
filename = fullfile(fileparts(which('metis')),'certification','metis','ITER_rampup_ECCD.mat');
metis_load(filename);
% post = load_metis_imas(filename); % to be used in a function
%
% not the same algorithm than zerodevolution (metis_run use waveform relaxation, 
% also time derivative are not the same than in zerodevolution as erodevolution use finite difference)
if true
    % metis_run; % already done when you load a METIS simulation
else
    zerodrunevolution;
    post.z0dinput = z0dstruct.z0dinput;
    post.zerod = z0dstruct.zs;
    post.profil0d = z0dstruct.profil;
end
%time = post.z0dinput.cons.temps(1); % starting time
time = post.z0dinput.cons.temps(3); % starting time
figure;
hold on
plot(post.zerod.temps,post.zerod.betaptot,'r');
[zs,profil,z0dstruct] = make_z0dstruct_from_post(post,time); 
time = zs.temps; % update time if change in make_z0dstruct_from_post
plot(zs.temps,zs.betaptot,'.');
drawnow
option = z0dstruct.z0dinput.option;
option.short_dt = 'off'; % optionnal (off, on or full: full is more stable)
option.mode_expo_inte = 1; % can be set to 0 with small dt, it is faster.
option.dwdt_method = 'implicit'; % method used to compute energy time derivative ('none' can be used if there is some convergence problem).
% loop on time
dt = 1; %s
while (time<  post.z0dinput.cons.temps(end))
    time = time + dt;
    % extract data for evolution (this is optionnal, you can use directly
    % post.z0dinput.cons, ...). It is to show how to use it with an
    % external code providing these data.
    cons1t  = zerod_get_one_time(post.z0dinput.cons,post.z0dinput.cons.temps,time,'linear','extrap');
	geo1t   = zerod_get_one_time(post.z0dinput.geo,post.z0dinput.cons.temps,time,'linear','extrap');
	exp0d1t = zerod_get_one_time(post.z0dinput.exp0d,post.z0dinput.cons.temps,time,'linear','extrap'); % for LCFS, other fields are not used.
    % call METIS in evoltuion time
    [zs,profil,z0dstruct] = zerodevolution(z0dstruct,option,time,cons1t,geo1t,post.z0dinput.exp0d,exp0d1t);
    % accumulated data are in z0dstruct zs and profil.
    plot(zs.temps,zs.betaptot,'.');
    
end