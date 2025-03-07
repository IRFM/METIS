function  jet_neutron_transp(post,revision)  

% test inputs 
if nargin < 2
    if isappdata(0,'TRANSP_RUN')
        revision = getappdata(0,'TRANSP_RUN');
    else
        revision = 0;
    end
end
if isappdata(0,'TRANSP_USER')
    altuser = getappdata(0,'TRANSP_USER');
else
    altuser = 'tranppf/trt0';
end

% read  TRANSP data
disp('Reading TRANSP data for neutron:')
% total
shot = post.z0dinput.shot;
disp('1/ total')
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/tot',shot,altuser),[],revision);
if ~isempty(rep)
    ntr_total       = rep.object.data.value.data;
    time_tr         = rep.object.dimensions.value.x0.value.data.value.data;
else
    ntr_total       = [];
    time_tr         = [];
end
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/todt',shot,altuser),[],revision);
if isempty(rep)
    rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/totn',shot,altuser),[],revision);
end 
if ~isempty(rep)
    ntr_total_dt    = rep.object.data.value.data;
    if isempty(time_tr)
        time_tr         = rep.object.dimensions.value.x0.value.data.value.data;
    end
else
    ntr_total_dt    = [];
    error('No data for neutron probuction found in this PPF TRANSP data');
end
%
% thermal
%
disp('2/ thermal')
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/th',shot,altuser),[],revision);
if isempty(rep)
     rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/thnt',shot,altuser),[],revision);
end    
if ~isempty(rep)
    ntr_th_total    = rep.object.data.value.data;
else
    ntr_th_total    = [];
end
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/thdd',shot,altuser),[],revision);
ntr_th_dd       = rep.object.data.value.data;
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/thdt',shot,altuser),[],revision);
ntr_th_dt       = rep.object.data.value.data;
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/thtt',shot,altuser),[],revision);
ntr_th_tt       = rep.object.data.value.data;
ntr_th_total = ntr_th_dd + ntr_th_dt + ntr_th_tt;

    
%
% beam plasma
% 
disp('3/ beam plasma')
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bp',shot,altuser),[],revision);
if isempty(rep)
  rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/btnt',shot,altuser),[],revision);
end  
if ~isempty(rep)
    ntr_bp_total    = rep.object.data.value.data;
else
    ntr_bp_total = [];
end
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bpdd',shot,altuser),[],revision);
if isempty(rep)
  rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/btdd',shot,altuser),[],revision);
end  
ntr_bp_dd       = rep.object.data.value.data;
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bpdt',shot,altuser),[],revision);
if isempty(rep)
  rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/btdt',shot,altuser),[],revision);
end  
ntr_bp_dt       = rep.object.data.value.data;
try
    rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bptd',shot,altuser),[],revision);
    ntr_bp_td       = rep.object.data.value.data;
catch
    try
        rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bttd',shot,altuser),[],revision);
        ntr_bp_td       = rep.object.data.value.data;
        
    catch
        disp('incomplete tree: no beam plasma d -> t');
        ntr_bp_td = zeros(size(time_tr));
    end
end
try
    rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bptt',shot,altuser),[],revision);
    ntr_bp_tt       = rep.object.data.value.data;
catch
    try
        rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bttt',shot,altuser),[],revision);
        ntr_bp_tt       = rep.object.data.value.data;
        
    catch
        disp('incomplete tree: no beam plasma t -> t');
        ntr_bp_tt = zeros(size(time_tr));
    end
end
if isempty(ntr_bp_total)
    ntr_bp_total = ntr_bp_dd + ntr_bp_dt;
    if all(isfinite(ntr_bp_td))
        ntr_bp_total = ntr_bp_total + ntr_bp_td;
    end
    if all(isfinite(ntr_bp_tt))
        ntr_bp_total = ntr_bp_total + ntr_bp_tt;
    end
end
%
% beam bream
%
disp('3/ beam beam')
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bb',shot,altuser),[],revision);
if isempty(rep)
    rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bbnt',shot,altuser),[],revision);   
end
if ~isempty(rep)
    ntr_bb_total    = rep.object.data.value.data;
else
    ntr_bb_total    = [];
end
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bbdd',shot,altuser),[],revision);
ntr_bb_dd       = rep.object.data.value.data;
try
    rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bbdt',shot,altuser),[],revision);
    ntr_bb_dt       = rep.object.data.value.data;
catch
    disp('incomplete tree: no beam plasma d -> t');
    ntr_bb_dt = zeros(size(time_tr));   
end
try
    rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/bbtt',shot,altuser),[],revision);
    ntr_bb_tt       = rep.object.data.value.data;
catch
    disp('incomplete tree: no beam plasma t -> t');
    ntr_bb_tt = zeros(size(time_tr));
end
if isempty(ntr_bb_total)
    ntr_bb_total = ntr_bb_dd + ntr_bb_dt + ntr_bb_tt;
end

if isempty(ntr_total)
    ntr_total = ntr_th_total + ntr_bp_total + ntr_bb_total;
end

% reading profiles
disp('4/ profiles');
data_te   = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/te',shot,altuser),[],revision);
data_ni   = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/ni',shot,altuser),[],revision);
data_nd   = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/nd',shot,altuser),[],revision);
data_nt   = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/nt',shot,altuser),[],revision);
data_ti   = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/ti',shot,altuser),[],revision);
data_dvol = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/dvol',shot,altuser),[],revision);
%
[time,x,te]   = decode_data(data_te);
[time,x,ni]   = decode_data(data_ni);
[time,x,nd]   = decode_data(data_nd);
[time,x,nt]   = decode_data(data_nt);
[time,x,ti]   = decode_data(data_ti);
[time,x,dvol] = decode_data(data_dvol);

%reading geometry
disp('5/ geometry')
% rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/rout',shot,altuser),[],revision);
% rout      = rep.object.data.value.data;
% rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/rin',shot,altuser),[],revision);
% rin       = rep.object.data.value.data;
rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/rmag',shot,altuser),[],revision);
if ~isempty(rep)
    rmag      = rep.object.data.value.data;
else
    rmag = [];
end
%
disp('end reading');


% compute volume avaraged et centre values
vt  = ones(length(time),1);
ve  = ones(1,length(x));
if ~isempty(rmag)
    mask = double((vt*x(:)') >= (rmag * ve));
else
    mask = vt * ve; 
end
tem_tr = sum(te .* dvol .* mask,2) ./  sum(dvol .* mask,2);
te0_tr = max(te,[],2);
tim_tr = sum(ti .* dvol .* mask,2) ./  sum(dvol .* mask,2);
ti0_tr = max(ti,[],2);
nim_tr = sum(ni .* dvol .* mask,2) ./  sum(dvol .* mask,2);
ni0_tr = max(ni,[],2);
ndm_tr = sum(nd .* dvol .* mask,2) ./  sum(dvol .* mask,2);
nd0_tr = max(nd,[],2);
ntm_tr = sum(nt .* dvol .* mask,2) ./  sum(dvol .* mask,2);
nt0_tr = max(nt,[],2);





disp('recomputing neutron contribution from METIS:')
if length(post.profil0d.temps) == length(post.zerod.temps)
    zs   = post.zerod;
    geo  = post.z0dinput.geo;
    cons = post.z0dinput.cons;
else
    noms = fieldnames(post.zerod);
    temps = post.zerod.temps;
    for k=1:length(noms)
        nomc = noms{k};
        var  = post.zerod.(nomc);
        if length(var) == length(temps)
            zs.(nomc) = interp1(temps,var,post.profil0d.temps,'nearest','extrap');
        else
            zs.(nomc) = var;
        end
    end
    zs.temps = post.profil0d.temps;
    noms = fieldnames(post.z0dinput.cons);
    temps = post.z0dinput.cons.temps;
    for k=1:length(noms)
        nomc = noms{k};
        var  = post.z0dinput.cons.(nomc);
        if length(var) == length(temps)
            cons.(nomc) = interp1(temps,var,post.profil0d.temps,'nearest','extrap');
        else
            cons.(nomc) = var;
        end
    end
    cons.temps = post.profil0d.temps;
    noms = fieldnames(post.z0dinput.geo);
    for k=1:length(noms)
        nomc = noms{k};
        var  = post.z0dinput.geo.(nomc);
        if length(var) == length(temps)
            geo.(nomc) = interp1(temps,var,post.profil0d.temps,'nearest','extrap');
        else
            geo.(nomc) = var;
        end
    end
end
[pfus,salpha,~,~,~,~,~,~,~,pfus_loss,~,~,~,splustd,splusdt,splusff,splusicrh] = ...
    zfus0tae(zs.nDm,zs.nTm,zs.tem,zs.nem,zs.zeff,zs.tite,geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,zs.vp,zs.sp, ...
    zs.pnbi_th,zs.taus_nbi,zs.ecrit_nbi,post.z0dinput.option.einj,post.z0dinput.option.einj2 ,cons.ftnbi, ...
    zs.pion_icrh,zs.taus_icrh,zs.ecrit_nbi,zs.einj_icrh,post.profil0d.temps,cons.pnbi,zs.d0,zs.qa,zs.qmin, ...
    zs.te0,zs.nebord,zs.tebord,zs.pped,zs.nim,zs.wth, ...
    post.z0dinput.option.tae,post.z0dinput.option.nb_nbi,post.z0dinput.option.fspot,post.z0dinput.option.e_shielding,post.profil0d,...
    post.z0dinput.option.fpolarized,post.z0dinput.option.forced_H_NBI);
% added losses for neutron total
salpha = salpha .* (pfus + pfus_loss) ./ pfus;

[neutron_total_tt,neutron_th_tt,neutron_nbi_th_tt,neutron_nbi_nbi_tt,pttfus,proton_tt,picrh_nbi_tt,einj_tt] = ...
    z0neutron_tt(post.z0dinput.option,cons,zs,post.profil0d);

[neutron_total_dd,neutron_th_dd,neutron_nbi_th_dd,neutron_nbi_nbi_dd,pddfus,proton_dd, picrh_nbi,einj] = ...
    z0neutron_dd(post.z0dinput.option,cons,zs,post.profil0d);

disp('end recomputing');

%
tm = post.profil0d.temps;

% display results
h = findobj(0,'type','figure','tag','neutron_jet_transp');
if isempty(h)
    h=figure('tag','neutron_jet_transp');
else
    figure(h);
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
    'defaultlinelinewidth',1,'color',[1 1 1])

%total (first line)
subplot(4,5,1)
plot(tm,salpha + neutron_total_tt + neutron_total_dd,'b',time_tr,ntr_total,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
title(sprintf('JET @ %d',shot));
ylabel('total')
subplot(4,5,2)
plot(tm,salpha,'b',time_tr,ntr_total_dt,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
title('Neutron flux');
ylabel('total DT');
subplot(4,5,3)
plot(tm,salpha,'b',time_tr,ntr_total_dt,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
legend('METIS','TRANSP');
% thermal (second line-
subplot(4,5,6)
nth_dt = salpha - splustd - splusdt - splusff -splusicrh;
plot(tm,nth_dt + neutron_th_tt + neutron_th_dd,'b',time_tr,ntr_th_total,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('total thermal');
subplot(4,5,7)
plot(tm,nth_dt,'b',time_tr,ntr_th_dt,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('thermal DT');
subplot(4,5,9)
plot(tm,neutron_th_dd,'b',time_tr,ntr_th_dd,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('thermal DD');
subplot(4,5,10)
plot(tm,neutron_th_tt,'b',time_tr,ntr_th_tt,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('thermal TT');
% beam plasma (thrid line)
subplot(4,5,11)
plot(tm,splustd + splusdt + neutron_nbi_th_tt + neutron_nbi_th_dd,'b',time_tr,ntr_bp_total,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('total beam plasma');
subplot(4,5,12)
plot(tm,splusdt ,'b',time_tr,ntr_bp_dt,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
ylabel('beam plasma DT');
subplot(4,5,13)
plot(tm,splustd ,'b',time_tr,ntr_bp_td,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('beam plasma TD');
subplot(4,5,14)
plot(tm,neutron_nbi_th_dd,'b',time_tr,ntr_bp_dd,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('beam plasma DD');
subplot(4,5,15)
plot(tm,neutron_nbi_th_tt,'b',time_tr,ntr_bp_tt,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('beam plasma TT');
% beam plasma (fourth line)
subplot(4,5,16)
plot(tm,splusff + splusicrh + neutron_nbi_nbi_dd + neutron_nbi_nbi_tt,'b',time_tr,ntr_bb_total,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('total beam beam');
xlabel('time (s)');
subplot(4,5,17)
plot(tm,splusff + splusicrh,'b',time_tr,ntr_bb_dt,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('beam beam DT');
xlabel('time (s)');
subplot(4,5,19)
plot(tm,neutron_nbi_nbi_dd,'b',time_tr,ntr_bb_dd,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('beam beam DD');
xlabel('time (s)');
subplot(4,5,20)
plot(tm,neutron_nbi_nbi_tt,'b',time_tr,ntr_bb_tt,'r');
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('beam beam TT');
xlabel('time (s)');

% other data
subplot(4,5,8)
plot(tm,post.profil0d.tip(:,1)/1e3,tm,zs.tem .* zs.tite/1e3,time,ti0_tr/1e3,'-.',time,tim_tr/1e3,'-.')
legend({'T_{i,0} METIS','<T_i> METIS','T_{i,0} TRANSP','<T_i> TRANSP'})
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('keV');
subplot(4,5,18)
plot(tm,post.profil0d.nip(:,1)/1e19,tm,zs.nim/1e19,time,ni0_tr/1e19,'-.',time,nim_tr/1e19,'-.')
legend({'n_{i,0} METIS','<n_i> METIS','n_{i,0} TRANSP','<n_i> TRANSP'})
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('10^{19} m^{-3}');
xlabel('time (s)');

vem  = ones(size(post.profil0d.xli));
nDp  = post.profil0d.n1p .* ((zs.nDm ./ zs.n1m) * vem);
nTp  = post.profil0d.n1p .* ((zs.nTm ./ zs.n1m) * vem);


subplot(4,5,4)
plot(tm,nDp(:,1)/1e19,tm,zs.nDm/1e19,time,nd0_tr/1e19,'-.',time,ndm_tr/1e19,'-.')
legend({'n_{D,0} METIS','<n_D> METIS','n_{D,0} TRANSP','<n_D> TRANSP'})
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('10^{19} m^{-3}');
title(sprintf('user = %s',altuser));
subplot(4,5,5)
plot(tm,nTp(:,1)/1e19,tm,zs.nTm/1e19,time,nt0_tr/1e19,'-.',time,ntm_tr/1e19,'-.')
legend({'n_{T,0} METIS','<n_T> METIS','n_{T,0} TRANSP','<n_T> TRANSP'})
set(gca,'xlim',[min(time_tr),max(time_tr)]);
set(gca,'ylim',[0,max(get(gca,'ylim'))]);
ylabel('10^{19} m^{-3}');
title(sprintf('revision # %d',revision));


function [time,x,data,description,data_units,time_units,time_desc,x_units,x_desc,certication] = decode_data(rep)


if isempty(rep)
    disp('no data')
    return
end
try
    description = rep.object.description.value;
catch
    description = '';
end
try
    data_units = rep.object.units.value;
catch
    data_units = '';
end
try
    data = rep.object.data.value.data;
catch
    data = [];
end
try
    time = rep.object.dimensions.value.x0.value.data.value.data(:);
catch
    time = [];
end
try
    time_units = rep.object.dimensions.value.x0.value.units.value;
catch
    time_units = 's';
end
try
    time_desc = rep.object.dimensions.value.x0.value.description.value;
catch
    time_desc = 'time';
end
try
    x_units = rep.object.dimensions.value.x1.value.units.value;
catch
    x_units = 'au';
end
try
    x_desc = rep.object.dimensions.value.x1.value.description.value;
catch
    x_desc = 'x';
end
try
    x = rep.object.dimensions.value.x1.value.data.value.data(:)';
catch
    x = []; 
end
try
    certication = rep.object.mask.value.status.value.data;
catch
    certication = NaN *data;
end





