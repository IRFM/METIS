function transp4metis(shot,revision,alternative_user)

% read data from public transp data
if (nargin < 2) || isempty(revision) || ~isfinite(revision)
    revision = 0;
end

% exemple of alternative user
%[s,t] = sal_get('data/pulse/99971/ppf/signal/zstancar/trau/te',[],1450);
% alternative_user = 'zstancar/trau' for testing
if (nargin <3) || isempty(alternative_user)
    alternative_user = 'tranppf/tra0';
    other_source = false;
else
    other_source = true;
end

data_te   = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/te',shot,alternative_user),[],revision);
if isempty(data_te)
    error('No TRANSP data available in SAL tree');
end
data_ne   = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/ne',shot,alternative_user),[],revision);
data_ti   = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/ti',shot,alternative_user),[],revision);
%data_dvol = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/dvol',shot,alternative_user),[],revision);
%
[time,x,te]   = decode_data(data_te);
[time,x,ne]   = decode_data(data_ne);
[time,x,ti]   = decode_data(data_ti);
%[time,x,dvol] = decode_data(data_dvol);

%
% if other_source
%    liste =  {'ppf/@shot/MAGN/IPLA'};
%    liste{1} = 'ppf/@shot/EFIT/RBND';
%    liste{2} = 'ppf/@shot/EFIT/RMAG';
%    data_EFIT  = cgcgetjet(shot,liste,'','');
%    rmag = interp1(data_EFIT.ppf.EFIT.RMAG.t,data_EFIT.ppf.EFIT.RMAG.data,time,'linear','extrap');
%    rin  = interp1(data_EFIT.ppf.EFIT.RBND.t, min(data_EFIT.ppf.EFIT.RBND.data,[],2),time,'linear','extrap');
%    rout = interp1(data_EFIT.ppf.EFIT.RBND.t, min(data_EFIT.ppf.EFIT.RBND.data,[],2),time,'linear','extrap');

if ~other_source
    rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/rout',shot),[],revision);
    rout      = rep.object.data.value.data;
    rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/rin',shot),[],revision);
    rin       = rep.object.data.value.data;
    rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/rmag',shot),[],revision);
    rmag      = rep.object.data.value.data;
else
    rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/rout',shot,alternative_user),[],revision);
    if~isempty(rep)
        rout      = rep.object.data.value.data;
    else
        rout = [];
    end
    rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/rin',shot,alternative_user),[],revision);
    if ~isempty(rep)
        rin       = rep.object.data.value.data;
    else
        rin = [];
    end
    rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/rmag',shot,alternative_user),[],revision);
    if ~isempty(rep)
        rmag      = rep.object.data.value.data;
    else
        rmag = [];
    end
    if ~isempty(rout) && ~isempty(rin) && ~isempty(rmag)
        other_source = false;
    else
        disp('some data are missing: trying alternative data');
    end
end

%search for metis data
% try
%     Raxe = evalin('base','post.profil0d.Raxe');
%     epsi = evalin('base','post.profil0d.epsi');
%     %rmx  = evalin('base','post.profil0d.epsi');
%     xli  = evalin('base','post.profil0d.xli');
%     temps  = evalin('base','post.profil0d.temps');
%     %
%     Rloc = cat(2,Raxe(:,end:-1:2) .* (1 - epsi(:,end:-1:2)),Raxe .* (1 + epsi));
%     xxx   = cat(2,-xli(end:-1:2),xli);
%     %rmx  = cat(2,-rmx(:,end:-1:2),rmx);
%     if any(size(x) == 1)
%         Rtr  = ones(length(time),1) * x(:)';
%     else
%         Rtr  = x;
%     end
% 
% catch
%     disp('no METIS data, use proxy for coordinate change');
%     Raxe = [];
% end


xx = linspace(0,1,21);
te_out = NaN * ones(length(time),length(xx));
ti_out = NaN * ones(length(time),length(xx));
ne_out = NaN * ones(length(time),length(xx));
nbar_out  = NaN * time;
% loop on time slices
for k=1:length(time)
%     if ~isempty(Raxe)
%         indt = find(temps >= time(k),1);
%         if isempty(indt)
%             indt = length(temps);
%         end
%         xloc = interp1(Rloc(indt,:),xxx,Rtr(k,:),'pchip','extrap');
%     else
%         pl     = ((te(k,:) + ti(k,:)) .* ne(k,:));
%         %rloc   = x(dvol(k,:) > 0);
%         rgeo   = sum(x.* double(pl == max(pl))) ./ sum(double(pl == max(pl)));
%         %aminor = max(rloc) - rgeo;
%         aminor = max(x) - rgeo;
%         xloc   = (x - rgeo) ./ aminor;
%     end
    if other_source
        xloc = x;
    else
        xloc        = (x - rmag(k)) ./ (rout(k) - rmag(k)) .* (x >= rmag(k)) + (x - rmag(k)) ./ (rmag(k) - rin(k)).* (x < rmag(k));
    end
    te_out(k,:) = interp1(xloc,max(30,te(k,:)),xx,'pchip','extrap');
    ti_out(k,:) = interp1(xloc,max(30,ti(k,:)),xx,'pchip','extrap');
    ne_out(k,:) = interp1(xloc,max(1e17,ne(k,:)),xx,'pchip','extrap');
    nbar_out(k) = trapz(x,max(1e17,ne(k,:))) ./ trapz(x,ones(size(x)));
end
if length(time) == length(xx)
    te_out(end+1,:) = te_out(end,:);
    ti_out(end+1,:) = ti_out(end,:);
    ne_out(end+1,:) = ne_out(end,:);
    time(end+1)     = time(end) + 0.01;
end

% make extrnal data
NE_EXP.temps =  time;
NE_EXP.x     =  xx;
NE_EXP.ne    =  ne_out;
setappdata(0,'NE_EXP',NE_EXP);
%
TE_EXP.temps =  time;
TE_EXP.x     =  xx;
TE_EXP.te    =  te_out;
setappdata(0,'TE_EXP',TE_EXP);
TI_EXP.temps =  time;
TI_EXP.x     =  xx;
TI_EXP.ti    =  ti_out;
setappdata(0,'TI_EXP',TI_EXP);

z0dinput = evalin('base','z0dinput');
nbar_new = interp1(time,nbar_out,z0dinput.cons.temps,'linear',NaN);
nbar_new(~isfinite(nbar_new)) = z0dinput.cons.nbar(~isfinite(nbar_new));
zassignin('base','z0dinput.cons.nbar',nbar_new);
disp('Interpretative data for N_e, T_e and T_i hab been read from TRANSP results');
 
 
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


