function test_transp_dt(shot,revision)

% read data from public transp data
if nargin < 2
    revision = 0;
end
data_te   = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/tra0/te',shot),[],revision);
data_ne   = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/tra0/ne',shot),[],revision);
data_nd   = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/tra0/nd',shot),[],revision);
data_nt   = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/tra0/nt',shot),[],revision);
data_ti   = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/tra0/ti',shot),[],revision);
data_dvol = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/tra0/dvol',shot),[],revision);
%

[time,x,te]   = decode_data(data_te);
[time,x,ne]   = decode_data(data_ne);
[time,x,nd]   = decode_data(data_nd);
[time,x,nt]   = decode_data(data_nt);
[time,x,ti]   = decode_data(data_ti);
[time,x,dvol] = decode_data(data_dvol);

rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/rout',shot),[],revision);
rout      = rep.object.data.value.data;
rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/rin',shot),[],revision);
rin       = rep.object.data.value.data;
rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/rmag',shot),[],revision);
rmag      = rep.object.data.value.data;


%
rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/thdt',shot),[],revision);
ntr_th_dt       = rep.object.data.value.data;

% compute sn_dt
[dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion(max(ti,30));
sn_dt = nd .* nt .* dt.sv;
vt  = ones(length(time),1);
ve  = ones(1,length(x));
%mask = double(((vt*x(:)') >= (rin * ve)) &( (vt*x(:)') <= (rout * ve)));
mask = double((vt*x(:)') >= (rmag * ve));
%ntest = sum(sn_dt .* dvol,2) ./ (rout - rin);
ntest = sum(sn_dt .* dvol .* mask,2);


figure;plot(time,ntr_th_dt,'r',time,ntest,'b');

 
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


