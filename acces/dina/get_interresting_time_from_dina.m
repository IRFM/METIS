function [selected_time,indice_selected] = get_interresting_time_from_dina(sm,nbt,ip_min,R_min,plotonoff)


% for testing
if nargin == 0
    idx = imas_open_env_backend(105092,1,getenv('USER'),'iter','3',13);
    sm = ids_get(idx,'summary');
    imas_close(idx);
end
if (nargin < 2) || isempty(nbt)
    nbt = 51;
end
if (nargin < 3) || isempty(ip_min)
    ip_min = 50e3;
end
if (nargin < 4) || isempty(R_min)
R_min = 1;
end
if (nargin < 5) || isempty(plotonoff)
    plotonoff = true;
end

% time
t = sm.time;

% energy signal
wth = sm.global_quantities.energy_thermal.value;
wbp = sm.global_quantities.energy_b_field_pol.value;

% for Bv
ip = sm.global_quantities.ip.value;
if isfield(sm.global_quantities,'li')
    li = sm.global_quantities.li.value;
else
    li = sm.global_quantities.li_3.value;
end
betap = sm.global_quantities.beta_pol.value;

R = sm.boundary.magnetic_axis_r.value;
a = sm.boundary.minor_radius.value;
K = sm.boundary.elongation.value;

indice_valid = (R>R_min) & (abs(ip) > ip_min);
R = max(R,R_min);



% constante 
mu0 = 4 .* pi .* 1e-7;

% proxy for vertical magneitic field
denom = 4 .* pi .* R  ./ (8 .* R ./ a ./ sqrt(K) + betap + li /2 - 3/2);
bv = mu0 .* ip ./ denom;

% time derivative
dwthdt = gradient(wth,t);
dwbpdt = gradient(wbp,t);
dbvdt  = gradient(bv,t);
dipdt  = gradient(ip,t);

% control variable
fwi  = cumtrapz(t,abs(dwthdt) + abs(dwbpdt));
fwi  = (fwi - min(fwi)) ./ (max(fwi) - min(fwi));
fbvi = cumtrapz(t,abs(dbvdt));
fbvi = (fbvi - min(fbvi)) ./ (max(fbvi)-min(fbvi));
fip  = cumtrapz(t,abs(dipdt));
fip  = (fip - min(fip)) ./ (max(fip) - min(fip));
% added to be strictely monotonic and have some points in flattop
ft   = (t - min(t)) ./ (max(t)-min(t));
%
f = (fwi + fbvi + fip + ft) ./ 4;

% juste to take into account validity
fv = f(indice_valid);
tv = t(indice_valid);
fvn = (fv - min(fv)) ./ (max(fv - min(fv)));

% time selection
list_indice     = find(indice_valid);
indice_selected = list_indice(fix(interp1(fvn,1:length(fvn),linspace(0,1,nbt),'nearest')));
selected_time   = t(indice_selected);

if ~plotonoff
    return
end

figure;
plot(t,f,t(indice_selected),f(indice_selected),'or')
title('time slices selection')
ylabel('control variable')
xlabel('time (s)')

figure;
subplot(3,3,1)
plot(t,wth,t(indice_selected),wth(indice_selected),'ro');
ylabel('wth (J)');
subplot(3,3,2)
plot(t,wbp,t(indice_selected),wbp(indice_selected),'ro');
ylabel('wbp (J)')
subplot(3,3,3)
plot(t,ip,t(indice_selected),ip(indice_selected),'ro');
ylabel('ip (A)')
subplot(3,3,4)
plot(t,li,t(indice_selected),li(indice_selected),'ro');
ylabel('li')
subplot(3,3,5)
plot(t,betap,t(indice_selected),betap(indice_selected),'ro')
ylabel('betap')
subplot(3,3,6);
plot(t,R,t(indice_selected),R(indice_selected),'ro');
ylabel('R (m)')
subplot(3,3,7);
plot(t,a,t(indice_selected),a(indice_selected),'ro');
ylabel('a (m)')
subplot(3,3,8);
plot(t,K,t(indice_selected),K(indice_selected),'ro');
ylabel('K')
