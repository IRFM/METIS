function input = metis4muscle3_update_metis_input_on_equilibrium(input,equi,LCFS_nbpts,only_ip)


% backward compatibility
%  extrapolation for external data
if isfield(input.option,'extrapolation_ext_data')
    setappdata(0,'EXTRAPOLATION_NEAREST_EX',input.option.extrapolation_ext_data)
end
if isfield(input.option,'nonan_ext_data')
    setappdata(0,'NONAN_EX',input.option.nonan_ext_data);
end

% collect data
equi_time    = equi.time(:);
equi_ip      = NaN * ones(size(equi_time));
equi_vloop   = NaN * ones(size(equi_time));
equi_flux    = NaN * ones(size(equi_time));
equi_R       = NaN * ones(size(equi_time));
equi_z0      = NaN * ones(size(equi_time));
equi_a       = NaN * ones(size(equi_time));
equi_K       = NaN * ones(size(equi_time));
equi_d       = NaN * ones(size(equi_time));
equi_b0      = NaN * ones(size(equi_time));
equi_vp      = NaN * ones(size(equi_time));
equi_sp      = NaN * ones(size(equi_time));
equi_sext    = NaN * ones(size(equi_time));
fbe_lcfs     = {};
%
for k=1:length(equi_time)
    equi_ip(k)    = abs(equi.time_slice{k}.global_quantities.ip);
    equi_flux(k)  = - equi.time_slice{k}.global_quantities.psi_boundary / (2*pi); % data are in COCOS 11 and METIS natively COCOS 7
    % geometry
    [~,~,~,~,geo_in] = zgeo0(input.geo1t, ...
                   equi.time_slice{k}.boundary.outline.r(:).', ...
                   equi.time_slice{k}.boundary.outline.z(:).',1);
    equi_R(k)     = (max(equi.time_slice{k}.boundary.outline.r) + min(equi.time_slice{k}.boundary.outline.r)) / 2;
    equi_z0(k)    = (max(equi.time_slice{k}.boundary.outline.z) + min(equi.time_slice{k}.boundary.outline.z)) / 2;
    equi_a(k)     = geo_in.a;
    equi_K(k)     = geo_in.K;
    equi_d(k)     = geo_in.d;
    equi_b0(k)    = abs(equi.vacuum_toroidal_field.b0(k) .*  equi.vacuum_toroidal_field.r0) ./equi_R(k);
    equi_vp(k)    = equi.time_slice{k}.global_quantities.volume;
    equi_sp(k)    = equi.time_slice{k}.global_quantities.area;
    equi_sext(k)  = equi.time_slice{k}.global_quantities.surface;
    %
    fbe_lcfs{k}.R = equi.time_slice{k}.boundary.outline.r(:).';
    fbe_lcfs{k}.Z = equi.time_slice{k}.boundary.outline.z(:).';
end

% update values
% cons1t
input.cons1t.ip = interp1_ex(equi_time,equi_ip,input.cons1t.temps,'linear','extrap');
if (nargin > 3) && ~isempty(only_ip)
    return;
end
% input.cons1t.vloop not available from equilibrium but not used here
input.cons1t.flux = interp1_ex(equi_time,equi_flux,input.cons1t.temps,'linear','extrap');
% geo1t
input.geo1t.R     = interp1_ex(equi_time,equi_R,input.cons1t.temps,'linear','extrap'); 
input.geo1t.z0    = interp1_ex(equi_time,equi_z0,input.cons1t.temps,'linear','extrap');
input.geo1t.a     = interp1_ex(equi_time,equi_a,input.cons1t.temps,'linear','extrap');
input.geo1t.K     = interp1_ex(equi_time,equi_K,input.cons1t.temps,'linear','extrap');
input.geo1t.d     = interp1_ex(equi_time,equi_d,input.cons1t.temps,'linear','extrap');
input.geo1t.b0    = interp1_ex(equi_time,equi_b0,input.cons1t.temps,'linear','extrap');
input.geo1t.vp    = interp1_ex(equi_time,equi_vp,input.cons1t.temps,'linear','extrap');
input.geo1t.sp    = interp1_ex(equi_time,equi_sp,input.cons1t.temps,'linear','extrap');
input.geo1t.sext  = interp1_ex(equi_time,equi_sext,input.cons1t.temps,'linear','extrap');


% sepa1t if exist
if isfield(input,'sepa1t') && isfield(input.sepa1t,'Rsepa')
    input.sepa1t = interpol_lcfs(equi_time,fbe_lcfs,LCFS_nbpts,input.cons1t.temps);
end


function [sepa,z0_out] = interpol_lcfs(equi_time,fbe_lcfs,LCFS_nbpts,time_out)


% initialisation
uni =  linspace(0, 2.* pi,LCFS_nbpts)';
r0  = NaN * ones(1,length(equi_time));
z0  = NaN * ones(1,length(equi_time));
rho = NaN * ones(length(uni),length(equi_time));
equi_time = equi_time(:)';

for k= 1:length(equi_time)

	% extraction de la sepa
	r = fbe_lcfs{k}.R;
	z = fbe_lcfs{k}.Z;
	r = r(:);
	z = z(:);
    indbad = find(r <= 0);
    if ~isempty(indbad)
        r(indbad) = [];
        z(indbad) = [];
    end
% 	KH = sort(unique(convhull(r,z)));
% 	if (length(KH) ~= length(r))
% 	    index_full = 1:length(r);
% 	    r = r(KH);
% 	    z = z(KH);
% 	    r = interp1(KH,r,index_full,'linear');
% 	    z = interp1(KH,z,index_full,'linear');
% 	    indbad_lcfs = find(~isfinite(r) | ~isfinite(z));
% 	    if ~isempty(indbad_lcfs)
% 		r(indbad_lcfs) = [];
% 		z(indbad_lcfs) = [];
% 	    end
% 	    r = r(:);
% 	    z = z(:);
% 	end

	r0(k)   = (min(r) + max(r)) ./ 2;
	%mask    = (r == max(r));
	%z0(k)   = sum(z .* mask) ./ max(1,sum(mask));
	z0(k)   = (min(z) + max(z)) ./ 2;
	cc   = (r - r0(k)) + sqrt(-1) .* (z - z0(k));
	thc  = unwrap(angle(cc));
	thc(thc <0) = thc(thc<0) + 2 .* pi;
	rhoc = abs(cc);
	[thc,indc] = sort(thc);
	rhoc       = rhoc(indc);
	rhoc = cat(1,rhoc,rhoc,rhoc);
	thc = cat(1,thc -2.*pi,thc,thc+2.*pi);
	indnok = find(diff(thc)<=0);
	thc(indnok) =[];
	rhoc(indnok)  = [];
	rho(:,k) = spline(thc,rhoc,uni);     
end

if time_out < min(equi_time)
    % adding LCFS for first time slice
    time_ini  = equi_time(1) - max(1e-6,mean(diff(equi_time)));
    r0        = cat(2,r0(1),r0);
    z0        = cat(2,z0(1),z0);
    equi_time  = cat(2,time_ini,equi_time);
    rho       = cat(2,rho(:,1),rho);
end
if time_out > max(equi_time)
    % adding LCFS for last time slice
    time_ini  = equi_time(end) + max(1e-6,mean(diff(equi_time)));
    r0        = cat(2,r0,r0(end));
    z0        = cat(2,z0,z0(end));
    equi_time  = cat(2,equi_time,time_ini);
    rho       = cat(2,rho,rho(:,end));
end
% resampling in time
rho_out  = interp1_ex(equi_time',rho',time_out,'linear','extrap')';
r0_out   = interp1_ex(equi_time',r0',time_out,'linear','extrap')';
z0_out   = interp1_ex(equi_time',z0',time_out,'linear','extrap')';

% new LCFS for METIS
sepa.Rsepa  = (ones(size(uni)) * r0_out + rho_out .* cos(uni * ones(1,size(rho_out,2))))';
sepa.Zsepa  = (rho_out .* sin(uni * ones(1,size(rho_out,2))))';





