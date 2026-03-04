function input = metis4muscle3_update_metis_input_on_equilibrium(input,equi,only_ip)

keyboard

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
%
for k=1:length(equi_time)
    equi_ip(k)    = abs(equi.time_slice{k}.global_quantitie.ip);
    equi_flux(k)  = - equi.time_slice{k}.global_quantitie.psi_boundary / (2*pi); % data are in COCOS 11 and METIS natively COCOS 7
    % geometry
    [vps,sps,sexts,peris,geo_in,xpoint,R,Z] = zgeo0(input.geo1t, ...
                   equi.time_slice{k}.boundary.ouline.r, ...
                   equi.time_slice{k}.boundary.ouline.z,1)
    equi_R(k)     = geo_in.R;
    equi_z0(k)    = geo_in.z0;
    equi_a(k)     = geo_in.a;
    equi_K(k)     = geo_in.K;
    equi_d(k)     = geo_in.d;
    equi_b0(k)    = abs(equi.vacuum_toroidal_field.b0(k) .*  equi.vacuum_toroidal_field.r0) ./equi_R(k);
    equi_vp(k)    = equi.time_slice{k}.global_quantitie.volume;
    equi_sp(k)    = equi.time_slice{k}.global_quantitie.area;
    equi_sext(k)  = equi.time_slice{k}.global_quantitie.surface;
end

% update value
% cons1t
input.cons1t.ip
if nargin > 2
    return;
end
input.cons1t.vloop
input.cons1t.flux
% geo1t
input.geo1t

ans = 

  struct with fields:

       R: 6.6800
      z0: 0.6550
       a: 1.6000
       K: 0.9962
       d: 1.3323e-15
      b0: 4.9192
      vp: NaN
      sp: NaN
    sext: NaN


% sepa1t if exist
if isfield(input,'sepa1t') && isfield(input.sepa1t,'Rsepa')
    
    
end

