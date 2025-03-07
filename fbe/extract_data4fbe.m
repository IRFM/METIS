% extract data for FBE in inverse evolutive mode
function metis_data4fbe = extract_data4fbe(post,nbtime,plotonoff,tmin,tmax)

% test entry
if nargin < 3
  plotonoff = 0;
end

% using energy content (thermal, magnetic, fast and rotation)
dwidt = abs(z0dxdt(post.zerod.w,post.zerod.temps)) + abs(z0dxdt(post.zerod.wbp,post.zerod.temps)) + ...
        abs(z0dxdt(post.zerod.wrot,post.zerod.temps)) + abs(z0dxdt(post.zerod.wmagtor,post.zerod.temps));
wi = cumtrapz(post.zerod.temps,dwidt);
% vertical field
mu0 = 4 * pi * 1e-7;
bv = mu0 .* post.zerod.ip ./ (4 .* pi .* post.z0dinput.geo.R) .* (log(8 .*  post.z0dinput.geo.R ./ ...
((post.zerod.sp/pi).^0.5)) + post.zerod.betaptot + post.zerod.li ./ 2 - 3/2);
bv = cumtrapz(post.zerod.temps,abs(z0dxdt(post.zerod.w,post.zerod.temps)));
%  normalisation
wi = (wi - min(wi)) ./ (max(wi) - min(wi)) +  ...
     (post.zerod.temps - min(post.zerod.temps)) ./ (max(post.zerod.temps) - min(post.zerod.temps)) + ...
     (bv - min(bv)) ./ (max(bv) - min(bv));
% selection of time slices
% loop to have at the end the requested number of time slices
for k=0:100
    % selection of time slices
    if isempty(nbtime) || ~isfinite(nbtime)
	wit = wi; 
	time_slices = post.zerod.temps;
	nbtime = length(time_slices);
	nbtime_loc = nbtime;
	break;
    else
	nbtime_loc = nbtime + k;    
    end 
    if nargin == 5
        wi_min = interp1(post.zerod.temps,wi,max(min(post.zerod.temps),tmin),'nearest',min(wi));
        wi_max = interp1(post.zerod.temps,wi,min(max(post.zerod.temps),tmax),'nearest',max(wi));
        wit    = linspace(wi_min,wi_max,nbtime_loc)';
    elseif nargin == 4
        wit = linspace(0,tmin,nbtime_loc)'; %HH
    else
        wit = linspace(0,max(wi),nbtime_loc)';
    end
    time_slices = unique(interp1(wi,post.zerod.temps,wit,'nearest','extrap'));
    %disp([nbtime,length(time_slices)])
    %time_slices = unique(union(time_slices,cat(1,min(post.zerod.temps),max(post.zerod.temps))));
    if length(time_slices) >= nbtime
        break;
    end
end
if length(time_slices) ~= nbtime
    fprintf('Requested number of time slices (%d) differs from obtained number of time slices (%d)\n',nbtime,length(time_slices));
end
%time_slices = interp1(wi,post.zerod.temps,wit,'nearest','extrap');
%time_slices = unique(union(time_slices,cat(1,min(post.zerod.temps),max(post.zerod.temps))));
wic         = interp1(post.zerod.temps,wi,time_slices,'nearest','extrap');

% time dependant data
metis_data4fbe.time       = time_slices;
metis_data4fbe.ip         = interp1_imas(post.zerod.temps,post.zerod.ip,time_slices,'nearest','extrap');
metis_data4fbe.vloop      = interp1_imas(post.zerod.temps,post.zerod.vloop,time_slices,'nearest','extrap');
metis_data4fbe.li         = interp1_imas(post.zerod.temps,post.zerod.li,time_slices,'nearest','extrap');
metis_data4fbe.betap      = interp1_imas(post.zerod.temps,post.zerod.betaptot,time_slices,'nearest','extrap');
metis_data4fbe.psi        = interp1_imas(post.profil0d.temps,post.profil0d.psi,time_slices,'nearest','extrap');
metis_data4fbe.q          = interp1_imas(post.profil0d.temps,post.profil0d.qjli,time_slices,'nearest','extrap');
metis_data4fbe.rho        = interp1_imas(post.profil0d.temps,post.profil0d.rmx,time_slices,'nearest','extrap');
metis_data4fbe.phi        = interp1_imas(post.profil0d.temps,post.profil0d.phi,time_slices,'nearest','extrap');
metis_data4fbe.epar       = interp1_imas(post.profil0d.temps,post.profil0d.epar,time_slices,'nearest','extrap');
metis_data4fbe.ptot       = interp1_imas(post.profil0d.temps,post.profil0d.ptot,time_slices,'nearest','extrap');
metis_data4fbe.dptotdpsi  = interp1_imas(post.profil0d.temps,post.profil0d.dptotdpsi,time_slices,'nearest','extrap');
metis_data4fbe.fdia       = interp1_imas(post.profil0d.temps,post.profil0d.fdia,time_slices,'nearest','extrap');
metis_data4fbe.df2dpsi    = interp1_imas(post.profil0d.temps,post.profil0d.df2dpsi,time_slices,'nearest','extrap');
metis_data4fbe.jphi       = interp1_imas(post.profil0d.temps,post.profil0d.jli,time_slices,'nearest','extrap');
metis_data4fbe.rmx        = interp1_imas(post.profil0d.temps,post.profil0d.rmx,time_slices,'nearest','extrap');
metis_data4fbe.C2         = interp1_imas(post.profil0d.temps,post.profil0d.C2,time_slices,'nearest','extrap');
metis_data4fbe.Raxe       = interp1_imas(post.profil0d.temps,post.profil0d.Raxe,time_slices,'nearest','extrap');
metis_data4fbe.ri         = interp1_imas(post.profil0d.temps,post.profil0d.ri,time_slices,'nearest','extrap');
metis_data4fbe.r2i        = interp1_imas(post.profil0d.temps,post.profil0d.r2i,time_slices,'nearest','extrap');
metis_data4fbe.psin       = (metis_data4fbe.psi - metis_data4fbe.psi(:,1) * ones(1,size(metis_data4fbe.psi,2))) ./ ...
                              ((metis_data4fbe.psi(:,end) - metis_data4fbe.psi(:,1)) * ones(1,size(metis_data4fbe.psi,2)));

% alternative derivation of dPtotdPsi and df2dpsi
% improved value on magnetic axis
psid1             = pdederive(post.profil0d.xli,metis_data4fbe.psi,0,2,2,1);
psid1(:,end)      = -(2*pi) .* mu0 .* metis_data4fbe.rmx(:,end) .* metis_data4fbe.ip ./ metis_data4fbe.C2(:,end);
% dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
psid2    = pdederive(post.profil0d.xli,metis_data4fbe.psi,1,0,2,2);
dpdpsi   = pdederive(post.profil0d.xli,metis_data4fbe.ptot,0,1,2,1) ./ min(-eps,psid1);
dpdpsi(psid1 > -eps) = 0;
dpdpsi(:,1) = max(0,dpdpsi(:,1)); 
% estimation simple de la valeur sur l'axe : (dP/dpsi ~= 0 sur l'axe meme si dP/drho = 0)
% peu d'incidence sur l'integrale
dpdx_1 = (metis_data4fbe.ptot(:,2) - metis_data4fbe.ptot(:,1)) ./ post.profil0d.xli(2) .^ 2;
dpdpsi_1 = - 2 .* dpdx_1 ./ (metis_data4fbe.fdia(:,1) ./ metis_data4fbe.Raxe(:,1)) .* metis_data4fbe.q(:,1) ./ metis_data4fbe.rmx(:,end) .^ 2;
dpdpsi(:,1) = max(dpdpsi(:,1),dpdpsi_1);
df2dpsi  = 2 .* mu0 .* (max(0,metis_data4fbe.jphi) .* metis_data4fbe.ri - dpdpsi)./ metis_data4fbe.r2i;
% direct derivation
for k=1:size(metis_data4fbe.psi,1)
    dpdpsi_direct(k,:) =   pdederive(metis_data4fbe.psi(k,:),metis_data4fbe.ptot(k,:),2,2,2,1);
    df2dpsi_direct(k,:) =   pdederive(metis_data4fbe.psi(k,:),metis_data4fbe.fdia(k,:) .^ 2,2,2,2,1);
end
jphi_rebuilt = (df2dpsi .* metis_data4fbe.r2i ./ (2 .* mu0)  + dpdpsi) ./ metis_data4fbe.ri;
if plotonoff ~= 0    
    figure;
    subplot(2,2,1)
    plot(metis_data4fbe.psin',metis_data4fbe.dptotdpsi','r',metis_data4fbe.psin',dpdpsi','b',metis_data4fbe.psin',dpdpsi_direct','g');
    xlabel('\psi_n');
    ylabel('dPd\psi');
    title('red = METIS; blue = For FEEQS; green =Direct');
    hold on
    plot(metis_data4fbe.psin',metis_data4fbe.dptotdpsi','.r',metis_data4fbe.psin',dpdpsi','.b',metis_data4fbe.psin',dpdpsi_direct','.g');
    subplot(2,2,2)
    plot(metis_data4fbe.psin',metis_data4fbe.df2dpsi','r',metis_data4fbe.psin',df2dpsi','b',metis_data4fbe.psin',df2dpsi_direct','g');
    xlabel('\psi_n');
    ylabel('dF^2d\psi');
    hold on
    plot(metis_data4fbe.psin',metis_data4fbe.df2dpsi','.r',metis_data4fbe.psin',df2dpsi','.b',metis_data4fbe.psin',df2dpsi_direct','.g');
    subplot(2,2,3);
    plot(metis_data4fbe.psin',metis_data4fbe.jphi','r',metis_data4fbe.psin',jphi_rebuilt','b');
    xlabel('\psi_n');
    ylabel('J_\phi');
    hold on
    plot(metis_data4fbe.psin',metis_data4fbe.jphi','.r',metis_data4fbe.psin',jphi_rebuilt','.b');
    

end
metis_data4fbe.dptotdpsi = dpdpsi;
metis_data4fbe.df2dpsi   = df2dpsi;

% LCFS
if isfield(post.z0dinput.exp0d,'Rsepa') && ~isempty(post.z0dinput.exp0d.Rsepa)
    metis_data4fbe.LCFS.R = interp1_imas(post.zerod.temps,post.z0dinput.exp0d.Rsepa,time_slices,'nearest','extrap');
    metis_data4fbe.LCFS.Z = interp1_imas(post.zerod.temps,post.z0dinput.exp0d.Zsepa,time_slices,'nearest','extrap');
else
    metis_data4fbe.LCFS.R = NaN * ones(size(time_slices),201);
    metis_data4fbe.LCFS.Z = NaN * ones(size(time_slices),201);
    for k= 1:length(time_slices)
            index        = find(post.z0dinput.cons.temps >= time_slices(1),1);
	    geo_metis    = zerod_get1t(post.z0dinput.geo,index);
	    [vps_next,void_sps,void_sexts,void_peris,geo_metis,void_xpoint,sepa_metis.Rsepa,sepa_metis.Zsepa] =  ...
		    zgeo0(geo_metis,[],[],1);  
	    sepa_metis.Zsepa = sepa_metis.Zsepa  + geo_metis.z0;  
	    metis_data4fbe.LCFS.R(k,:) = sepa_metis.Rsepa(:);
	    metis_data4fbe.LCFS.Z(k,:) = sepa_metis.Zsepa(:); 
    end
end

% flux managment 
flux_data_output_ = poloidal_flux(post,0.33);
% break down duration
metis_data4fbe.breakdown.dt_break = flux_data_output_.dt_break;
if ~isempty(post.z0dinput.option.available_flux) && isfinite(post.z0dinput.option.available_flux)
    available_flux = post.z0dinput.option.available_flux;
else
    available_flux = 2 .* pi .*  (metis_data4fbe.psi(end,end) - metis_data4fbe.psi(1,end)) + flux_data_output_.breakdown_flux + max(flux_data_output_.flux_bv); 
   % available_flux = 41; % default JT-60SA
end
breakdown_flux = flux_data_output_.breakdown_flux;
unsymetry  = max(flux_data_output_.flux_bv);
flux_plus  = (available_flux - breakdown_flux - unsymetry) ./ 2;
flux_break = flux_plus + breakdown_flux;
flux_moins = (available_flux - breakdown_flux + unsymetry) ./ 2;
flux_induc = interp1_imas(post.zerod.temps,flux_data_output_.flux_induc,time_slices(1),'nearest','extrap');
flux_res   = interp1_imas(post.zerod.temps,flux_data_output_.flux_res,time_slices(1),'nearest','extrap');
flux_start = flux_plus - flux_induc - flux_res;
%
if length(flux_start) == 1
    metis_data4fbe.psi = metis_data4fbe.psi - metis_data4fbe.psi(1,end) + flux_start ./ 2 ./ pi;
else
    metis_data4fbe.psi = metis_data4fbe.psi - metis_data4fbe.psi(1,end) + (flux_start * ones(1,size(metis_data4fbe.psi,2))) ./ 2 ./ pi;
end

% from Optimization of Plasma Initiation Scenarios in JT-60SA Makoto MATSUKAWA et al, J. Plasma Fusion Res. SERIES, Vol. 9 (2010)
if ~isempty(strfind(upper(post.z0dinput.machine),'JT-60SA'))
    % test of compatibility
    if flux_break > 17.8
        flux_avail_real = 2 .* 17.8 + unsymetry;
        fprintf('Flux requested in METIS simulation is greater than real available flux: METIS = %g & JT-60SA = %g\n',available_flux,flux_avail_real);
    end    
    metis_data4fbe.init_current.CS1 = 19.58e3 ./  17.8 .* flux_start;
    metis_data4fbe.init_current.CS2 = 19.07e3 ./  17.8 .* flux_start; 
    metis_data4fbe.init_current.CS3 = 19.02e3 ./  17.8 .* flux_start;
    metis_data4fbe.init_current.CS4 = 18.09e3 ./  17.8 .* flux_start;
    metis_data4fbe.init_current.EF1 = 1.367e3 ./  17.8 .* flux_start;
    metis_data4fbe.init_current.EF2 = -0.1112e3 ./  17.8 .* flux_start;
    metis_data4fbe.init_current.EF3 = 16.17e3 ./  17.8 .* flux_start;
    metis_data4fbe.init_current.EF4 = 11.38e3 ./  17.8 .* flux_start;
    metis_data4fbe.init_current.EF5 = -0.8427e3 ./  17.8 .* flux_start;
    metis_data4fbe.init_current.EF6 = 1.298e3 ./  17.8 .* flux_start;
    % breakdown current
    metis_data4fbe.breakdown.init_current.CS1 = 19.58e3 ./  17.8 .* flux_break;
    metis_data4fbe.breakdown.init_current.CS2 = 19.07e3 ./  17.8 .* flux_break; 
    metis_data4fbe.breakdown.init_current.CS3 = 19.02e3 ./  17.8 .* flux_break;
    metis_data4fbe.breakdown.init_current.CS4 = 18.09e3 ./  17.8 .* flux_break;
    metis_data4fbe.breakdown.init_current.EF1 = 1.367e3 ./  17.8 .* flux_break;
    metis_data4fbe.breakdown.init_current.EF2 = -0.1112e3 ./  17.8 .* flux_break;
    metis_data4fbe.breakdown.init_current.EF3 = 16.17e3 ./  17.8 .* flux_break;
    metis_data4fbe.breakdown.init_current.EF4 = 11.38e3 ./  17.8 .* flux_break;
    metis_data4fbe.breakdown.init_current.EF5 = -0.8427e3 ./  17.8 .* flux_break;
    metis_data4fbe.breakdown.init_current.EF6 = 1.298e3 ./  17.8 .* flux_break;
end
% intsert data for other tokamake here

% break down duration
metis_data4fbe.dt_break = flux_data_output_.dt_break;

% plot of time selection
if plotonoff == 1
  figure;
  subplot(2,2,1);
  R = (max(metis_data4fbe.LCFS.R,[],2) + min(metis_data4fbe.LCFS.R,[],2)) /2;
  a = (max(metis_data4fbe.LCFS.R,[],2) - min(metis_data4fbe.LCFS.R,[],2)) /2;  
  plot(post.zerod.temps,wi,'b',time_slices,wic,'.r', ...
  post.z0dinput.cons.temps,post.z0dinput.geo.R,'c',time_slices,R,'.m', ...
  post.z0dinput.cons.temps,post.z0dinput.geo.a,'g',time_slices,a,'.k');
  title('extraction data from METIS for FEEQS in iverse evolutime mode');
  ylabel('Varition index, R(m) & a(m)');
  subplot(2,2,2);
  plot(metis_data4fbe.LCFS.R',metis_data4fbe.LCFS.Z');
  xlabel('R (m)');
  ylabel('R (m)');
  title('LCFS sequence');
  axis('equal');
  subplot(2,2,3)
  plot(post.zerod.temps,post.zerod.ip ./1e6,'b',metis_data4fbe.time,metis_data4fbe.ip ./1e6,'.k'); 
  xlabel('time (s)');
  ylabel('Ip (MA)');
  subplot(2,2,4)
  plot(post.zerod.temps,post.zerod.li,'c',metis_data4fbe.time,metis_data4fbe.li,'.k'); 
  hold on
  plot(post.zerod.temps,post.zerod.betaptot,'m',metis_data4fbe.time,metis_data4fbe.betap,'.k'); 
  xlabel('time (s)');
  ylabel('l_i & \beta_P');
  edition2
end

function flux_data_output_ = poloidal_flux(post,couplage)
        
        % z0plotflux.m
        NOPLUTFLUX_HELIOS = 1;
        z0plotflux;
        
        % create data structure
        liste_of_fields_ = who;
        for kkk_ = 1:length(liste_of_fields_)
            flux_data_output_.(liste_of_fields_{kkk_}) = eval(liste_of_fields_{kkk_});
        end
 
