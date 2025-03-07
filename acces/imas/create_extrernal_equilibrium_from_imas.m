function [data,equi] = create_extrernal_equilibrium_from_imas(time_out,time_offset,shot, imas_run, imas_occurrence, imas_user, imas_machine,imas_backend,flag_extrap_data_flag)

% initialisation
data = [];
% constant
mu0 = 4 * pi * 1e-7;

if nargin < 9
    flag_extrap_data_flag = false;
end
%
if nargin < 3
    error('syntax: equi_ext = create_extrernal_equilibrium_from_imas(time_out,time_offset,shot, imas_run, imas_occurrence, imas_user, imas_machine)');
elseif isstruct(shot)
    equi = shot;
    if nargin > 3
        flag_extrap_data_flag = imas_run;
    end
else
    if (nargin < 4) || isempty(imas_run)
        imas_run = 0;
    end
    if (nargin < 5) || isempty(imas_occurrence)
        imas_occurrence = 0;
    end
    if (nargin < 6) || isempty(imas_user)
        imas_user = '';
    end
    if (nargin < 7) || isempty(imas_machine)
        imas_machine = '';
    end
    if (nargin < 8) || isempty(imas_backend)
        imas_backend  = '';
    end
    imas_version = '3'; 
    [imas_machine,imas_user,imas_version,imas_backend] = imasdb(imas_machine,imas_user,imas_version,imas_backend);
    id_backend = get_imas_backend_id(imas_backend)
    if ~isempty(id_backend)
        try
            idx = imas_open_env_backend(shot, imas_run, imas_user, imas_machine, imas_version,imas_backend);
        catch
            % go back to MDS+ if not available
            idx = imas_open_env('ids', shot, imas_run, imas_user, imas_machine, imas_version);
            fprintf('selected backend (%d)is nos available switching back to default backend: %d\n',id_backend,imas_get_backendID(idx));            
        end
   else
        idx = imas_open_env('ids', shot, imas_run, imas_user, imas_machine, imas_version);
    end
    if imas_occurrence == 0
        equi = ids_get(idx,'equilibrium');
    else
        equi = ids_get(idx,sprintf('equilibrium/%d',imas_occurrence));
    end
    imas_close(idx);
end

if isempty(equi)
    error('No data for equilibrium');
end

if isempty(time_offset)
    time_offset = 0;
end

% get equilibrium time 
if equi.ids_properties.homogeneous_time == 0
    time_equi = NaN * ones(length(equi.time_slice),1);
    for k = 1:length(equi.time_slice)
        time_equi(k) =  equi.time_slice{k}.time;
    end
else
    time_equi = equi.time;
end
if length(equi.code.output_flag) == length(time_equi)
    indok = find(equi.code.output_flag == 0);
else
    indok = 1:length(time_equi);
end

%lookup table
if isempty(time_out)
    time_out = time_equi(indok) - time_offset;
    ind_out  = indok;
else
    ind_out = round(interp1(time_equi(indok) - time_offset,indok,time_out,'nearest','extrap'));
end

% data fields initialisation
nbt                = length(ind_out);
nbx                = length(equi.time_slice{indok(1)}.profiles_1d.rho_tor_norm);
v0d                = NaN * ones(nbt,1);
v1d                = NaN * ones(nbt,nbx);
data.x             = v1d;
data.time          = time_out;
data.psi           = v1d ;
data.phi           = v1d;
data.dphidt        = v1d;
data.psid1         = v1d;
data.psid2         = v1d;
data.dpsidt        = v1d;
data.qjli          = v1d;
data.jli           = v1d;
data.epar          = v1d;
data.bpol          = v1d;
data.rmx           = v1d;
data.rm            = v0d;
data.drmdt         = v0d;
data.grho2r2       = v1d;
data.r2i           = v1d;
data.ri            = v1d;
data.grho          = v1d;
data.grho2         = v1d;
data.spr_tor       = v1d;
data.vpr_tor       = v1d; 
data.volume        = v1d;
data.surface_pol   = v1d;
data.surface_flux  = v1d;
data.dphidx        = v1d;
data.vpr           = v1d;
data.spr           = v1d;
data.C2            = v1d;
data.C3            = v1d;
data.ej            = v1d;
data.jeff          = v1d;
data.jres          = v1d;
data.Raxe          = v1d;
data.Zaxe          = v1d;
data.epsi          = v1d;
data.a             = v1d;
data.difcurconv    = v0d;
data.df2dpsi       = v1d;
data.dptotdpsi     = v1d;
data.ptot          = v1d;
data.fdia          = v1d;
data.wbp           = v0d;
data.dwbpdt        = v0d;
data.kx            = v1d;
data.dx            = v1d;
data.ipout         = v0d;
data.lif           = v0d;
data.betap         = v0d;
data.pohm          = v0d;
data.q95           = v0d;
data.qmin          = v0d;
data.q0            = v0d;
data.psin          = v1d;
data.d95           = v0d;                
data.K95           = v0d;   
data.piqj          = v0d;
data.geo.a         = v0d;
data.geo.R         = v0d;
data.geo.z0        = v0d;
data.geo.K         = v0d;
data.geo.d         = v0d;
data.geo.b0        = v0d;
data.geo.vp        = v0d;
data.geo.sp        = v0d;
data.geo.peri      = v0d;
data.geo.sext      = v0d;
data.geo.xpoint    = v0d;
data.geo.Rsepa     = cell(nbt,1);
data.geo.Zsepa     = cell(nbt,1);
data.geo.peri      = v0d;
data.indice_inv    = v0d;
data.poynting      = v0d;
data.ftrap         = v1d;

% use all valid data to compute time derivative
psi_time  = NaN * ones(length(indok),nbx);
phi_time  = NaN * ones(length(indok),nbx);
rho_time  = NaN * ones(length(indok),nbx);
raxe_time = NaN * ones(length(indok),nbx);
zaxe_time = NaN * ones(length(indok),nbx);
bpol_time = NaN * ones(length(indok),nbx);
wbp_time  = NaN * ones(length(indok),1);
for k = 1:length(indok)
    psi_time(k,:)  = - equi.time_slice{indok(k)}.profiles_1d.psi ./ (2*pi);
    phi_time(k,:)  =   equi.time_slice{indok(k)}.profiles_1d.phi;
    rho_time(k,:)  =   equi.time_slice{indok(k)}.profiles_1d.rho_tor;
    raxe_time(k,:) =   (equi.time_slice{indok(k)}.profiles_1d.r_outboard' + equi.time_slice{indok(k)}.profiles_1d.r_inboard') ./ 2;
    wbp_time(k)    =   mu0 .* equi.time_slice{indok(k)}.global_quantities.ip .^ 2 .* raxe_time(k,end) ./ 4 .* equi.time_slice{indok(k)}.global_quantities.li_3;
    bpol_time(k,:) =   sqrt(equi.time_slice{indok(k)}.profiles_1d.gm2 .* equi.time_slice{indok(k)}.profiles_1d.dpsi_drho_tor .^ 2) ./ 2 ./ pi;
end
if (length(indok) > 1)
    dpsidt_time = z0dxdt_c(psi_time,time_equi(indok));
    dphidt_time = z0dxdt_c(phi_time,time_equi(indok));
    drhodt_time = z0dxdt_c(rho_time,time_equi(indok));
    if ~flag_extrap_data_flag
        dpsidt_time(end,:) = dpsidt_time(end - 1,:);
        dphidt_time(end,:) = dphidt_time(end - 1,:);
        drhodt_time(end,:) = drhodt_time(end - 1,:);
    end
    %
    dwbpdt_time          = z0dxdt_c(wbp_time,time_equi(indok));
    if ~flag_extrap_data_flag
        dwbpdt_time(end)     = dwbpdt_time(end-1);
    end
    dEBpoldt_time        = z0dxdt_c(bpol_time .^ 2 ./ 2 ./ mu0,time_equi(indok));
    if ~flag_extrap_data_flag
        dEBpoldt_time(end,:) = dEBpoldt_time(end-1,:);
    end
else
    dpsidt_time   = zeros(size(psi_time));
    dphidt_time   = zeros(size(phi_time));
    drhodt_time   = zeros(size(rho_time));   
    dwbpdt_time   = zeros(size(wbp_time));
    dEBpoldt_time = zeros(size(bpol_time));
end
drmdt_time = drhodt_time(:,end);
%
time_time = time_equi(indok) - time_offset;
%lookup table
ind_time = round(interp1(time_time,1:length(indok),time_out,'nearest','extrap'));

    
% loop to fill the field
% time derivative in indexed on m
for k=1:length(ind_out)
    l   = ind_out(k);
    m   = ind_time(k);
    tsl = equi.time_slice{l}.profiles_1d;
    %
    jphi_1d = tsl.j_tor(:)';
    Rim     = tsl.gm9(:)';
    Rm      = tsl.gm8(:)';
    volume  = tsl.volume(:)';
    area    = tsl.area(:)';
    ephi    = - Rim .* dpsidt_time(m,:);
    ej      = jphi_1d .* ephi;
    pohm_equi1d = sum(diff(volume,1,2) .* (ej(1:end-1) + ej(2:end)) ./ 2,2);
    ip_equi1d   = sum(diff(area,1,2) .* (jphi_1d(1:end-1) + jphi_1d(2:end)) ./ 2,2);
    vres_equi1d = pohm_equi1d ./ ip_equi1d;
    %
    % computation of derived quantities
    aminor = (tsl.r_outboard' - tsl.r_inboard') ./ 2;
    raxe   = (tsl.r_outboard' + tsl.r_inboard') ./ 2;
    if isfield(tsl,'geometric_axis') && ~isempty(tsl.geometric_axis.z)
        zaxe   = tsl.geometric_axis.z';
    else
        zaxe   = zeros(size(raxe));
    end
    epsi   = aminor ./ raxe;
    x      = aminor ./ max(eps,aminor(end));
    %
    Bt     = abs(equi.vacuum_toroidal_field.r0 .* equi.vacuum_toroidal_field.b0(l) ./ raxe(end));
    rmx    = sqrt(abs(phi_time(m,:) ./ pi ./ Bt ));
    % dpsidx s'annule au centre
    psid1    = pdederive(x,psi_time(m,:),0,2,2,1);
    % dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
    psid2    = pdederive(x,psi_time(m,:),1,0,2,2);
    psid2(1) = - Bt ./ abs(tsl.q(1)) .* rmx(end) .^ 2;
    %
    bpol  = sqrt(tsl.gm2 .* tsl.dpsi_drho_tor .^ 2) ./ 2 ./ pi;
    
    % create external data structure 
    data.x(k,:)         = x;
    data.psi(k,:)       = psi_time(m,:);
    data.phi(k,:)       = abs(phi_time(m,:));
    data.dphidt(k,:)    = dphidt_time(m,:);
    data.psid1(k,:)     = psid1;
    data.psid2(k,:)     = psid2;
    data.dpsidt(k,:)    = dpsidt_time(m,:);
    data.qjli(k,:)      = abs(tsl.q);
    data.jli(k,:)       = jphi_1d;
    data.epar(k,:)      = ephi;
    data.bpol(k,:)      = bpol;
    data.rmx(k,:)       = rmx;
    data.rm(k)          = max(rmx);
    data.drmdt(k)       = drmdt_time(m);
    data.grho2r2(k,:)   = tsl.gm2;
    data.r2i(k,:)       = tsl.gm1;
    data.ri(k,:)        = tsl.gm9;
    data.grho(k,:)      = abs(tsl.gm7);
    data.grho2(k,:)     = abs(tsl.gm3);
    data.spr_tor(k,:)   = tsl.darea_drho_tor;
    if any(tsl.darea_drho_tor == 0)
        mask0 = (tsl.darea_drho_tor == 0);
        data.spr_tor(k,mask0) =  tsl.gm9(mask0) .* tsl.dvolume_drho_tor(mask0) ./ (2*pi);
    end
    data.vpr_tor(k,:)       = tsl.dvolume_drho_tor;
    data.volume(k,:)        = tsl.volume;
    data.surface_pol(k,:)   = tsl.area;
    data.surface_flux(k,:)  = tsl.surface;
    data.dphidx(k,:)        = pdederive(x,data.phi(k,:),0,2,2,1);
    data.vpr(k,:)           = pdederive(x,data.volume(k,:),0,2,2,1);
    data.spr(k,:)           = pdederive(x,data.surface_pol(k,:),0,2,2,1);
    data.C2(k,:)            = tsl.dvolume_drho_tor(:)' .* data.grho2r2(k,:);
    data.C3(k,:)            = tsl.dvolume_drho_tor(:)' .* data.r2i(k,:);
    data.ej(k,:)            = ej;
    data.jeff(k,:)          = tsl.j_parallel;
    data.jres(k,:)          = NaN * data.jeff(k,:); % must be computed with METIS data
    data.Raxe(k,:)          = raxe;
    data.Zaxe(k,:)          = zaxe;
    data.epsi(k,:)          = epsi;
    data.a(k,:)             = aminor;
    data.difcurconv(k)      = equi.time_slice{l}.convergence.iterations_n;
    data.df2dpsi(k,:)       = - (4*pi) .* tsl.f_df_dpsi;
    data.dptotdpsi(k,:)     = - (2*pi) .* tsl.dpressure_dpsi;
    data.ptot(k,:)          = abs(tsl.pressure);
    data.fdia(k,:)          = abs(tsl.f);
    data.wbp (k)            = wbp_time(m);
    data.dwbpdt(k)          = dwbpdt_time(m);
    data.kx(k,:)            = tsl.elongation(:)';
    data.dx(k,:)            = (tsl.triangularity_lower(:)' + tsl.triangularity_upper(:)') ./ 2;
    data.ipout(k)           = abs(equi.time_slice{l}.global_quantities.ip);
    data.lif(k)             = equi.time_slice{l}.global_quantities.li_3;
    data.betap(k)           = equi.time_slice{l}.global_quantities.beta_pol;
    data.pohm(k)            = pohm_equi1d;
    data.q95(k)             = abs(equi.time_slice{l}.global_quantities.q_95);
    data.qmin(k)            = min(data.qjli(k,:),[],2);
    data.q0(k)              = abs(equi.time_slice{l}.global_quantities.q_axis);
    data.psin(k,:)          = (data.psi(k,:) - data.psi(k,1)) ./ (data.psi(k,end) - data.psi(k,1));
    dd95                    = abs(data.psin(k,:) - 0.95);
    mask95                  = double(dd95 == min(dd95,[],2));
    data.d95(k)             = sum(data.dx(k,:) .* mask95,2) ./ max(1, sum(mask95,2));
    data.K95(k)             = sum(data.kx(k,:) .* mask95,2) ./ max(1, sum(mask95,2));
    data.piqj(k)            = max(0.1,min(10,data.jli(k,1) ./ (sum(diff(data.surface_pol(k,:),1,2) .*  ...
        (data.jli(k,1:end-1) + data.jli(k,2:end)) ./ 2,2) ./ data.surface_pol(k,end))));
    
    % LCFS
    bnd = equi.time_slice{l}.boundary; 
    data.geo.a(k)         = bnd.minor_radius;
    data.geo.R(k)         = bnd.geometric_axis.r;
    data.geo.z0(k)        = bnd.geometric_axis.z;
    data.geo.K(k)         = bnd.elongation;
    data.geo.d(k)         = bnd.triangularity;
    data.geo.b0(k)        = abs(equi.vacuum_toroidal_field.r0 .* equi.vacuum_toroidal_field.b0(l) ./ data.geo.R(k));
    data.geo.vp(k)        = data.volume(k,end);
    data.geo.sp(k)        = data.surface_pol(k,end);
    data.geo.sext(k)      = data.surface_flux(k,end);
    xpoint = false;
    for kl = 1:length(bnd.x_point)
        if ~isempty(bnd.x_point{kl})
            if ~isempty(bnd.x_point{kl}.r) && isfinite(bnd.x_point{kl}.r) && ...
               ~isempty(bnd.x_point{kl}.z) && isfinite(bnd.x_point{kl}.z) && ...
                (bnd.x_point{kl}.r > 0) 
                xpoint = true;
            end
        end
    end
    if xpoint
        data.geo.xpoint(k)    = 1;
    else
        data.geo.xpoint(k)    = 0;        
    end
    data.geo.Rsepa{k}     = bnd.outline.r;
    data.geo.Zsepa{k}    =  bnd.outline.z;
    rep = z0polygeom(data.geo.Rsepa{k},data.geo.Zsepa{k});
    data.geo.peri(k)      = rep(4);
    
    % inversion radius (no prescription on trigger value
    if min(data.qjli(k,:)) < 0.9
        mask = (data.qjli(k,:) <= 1);
        ggr  = (1:size(data.qjli,2));
        data.indice_inv(k) = min(length(data.qjli(k,:)),max(mask .* ggr,[],2));
    else
        data.indice_inv(k) = 0;
    end
    % Poynting flux
    % ref :Ohmic ???ux consumption during initial operation of the NSTX spherical torus
    % J.E. Menard et al, NF 2001 41 1197
    % ici on prend en compte que la composante est utile pour la consommation de flux poloidal.
    data.poynting(k) = trapz(data.rmx(k,:),data.vpr_tor(k,:) .* (dEBpoldt_time(m,:) + data.ej(k,:)),2);
    
    % effective trapped fraction
    if  ~isempty(tsl.trapped_fraction)
        data.ftrap(k,:) = tsl.trapped_fraction;
    else
        data.ftrap(k,:) = ftrap(data.x(k,:),data.epsi(k,:),data.psi(k,:),data.bpol(k,:),data.fdia(k,:), ...
                          data.Raxe(k,:),data.r2i(k,:),data.ri(k,:),data.grho2r2(k,:),data.grho2(k,:),data.a(k,:));
    end
    
end

% set same stucture for current diffusion and eqauilibrium
setappdata(0,'CURDIFF_EXP',data);
setappdata(0,'EQUILIBRIUM_EXP',data);


function ft = ftrap(x,epsi,psi,bpol,fdia,Raxe,r2i,ri,grho2r2,grho2,a)

b2m       = bpol .^ 2  + fdia .^ 2 .* r2i;
bm        = fdia .* ri .* (1 + 0.5 .* bpol ./ sqrt(grho2r2) ./ fdia .* sqrt(grho2));
% calcul de bmax
btor         = (fdia ./ (Raxe - a));
grho         = abs(1 ./ max(eps,abs(pdederive(x,Raxe - a,0,2,2,1))));
grho(1)      = 2 .* grho(2) - grho(3);
bpol         = abs(pdederive(x,psi,0,2,2,1))./ (Raxe - a) .* grho;
bmax         = sqrt(btor .^ 2 + bpol .^ 2);

% variable
h  = min(1,bm ./ bmax);
h2 = min(1,b2m ./ bmax .^ 2);

% expression de ftrap  Lin-Liu and Miller Phys. of Plasmas 2 (5) 1995
ftu          = 1 - h2 ./ h .^ 2 .* ( 1 - sqrt(1 - h) .* (1 + 0.5 .* h));
ftl          = 1 - h2 .* (1./ h2 - sqrt(1-sqrt(h2)) ./ h2  - sqrt(1 - sqrt(h2)) ./ 2 ./h);
ftc          = 1 - (1-epsi) .^ 2 ./ sqrt(1 - epsi .^ 2) ./ (1 + 1.46 .* sqrt(epsi));  % Wesson , Tokamak
fte          = 0.75 .* ftu + 0.25 .* ftl;
ft           = x .* fte + (1 - x) .* ftc;

ft_alt     = 1 - (1-epsi) .^ 2 ./ sqrt(1 - epsi .^ 2) ./ (1 + 1.46 .* sqrt(epsi));
ft(~isfinite(ft)) = ft_alt(~isfinite(ft));

function xp = z0dxdt_c(x,t)

% dimension des entrees
sx  = size(x);
nd  = length(sx);
t   = t(:);
ve = ones(1,size(x,2));
xp = diff(x,1,1) ./ diff(t*ve,1,1);
xp  = cat(1,xp(1,:),xp);