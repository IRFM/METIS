function equi_ext = prepare_equi_ext(option,cons,geo,exp0d,profil0d)

persistent psi_offset
if isempty(psi_offset)
    psi_offset = 0;
end

% data from METIS
temps = cons.temps;
x     = profil0d.xli;

% clear persistent interpollant
resample_time_space_fbe;
if isappdata(0,'EQUILIBRIUM_EXP') 
    cur_diff_ext = getappdata(0,'EQUILIBRIUM_EXP');
elseif isappdata(0,'CURDIFF_EXP')
    cur_diff_ext = getappdata(0,'CURDIFF_EXP');
else
    error('No available external data either for equilibrium or current diffusion');
end
x_in = cur_diff_ext.x;
time_in = cur_diff_ext.time;
% update psi_offset
if option.evolution == 0
    psi_LCFS = resample_time_fbe(time_in,cur_diff_ext.psi(:,end),temps);
    psi_offset = psi_LCFS(1);
end

%
% assigin field to data with resampling
equi_ext.time         = temps;
equi_ext.x            = x;
equi_ext.psi_offset   = psi_offset;
equi_ext.psi          = resample_time_space_fbe(time_in,x_in,cur_diff_ext.psi,temps,x) - psi_offset;
equi_ext.dpsidt       = resample_time_space_fbe(time_in,x_in,cur_diff_ext.dpsidt,temps,x);
equi_ext.phi_tor      = resample_time_space_fbe(time_in,x_in,cur_diff_ext.phi,temps,x);
equi_ext.dphidx_tor   = resample_time_space_fbe(time_in,x_in,cur_diff_ext.dphidx,temps,x);
equi_ext.qjli         = resample_time_space_fbe(time_in,x_in,cur_diff_ext.qjli,temps,x);
equi_ext.jli          = resample_time_space_fbe(time_in,x_in,cur_diff_ext.jli,temps,x);
equi_ext.jeff         = resample_time_space_fbe(time_in,x_in,cur_diff_ext.jeff,temps,x);
equi_ext.epar         = resample_time_space_fbe(time_in,x_in,cur_diff_ext.epar,temps,x);
equi_ext.ej           = resample_time_space_fbe(time_in,x_in,cur_diff_ext.ej,temps,x);
equi_ext.fdia         = resample_time_space_fbe(time_in,x_in,cur_diff_ext.fdia,temps,x);
equi_ext.bpol         = resample_time_space_fbe(time_in,x_in,cur_diff_ext.bpol,temps,x);
%
equi_ext.rmx          = resample_time_space_fbe(time_in,x_in,cur_diff_ext.rmx,temps,x);
equi_ext.Raxe         = resample_time_space_fbe(time_in,x_in,cur_diff_ext.Raxe,temps,x);
if isfield(cur_diff_ext,'Zaxe')
    equi_ext.Zaxe     = resample_time_space_fbe(time_in,x_in,cur_diff_ext.Zaxe,temps,x); 
    % internaly METIS use z0 = 0.
    equi_ext.Zaxe     = equi_ext.Zaxe - equi_ext.Zaxe(:,end) * ones(1,size(equi_ext.Zaxe,2));
else
    warning('In external equilibrium, the field Zaxe is not defined');
    equi_ext.Zaxe     = zeros(size(equi_ext.Raxe));
end
equi_ext.epsi         = resample_time_space_fbe(time_in,x_in,cur_diff_ext.epsi,temps,x);
equi_ext.a            = resample_time_space_fbe(time_in,x_in,cur_diff_ext.a,temps,x);
equi_ext.vpr_tor      = resample_time_space_fbe(time_in,x_in,cur_diff_ext.vpr_tor,temps,x);
equi_ext.spr_tor      = resample_time_space_fbe(time_in,x_in,cur_diff_ext.spr_tor,temps,x);
equi_ext.volume       = resample_time_space_fbe(time_in,x_in,cur_diff_ext.volume,temps,x);
equi_ext.vpr          = pdederive(x,equi_ext.volume,0,2,2,1);
equi_ext.vpr          = equi_ext.vpr .* ((equi_ext.volume(:,end) ./ max(eps,trapz(x,equi_ext.vpr,2))) * ones(size(x)));
equi_ext.surface_pol  = resample_time_space_fbe(time_in,x_in,cur_diff_ext.surface_pol,temps,x);
equi_ext.spr          = pdederive(x,equi_ext.surface_pol,0,2,2,1);
equi_ext.spr          = equi_ext.spr .* ((equi_ext.surface_pol(:,end) ./ max(eps,trapz(x,equi_ext.spr,2))) * ones(size(x)));
%
equi_ext.grho2r2      = resample_time_space_fbe(time_in,x_in,cur_diff_ext.grho2r2,temps,x);
equi_ext.r2i          = resample_time_space_fbe(time_in,x_in,cur_diff_ext.r2i,temps,x);
equi_ext.ri           = resample_time_space_fbe(time_in,x_in,cur_diff_ext.ri,temps,x);
equi_ext.C2           = resample_time_space_fbe(time_in,x_in,cur_diff_ext.C2,temps,x);
equi_ext.C3           = resample_time_space_fbe(time_in,x_in,cur_diff_ext.C3,temps,x);
equi_ext.grho         = resample_time_space_fbe(time_in,x_in,cur_diff_ext.grho,temps,x);
equi_ext.grho2        = resample_time_space_fbe(time_in,x_in,cur_diff_ext.grho2,temps,x);
%
equi_ext.jgs          = equi_ext.jli;
equi_ext.df2dpsi      = resample_time_space_fbe(time_in,x_in,cur_diff_ext.df2dpsi,temps,x);
equi_ext.dptotdpsi    = resample_time_space_fbe(time_in,x_in,cur_diff_ext.dptotdpsi,temps,x);
%
equi_ext.lif          = resample_time_fbe(time_in,cur_diff_ext.lif,temps);
equi_ext.betap        = resample_time_fbe(time_in,cur_diff_ext.betap,temps);
equi_ext.ipout        = resample_time_fbe(time_in,cur_diff_ext.ipout,temps);
equi_ext.difcurconv   = resample_time_fbe(time_in,cur_diff_ext.difcurconv,temps);
equi_ext.volume       = resample_time_space_fbe(time_in,x_in,cur_diff_ext.volume,temps,x);
equi_ext.dvdx         = pdederive(x,equi_ext.volume,0,2,2,1);
equi_ext.dvdx         = ((equi_ext.volume(:,end) ./ trapz(x,equi_ext.dvdx,2)) * ones(size(x))) .* equi_ext.dvdx;
equi_ext.phiplasma    = trapz(x,equi_ext.dvdx .* equi_ext.r2i .* (equi_ext.fdia - equi_ext.fdia(:,end) * ones(size(x))),2) ./ (2 * pi);
equi_ext.indice_inv   = resample_time_fbe(time_in,cur_diff_ext.indice_inv ./ size(cur_diff_ext.psi,2),temps) .* size(equi_ext.psi,2);
equi_ext.poynting     = resample_time_fbe(time_in,cur_diff_ext.poynting,temps);
equi_ext.wbp          = resample_time_fbe(time_in,cur_diff_ext.wbp,temps);
equi_ext.dwbpdt       = resample_time_fbe(time_in,cur_diff_ext.dwbpdt,temps);
equi_ext.rm           = resample_time_fbe(time_in,cur_diff_ext.rm,temps);
equi_ext.drmdt        = resample_time_fbe(time_in,cur_diff_ext.drmdt,temps);
equi_ext.dphidt       = resample_time_fbe(time_in,cur_diff_ext.dphidt(:,end),temps);
equi_ext.psid1        = resample_time_space_fbe(time_in,x_in,cur_diff_ext.psid1,temps,x);
equi_ext.psid2        = resample_time_space_fbe(time_in,x_in,cur_diff_ext.psid2,temps,x);
equi_ext.pohm         = resample_time_fbe(time_in,cur_diff_ext.pohm,temps);
equi_ext.q0           = resample_time_fbe(time_in,cur_diff_ext.q0,temps);
equi_ext.qmin         = resample_time_fbe(time_in,cur_diff_ext.qmin,temps);
equi_ext.q95          = resample_time_fbe(time_in,cur_diff_ext.q95,temps);

% effective trapped fraction
equi_ext.ftrap        = resample_time_space_fbe(time_in,x_in,cur_diff_ext.ftrap,temps,x);

% moments
equi_ext.kx           = resample_time_space_fbe(time_in,x_in,cur_diff_ext.kx,temps,x);
equi_ext.dx           = resample_time_space_fbe(time_in,x_in,cur_diff_ext.dx,temps,x);
equi_ext.d0           = equi_ext.Raxe(:,1) - equi_ext.Raxe(:,end);
equi_ext.K95          = resample_time_fbe(time_in,cur_diff_ext.K95,temps);
equi_ext.d95          = resample_time_fbe(time_in,cur_diff_ext.d95,temps);
equi_ext.piqj         = resample_time_fbe(time_in,cur_diff_ext.piqj,temps);


% LCFS
fbe_lcfs = cell([length(time_in),1]);
for k = 1:length(time_in)
    fbe_lcfs{k}.time = time_in(k);
    fbe_lcfs{k}.ip   = cur_diff_ext.ipout(k);
    fbe_lcfs{k}.R    = cur_diff_ext.geo.Rsepa{k};
    fbe_lcfs{k}.Z    = cur_diff_ext.geo.Zsepa{k};
end
z0dinput.cons     = cons;
z0dinput.geo      = geo;
z0dinput.option   = option;
z0dinput.exp0d    = exp0d;
z0dinput          = z0separatrix_fromfbetometis(z0dinput,fbe_lcfs);
equi_ext.geo      = z0dinput.geo;
equi_ext.cons     = z0dinput.cons;
rb0               = resample_time_fbe(time_in,cur_diff_ext.geo.b0 .* cur_diff_ext.geo.R,temps);
equi_ext.geo.b0   = rb0 ./ equi_ext.geo.R;
equi_ext.vp       = resample_time_fbe(time_in,cur_diff_ext.geo.vp,temps);
equi_ext.geo.vp   = equi_ext.vp;
equi_ext.sp       = resample_time_fbe(time_in,cur_diff_ext.geo.sp,temps);
equi_ext.geo.sp   = equi_ext.sp;
equi_ext.peri     = resample_time_fbe(time_in,cur_diff_ext.geo.peri,temps);
equi_ext.geo.peri = equi_ext.peri;
equi_ext.sext     = resample_time_fbe(time_in,cur_diff_ext.geo.sext,temps);
equi_ext.geo.sext = equi_ext.sext;
equi_ext.xpoint   = double(resample_time_fbe(time_in,cur_diff_ext.geo.xpoint,temps) > 0);
equi_ext.Rsepa    = z0dinput.exp0d.Rsepa;   
equi_ext.Zsepa    = z0dinput.exp0d.Zsepa;
















