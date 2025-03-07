% Make external equilibrium data structure from METIS data using FEEQS fast
% fixed equilibrium mode of FEEQS.M to recompute precise solutionof G-S
% equation. FEEQS.M must be installed separately on the same computer.
%
% syntax:
%
%  make_external_equi_with_feeqs(post,plotonoff);
%
% input:
%
%   post = METIS data structure
%
% optional inputs:
%
%   plotonoff   = if true, diplay graphs
%
% output:
%
%  External equilibrium for METIS
%
% J-F Artaud 2021
function [data,profiles_2d,FPSI,FBR,FBZ,FBPHI,plist_out] = make_external_equi_with_feeqs(post,plotonoff)

% initialisation outputs
profiles_2d = [];
FPSI = [];
FBR = [];
FBZ = [];
FBPHI = [];
plist_out = [];
data = [];

if nargin < 1
    error('syntax:  make_eqdsk_with_feeqs(post,[precision,plotonoff]);');
elseif isempty(post)
    error('post empty')
elseif ~isfield(post,'zerod') || ~isfield(post,'profil0d') || ~isfield(post,'z0dinput')
    error('post is not a METIS result data structure');
end

if (nargin < 2) || isempty(plotonoff)
    plotonoff = false;
end

% constant
mu0 = 4 * pi * 1e-7;

% make path for FEEQS
set_reactor_feeqs_path

% extract data
% time dependant data
time_slices               = post.profil0d.temps;
%
metis_data4fbe.time                 = time_slices;
metis_data4fbe.index                = round(interp1_imas(post.zerod.temps,1:length(post.zerod.temps),time_slices,'nearest','extrap'));
metis_data4fbe.ip                   = interp1_imas(post.zerod.temps,post.zerod.ip,time_slices,'nearest','extrap');
metis_data4fbe.vloop                = interp1_imas(post.zerod.temps,post.zerod.vloop,time_slices,'nearest','extrap');
metis_data4fbe.li                   = interp1_imas(post.zerod.temps,post.zerod.li,time_slices,'nearest','extrap');
metis_data4fbe.betap                = interp1_imas(post.zerod.temps,post.zerod.betaptot,time_slices,'nearest','extrap');
metis_data4fbe.vp                   = interp1_imas(post.zerod.temps,post.zerod.vp,time_slices,'nearest','extrap');
metis_data4fbe.psi                  = interp1_imas(post.profil0d.temps,post.profil0d.psi,time_slices,'nearest','extrap');
metis_data4fbe.q                    = interp1_imas(post.profil0d.temps,post.profil0d.qjli,time_slices,'nearest','extrap');
metis_data4fbe.rho                  = interp1_imas(post.profil0d.temps,post.profil0d.rmx,time_slices,'nearest','extrap');
metis_data4fbe.phi                  = interp1_imas(post.profil0d.temps,post.profil0d.phi,time_slices,'nearest','extrap');
metis_data4fbe.epar                 = interp1_imas(post.profil0d.temps,post.profil0d.epar,time_slices,'nearest','extrap');
metis_data4fbe.ptot                 = interp1_imas(post.profil0d.temps,post.profil0d.ptot,time_slices,'nearest','extrap');
metis_data4fbe.dptotdpsi            = interp1_imas(post.profil0d.temps,post.profil0d.dptotdpsi,time_slices,'nearest','extrap');
metis_data4fbe.fdia                 = interp1_imas(post.profil0d.temps,post.profil0d.fdia,time_slices,'nearest','extrap');
metis_data4fbe.df2dpsi              = interp1_imas(post.profil0d.temps,post.profil0d.df2dpsi,time_slices,'nearest','extrap');
metis_data4fbe.jphi                 = interp1_imas(post.profil0d.temps,post.profil0d.jli,time_slices,'nearest','extrap');
metis_data4fbe.rmx                  = interp1_imas(post.profil0d.temps,post.profil0d.rmx,time_slices,'nearest','extrap');
metis_data4fbe.C2                   = interp1_imas(post.profil0d.temps,post.profil0d.C2,time_slices,'nearest','extrap');
metis_data4fbe.grad_rho2or2         = interp1_imas(post.profil0d.temps,post.profil0d.grho2r2,time_slices,'nearest','extrap');
metis_data4fbe.Raxe                 = interp1_imas(post.profil0d.temps,post.profil0d.Raxe,time_slices,'nearest','extrap');
metis_data4fbe.kx                   = interp1_imas(post.profil0d.temps,post.profil0d.kx,time_slices,'nearest','extrap');
metis_data4fbe.dx                   = interp1_imas(post.profil0d.temps,post.profil0d.dx,time_slices,'nearest','extrap');
metis_data4fbe.ri                   = interp1_imas(post.profil0d.temps,post.profil0d.ri,time_slices,'nearest','extrap');
metis_data4fbe.r2i                  = interp1_imas(post.profil0d.temps,post.profil0d.r2i,time_slices,'nearest','extrap');
metis_data4fbe.peri                 = interp1_imas(post.zerod.temps,post.zerod.peri,time_slices,'nearest','extrap');
metis_data4fbe.psin                 = (metis_data4fbe.psi - metis_data4fbe.psi(:,1) * ones(1,size(metis_data4fbe.psi,2))) ./ ...
    ((metis_data4fbe.psi(:,end) - metis_data4fbe.psi(:,1)) * ones(1,size(metis_data4fbe.psi,2)));
metis_data4fbe.betap_cir            = 8/3 .* interp1_imas(post.zerod.temps,post.zerod.w,time_slices,'nearest','extrap')  ./  ...
    (mu0 .* metis_data4fbe.ip .^ 2 .*  metis_data4fbe.Raxe(:,end));
metis_data4fbe.betap_loc            = (2/3) * interp1_imas(post.zerod.temps,post.zerod.w,time_slices,'nearest','extrap') ./ ...
    metis_data4fbe.vp ./ (mu0 * metis_data4fbe.ip .^ 2 ./  metis_data4fbe.peri .^ 2 / 2);
%
metis_data4fbe.xpoint      = interp1_imas(post.zerod.temps,post.zerod.xpoint,time_slices,'nearest','extrap');


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
if plotonoff
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
    for k= 1:length(time_slices)
        metis_data4fbe.z0(k) = (max(metis_data4fbe.LCFS.Z(k,:)) + min(metis_data4fbe.LCFS.Z(k,:))) ./ 2;
        metis_data4fbe.LCFS.Z(k,:) = metis_data4fbe.LCFS.Z(k,:) -  metis_data4fbe.z0(k);
        metis_data4fbe.z0(k) = metis_data4fbe.z0(k) + post.z0dinput.geo.z0(k);
    end
else
    metis_data4fbe.LCFS.R = NaN * ones(length(time_slices),201);
    metis_data4fbe.LCFS.Z = NaN * ones(length(time_slices),201);
    for k= 1:length(time_slices)
        index        = find(post.z0dinput.cons.temps >= time_slices(k),1);
        geo_metis    = zerod_get1t(post.z0dinput.geo,index);
        [vps_next,void_sps,void_sexts,void_peris,geo_metis,void_xpoint,sepa_metis.Rsepa,sepa_metis.Zsepa] =  ...
            zgeo0(geo_metis,[],[],1);
        %sepa_metis.Zsepa = sepa_metis.Zsepa
        metis_data4fbe.LCFS.R(k,:) = sepa_metis.Rsepa(:);
        metis_data4fbe.LCFS.Z(k,:) = sepa_metis.Zsepa(:); 
        metis_data4fbe.z0(k) = geo_metis.z0;
    end
end

% initialise data for external equilibrium
% data fields initialisation
nbt                = length(time_slices);
nbx                = 21;
v0d                = NaN * ones(nbt,1);
v1d                = NaN * ones(nbt,nbx);
data.x             = v1d;
data.time          = time_slices;
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

% void figure for parasitic graph coming from FEEQS
hf_void = figure;

% loop on profiles time
for k=1:length(time_slices)
    fprintf('Computation of precise equilibrium with FEEQS fast mode @ t = %g (s); index = %d/%d (starting with p=7):\n',time_slices(k),metis_data4fbe.index(k),length(time_slices));
    equiok = false;
    plist = [];
    % plasma current from Metis
    plist.Ip = metis_data4fbe.ip(k);   % total plasma current
    
    % LCFS
    LCFS = cat(2,metis_data4fbe.LCFS.R(k,:)',metis_data4fbe.LCFS.Z(k,:)');
    
    % value of psi on LCFS
    psi_ref = metis_data4fbe.psi(k,end);
    % Edge flux is unknown
    plist.psi_bd     = psi_ref;
    %  magnetic reigidity
    plist.RB0        = metis_data4fbe.fdia(k,end);
    
    % turn on post processing and tune options
    %    plist.iter_max         = maximal number of iterations (on Newton convergence loop, 100)
    %    plist.NewtonStop       = Newton loop stoping criterium (1e-10)
    plist.verbose          = double(plotonoff);
    plist.plotonoff        = plotonoff;
    %    plist.p                = hpFEM order  (7 is the best trend off for JT-60SA ; 3 provide smoother 2D grid for inverse grid and are faster)
    plist.twoD_data        = false; %set on or off computation of 2D data Psi, BR, BZ and Bphi needed for graphical display and inverse grid computation (false, quite slow, can be computed later)
    plist.other            = true; %set on or off false computation of other quantities like q, betap,li and flux averaged quantities (gm1..gm8) (false, can be computed later)
    plist.nbp_plot         = 11;   %number of radial points in plot and number of points in 1D output profiles (101, don't change by mutch computation time)
    plist.nbp_grid2d       = 21; %number of points in each direction taken in sampling of hpFEM for graphical display and inverse grid computation (101, computation time go like p^2 *  nbp_grid2d^2)
    %    plist.match_li_beta    = when profile current are prescribed by alpha, beta and gamma; assume beta been beta_p and gamma been l_i(3) and perform convergence to match beta_p and l_i(3) (false)
    %    plist.tol_li_beta      = convergence criterium for the convergence to beta_p and l_i(3) (1e-3)
    plist.limiter_name     = post.z0dinput.option.first_wall; %optional polodidal limiter description (empty string)
    plist.space_filter = 3;
    %    poly_1d_order          = polynomial order use in flux surface averaged quatities computation (30)
    plist.newton_convergence_mode = 1;            
    plist.match_li_beta    = false;

    % high
    plist.p = 7;
    try
        [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS, metis_data4fbe.ip(k),...
            1/2*metis_data4fbe.df2dpsi(k,:),metis_data4fbe.dptotdpsi(k,:),metis_data4fbe.psi(k,:),plist);
        if isempty(psi_ax) || ~isfinite(psi_ax) || (lambda < sqrt(eps))
            error('just to trap it');            
        elseif (plist_out.residual  <= plist_out.NewtonStop) &&  ...
                   (plist_out.fact_ip_inte > 0.5) && ...
                   (plist_out.fact_ip_inte < 1.5)
            equiok = true;
        else
            error('just to trap it');
        end
    catch
        equiok = false;
        fprintf('retry with p = 5\n');
    end
    if ~equiok
        plist.p = 5;
        try
            [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS, metis_data4fbe.ip(k),...
                1/2*metis_data4fbe.df2dpsi(k,:),metis_data4fbe.dptotdpsi(k,:),metis_data4fbe.psi(k,:),plist);          
            if isempty(psi_ax) || ~isfinite(psi_ax) || (lambda < sqrt(eps))
                error('just to trap it');
            elseif (plist_out.residual  <= plist_out.NewtonStop) &&  ...
                   (plist_out.fact_ip_inte > 0.5) && ...
                   (plist_out.fact_ip_inte < 1.5)
                equiok = true;
            else
                error('just to trap it');
            end
        catch
            equiok = false;
            fprintf('retry with p = 3\n');
        end
    end
    if ~equiok
        plist.p = 3;
        try
            [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS, metis_data4fbe.ip(k),...
                1/2*metis_data4fbe.df2dpsi(k,:),metis_data4fbe.dptotdpsi(k,:),metis_data4fbe.psi(k,:),plist);
            if isempty(psi_ax) || ~isfinite(psi_ax) || (lambda < sqrt(eps))
                error('just to trap it');
            elseif (plist_out.residual  <= plist_out.NewtonStop) &&  ...
                   (plist_out.fact_ip_inte > 0.5) && ...
                   (plist_out.fact_ip_inte < 1.5)
                equiok = true;
            else
                error('just to trap it');
            end
        catch
            fprintf('Unable to converge with profiles, try with moments\n');
            try
                plist.match_li_beta    = true;
                [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS, metis_data4fbe.ip(k),...
                    metis_data4fbe.li(k),metis_data4fbe.betap_loc(k),1,plist);
                if ~isempty(psi_ax)  && isfinite(psi_ax)
                    if (plist_out.residual  <= plist_out.NewtonStop) &&  ...
                            (plist_out.fact_ip_inte > 0.5) && ...
                            (plist_out.fact_ip_inte < 1.5)
                        equiok = true;
                    else
                        disp('Unable to converge also with moments');
                    end
                else
                    disp('Unable to converge also with moments');
                end
            catch
                disp('Unable to converge also with moments');
            end
        end
    end
    if plotonoff
        drawnow
    end
       
    if equiok
        disp('equilibrium found, starting post processing');
        if isfield(plist_out,'hf')
            plist.hf = plist_out.hf;
        end
        if isfield(plist_out,'hf_q')
            plist.hf_q = plist_out.hf_q;
            figure(plist.hf_q)
            hold on
            plot(metis_data4fbe.psin(k,:),metis_data4fbe.q(k,:),'r');
            legend('FEEQS','METIS');
            drawnow
        end
        % detail comparison with METIS
        if plotonoff
            hfm = findobj(0,'type','figure','tag','metis_feeqs_cp');
            if isempty(hfm)
                hfm = figure('tag','metis_feeqs_cp');
            else
                figure(hfm);
                clf;
            end
            subplot(2,4,1)
            plot(plist_out.psin,plist_out.q,'b',metis_data4fbe.psin(k,:),metis_data4fbe.q(k,:),'r');
            %xlabel('Psi normalized');
            ylabel('q');
            legend('FEEQS','METIS');
            subplot(2,4,2)
            plot(plist_out.psin,plist_out.jphi,'b',metis_data4fbe.psin(k,:),metis_data4fbe.jphi(k,:),'r');
            %xlabel('Psi normalized');
            ylabel('<J_{phi}> (A/m^2)');
            %legend('FEEQS','METIS');
            subplot(2,4,3)
            plot(plist_out.psin,plist_out.ptot,'b',metis_data4fbe.psin(k,:),metis_data4fbe.ptot(k,:),'r');
            %xlabel('Psi normalized');
            ylabel('P_{tot} (Pa)');
            %legend('FEEQS','METIS');
            subplot(2,4,4)
            plot(plist_out.psin,plist_out.fdia,'b',metis_data4fbe.psin(k,:),metis_data4fbe.fdia(k,:),'r');
            %xlabel('Psi normalized');
            ylabel('F (T.m)');
            %legend('FEEQS','METIS');
            subplot(2,4,5)
            plot(plist_out.psin,plist_out.phi,'b',metis_data4fbe.psin(k,:),metis_data4fbe.phi(k,:),'r');
            xlabel('Psi normalized');
            ylabel('Phi (Wb)');
            %legend('FEEQS','METIS');
            subplot(2,4,6)
            plot(plist_out.psin,plist_out.grad_rho2or2,'b',metis_data4fbe.psin(k,:),metis_data4fbe.grad_rho2or2(k,:),'r');
            xlabel('Psi normalized');
            ylabel('<grad(rho)^2/R^2>');
            %legend('FEEQS','METIS');
            subplot(2,4,7)
            plot(plist_out.psin,plist_out.profileHandle.Sffp(plist_out.psin,1) .* plist_out.lambda,'b',metis_data4fbe.psin(k,:),metis_data4fbe.df2dpsi(k,:)/2,'r', ...
                metis_data4fbe.psin(k,:),metis_data4fbe.df2dpsi(k,:)/2,'r.');
            xlabel('Psi normalized');
            ylabel('F*dF/dPsi');
            subplot(2,4,8)
            plot(plist_out.psin,plist_out.profileHandle.Spp(plist_out.psin,1) .* plist_out.lambda,'b',metis_data4fbe.psin(k,:),metis_data4fbe.dptotdpsi(k,:),'r', ...
                metis_data4fbe.psin(k,:),metis_data4fbe.dptotdpsi(k,:),'r.');
            xlabel('Psi normalized');
            ylabel('dP/dPsi');
            drawnow
            fprintf('Ip FEEQS     = %g MA & Ip METIS   = %g MA\n',plist_out.ip_feeqs/1e6,metis_data4fbe.ip(k)/1e6);
            fprintf('li_3 FEEQS   = %g & li_3 METIS    = %g\n',plist_out.li_3,metis_data4fbe.li);
            fprintf('beta_p FEEQS = %g & beta_p METIS  = %g\n',plist_out.beta_p,metis_data4fbe.betap_loc);
            
        end
        
        % 2D grid
        [profiles_2d,FPSI,FBR,FBZ,FBPHI,hout,plist_out] = compute_inverse_grid_interp_hpFEM(plist_out,[],linspace(0,1,21),1/pi/21);
        plist_out.hout = hout;
        
        % fill data for external equilibrium
%          Rgeo_loc = (max(profiles_2d.r,[],2) +  min(profiles_2d.r,[],2)) ./ 2;
%          Rgeo_loc = Rgeo_loc(end);
%          Zgeo_loc = (max(profiles_2d.z,[],2) +  min(profiles_2d.z,[],2)) ./ 2;
%          Zgeo_loc = Zgeo_loc(end);
%          fprintf('Axis (GS and geo):');
%          disp([plist_out.Raxis,plist_out.Zaxis,Rgeo_loc,Zgeo_loc]);
        % data extraction
        Raxe   = (max(profiles_2d.r,[],2) +  min(profiles_2d.r,[],2)) ./ 2;
        Zaxe   = (max(profiles_2d.z,[],2) +  min(profiles_2d.z,[],2)) ./ 2;
        aminor = (max(profiles_2d.r,[],2) -  min(profiles_2d.r,[],2)) ./ 2;
        x      = aminor ./ max(aminor);
        kx     = (max(profiles_2d.z,[],2) -  min(profiles_2d.z,[],2)) ./ ...
                 (max(profiles_2d.r,[],2) -  min(profiles_2d.r,[],2));
        kx(1) = 2 .* kx(2) - kx(3);
        for l=1:size(profiles_2d.z,1)
            [~,indzmax] = max(profiles_2d.z(l,:),[],2);
            [~,indzmin] = min(profiles_2d.z(l,:),[],2);
            dx(l)          = abs((Raxe(l) - profiles_2d.r(l,indzmax)) ./ aminor(l) + ...
                                (Raxe(l) - profiles_2d.r(l,indzmin)) ./ aminor(l)) / 2;
        end
        dx(1) = 2 .* dx(2) - dx(3);
        dx    = abs(polyval(polyfit(x(:),dx(:),5),x(:)));
        %
        psid1             = pdederive(x(:)',plist_out.psi1D(:)',0,2,2,1);
        % dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
        psid2    = pdederive(x(:)',plist_out.psi1D(:)',1,0,2,2);

        
        
        % create external data structure
        data.x(k,:)         = x(:)';
        data.psi(k,:)       = plist_out.psi1D(:)';
        data.phi(k,:)       = abs(plist_out.phi);
        data.psid1(k,:)     = psid1;
        data.psid2(k,:)     = psid2;
        data.qjli(k,:)      = abs(plist_out.q(:)');
        data.jli(k,:)       = plist_out.jphi(:)';
        data.bpol(k,:)      = plist_out.bpol(:)';
        data.rmx(k,:)       = plist_out.rho(:)';
        data.rm(k)          = max(data.rmx(k,:));
        data.grho2r2(k,:)   = plist_out.grad_rho2or2(:)';
        data.r2i(k,:)       = plist_out.ri2(:)';
        data.ri(k,:)        = plist_out.ri(:)';
        data.grho(k,:)      = abs(plist_out.grad_rho(:));
        data.grho2(k,:)     = abs(plist_out.grad_rho2(:));
        data.spr_tor(k,:)   = abs(plist_out.dSdrho(:)');
        data.spr_tor(k,1)   = 2.* data.spr_tor(k,2) - data.spr_tor(k,3);
        data.vpr_tor(k,:)       = abs(plist_out.dVdrho(:)');
        data.vpr_tor(k,1)       = 2.* data.vpr_tor(k,2) - data.vpr_tor(k,3);
        data.volume(k,:)        = cumtrapz(plist_out.rho(:)',abs(plist_out.dVdrho(:)'));
        data.volume(k,:)        = data.volume(k,:) ./ data.volume(k,end) .*  plist_out.volume;        
        data.surface_pol(k,:)   = cumtrapz(plist_out.rho(:)',abs(plist_out.dSdrho(:)'));
        data.surface_pol(k,:)   = data.surface_pol(k,:) ./ data.surface_pol(k,end) .*  plist_out.section;        
        data.surface_flux(k,:)  = abs(plist_out.dVdrho(:)' .* abs(plist_out.grad_rho(:)'));
        data.dphidx(k,:)        = pdederive(x,data.phi(k,:),0,2,2,1);
        data.vpr(k,:)           = pdederive(x,data.volume(k,:),0,2,2,1);
        data.spr(k,:)           = pdederive(x,data.surface_pol(k,:),0,2,2,1);
        data.C2(k,:)            = abs(plist_out.dVdrho(:)' .* data.grho2r2(k,:));
        data.C3(k,:)            = abs(plist_out.dVdrho(:)' .* data.r2i(k,:));
        data.jeff(k,:)          = plist_out.j_par(:)';
        data.jres(k,:)          = NaN * data.jeff(k,:); % must be computed with METIS data
        data.Raxe(k,:)          = Raxe(:)';
        data.Zaxe(k,:)          = Zaxe(:)';
        data.epsi(k,:)          = plist_out.epsi(:)';
        data.a(k,:)             = aminor(:)';
        data.difcurconv(k)      = plist_out.nb_iter;
        data.df2dpsi(k,:)       = 2 .* plist_out.profileHandle.Sffp(plist_out.psin,1)';
        data.dptotdpsi(k,:)     = plist_out.profileHandle.Spp(plist_out.psin,1)';
        data.ptot(k,:)          = plist_out.ptot(:)' +  metis_data4fbe.ptot(1,end);
        data.fdia(k,:)          = plist_out.fdia(:)';
        data.wbp(k)             = plist_out.mu0 .* plist_out.ip_feeqs .^ 2 .* Raxe(end) ./ 4 .* plist_out.li_3;
        data.kx(k,:)            = kx(:)';
        data.dx(k,:)            = dx(:)';
        data.ipout(k)           = plist_out.ip_feeqs;
        data.lif(k)             = plist_out.li_3;
        data.betap(k)           = plist_out.beta_p; 
        data.qmin(k)            = min(data.qjli(k,:),[],2);
        data.q0(k)              = data.qjli(k,1);
        data.psin(k,:)          = (data.psi(k,:) - data.psi(k,1)) ./ (data.psi(k,end) - data.psi(k,1));
        dd95                    = abs(data.psin(k,:) - 0.95);
        mask95                  = double(dd95 == min(dd95,[],2));
        data.q95(k)             = sum(data.qjli(k,:) .* mask95,2) ./ max(1, sum(mask95,2));
        data.d95(k)             = sum(data.dx(k,:) .* mask95,2) ./ max(1, sum(mask95,2));
        data.K95(k)             = sum(data.kx(k,:) .* mask95,2) ./ max(1, sum(mask95,2));
        data.piqj(k)            = max(0.1,min(10,data.jli(k,1) ./ (sum(diff(data.surface_pol(k,:),1,2) .*  ...
                                  (data.jli(k,1:end-1) + data.jli(k,2:end)) ./ 2,2) ./ data.surface_pol(k,end))));
        
        % LCFS
        data.geo.a(k)         = aminor(end);
        data.geo.R(k)         = Raxe(end);
        data.geo.z0(k)        = metis_data4fbe.z0(k);
        data.geo.K(k)         = kx(end);
        data.geo.d(k)         = dx(end);
        data.geo.b0(k)        = plist_out.RB0 ./ Raxe(end);
        data.geo.vp(k)        = data.volume(k,end);
        data.geo.sp(k)        = data.surface_pol(k,end);
        data.geo.sext(k)      = data.surface_flux(k,end);
        % xpoint is inherited from METIS
        data.geo.xpoint(k)    = metis_data4fbe.xpoint(k);
        data.geo.Rsepa{k}     = profiles_2d.r(end,:);
        data.geo.Zsepa{k}    =  profiles_2d.z(end,:);
        data.geo.peri(k)      = plist_out.len_LCFS;
        
        % inversion radius (no prescription on trigger value
        if min(data.qjli(k,:)) < 0.9
            mask = (data.qjli(k,:) <= 1);
            ggr  = (1:size(data.qjli,2));
            data.indice_inv(k) = min(length(data.qjli(k,:)),max(mask .* ggr,[],2));
        else
            data.indice_inv(k) = 0;
        end
        
        % effective trapped fraction
        data.ftrap(k,:) = ftrap(data.x(k,:),data.epsi(k,:),data.psi(k,:),data.bpol(k,:),data.fdia(k,:), ...
                                data.Raxe(k,:),data.r2i(k,:),data.ri(k,:),data.grho2r2(k,:),data.grho2(k,:),data.a(k,:));

        
    end
    
end

% close void figure
if ishandle(hf_void)
    delete(hf_void);
end

% remove invalid time
if any(~isfinite(data.ipout))
    indok = isfinite(data.ipout);
    noms = fieldnames(data);
    for k=1:length(noms)
        if isstruct(data.(noms{k}))
            noms2 = fieldnames(data.(noms{k}));
            for l=1:length(noms2)
                 data.(noms{k}).(noms2{l}) = data.(noms{k}).(noms2{l})(indok,:);
            end
        else
            data.(noms{k}) = data.(noms{k})(indok,:);
        end
    end
end
if length(data.ipout) >= 3
    
    % compute time derivative
    data.dphidt       = z0dxdt_c(data.phi,data.time);
    data.dpsidt       = z0dxdt_c(data.psi,data.time);
    % approximation: neglecting ohmic power due to diamagnetic change
    ephi    = - data.ri .* data.dpsidt;
    ej      =   data.jli.* ephi;
    pohm_equi1d = sum(diff(data.volume,1,2) .* (ej(:,1:end-1) + ej(:,2:end)) ./ 2,2);
    ip_equi1d   = sum(diff(data.surface_pol,1,2) .* (data.jli(:,1:end-1) + data.jli(:,2:end)) ./ 2,2);
    vres_equi1d = pohm_equi1d ./ ip_equi1d;
    data.epar = ephi;
    data.ej   = ej;
    data.dwbpdt          = z0dxdt_c(data.wbp,data.time);
    data.pohm            = pohm_equi1d;
    data.drmdt           = z0dxdt_c(data.rm,data.time);
    % Poynting flux
    % ref :Ohmic ???ux consumption during initial operation of the NSTX spherical torus
    % J.E. Menard et al, NF 2001 41 1197
    % ici on prend en compte que la composante est utile pour la consommation de flux poloidal.
    dEBpoldt_time        = z0dxdt_c(data.bpol .^ 2 ./ 2 ./ plist_out.mu0,data.time);
    for k =1:length(data.time)
        data.poynting(k) = trapz(data.rmx(k,:),data.vpr_tor(k,:) .* (dEBpoldt_time(k,:) + data.ej(k,:)),2);
    end
    % set on external equilibrium if nbt > 3
    % set same stucture for current diffusion and equilibrium
    %%%setappdata(0,'CURDIFF_EXP',data);
    setappdata(0,'EQUILIBRIUM_EXP',data);
    
end

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



