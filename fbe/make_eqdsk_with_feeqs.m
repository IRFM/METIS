% Make eqdsk file from METIS data using FEEQS fast mode
% FEEQS must be installed separately on the same computer
% syntax:
%
%  make_eqdsk_with_feeqs(post,time_out,ngrid,filename_in,precision,plotonoff);
%
% input:
%
%   post = METIS data structure
%
% optional inputs:
%
%   time_out     = vector of time (s), if empty generate one eqdsk file for
%                  each time slice in the METIS simulation.
%
%   ngrid        = number of point in outputed grid, if vector of two
%                  element = [nR,nZ], if empty selected from plasma size 
%                  and requested, precision
%
%   filename_in = optional root name for output file 
%                (must also contains the path)
%
%   precision   = precision and size of RZ grid (low, medium, high)
%
%   plotonoff   = if true, diplay graphs
%
% output:
%
%  one geqdsk file for each element of time_out or for each time slice in
%  the METIS simulaton
%
% J-F Artaud 2021
function [profiles_2d,FPSI,FBR,FBZ,FBPHI,plist_out] = make_eqdsk_with_feeqs(post,time_out,ngrid,filename_in,precision,plotonoff,alternative_extrapolation)

% initialisation outputs
profiles_2d = [];
FPSI = [];
FBR = [];
FBZ = [];
FBPHI = [];
plist_out = [];

if nargin < 1
    error('syntax:  make_eqdsk_with_feeqs(post,[time_out,ngird,filename_in,precision,plotonoff]);');
elseif isempty(post)
    error('post empty')
elseif ~isfield(post,'zerod') || ~isfield(post,'profil0d') || ~isfield(post,'z0dinput')
    error('post is not a METIS result data structure');
end
if (nargin < 2)
    time_out = [];
end
if (nargin < 3)
    ngrid = [];
end
if (nargin < 4) || isempty(filename_in)
    filename_in = sprintf('METIS_FEEQS_FAST_%s_%d',post.z0dinput.machine,post.z0dinput.shot);
end
if (nargin < 5) || isempty(precision)
    precision = 'high';   % low, medium or high
 end

if (nargin < 6) || isempty(plotonoff)
    if length(time_out) == 1
         plotonoff = true;
    else
         plotonoff = false;
    end
end
if (nargin < 7) || isempty(alternative_extrapolation)
    alternative_extrapolation = 'Interpolation';
end

% constant
mu0 = 4 * pi * 1e-7;

% make path for FEEQS
set_reactor_feeqs_path

% extract data
% time dependant data
if ~isempty(time_out)
    ind_time                  = interp1(post.profil0d.temps,1:length(post.profil0d.temps),time_out,'nearest',NaN);
    ind_time(~isfinite(ind_time)) = [];
    time_slices               = post.profil0d.temps(unique(ind_time));
else
    time_slices               = post.profil0d.temps;
end
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
else
    metis_data4fbe.LCFS.R = NaN * ones(length(time_slices),201);
    metis_data4fbe.LCFS.Z = NaN * ones(length(time_slices),201);
    for k= 1:length(time_slices)
        index        = find(post.z0dinput.cons.temps >= time_slices(k),1);
        geo_metis    = zerod_get1t(post.z0dinput.geo,index);
        [vps_next,void_sps,void_sexts,void_peris,geo_metis,void_xpoint,sepa_metis.Rsepa,sepa_metis.Zsepa] =  ...
            zgeo0(geo_metis,[],[],1);
        sepa_metis.Zsepa = sepa_metis.Zsepa  + geo_metis.z0;
        metis_data4fbe.LCFS.R(k,:) = sepa_metis.Rsepa(:);
        metis_data4fbe.LCFS.Z(k,:) = sepa_metis.Zsepa(:);
    end
end


% loop on profiles time
for k=1:length(time_slices)
      fprintf('Computation of precise equilibrium with FEEQS fast mode @ t = %g (s); index = %d:\n',time_slices(k),metis_data4fbe.index(k));
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
      %    plist.nbp_plot         = number of radial points in plot and number of points in 1D output profiles (101, don't change by mutch computation time)
      %    plist.nbp_grid2d       = number of points in each direction taken in sampling of hpFEM for graphical display and inverse grid computation (101, computation time go like p^2 *  nbp_grid2d^2)
      %    plist.match_li_beta    = when profile current are prescribed by alpha, beta and gamma; assume beta been beta_p and gamma been l_i(3) and perform convergence to match beta_p and l_i(3) (false)
      %    plist.tol_li_beta      = convergence criterium for the convergence to beta_p and l_i(3) (1e-3)
      plist.limiter_name     = post.z0dinput.option.first_wall; %optional polodidal limiter description (empty string)
      plist.space_filter = 7;
      %    poly_1d_order          = polynomial order use in flux surface averaged quatities computation (30)
      switch precision
          case 'high'
              plist.p = 7;
              plist.nbp_grid2d = 101;                          
          case 'medium'
              plist.p = 5;
              plist.nbp_grid2d = 51;            
          otherwise
              plist.p = 3;
              plist.nbp_grid2d = 21;
      end
      plist.newton_convergence_mode = 1;    
      plist.match_li_beta           = false;
   
      try
          [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS, metis_data4fbe.ip(k),...
              1/2*metis_data4fbe.df2dpsi(k,:),metis_data4fbe.dptotdpsi(k,:),metis_data4fbe.psi(k,:),plist);
          if isempty(psi_ax) || ~isfinite(psi_ax)
              error('just to trap it');
          elseif (plist_out.residual  <= plist_out.NewtonStop) &&  ...
                  (plist_out.fact_ip_inte > 0.5) && ...
                  (plist_out.fact_ip_inte < 1.5)
              equiok = true;
          else
                error('just to trap it');
          end
      catch
          fprintf('Unable to converge with profiles, try with moments; ');
          try
              plist.match_li_beta    = true;
              [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS, metis_data4fbe.ip(k),...
                  metis_data4fbe.li(k),metis_data4fbe.betap_loc(k),1,plist);
              if ~isempty(psi_ax) && isfinite(psi_ax)
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
              fprintf('li_3 FEEQS   = %g & li_3 METIS    = %g\n',plist_out.li_3,metis_data4fbe.li(k));
              fprintf('beta_p FEEQS = %g & beta_p METIS  = %g\n',plist_out.beta_p,metis_data4fbe.betap_loc(k));   
                            
          end
          
          [profiles_2d,FPSI,FBR,FBZ,FBPHI,hout,plist_out] = compute_inverse_grid_interp_hpFEM(plist_out,[],'',1/pi/plist.nbp_grid2d);
          if plotonoff
                % plot equilibrium moment for comparison
                hfm2 = findobj(0,'type','figure','tag','metis_feeqs_cp2');
                if isempty(hfm2)
                    hfm2 = figure('tag','metis_feeqs_cp2');
                else
                    figure(hfm2);
                    clf;
                end

                psin_feeqs = profiles_2d.psi(:,1);
                psin_feeqs  = (psin_feeqs - psin_feeqs(1)) ./ (psin_feeqs(end) - psin_feeqs(1));
                Raxe_feeqs  = (max(profiles_2d.r,[],2) + min(profiles_2d.r,[],2)) ./ 2;
                kx_feeqs    = (max(profiles_2d.z,[],2) - min(profiles_2d.z,[],2)) ./ (max(profiles_2d.r,[],2) - min(profiles_2d.r,[],2));
                kx_feeqs(1) = kx_feeqs(2);
                dzmax       = abs(profiles_2d.z - (max(profiles_2d.z,[],2) * ones(1,size(profiles_2d.z,2)))) ./ (max(profiles_2d.z,[],2) - min(profiles_2d.z,[],2));
                dzmax       = exp(- 30 .* dzmax);
                %mask_max    = double(profiles_2d.z == (max(profiles_2d.z,[],2) * ones(1,size(profiles_2d.z,2))));
                rzmax       = sum(profiles_2d.r .* dzmax,2) ./ max(1,sum(dzmax,2));
                feeqs_dup   = 2 .* (Raxe_feeqs - rzmax) ./ (max(profiles_2d.r,[],2) - min(profiles_2d.r,[],2));
                dzmin       = abs(profiles_2d.z - (min(profiles_2d.z,[],2) * ones(1,size(profiles_2d.z,2)))) ./ (max(profiles_2d.z,[],2) - min(profiles_2d.z,[],2));
                dzmin       = exp(- 30 .* dzmin);
                %mask_min    = double(profiles_2d.z == (min(profiles_2d.z,[],2) * ones(1,size(profiles_2d.z,2))));
                rzmin       = sum(profiles_2d.r .* dzmin,2) ./ max(1,sum(dzmin,2));
                feeqs_dlow  = 2 .* (Raxe_feeqs - rzmin) ./ (max(profiles_2d.r,[],2) - min(profiles_2d.r,[],2));
                dx_feeqs    = (feeqs_dup + feeqs_dlow) / 2;
                subplot(2,2,2)
                plot(psin_feeqs,Raxe_feeqs - Raxe_feeqs(end),metis_data4fbe.psin(k,:),metis_data4fbe.Raxe(k,:) - metis_data4fbe.Raxe(k,end))
                xlabel('Psi_n')
                ylabel('delta (m)')
                subplot(2,2,1)
                plot(psin_feeqs,Raxe_feeqs,metis_data4fbe.psin(k,:),metis_data4fbe.Raxe(k,:))
                xlabel('Psi_n')
                ylabel('Raxe (m)')
                legend('FEEQS','METIS')
                subplot(2,2,3)
                plot(psin_feeqs,kx_feeqs,metis_data4fbe.psin(k,:),metis_data4fbe.kx(k,:))
                xlabel('Psi_n')
                ylabel('kx')
                subplot(2,2,4)
                plot(psin_feeqs,dx_feeqs,metis_data4fbe.psin(k,:),metis_data4fbe.dx(k,:))
                xlabel('Psi_n')
                ylabel('dx')
         end


          plist.hout = hout;
          % EQDSK file
          if isempty(ngrid)
            nw = 2 * plist.nbp_grid2d +1;
            nh = ceil(nw * sum(jj_LCFS.lcfs_kappa) /2);
          else
            nw = ngrid(1);
            nh = ngrid(end);
          end
          %
          rmaxis = plist_out.Raxis;
          zmaxis = plist_out.Zaxis;
          simag  = plist_out.psi_ax;
          sibry  = plist_out.psi_bd;
          rcentr = jj_LCFS.R0;
          bcentr = plist_out.RB0 ./ rcentr;
          current = plist_out.ip_feeqs;
          %
          fpol = interp1(plist_out.psin,plist_out.fdia,linspace(0,1,nw));
          pres = interp1(plist_out.psin,plist_out.ptot,linspace(0,1,nw));
          ffprim = plist_out.profileHandle.Sffp(linspace(0,1,nw),1) .* plist_out.lambda;
          pprime = plist_out.profileHandle.Spp(linspace(0,1,nw),1) .* plist_out.lambda;
          qpsi = interp1(plist_out.psin,plist_out.q,linspace(0,1,nw));
          %
          if isempty(plist_out.limiter_contour_RZ)
              rin  = jj_LCFS.lcfs_C(1) - 1.25 * jj_LCFS.lcfs_a;
              rout = jj_LCFS.lcfs_C(1) + 1.25 * jj_LCFS.lcfs_a;
              zup  = jj_LCFS.lcfs_C(2) + 1.25 * jj_LCFS.lcfs_a * jj_LCFS.lcfs_kappa(1);
              zlow = jj_LCFS.lcfs_C(2) - 1.25 * jj_LCFS.lcfs_a * jj_LCFS.lcfs_kappa(2);
              rlim = [rin,rout,rout,rin,rin];
              zlim = [zup,zup,zlow,zlow,zup];
          elseif size(plist_out.limiter_contour_RZ,2) == 2
             rlim = plist_out.limiter_contour_RZ(:,1);
             zlim = plist_out.limiter_contour_RZ(:,2);
         else    
            rlim = plist_out.limiter_contour_RZ(1,:);
            zlim = plist_out.limiter_contour_RZ(2,:);
          end 
          if (rlim(1) ~= rlim(end))||(zlim(1)~=zlim(end))
              rlim(end+1) = rlim(1);
              zlim(end+1) = zlim(1);
          end
          limitr = length(rlim);
          lim    = NaN * ones(1,2*limitr);         
          lim(1:2:2*limitr) = rlim;
          lim(2:2:2*limitr) = zlim;
          rdim = max(rlim)-min(rlim)+0.05;
          zdim = max(zlim)-min(zlim)+0.05;
          rleft = min(rlim)-0.025;
          zmid = (max(zlim) + min(zlim)) / 2;
          
          [rgrid, zgrid] = meshgrid(linspace(rleft,rleft+rdim,nw),linspace(zmid-zdim/2,zmid+zdim/2,nh));
          psirz = FPSI(rgrid, zgrid);
          
          switch alternative_extrapolation
              case 'G-S polynomial'
                  R_LCFS = metis_data4fbe.LCFS.R(k,:);
                  Z_LCFS = metis_data4fbe.LCFS.Z(k,:);
                  PSI_LCFS = metis_data4fbe.psi(k,end) * ones(size(R_LCFS));
                  BR_LCFS = FBR(R_LCFS,Z_LCFS);
                  BZ_LCFS = FBZ(R_LCFS,Z_LCFS);
                  
                  %               figure(21)
                  %               lcl = linspace(metis_data4fbe.psi(k,1), 1.3 .* metis_data4fbe.psi(k,end) - 0.3 .* metis_data4fbe.psi(k,1),131);
                  %               contour(rgrid,zgrid,psirz,lcl,'color','r');
                  %               hold on
                  
                  [FPSI_ext,FBR_ext,FBZ_ext] = extrapolate_from_LCFS(R_LCFS,Z_LCFS,PSI_LCFS,BR_LCFS,BZ_LCFS);
                  mask_out = ~zinout(R_LCFS,Z_LCFS,rgrid,zgrid);
                  psirz(mask_out) = FPSI_ext(rgrid(mask_out),zgrid(mask_out));
                  
                  %               contour(rgrid,zgrid,psirz,lcl,'color','b');
                  %               plot(R_LCFS,Z_LCFS,'k',R_LCFS,Z_LCFS,'.k');
                  %               keyboard
          end
          
          rbbbs = metis_data4fbe.LCFS.R(k,:);
          zbbbs = metis_data4fbe.LCFS.Z(k,:);
          if (rbbbs(1) ~= rbbbs(end)) ||(zbbbs(1) ~= zbbbs(end))
              rbbbs(end+1) = rbbbs(1);
              zbbbs(end+1) = zbbbs(1);
          end
          nbbbs = length(rbbbs);
          bbbs  = NaN * ones(1,2*nbbbs);
          bbbs(1:2:2*nbbbs) = rbbbs;
          bbbs(2:2:2*nbbbs) = zbbbs;
          xdum = 0;
          
          if plotonoff

              figure
              subplot(2,2,1)
              contour(rgrid,zgrid,psirz,50)
              colorbar
              hold on
              plot(rbbbs,zbbbs,'r','LineWidth',2)
              plot(rlim,zlim,'k')
              xlabel('R [m]','FontSize',12)
              ylabel('Z [m]','FontSize',12)
              axis equal
              axis tight
              hold off             
              
              
              subplot(2,2,2)
              plot(qpsi)
              ylabel('qpsi')
              subplot(2,2,3)
              plot(fpol)
              ylabel('fpol')
              subplot(2,2,4)
              plot(pres)
              ylabel('pres')
              
              drawnow

          end
          
          % Handle NaN

          %%
          time_eqdsk = time_slices(k);
          shot       = post.z0dinput.shot;
          comment    = 'FEEQS.M from METIS';
          filename   = sprintf('%s_%g_%dx%d.geqdsk',filename_in,time_eqdsk,nw,nh);
          fid        = fopen(filename,'w');
          fprintf(fid,'%8.8d %16.9e %21.21s %4.4d %4.4d %4.4d\n',shot,time_eqdsk,comment(1:min(21,length(comment))),1,nw,nh);
          fclose(fid);
          %dlmwrite(filename,[shot time_eqdsk],'delimiter',' ')
          %dlmwrite(filename,[1 nw nh],'delimiter',' ','coffset',51,'roffset',0,'-append')
          dlmwrite(filename,[rdim zdim rcentr rleft zmid],'-append','precision','%16.9e','delimiter','')
          dlmwrite(filename,[rmaxis zmaxis simag sibry bcentr],'-append','precision','%16.9e','delimiter','')
          dlmwrite(filename,[current simag xdum rmaxis xdum],'-append','precision','%16.9e','delimiter','')
          dlmwrite(filename,[zmaxis xdum sibry xdum xdum],'-append','precision','%16.9e','delimiter','')
          for ind=1:floor(nw/5)
              dlmwrite(filename,fpol(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
          end
          dlmwrite(filename,fpol(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
          for ind=1:floor(nw/5)
              dlmwrite(filename,pres(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
          end
          dlmwrite(filename,pres(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
          for ind=1:floor(nw/5)
              dlmwrite(filename,ffprim(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
          end
          dlmwrite(filename,ffprim(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
          for ind=1:floor(nw/5)
              dlmwrite(filename,pprime(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
          end
          dlmwrite(filename,pprime(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
          psirz_reshape=reshape(psirz',1,nw*nh);
          for ind=1:floor(nw*nh/5)
              dlmwrite(filename,psirz_reshape(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
          end
          dlmwrite(filename,psirz_reshape(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
          for ind=1:floor(nw/5)
              dlmwrite(filename,qpsi(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
          end
          dlmwrite(filename,qpsi(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
          dlmwrite(filename,[nbbbs limitr],'-append','precision','  %i','delimiter','')
          for ind=1:floor(nbbbs*2/5)
              dlmwrite(filename,bbbs(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
          end
          dlmwrite(filename,bbbs(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
          for ind=1:floor(limitr*2/5)
              dlmwrite(filename,lim(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
          end
          dlmwrite(filename,lim(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
          
          disp([filename ' saved'])
          %%
          % dlmwrite('wall_rz',rbbbs')
          % dlmwrite('wall_rz',zbbbs','-append')

          %type(filename)
           
      end
      
end
