function ids_equilibrium_out = improve_equilibrium_with_feeqs(ids_equilibrium_in,mode_debug,alternative_extrapolation)

% set path to FEEQS.M
% get path to FEEQS
set_reactor_feeqs_path;

% mode debug
if (nargin < 2) || isempty(mode_debug)
    mode_debug = false;
end

if (nargin < 3) || isempty(alternative_extrapolation)
    alternative_extrapolation = 'Interpolation';
end


% get time
if ids_equilibrium_in.ids_properties.homogeneous_time == 1
   time = ids_equilibrium_in.time; 
else
   time  = NaN * ones(legnth(ids_equilibrium_in.time_slice),1);
   for k = 1:length(time)
      time(k) =  ids_equilibrium_in.time_slice{k}.time;
   end
end
    
% allocate output
ids_equilibrium_out = ids_init('equilibrium');
ids_equilibrium_out.time_slice = ids_allocate('equilibrium','time_slice',length(time));
ids_equilibrium_out.time       = time;
for k=1:length(time)
    ids_equilibrium_out.time_slice{k}.time = time(k);
end
ids_equilibrium_out.vacuum_toroidal_field.r0  = ids_equilibrium_in.vacuum_toroidal_field.r0;
ids_equilibrium_out.vacuum_toroidal_field.b0  = ids_equilibrium_in.vacuum_toroidal_field.b0;

% fill ids_prorperties and code
noms = fieldnames(ids_equilibrium_in.code);
for k=1:length(noms)
    ids_equilibrium_out.code.(noms{k}) = ids_equilibrium_in.code.(noms{k});
end
ids_equilibrium_out.code.name = 'Fixed boundary FEEQS.M postprocessing of METIS';
ids_equilibrium_out.code.output_flag = NaN .* ids_equilibrium_out.code.output_flag;
%
noms = fieldnames(ids_equilibrium_in.ids_properties);
for k=1:length(noms)
    ids_equilibrium_out.ids_properties.(noms{k}) = ids_equilibrium_in.ids_properties.(noms{k});
end
ids_equilibrium_out.ids_properties.comment = 'FEEQS.M fixed equilibrium from METIS one';
ids_equilibrium_out.ids_properties.source  = 'FEEQS.M call from METIS';  


% loop on time slices
for k=1:length(time)
   % prepare call of FEEQS.M
   q            = ids_equilibrium_in.time_slice{k}.profiles_1d.q;
   sign_q       = sign(q(end));
   F            = ids_equilibrium_in.time_slice{k}.profiles_1d.f;
   Ptot         = ids_equilibrium_in.time_slice{k}.profiles_1d.pressure;
   sign_F       = sign(F(end));      
   psi_eq       = ids_equilibrium_in.time_slice{k}.profiles_1d.psi;
   if psi_eq(1) > psi_eq(end)
       sign_psi = 1;
   else
       sign_psi = -1;
   end
   psin         = (psi_eq - psi_eq(1)) ./ (psi_eq(end) - psi_eq(1));
   phi          = ids_equilibrium_in.time_slice{k}.profiles_1d.phi;   
   % 2*pi factor detection
   qloc = pdederive(psi_eq(:),phi(:),2,2,1,1);
   % psi_feeqs = psi_imas / factor_2pi * sign_psi
   if abs(q(end)./qloc(end)) > 0.8
      factor_2pi = 2 * pi;
   else
      factor_2pi = 1;
   end
   fdfdpsi      = ids_equilibrium_in.time_slice{k}.profiles_1d.f_df_dpsi;
   dptotdpsi    = ids_equilibrium_in.time_slice{k}.profiles_1d.dpressure_dpsi;
   j_tor        = ids_equilibrium_in.time_slice{k}.profiles_1d.j_tor;   
   gm2          = ids_equilibrium_in.time_slice{k}.profiles_1d.gm2;     
   li           = ids_equilibrium_in.time_slice{k}.global_quantities.li_3;
   betap        = ids_equilibrium_in.time_slice{k}.global_quantities.beta_pol;
   ip           = ids_equilibrium_in.time_slice{k}.global_quantities.ip;
   sign_ip      = sign(ip);
   %
   
   % handle signs and 2*pi factor
   ip       = sign_ip .* ip;
   j_tor    = sign_ip * j_tor;
   psi_eq   = sign_psi .* psi_eq ./ factor_2pi;
   q        = sign_q * q;
   F        = sign_F .* F;
   phi      = sign_F .* phi;
   %F(:) .* pdederive(psi_eq(:),F(:),2,2,1,1)
   fdfdpsi = sign_psi .* factor_2pi .* fdfdpsi;
   % security on sign for FEEQS
   fdfdpsi = sign(mean(F .* pdederive(psi_eq,F,0,2,2,1)) .* fdfdpsi) .* fdfdpsi;
   %pdederive(psi_eq(:),Ptot(:),2,2,1,1)
   dptotdpsi = sign_psi .* factor_2pi .* dptotdpsi;
   % security on sign for FEEQS
   dptotdpsi = sign(mean(pdederive(psi_eq,Ptot,0,2,2,1)) .* dptotdpsi) .* dptotdpsi;

   % boundary
   R_LCFS  = ids_equilibrium_in.time_slice{k}.boundary.outline.r;
   Z_LCFS  = ids_equilibrium_in.time_slice{k}.boundary.outline.z;
   if (R_LCFS(1) ~= R_LCFS(end)) || (Z_LCFS(1) ~= Z_LCFS(end))
        R_LCFS(end + 1) = R_LCFS(1);
        Z_LCFS(end + 1) = Z_LCFS(1);
   end
   
   % call of fixed boundary solver of FEEQS
   % prepare data for FEEQS fast fixed boundary equilibrium
   fprintf('Computation of precise equilibrium with FEEQS fast mode @ t = %g (s); index = %d:\n',time(k),k);
   equiok = false;
   plist = [];
   % plasma current from Metis
   plist.Ip = abs(ip);   % total plasma current
   
   % LCFS
   LCFS = cat(2,R_LCFS(:),Z_LCFS(:));
   
   % value of psi on LCFS
   psi_ref = psi_eq(end);
   % Edge flux is unknown
   plist.psi_bd     = psi_ref;
   %  magnetic reigidity
   plist.RB0        = abs(F(end));
   
   % turn on post processing and tune options
   %    plist.iter_max         = maximal number of iterations (on Newton convergence loop, 100)
   %    plist.NewtonStop       = Newton loop stoping criterium (1e-10)
   if mode_debug
       plist.verbose          = true; % for testing
       plist.plotonoff        = true; % for testing
   else
       plist.verbose          = false; % for testing
       plist.plotonoff        = false; % for testing      
   end
   %    plist.p                = hpFEM order  (7 is the best trend off for JT-60SA ; 3 provide smoother 2D grid for inverse grid and are faster)
   plist.twoD_data        = false; %set on or off computation of 2D data Psi, BR, BZ and Bphi needed for graphical display and inverse grid computation (false, quite slow, can be computed later)
   plist.other            = true; %set on or off false computation of other quantities like q, betap,li and flux averaged quantities (gm1..gm8) (false, can be computed later)
   %    plist.nbp_plot         = number of radial points in plot and number of points in 1D output profiles (101, don't change by mutch computation time)
   %    plist.nbp_grid2d       = number of points in each direction taken in sampling of hpFEM for graphical display and inverse grid computation (101, computation time go like p^2 *  nbp_grid2d^2)
   %    plist.match_li_beta    = when profile current are prescribed by alpha, beta and gamma; assume beta been beta_p and gamma been l_i(3) and perform convergence to match beta_p and l_i(3) (false)
   %    plist.tol_li_beta      = convergence criterium for the convergence to beta_p and l_i(3) (1e-3)
   %plist.limiter_name     = ''; %optional polodidal limiter description (empty string)
   %plist.limiter_contour_RZ(1,:) = Rw;
   %plist.limiter_contour_RZ(2,:) = Zw;
   %
   plist.space_filter = 7;
   %    poly_1d_order          = polynomial order use in flux surface averaged quatities computation (30)
   plist.p = 7;
   plist.nbp_grid2d = 101;
   plist.newton_convergence_mode = 1;
   plist.match_li_beta           = false;
   
   % mode p = 7
   disp('try p = 7');
   try
       [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS,plist.Ip,fdfdpsi,dptotdpsi,psi_eq,plist);
       if isempty(psi_ax) || ~isfinite(psi_ax)
           error('just to trap it');
       elseif lambda < sqrt(eps)
           error('just to trap it');
       elseif (plist_out.residual  <= plist_out.NewtonStop) &&  ...
               (plist_out.fact_ip_inte > 0.5) && ...
               (plist_out.fact_ip_inte < 1.5)
           equiok = true;
       else
           error('just to trap it');
       end
   catch
           disp('non convergence, switch to p = 5');
           plist.p = 5;
           equiok = false;
   end
   if ~equiok
       try
           [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS,plist.Ip,fdfdpsi,dptotdpsi,psi_eq,plist);
           if isempty(psi_ax) || ~isfinite(psi_ax)
               error('just to trap it');
           elseif lambda < sqrt(eps)
               error('just to trap it');
           elseif (plist_out.residual  <= plist_out.NewtonStop) &&  ...
                   (plist_out.fact_ip_inte > 0.5) && ...
                   (plist_out.fact_ip_inte < 1.5)
               equiok = true;
           else
               error('just to trap it');
           end
       catch
           disp('non convergence, switch to p = 3');
           plist.p = 3;
           equiok = false;
       end
   end
   if ~equiok
       try
           [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS,plist.Ip,fdfdpsi,dptotdpsi,psi_eq,plist);
           if isempty(psi_ax) || ~isfinite(psi_ax)
               error('just to trap it');
           elseif lambda < sqrt(eps)
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
               [jj_LCFS,psi,DMaps,psi_ax,lambda,plist_out] = FixedBoundaryEquilibrium_hpFEM(LCFS,plist.Ip,li,betap,1,plist);
               if lambda < sqrt(eps)
                   disp('Unable to converge also with moments');
               elseif ~isempty(psi_ax) || ~isfinite(psi_ax)
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
  
   if plist.plotonoff
       drawnow
   end
   if equiok
       disp('equilibrium found, starting post processing');
       ids_equilibrium_out.code.output_flag(k) = 0;
       if isfield(plist_out,'hf')
           plist.hf = plist_out.hf;
       end
       if isfield(plist_out,'hf_q')
           plist.hf_q = plist_out.hf_q;
           figure(plist.hf_q)
           hold on
           plot(psin,q,'r');
           legend('FEEQS','METIS');
           drawnow
       end
       % detail comparison with METIS
       if plist.plotonoff
           hfm = findobj(0,'type','figure','tag','nice_feeqs_cp');
           if isempty(hfm)
               hfm = figure('tag','nice_feeqs_cp');
           else
               figure(hfm);
               clf;
           end
           subplot(2,4,1)
           plot(plist_out.psin,plist_out.q,'b',psin,q,'r');
           %xlabel('Psi normalized');
           ylabel('q');
           legend('FEEQS','METIS');
           subplot(2,4,2)
           plot(plist_out.psin,plist_out.jphi,'b',psin,j_tor,'r');
           %xlabel('Psi normalized');
           ylabel('<J_{phi}> (A/m^2)');
           %legend('FEEQS','METIS');
           subplot(2,4,3)
           plot(plist_out.psin,plist_out.ptot,'b',psin,Ptot,'r');
           %xlabel('Psi normalized');
           ylabel('P_{tot} (Pa)');
           %legend('FEEQS','METIS');
           subplot(2,4,4)
           plot(plist_out.psin,plist_out.fdia,'b',psin,F,'r');
           %xlabel('Psi normalized');
           ylabel('F (T.m)');
           %legend('FEEQS','METIS');
           subplot(2,4,5)
           plot(plist_out.psin,plist_out.phi,'b',psin,phi,'r');
           xlabel('Psi normalized');
           ylabel('Phi (Wb)');
           %legend('FEEQS','METIS');
           subplot(2,4,6)
           plot(plist_out.psin,plist_out.grad_rho2or2,'b',psin,gm2,'r');
           xlabel('Psi normalized');
           ylabel('<grad(rho)^2/R^2>');
           %legend('FEEQS','METIS');
           subplot(2,4,7)
           plot(plist_out.psin,plist_out.profileHandle.Sffp(plist_out.psin,1) .* plist_out.lambda,'b',psin,fdfdpsi,'r', ...
               psin,fdfdpsi,'r.');
           xlabel('Psi normalized');
           ylabel('F*dF/dPsi');
           subplot(2,4,8)
           plot(plist_out.psin,plist_out.profileHandle.Spp(plist_out.psin,1) .* plist_out.lambda,'b',psin,dptotdpsi,'r', ...
               psin ,dptotdpsi,'r.');
           xlabel('Psi normalized');
           ylabel('dP/dPsi');
           drawnow
           fprintf('Ip FEEQS     = %g MA & Ip METIS   = %g MA\n',plist_out.ip_feeqs/1e6,abs(ip)/1e6);
           fprintf('li_3 FEEQS   = %g & li_3 METIS    = %g\n',plist_out.li_3,li);
           fprintf('beta_p FEEQS = %g & beta_p METIS  = %g\n',plist_out.beta_p,betap);
           
       end
       % get interpollant
       %[profiles_2d,FPSI,FBR,FBZ,FBPHI,hout,plist_out] = compute_inverse_grid_interp_hpFEM(plist_out,[],'rho_tor',1/pi/plist.nbp_grid2d);
       %[profiles_2d,FPSI,FBR,FBZ,FBPHI,hout,plist_out] = compute_inverse_grid_interp_hpFEM(plist_out,[],'rho_tor',0);
       [profiles_2d,FPSI,FBR,FBZ,FBPHI,hout,plist_out] = compute_inverse_grid_interp_hpFEM(plist_out,[],linspace(0,1,201),0);

       if plist.plotonoff
           % plot equilibrium moment for comparison
           hfm2 = findobj(0,'type','figure','tag','nice_feeqs_cp2');
           if isempty(hfm2)
               hfm2 = figure('tag','nice_feeqs_cp2');
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
           
           % NICE data
           Raxe = (ids_equilibrium_in.time_slice{k}.profiles_1d.r_outboard + ids_equilibrium_in.time_slice{k}.profiles_1d.r_inboard) / 2;           
           subplot(2,2,2)
           plot(psin_feeqs,Raxe_feeqs - Raxe_feeqs(end),psin,Raxe - Raxe(end))
           xlabel('Psi_n')
           ylabel('delta (m)')
           subplot(2,2,1)
           plot(psin_feeqs,Raxe_feeqs,psin,Raxe)
           xlabel('Psi_n')
           ylabel('Raxe (m)')
           legend('FEEQS','METIS')
           subplot(2,2,3)
           plot(psin_feeqs,kx_feeqs,psin,ids_equilibrium_in.time_slice{k}.profiles_1d.elongation)
           xlabel('Psi_n')
           ylabel('kx')
           subplot(2,2,4)
           plot(psin_feeqs,dx_feeqs,psin,(ids_equilibrium_in.time_slice{k}.profiles_1d.triangularity_upper + ids_equilibrium_in.time_slice{k}.profiles_1d.triangularity_lower)/2)
           xlabel('Psi_n')
           ylabel('dx')
       end
       
       % fill equilibrium_out
       % vaccum field
       
       % boundary by definition is the same
       ids_equilibrium_out.time_slice{k}.boundary = ids_equilibrium_in.time_slice{k}.boundary;
       ids_equilibrium_out.time_slice{k}.boundary.elongation = mean(jj_LCFS.lcfs_kappa);
       ids_equilibrium_out.time_slice{k}.boundary.elongation_upper = jj_LCFS.lcfs_kappa(1);
       ids_equilibrium_out.time_slice{k}.boundary.elongation_lower = jj_LCFS.lcfs_kappa(2);
       ids_equilibrium_out.time_slice{k}.boundary.triangularity    = mean(jj_LCFS.lcfs_delta);
       ids_equilibrium_out.time_slice{k}.boundary.triangularity_upper = jj_LCFS.lcfs_delta(1);
       ids_equilibrium_out.time_slice{k}.boundary.triangularity_lower = jj_LCFS.lcfs_delta(2);
       ids_equilibrium_out.time_slice{k}.boundary.geometric_axis.r    = jj_LCFS.lcfs_C(1);
       ids_equilibrium_out.time_slice{k}.boundary.geometric_axis.z    = jj_LCFS.lcfs_C(2);
       ids_equilibrium_out.time_slice{k}.boundary.minor_radius        = jj_LCFS.lcfs_a;
       
       % global quantities
       ids_equilibrium_out.time_slice{k}.global_quantities.beta_pol = plist_out.beta_p;
       betan  = plist_out.W_MHD .* (1.6.*pi./3) .* jj_LCFS.lcfs_a ./ plist_out.volume ./ (plist_out.fdia(end) ./ jj_LCFS.R0) ./ plist_out.ip_feeqs;
       ids_equilibrium_out.time_slice{k}.global_quantities.beta_tor = (2/3) * plist_out.W_MHD./ plist_out.volume ./ (plist_out.fdia(end) ./ jj_LCFS.R0) .^ 2 .* 2 .* plist_out.mu0; 
       ids_equilibrium_out.time_slice{k}.global_quantities.beta_normal = 100 .* betan;
       ids_equilibrium_out.time_slice{k}.global_quantities.ip = sign_ip .* plist_out.ip_feeqs;
       ids_equilibrium_out.time_slice{k}.global_quantities.li_3 = plist_out.li_3;
       ids_equilibrium_out.time_slice{k}.global_quantities.volume = plist_out.volume;
       ids_equilibrium_out.time_slice{k}.global_quantities.area =   plist_out.section;
       ids_equilibrium_out.time_slice{k}.global_quantities.surface = plist_out.len_LCFS .* 2 .* pi .* ids_equilibrium_out.time_slice{k}.boundary.geometric_axis.r;
       ids_equilibrium_out.time_slice{k}.global_quantities.length_pol = plist_out.len_LCFS;
       ids_equilibrium_out.time_slice{k}.global_quantities.psi_axis = sign_psi .*  factor_2pi .* plist_out.psi_ax;
       ids_equilibrium_out.time_slice{k}.global_quantities.psi_boundary = sign_psi .*  factor_2pi .* plist_out.psi_bd;
       ids_equilibrium_out.time_slice{k}.global_quantities.magnetic_axis.r  = profiles_2d.r(1,1);
       ids_equilibrium_out.time_slice{k}.global_quantities.magnetic_axis.z  = profiles_2d.z(1,1);
       ids_equilibrium_out.time_slice{k}.global_quantities.magnetic_axis.b_field_tor = sign_F .* profiles_2d.b_field_tor(1,1);
       ids_equilibrium_out.time_slice{k}.global_quantities.magnetic_axis.b_tor  = sign_F .* profiles_2d.b_tor(1,1);
       %ids_equilibrium_out.time_slice{k}.global_quantities.current_centre: [1x1 struct]
       ids_equilibrium_out.time_slice{k}.global_quantities.q_axis = sign_q .* plist_out.q(1);
       ids_equilibrium_out.time_slice{k}.global_quantities.q_95 = sign_q .* interp1(plist_out.psin,plist_out.q,0.95,'pchip','extrap');
       [ids_equilibrium_out.time_slice{k}.global_quantities.q_min.value,iqmin] = min(plist_out.q);
       ids_equilibrium_out.time_slice{k}.global_quantities.q_min.value = sign_q .* ids_equilibrium_out.time_slice{k}.global_quantities.q_min.value; 
       ids_equilibrium_out.time_slice{k}.global_quantities.q_min.rho_tor_norm  = mean(plist_out.rho_tor_norm(iqmin));
       ids_equilibrium_out.time_slice{k}.global_quantities.energy_mhd =  plist_out.W_MHD;
       ids_equilibrium_out.time_slice{k}.global_quantities.w_mhd = plist_out.W_MHD;
       %ids_equilibrium_out.time_slice{k}.global_quantities.psi_external_average = 
       %ids_equilibrium_out.time_slice{k}.global_quantities.v_external = 
       %ids_equilibrium_out.time_slice{k}.global_quantities.plasma_inductance = 
       %ids_equilibrium_out.time_slice{k}.global_quantities.plasma_resistance = 
       
       if mode_debug
           % for testing
           disp('global_quantities in & out');
           noms = fieldnames(ids_equilibrium_out.time_slice{k}.global_quantities);
           for ln = 1:length(noms)
               if ~isstruct(ids_equilibrium_out.time_slice{k}.global_quantities.(noms{ln}))
                   fprintf('%s = %g & %g\n',noms{ln},ids_equilibrium_in.time_slice{k}.global_quantities.(noms{ln}),ids_equilibrium_out.time_slice{k}.global_quantities.(noms{ln}));
               else
                   subnoms = fieldnames(ids_equilibrium_out.time_slice{k}.global_quantities.(noms{ln}));
                   for lm =1:length(subnoms)
                       fprintf('%s.%s = %g & %g\n',noms{ln},subnoms{lm}, ...
                           ids_equilibrium_in.time_slice{k}.global_quantities.(noms{ln}).(subnoms{lm}),ids_equilibrium_out.time_slice{k}.global_quantities.(noms{ln}).(subnoms{lm}));
                       
                   end
               end
           end
       end
       
       % profiles 1D
       ids_equilibrium_out.time_slice{k}.profiles_1d.psi = sign_psi .*  factor_2pi .* plist_out.psi1D;
       ids_equilibrium_out.time_slice{k}.profiles_1d.phi = sign_F .* plist_out.phi;
       ids_equilibrium_out.time_slice{k}.profiles_1d.pressure = plist_out.ptot + ids_equilibrium_in.time_slice{k}.profiles_1d.pressure(end); % as it is just integral of dP/dpsi
       ids_equilibrium_out.time_slice{k}.profiles_1d.f = sign_F .* plist_out.fdia;
       ids_equilibrium_out.time_slice{k}.profiles_1d.dpressure_dpsi = sign_psi ./ factor_2pi .* plist_out.profileHandle.Spp(plist_out.psin,1) .* plist_out.lambda;
       ids_equilibrium_out.time_slice{k}.profiles_1d.f_df_dpsi = sign_F .*  sign_psi ./ factor_2pi .* plist_out.profileHandle.Sffp(plist_out.psin,1) .* plist_out.lambda;
       ids_equilibrium_out.time_slice{k}.profiles_1d.j_tor = sign_ip .* plist_out.jphi;
       ids_equilibrium_out.time_slice{k}.profiles_1d.j_parallel = sign_ip .* plist_out.j_par;
       ids_equilibrium_out.time_slice{k}.profiles_1d.q = sign_q .* plist_out.q;
       %pp = polyfit(plist_out.rho,plist_out.q,7);
       %sq = polyval(polyder(pp),plist_out.rho)./ polyval(pp,plist_out.rho) .* plist_out.rho; 
       ids_equilibrium_out.time_slice{k}.profiles_1d.magnetic_shear = pdederive(plist_out.rho,plist_out.q,0,2,1,1) ./ plist_out.q .* plist_out.rho;
       ids_equilibrium_out.time_slice{k}.profiles_1d.r_inboard  = min(profiles_2d.r,[],2);
       ids_equilibrium_out.time_slice{k}.profiles_1d.r_outboard = max(profiles_2d.r,[],2);
       ids_equilibrium_out.time_slice{k}.profiles_1d.rho_tor = plist_out.rho;
       ids_equilibrium_out.time_slice{k}.profiles_1d.rho_tor_norm = plist_out.rho_tor_norm;
       ids_equilibrium_out.time_slice{k}.profiles_1d.dpsi_drho_tor = sign_psi .*  factor_2pi .* pdederive(plist_out.rho,plist_out.psi1D,0,2,1,1);
       ids_equilibrium_out.time_slice{k}.profiles_1d.dpsi_drho_tor = -sign_psi .*  factor_2pi ./ plist_out.drhodpsi;
       ids_equilibrium_out.time_slice{k}.profiles_1d.dpsi_drho_tor(1) = 0;
       ids_equilibrium_out.time_slice{k}.profiles_1d.geometric_axis.r = (max(profiles_2d.r,[],2) + min(profiles_2d.r,[],2)) ./ 2;
       ids_equilibrium_out.time_slice{k}.profiles_1d.geometric_axis.z = (max(profiles_2d.z,[],2) + min(profiles_2d.z,[],2)) ./ 2;
       %
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
        %
       ids_equilibrium_out.time_slice{k}.profiles_1d.elongation = medfilt1(kx_feeqs,3);
       feeqs_dup(1) = feeqs_dup(2);
       ids_equilibrium_out.time_slice{k}.profiles_1d.triangularity_upper = medfilt1(feeqs_dup,3);
       feeqs_dlow(1) = feeqs_dlow(2);
       ids_equilibrium_out.time_slice{k}.profiles_1d.triangularity_lower = medfilt1(feeqs_dlow,3);
%        ids_equilibrium_out.time_slice{k}.profiles_1d.squareness_upper_inner = 
%        ids_equilibrium_out.time_slice{k}.profiles_1d.squareness_upper_outer = 
%        ids_equilibrium_out.time_slice{k}.profiles_1d.squareness_lower_inner = 
%        ids_equilibrium_out.time_slice{k}.profiles_1d.squareness_lower_outer = 
       ids_equilibrium_out.time_slice{k}.profiles_1d.volume = cumtrapz(plist_out.psin,plist_out.dVdpsin);
       ids_equilibrium_out.time_slice{k}.profiles_1d.volume = ids_equilibrium_out.time_slice{k}.profiles_1d.volume ./ ...
                  ids_equilibrium_out.time_slice{k}.profiles_1d.volume(end) .* plist_out.volume;     
       ids_equilibrium_out.time_slice{k}.profiles_1d.rho_volume_norm = ids_equilibrium_out.time_slice{k}.profiles_1d.volume ./ ...
                  ids_equilibrium_out.time_slice{k}.profiles_1d.volume(end);              
       ids_equilibrium_out.time_slice{k}.profiles_1d.dvolume_dpsi = plist_out.dVdpsin ./  ...
                  (ids_equilibrium_out.time_slice{k}.profiles_1d.psi(end) -ids_equilibrium_out.time_slice{k}.profiles_1d.psi(1));
       ids_equilibrium_out.time_slice{k}.profiles_1d.dvolume_drho_tor = plist_out.dVdrho;
       ids_equilibrium_out.time_slice{k}.profiles_1d.area = cumtrapz(plist_out.psin,plist_out.dSdpsin);
       ids_equilibrium_out.time_slice{k}.profiles_1d.area = ids_equilibrium_out.time_slice{k}.profiles_1d.area ./ ...
                  ids_equilibrium_out.time_slice{k}.profiles_1d.area(end) .* plist_out.section;
       ids_equilibrium_out.time_slice{k}.profiles_1d.darea_dpsi = plist_out.dSdpsin ./  ...
                  (ids_equilibrium_out.time_slice{k}.profiles_1d.psi(end) -ids_equilibrium_out.time_slice{k}.profiles_1d.psi(1));
       ids_equilibrium_out.time_slice{k}.profiles_1d.darea_dpsi(1) = ids_equilibrium_out.time_slice{k}.profiles_1d.darea_dpsi(2);     
       ids_equilibrium_out.time_slice{k}.profiles_1d.darea_drho_tor = plist_out.dSdrho;
       %
       r2d = profiles_2d.r;
       z2d = profiles_2d.z;
       r2d = cat(2,r2d,r2d(:,1));
       z2d = cat(2,z2d,z2d(:,1));
       dl  = sqrt(diff(r2d,1,2) .^ 2 + diff(z2d,1,2) .^ 2);
       len_surf = sum(dl,2);
       len_surf = len_surf ./ len_surf(end) .* plist_out.len_LCFS;     
       ids_equilibrium_out.time_slice{k}.profiles_1d.surface = len_surf .*  ...
              ids_equilibrium_out.time_slice{k}.profiles_1d.geometric_axis.r .* 2 .* pi;
       ids_equilibrium_out.time_slice{k}.profiles_1d.trapped_fraction = plist_out.ftrap;
       ids_equilibrium_out.time_slice{k}.profiles_1d.gm1 = plist_out.ri2;
       ids_equilibrium_out.time_slice{k}.profiles_1d.gm2 = plist_out.grad_rho2or2;
       ids_equilibrium_out.time_slice{k}.profiles_1d.gm3 = plist_out.grad_rho2;
       %
       % rougth approximation
       BNORM2 = profiles_2d.b_field_r .^2 + profiles_2d.b_field_z .^2 +  profiles_2d.b_field_tor .^ 2;
       b2ave  = sum(BNORM2 .* dl,2) ./ sum(dl,2);
       b2ave(1) = b2ave(2);
       bave  = sum(sqrt(BNORM2) .* dl,2) ./ sum(dl,2);
       bave(1) = bave(2);
       b2iave  = sum(1./ BNORM2 .* dl,2) ./ sum(dl,2);
       b2iave(1) = b2iave(2);
       grho2_2d  = profiles_2d.r .^ 2 .* plist_out.drhodpsi .^ 2 .*  ...
                        (profiles_2d.b_field_r .^2 + profiles_2d.b_field_z .^2);
       grho2b2   =  sum(grho2_2d./ BNORM2 .* dl,2) ./ sum(dl,2);      
       grho2b2(1) = grho2b2(2);
       %
       ids_equilibrium_out.time_slice{k}.profiles_1d.gm4 = b2iave;
       ids_equilibrium_out.time_slice{k}.profiles_1d.gm5 = b2ave;
       ids_equilibrium_out.time_slice{k}.profiles_1d.gm6 = grho2b2;
       ids_equilibrium_out.time_slice{k}.profiles_1d.gm7 = plist_out.grad_rho;
       ids_equilibrium_out.time_slice{k}.profiles_1d.gm8 = plist_out.rave;
       ids_equilibrium_out.time_slice{k}.profiles_1d.gm9 = plist_out.ri;
       %ids_equilibrium_out.time_slice{k}.profiles_1d.b_average = sqrt(plist_out.bpol2 + plist_out.fdia .^ 2  .* plist_out.ri2); % Metis estimation
       ids_equilibrium_out.time_slice{k}.profiles_1d.b_average = bave; % direct estimation
       ids_equilibrium_out.time_slice{k}.profiles_1d.b_field_average = bave; % direct estimation
       ids_equilibrium_out.time_slice{k}.profiles_1d.b_min = min(sqrt(BNORM2),[],2);
       ids_equilibrium_out.time_slice{k}.profiles_1d.b_field_min = min(sqrt(BNORM2),[],2);
       ids_equilibrium_out.time_slice{k}.profiles_1d.b_max = max(sqrt(BNORM2),[],2);
       ids_equilibrium_out.time_slice{k}.profiles_1d.b_field_max = max(sqrt(BNORM2),[],2);
       ids_equilibrium_out.time_slice{k}.profiles_1d.beta_pol = cumtrapz(ids_equilibrium_out.time_slice{k}.profiles_1d.volume, ...
                ids_equilibrium_out.time_slice{k}.profiles_1d.pressure) ./ ids_equilibrium_out.time_slice{k}.profiles_1d.volume ./ ...
                plist_out.bpol2 .* 2 .* plist_out.mu0; 
       ids_equilibrium_out.time_slice{k}.profiles_1d.beta_pol(1) = ids_equilibrium_out.time_slice{k}.profiles_1d.beta_pol(2); 
       ids_equilibrium_out.time_slice{k}.profiles_1d.beta_pol = ids_equilibrium_out.time_slice{k}.profiles_1d.beta_pol ./ ...
                ids_equilibrium_out.time_slice{k}.profiles_1d.beta_pol(end) .* ids_equilibrium_out.time_slice{k}.global_quantities.beta_pol;
       %ids_equilibrium_out.time_slice{k}.profiles_1d.mass_density = ?
       
       if mode_debug
           % for testing
           disp('profiles_1d in (blue) & out (red)');
           noms = fieldnames(ids_equilibrium_out.time_slice{k}.profiles_1d);
           figure;
           lsp = 1;
           for ln = 1:length(noms)
               if isempty(strfind(noms{ln},'_error_'))
                   if ~isstruct(ids_equilibrium_out.time_slice{k}.profiles_1d.(noms{ln}))
                       try
                           subplot(3,4,lsp);
                           try
                               plot(ids_equilibrium_out.time_slice{k}.profiles_1d.rho_tor_norm,ids_equilibrium_out.time_slice{k}.profiles_1d.(noms{ln}),'r');
                           end
                           hold on
                           try
                               plot(ids_equilibrium_in.time_slice{k}.profiles_1d.rho_tor_norm,ids_equilibrium_in.time_slice{k}.profiles_1d.(noms{ln}),'b');
                           end
                           xlabel('rho_tor_norm');
                           ylabel(noms{ln});
                           title('red -> FEEQS & blue -> METIS');
                           lsp = lsp +1;
                       catch
                           disp(noms{ln})
                       end
                   else
                       subnoms = fieldnames(ids_equilibrium_out.time_slice{k}.profiles_1d.(noms{ln}));
                       for lm =1:length(subnoms)
                           if isempty(strfind(subnoms{lm},'_error_'))
                               try
                                   subplot(3,4,lsp);
                                   try
                                       plot(ids_equilibrium_out.time_slice{k}.profiles_1d.rho_tor_norm,ids_equilibrium_out.time_slice{k}.profiles_1d.(noms{ln}).(subnoms{lm}),'r');
                                   end
                                   hold on
                                   try
                                       plot(ids_equilibrium_in.time_slice{k}.profiles_1d.rho_tor_norm,ids_equilibrium_in.time_slice{k}.profiles_1d.(noms{ln}).(subnoms{lm}),'b');
                                   end
                                   xlabel('rho_tor_norm');
                                   ylabel(sprintf('%s/%s',noms{ln},subnoms{lm}));
                                   title('red -> FEEQS & blue -> METIS');
                                   lsp = lsp +1;
                               catch
                                   fprintf('%s/%s\n',noms{ln},subnoms{lm})
                                   
                               end
                           end
                       end
                   end
                   if lsp > 12
                       figure;
                       lsp = 1;
                   end
               end
               drawnow
           end
       end

       
       % profiles 2D (2 grids polar and rectangular)       
       % search for grid size
       if ids_equilibrium_in.time_slice{k}.profiles_2d{1}.grid_type.index == 1
           index_rectangular = 1;
           index_polar = 2;
       else
           index_rectangular = 1;
           index_polar = 2;           
       end
       % size of EQDSK grid
       rr = ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.grid.dim1;
       zz = ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.grid.dim2;
       % improve resolution - > keep the original one controller in GUI
       %rr = linspace(min(rr),max(rr),min(151,length(rr)*3));
       %zz = linspace(min(zz),max(zz),min(271,length(zz)*3));
       % start 2D grid ala EQDSK
       ids_equilibrium_out.time_slice{k}.profiles_2d(index_rectangular) = ids_allocate('equilibrium','time_slice/profiles_2d',1);
       [rgrid, zgrid] = meshgrid(rr,zz);
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.type.name = 'total';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.type.index = 0;
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.type.description = 'from fixed boundary equilibrium with extrapolation';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.grid_type.name = 'rectangular';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.grid_type.index = 1;
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.grid_type.description = 'cylindrical R,Z ala eqdsk, within the corresponding COCOS convention (COCOS=11 is assumed in ITER)';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.grid.dim1   = rr;
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.grid.dim2   = zz;
       %ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.grid.volume_element: []
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.r           = rgrid';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.z           = zgrid';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.psi         = sign_psi .* factor_2pi .* FPSI(rgrid,zgrid)';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.theta       = unwrap(angle((rgrid - profiles_2d.r(1,1)) + sqrt(-1) .* (zgrid - profiles_2d.z(1,1))));
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.phi         = sign_F .* interp1(sign_psi .* factor_2pi .* plist_out.psi1D,plist_out.phi,ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.psi,'pchip',NaN);       
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_r         = -sign_ip .* FBR(rgrid,zgrid)';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_r   = ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_r;
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_z         = -sign_ip .* FBZ(rgrid,zgrid)';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_z   = ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_z;
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_tor       = sign_F .*  FBPHI(rgrid,zgrid)';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_tor = ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_tor;
       
       % alternative extrapolation ouside the LCFS
       switch alternative_extrapolation
           case 'G-S polynomial'
               PSI_LCFS = mean(FPSI(R_LCFS,Z_LCFS)) * ones(size(R_LCFS));
               BR_LCFS  = FBR(R_LCFS,Z_LCFS);
               BZ_LCFS  = FBZ(R_LCFS,Z_LCFS);  
               
%                figure
%                lcl = cat(2,linspace(profiles_2d.psi(1,1), 1.3 .* profiles_2d.psi(end,1),21),linspace(profiles_2d.psi(end,1), 2 .* profiles_2d.psi(end,1) - profiles_2d.psi(1,1),51));
%                contour(rgrid,zgrid,ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.psi',sign_psi .* factor_2pi .* lcl,'color','r');
%                hold on
               

               % 
               %dbstop if all error
               [FPSI_ext,FBR_ext,FBZ_ext] = extrapolate_from_LCFS(R_LCFS,Z_LCFS,PSI_LCFS,BR_LCFS,BZ_LCFS);
               mask_out = ~zinout(R_LCFS,Z_LCFS,rgrid,zgrid);
               %
               PSI_grid = FPSI(rgrid,zgrid);
               PSI_grid(mask_out) =  FPSI_ext(rgrid(mask_out),zgrid(mask_out));
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.psi = sign_psi .* factor_2pi .* PSI_grid';
               %
               BR_grid = FBR(rgrid,zgrid);
               BR_grid(mask_out) =  FBR_ext(rgrid(mask_out),zgrid(mask_out));
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_r                   = -sign_ip .* BR_grid';
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_r             = ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_r;
               %
               BZ_grid = FBZ(rgrid,zgrid);
               BZ_grid(mask_out) =  FBZ_ext(rgrid(mask_out),zgrid(mask_out));
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_z                   = -sign_ip .* BZ_grid';
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_z             = ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_z;
               
%                contour(rgrid,zgrid,ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.psi',sign_psi .* factor_2pi .* lcl,'color','b');
%                hold on
%                plot(R_LCFS,Z_LCFS,'k',R_LCFS,Z_LCFS,'k.');
%                drawnow
%                %keyboard
               
               
       end
       
       %
       pprim  = interp1(sign_psi .* factor_2pi .* plist_out.psi1D, plist_out.profileHandle.Spp(plist_out.psin,1) .* plist_out.lambda, ...
                       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.psi,'pchip',0);
       ffprim = interp1(sign_psi .* factor_2pi .* plist_out.psi1D, plist_out.profileHandle.Sffp(plist_out.psin,1) .* plist_out.lambda, ...
                       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.psi,'pchip',0);                   
       jphi  = pprim .* rgrid' + ffprim ./ rgrid' ./ plist_out.mu0; 
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.j_tor  = sign_ip .* jphi;
       jpar  = jphi .* ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_tor + ffprim ./ plist_out.mu0 ./  ...
                     (ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_tor .* rgrid') .*  ...
 	                  (ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_r.^ 2 +  ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_z.^ 2);
       jpar(~isfinite(jpar)) = 0;          
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.j_parallel = sign_ip .* sign_ip .* jpar ./ (plist_out.RB0 ./ jj_LCFS.R0);
       
       % test
       if mode_debug
           figure;
           subplot(3,3,1)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.psi,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.psi,'or');
           set(gca,'xlim',[min(R_LCFS(:)),max(R_LCFS(:))],'ylim',[min(Z_LCFS(:)),max(Z_LCFS(:))]);
           title('psi');
           subplot(3,3,2)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.phi,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.phi,'or');
           set(gca,'xlim',[min(R_LCFS(:)),max(R_LCFS(:))],'ylim',[min(Z_LCFS(:)),max(Z_LCFS(:))]);
           title('phi');
           subplot(3,3,3)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.j_tor,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.j_tor,'or');
           set(gca,'xlim',[min(R_LCFS(:)),max(R_LCFS(:))],'ylim',[min(Z_LCFS(:)),max(Z_LCFS(:))]);
           title('jtor');
           subplot(3,3,4)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.j_parallel,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.j_parallel,'or');
           set(gca,'xlim',[min(R_LCFS(:)),max(R_LCFS(:))],'ylim',[min(Z_LCFS(:)),max(Z_LCFS(:))]);
           title('j_par');         
           subplot(3,3,5)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.b_field_tor,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_tor,'or');
           set(gca,'xlim',[min(R_LCFS(:)),max(R_LCFS(:))],'ylim',[min(Z_LCFS(:)),max(Z_LCFS(:))]);
           title('Bphi');         
           subplot(3,3,6)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.b_field_r,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_r,'or');
           set(gca,'xlim',[min(R_LCFS(:)),max(R_LCFS(:))],'ylim',[min(Z_LCFS(:)),max(Z_LCFS(:))]);
           title('BR');         
           subplot(3,3,7)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.b_field_z,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_z,'or');
           set(gca,'xlim',[min(R_LCFS(:)),max(R_LCFS(:))],'ylim',[min(Z_LCFS(:)),max(Z_LCFS(:))]);
           title('BZ');         
           figure
           quiver(ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.b_field_r, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_rectangular}.b_field_z,'b');
           hold on
           quiver(ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_r, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular}.b_field_z,'b');
           set(gca,'xlim',[min(R_LCFS(:)),max(R_LCFS(:))],'ylim',[min(Z_LCFS(:)),max(Z_LCFS(:))]);
        
       end
       
       % polar grid
       % cocos 
       profiles_2d.psi               = sign_psi .* factor_2pi .* profiles_2d.psi;
       profiles_2d.psi_error_upper   = sign_psi .* factor_2pi .* profiles_2d.psi_error_upper;
       profiles_2d.psi_error_lower   = sign_psi .* factor_2pi .* profiles_2d.psi_error_lower;
       profiles_2d.phi               = sign_F   .* interp1(sign_psi .* factor_2pi .* plist_out.psi1D,plist_out.phi,profiles_2d.psi,'pchip',NaN);
       profiles_2d.j_tor             = sign_ip  .* profiles_2d.j_tor;
       profiles_2d.j_parallel        = sign_ip  .* sign_F   .* profiles_2d.j_parallel;
       profiles_2d.b_r               = -sign_ip  .* profiles_2d.b_r;
       profiles_2d.b_field_r         = -sign_ip  .* profiles_2d.b_field_r;
       profiles_2d.b_z               = -sign_ip  .* profiles_2d.b_z;
       profiles_2d.b_field_z         = -sign_ip  .* profiles_2d.b_field_z;
       profiles_2d.b_tor             = sign_F   .* profiles_2d.b_tor;
       profiles_2d.b_field_tor       = sign_F   .* profiles_2d.b_field_tor;

       %
       ids_equilibrium_out.time_slice{k}.profiles_2d(index_polar) = ids_allocate('equilibrium','time_slice/profiles_2d',1);
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.type.name = 'total';
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.type.index = 0;
       ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.type.description = 'from fixed boundary equilibrium with extrapolation';
       noms = fieldnames(profiles_2d);
       for ln=1:length(noms)
           if ~isstruct(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.(noms{ln}))
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.(noms{ln}) = profiles_2d.(noms{ln});
           else
               subnoms = fieldnames(profiles_2d.(noms{ln}));
               for lm = 1:length(subnoms)
                   ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.(noms{ln}).(subnoms{lm}) = profiles_2d.(noms{ln}).(subnoms{lm});
               end
           end
       end
        
        % test
       if mode_debug
           figure;
           subplot(3,3,1)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.psi,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.psi,'or');
           title('psi');
           subplot(3,3,2)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.phi,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.phi,'or');
           title('phi');
           subplot(3,3,3)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.j_tor,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.j_tor,'or');
           title('jtor');
           subplot(3,3,4)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.j_parallel,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.j_parallel,'or');
           title('j_par');         
           subplot(3,3,5)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.b_field_tor,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.b_field_tor,'or');
           title('Bphi');         
           subplot(3,3,6)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.b_field_r,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.b_field_r,'or');
           title('BR');         
           subplot(3,3,7)
           plot3(ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.b_field_z,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.b_field_z,'or');
           title('BZ');         
           figure
           quiver(ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.b_field_r, ...
               ids_equilibrium_in.time_slice{k}.profiles_2d{index_polar}.b_field_z,'b');
           hold on
           quiver(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.r,  ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.z, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.b_field_r, ...
               ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar}.b_field_z,'b');
        
       end
             
       % coordinate system
       % initilialisation substructure
       ids_equilibrium_out.time_slice{k}.coordinate_system.grid_type.name = profiles_2d.grid_type.name;
       ids_equilibrium_out.time_slice{k}.coordinate_system.grid_type.index = profiles_2d.grid_type.index;
       ids_equilibrium_out.time_slice{k}.coordinate_system.grid_type.description = profiles_2d.grid_type.description;
       noms = fieldnames(profiles_2d.grid);
       for ln=1:length(noms)
            ids_equilibrium_out.time_slice{k}.coordinate_system.grid.(noms{ln}) = profiles_2d.grid.(noms{ln});
       end
       dPSIdR = -sign_psi .* factor_2pi .* profiles_2d.r.* profiles_2d.b_field_z;
       dPSIdZ =  sign_psi .* factor_2pi .* profiles_2d.r .* profiles_2d.b_field_r;
       dRdth  = NaN * ones(size(profiles_2d.r));
       dZdth  = NaN * ones(size(profiles_2d.r));
       dRdPSI  = NaN * ones(size(profiles_2d.r));
       dZdPSI  = NaN * ones(size(profiles_2d.r));
       thloc = profiles_2d.theta + 16 .* pi;
       for ll =1:size(profiles_2d.r,2)
           lp = ll+1;
           if lp > size(profiles_2d.r,2)
               lp = 1;
           end
           lm = ll-1;
           if lm < 1
               lm = size(profiles_2d.r,2);
           end
           dRdth(:,ll) = (profiles_2d.r(:,lp) - profiles_2d.r(:,lm)) ./ ...
                         (thloc(:,lp) - thloc(:,lm));
           dZdth(:,ll) = (profiles_2d.z(:,lp) - profiles_2d.z(:,lm)) ./ ...
                         (thloc(:,lp) - thloc(:,lm));
           dRdPSI(:,ll) = pdederive(sign_psi .* factor_2pi .* profiles_2d.psi(:,ll),profiles_2d.r(:,ll),2,2,1,1);
           dZdPSI(:,ll) = pdederive(sign_psi .* factor_2pi .* profiles_2d.psi(:,ll),profiles_2d.z(:,ll),2,2,1,1);
                               
       end
       dRdth(:,1) = (dRdth(:,2) + dRdth(:,end-1)) / 2;
       dRdth(:,end) = (dRdth(:,2) + dRdth(:,end-1)) / 2;
       dZdth(:,1) = (dZdth(:,2) + dZdth(:,end-1)) / 2;
       dZdth(:,end) = (dZdth(:,2) + dZdth(:,end-1)) / 2;
     
       jacobian      = abs(2 .* pi .* profiles_2d.r .* (dRdPSI .* dZdth - dZdPSI .* dRdth));
       jacobian(1,:) = 0;
        % calcul de dthetadR et dthetadZ
       dthdR = sign_psi * dZdPSI ./ max(eps,jacobian ./ (2 .* pi .* profiles_2d.r));
       dthdR(1,:) = 0;
       dthdZ = - sign_psi * dRdPSI ./ max(eps,jacobian./ (2 .* pi .* profiles_2d.r));
       dthdZ(1,:) = 0;     
        
       ids_equilibrium_out.time_slice{k}.coordinate_system.jacobian = abs(jacobian);
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant = zeros(size(dPSIdR,1),size(dPSIdR,2),3,3);
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,1,1) = (dPSIdR .^ 2 +dPSIdZ .^ 2);
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,1,2) = dthdR .* dPSIdR +  dthdZ .* dPSIdZ;
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,1,3) = 0;
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,2,1) =  dthdR .* dPSIdR +  dthdZ .* dPSIdZ;
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,2,2) =  dthdR .^2 + dthdZ .^ 2;
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,2,3) = 0;
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,3,1) = 0;
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,3,2) = 0;
       ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,3,3) = 1./ profiles_2d.r .^2;
       ids_equilibrium_out.time_slice{k}.coordinate_system.r        = profiles_2d.r;
       ids_equilibrium_out.time_slice{k}.coordinate_system.z        = profiles_2d.z;
       
       
       if mode_debug 
           figure;
           subplot(2,3,1)
           plot3(ids_equilibrium_in.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_in.time_slice{k}.coordinate_system.z, ...
                 ids_equilibrium_in.time_slice{k}.coordinate_system.jacobian,'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_out.time_slice{k}.coordinate_system.z, ...
                 ids_equilibrium_out.time_slice{k}.coordinate_system.jacobian,'or');           
           title('jacobian');
          
             
           subplot(2,3,2)
           plot3(ids_equilibrium_in.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_in.time_slice{k}.coordinate_system.z, ...
                 squeeze(ids_equilibrium_in.time_slice{k}.coordinate_system.tensor_contravariant(:,:,1,1)),'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_out.time_slice{k}.coordinate_system.z, ...
                 squeeze(ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,1,1)),'or');           
           title('(1,1)');
       
           subplot(2,3,3)
           plot3(ids_equilibrium_in.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_in.time_slice{k}.coordinate_system.z, ...
                 squeeze(ids_equilibrium_in.time_slice{k}.coordinate_system.tensor_contravariant(:,:,1,2)),'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_out.time_slice{k}.coordinate_system.z, ...
                 squeeze(ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,1,2)),'or');           
           title('(1,2)');
           
           subplot(2,3,4)
           plot3(ids_equilibrium_in.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_in.time_slice{k}.coordinate_system.z, ...
                 squeeze(ids_equilibrium_in.time_slice{k}.coordinate_system.tensor_contravariant(:,:,2,2)),'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_out.time_slice{k}.coordinate_system.z, ...
                 squeeze(ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,2,2)),'or');           
           title('(2,2)');
           
           subplot(2,3,5)
           plot3(ids_equilibrium_in.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_in.time_slice{k}.coordinate_system.z, ...
                 squeeze(ids_equilibrium_in.time_slice{k}.coordinate_system.tensor_contravariant(:,:,3,3)),'.b');
           hold on
           plot3(ids_equilibrium_out.time_slice{k}.coordinate_system.r, ...
                 ids_equilibrium_out.time_slice{k}.coordinate_system.z, ...
                 squeeze(ids_equilibrium_out.time_slice{k}.coordinate_system.tensor_contravariant(:,:,3,3)),'or');           
           title('(3,3)');
           
           drawnow
      end    
           
       % ggd
       ids_void = ids_gen('equilibrium');
       model = ids_void.time_slice{1}.ggd;
       ggd = imas_write_regular_grid(ids_equilibrium_out.time_slice{k}.profiles_2d{index_rectangular},model{1});
       ids_equilibrium_out.time_slice{k}.ggd{1} = ggd;
       ggd = imas_write_regular_grid(ids_equilibrium_out.time_slice{k}.profiles_2d{index_polar},model{1});
       ids_equilibrium_out.time_slice{k}.ggd{2} = ggd;
       % new grid definition
       ids_equilibrium_out.grids_ggd{k}.grid{1} = ids_equilibrium_out.time_slice{k}.ggd{1}.grid;
       ids_equilibrium_out.grids_ggd{k}.grid{2} = ids_equilibrium_out.time_slice{k}.ggd{2}.grid;
       ids_equilibrium_out.grids_ggd{k}.time    = time(k);
       
   else
       ids_equilibrium_out.time_slice{k} = ids_equilibrium_in.time_slice{k}; % no improvement
       ids_equilibrium_out.grids_ggd{k}  = ids_equilibrium_in.grids_ggd{k};
       ids_equilibrium_out.code.output_flag(k) = -2;
   end
      
end



% function to map ggd from regular grid store in profiles_2d
% thank to Th. Aniel
function ggd = imas_write_regular_grid(p2d,ggd)

ggd.grid.identifier.name        = p2d.grid_type.name;
ggd.grid.identifier.index       = p2d.grid_type.index;
ggd.grid.identifier.description = p2d.grid_type.description; 
%
ggd.grid.space{1}.geometry_type.name        = [];
ggd.grid.space{1}.geometry_type.index       =  -999999999;
ggd.grid.space{1}.geometry_type.description = [];   
ggd.grid.space{1}.coordinates_type    = [];
%	
n1 = length(p2d.grid.dim1); 
n2 = length(p2d.grid.dim2); 

for k1 = [1:n1]
   % initialisation sustructure
   
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1} = ggd.grid.space{1}.objects_per_dimension{1}.object{1};
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.boundary{1}.index      = -999999999;
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.boundary{1}.neighbours = [];
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.geometry(1)            = p2d.grid.dim1(k1);
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.nodes                  = [];
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.measure                = -9.0e40;
   
end

% initialisation substructure
ggd.grid.space{2} = ggd.grid.space{1};
for k2 = [1:n2]
   % initialisation sustructure
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2} = ggd.grid.space{2}.objects_per_dimension{1}.object{1};		
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.boundary{1}.index      = -999999999;
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.boundary{1}.neighbours = [];
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.geometry(1)            = p2d.grid.dim2(k2);
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.nodes                  = [];
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.measure                = -9.0e40;
   
end
   
ggd.grid.grid_subset{1}.identifier.name                = '';
ggd.grid.grid_subset{1}.identifier.index               = -999999999;
ggd.grid.grid_subset{1}.identifier.description         = '';
ggd.grid.grid_subset{1}.dimension                      = -999999999;
ggd.grid.grid_subset{1}.element{1}.object{1}.space     = -999999999;
ggd.grid.grid_subset{1}.element{1}.object{1}.dimension = -999999999;
ggd.grid.grid_subset{1}.element{1}.object{1}.index     = -999999999;
ggd.grid.grid_subset{1}.base{1}.jacobian               = [];
ggd.grid.grid_subset{1}.base{1}.g11_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g12_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g13_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g21_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g22_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g23_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g31_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g32_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g33_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g11_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g12_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g13_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g21_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g22_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g23_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g31_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g32_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g33_contravariant      = [];

ggd.grid.grid_subset{1}.metric.jacobian          = [];
ggd.grid.grid_subset{1}.metric.g11_covariant     = [];
ggd.grid.grid_subset{1}.metric.g12_covariant     = [];
ggd.grid.grid_subset{1}.metric.g13_covariant     = [];
ggd.grid.grid_subset{1}.metric.g21_covariant     = [];
ggd.grid.grid_subset{1}.metric.g22_covariant     = [];
ggd.grid.grid_subset{1}.metric.g23_covariant     = [];
ggd.grid.grid_subset{1}.metric.g31_covariant     = [];
ggd.grid.grid_subset{1}.metric.g32_covariant     = [];
ggd.grid.grid_subset{1}.metric.g33_covariant     = [];
ggd.grid.grid_subset{1}.metric.g11_contravariant = [];
ggd.grid.grid_subset{1}.metric.g12_contravariant = [];
ggd.grid.grid_subset{1}.metric.g13_contravariant = [];
ggd.grid.grid_subset{1}.metric.g21_contravariant = [];
ggd.grid.grid_subset{1}.metric.g22_contravariant = [];
ggd.grid.grid_subset{1}.metric.g23_contravariant = [];
ggd.grid.grid_subset{1}.metric.g31_contravariant = [];
ggd.grid.grid_subset{1}.metric.g32_contravariant = [];
ggd.grid.grid_subset{1}.metric.g33_contravariant = [];

ggd.r{1}.grid_index                     = -999999999;
ggd.r{1}.grid_subset_index              = -999999999;
ggd.r{1}.values                         = reshape(p2d.r.',n1 * n2,1);
ggd.z{1}.grid_index                     = -999999999;
ggd.z{1}.grid_subset_index              = -999999999;
ggd.z{1}.values                         = reshape(p2d.z.',n1 * n2,1);
ggd.theta{1}.grid_index                 = -999999999;
ggd.theta{1}.grid_subset_index          = -999999999;
if isfield(p2d,'theta')
  ggd.theta{1}.values                   = reshape(p2d.theta.',n1 * n2,1);
else
  ggd.theta{1}.values                   = [];
end
ggd.psi{1}.grid_index                   = -999999999;
ggd.psi{1}.grid_subset_index            = -999999999;
ggd.psi{1}.values                       = reshape(p2d.psi.',n1 * n2,1);
ggd.phi{1}.grid_index                   = -999999999;
ggd.phi{1}.grid_subset_index            = -999999999;
if isfield(p2d,'phi')
  ggd.phi{1}.values                     = reshape(p2d.phi.',n1 * n2,1);;
else
  ggd.phi{1}.values                     = [];
end
ggd.j_tor{1}.grid_index                 = -999999999;
ggd.j_tor{1}.grid_subset_index          = -999999999;
if isfield(p2d,'j_tor')
  ggd.j_tor{1}.values                   = reshape(p2d.j_tor.',n1 * n2,1);
else
  ggd.j_tor{1}.values                   = [];
end
ggd.j_parallel{1}.grid_index            = -999999999;
ggd.j_parallel{1}.grid_subset_index     = -999999999;
if isfield(p2d,'j_parallel')
  ggd.j_parallel{1} .values             = reshape(p2d.j_parallel.',n1 * n2,1);
else
  ggd.j_parallel{1}.values              = [];
end
ggd.b_field_r{1}.grid_index             = -999999999;
ggd.b_field_r{1}.grid_subset_index      = -999999999;
ggd.b_field_r{1}.values                 = reshape(p2d.b_field_r.',n1 * n2,1);
ggd.b_field_z{1}.grid_index             = -999999999;
ggd.b_field_z{1}.grid_subset_index      = -999999999;
ggd.b_field_z{1}.values                 = reshape(p2d.b_field_z.',n1 * n2,1);
ggd.b_field_tor{1}.grid_index           = -999999999;
ggd.b_field_tor{1}.grid_subset_index    = -999999999;
ggd.b_field_tor{1}.values               = reshape(p2d.b_field_tor.',n1 * n2,1);

