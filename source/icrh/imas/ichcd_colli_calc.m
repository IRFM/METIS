%function [Pe,PM,pabs,R,Wperp,Wpar,Tpar_out,P_loss_conv,flag_non_conv,Xi,next] = ichcd_colli_calc(core_profiles,equilibrium,pulse_schedule,itime,ichcd_params)
 function [core_sources] = ichcd_colli_calc_v2(core_profiles,equilibrium,pulse_schedule,itime,ichcd_params)
%#codegen
% Function to derive ICRH power
% persisten,t remove to ensure Matlab automatic parallelisation will be
% working, replace by last
% persistent te_prev pabs_prev Pe_prev PM_prev U_prev R_prev Wperp_prev Wpar_prev Tpar_prev Xi

% Index for minority species   
ind_mino = 4;
   
nions = numel(core_profiles.profiles_1d{1}.ion);
% Find main ions
for jj=1:nions
   nimax(jj) = max(core_profiles.profiles_1d{itime}.ion{jj}.density);
end

% Find main ions
ind_Mions = find(nimax>1e19); 

% Map all core profiles to rho_tor_norm from equilibrium
q = interp1(core_profiles.profiles_1d{itime}.grid.rho_tor_norm,core_profiles.profiles_1d{itime}.q,equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm,'spline');
te = interp1(core_profiles.profiles_1d{itime}.grid.rho_tor_norm,core_profiles.profiles_1d{itime}.electrons.temperature,equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm,'spline');
ne0 = interp1(core_profiles.profiles_1d{itime}.grid.rho_tor_norm,core_profiles.profiles_1d{itime}.electrons.density,equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm,'spline');
ti = interp1(core_profiles.profiles_1d{itime}.grid.rho_tor_norm,core_profiles.profiles_1d{itime}.ion{ind_Mions(1)}.temperature,equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm,'spline');

for ii=1:ind_Mions
	 ni1(ii,:) = interp1(core_profiles.profiles_1d{itime}.grid.rho_tor_norm,core_profiles.profiles_1d{itime}.ion{ind_Mions(ii)}.density,equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm,'spline');
         AM(ii) = core_profiles.profiles_1d{itime}.ion{ind_Mions(ii)}.element{1}.a;
         ZM(ii) = core_profiles.profiles_1d{itime}.ion{ind_Mions(ii)}.element{1}.z_n;
end

% Adds minority to main ions
ni0 = sum(ni1);
ni0=ni0';
ni0 = ni0 + interp1(core_profiles.profiles_1d{itime}.grid.rho_tor_norm,core_profiles.profiles_1d{itime}.ion{ind_mino}.density,equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm,'spline');

ichcd_params.AM = mean(sum(ni1./ni0'.*AM'));
ichcd_params.ZM = mean(sum(ni1./ni0'.*ZM'));

nm0 = interp1(core_profiles.profiles_1d{itime}.grid.rho_tor_norm,core_profiles.profiles_1d{itime}.ion{ind_mino}.density,equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm,'spline');
ichcd_params.Am = core_profiles.profiles_1d{itime}.ion{ind_mino}.element{1}.a;
ichcd_params.Zm = core_profiles.profiles_1d{itime}.ion{ind_mino}.element{1}.z_n;

model.equi.R0 = equilibrium.vacuum_toroidal_field.r0;
model.equi.B0 = equilibrium.vacuum_toroidal_field.b0;
geop.Vp = equilibrium.time_slice{itime}.profiles_1d.dvolume_drho_tor;
model.equi.epsilon = (equilibrium.time_slice{itime}.profiles_1d.r_outboard(end)-equilibrium.time_slice{itime}.profiles_1d.r_inboard(end))/2./model.equi.R0;
model.equi.Raxe = (equilibrium.time_slice{itime}.profiles_1d.r_outboard(:)+equilibrium.time_slice{itime}.profiles_1d.r_inboard(:))/2;
model.equi.xRaxe = (equilibrium.time_slice{itime}.profiles_1d.r_outboard(:)-equilibrium.time_slice{itime}.profiles_1d.r_inboard(:))./(equilibrium.time_slice{itime}.profiles_1d.r_outboard(end)-equilibrium.time_slice{itime}.profiles_1d.r_inboard(end));
model.geom.rhogauss=model.equi.xRaxe;
model.equi.kappa = equilibrium.time_slice{itime}.profiles_1d.elongation(end);
model.rgrid.nrhogauss = numel(equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm);
model.rgrid.nrho = numel(equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm);

nrho = length(model.geom.rhogauss);

if isfield(ichcd_params,'Xi');
    Xi = ichcd_params.Xi;
else
    Xi = 0.02 .* ones(nrho,1);
end


% convergence flag
flag_non_conv = 0;
  
    phys = icrh_phys;
    ipa = ichcd_params;
    R0 = model.equi.R0;
    B0 = model.equi.B0;
    a  = R0*model.equi.epsilon;
    eps1 = model.equi.kappa;
    
    %%% ATTENTION
    % m for minority spicies (heated test ion)
    % M for majority species
    tie0 = ti*phys.e; % while ti [eV]
    tee0 = te*phys.e;
    tim0 = tie0; tiM0 = tie0; % Suppose tim = tiM = ti; temperatures are equal
    
    nM0 = ni0-nm0;
    mi = ipa.Am*phys.mp;
    % Given heating Power
    P_RF = pulse_schedule.ic.antenna{1}.power.reference.data; % U(ipa.uindices);
    P_loss = ipa.Ploss;
    
    %% Cal Process %%%%%%%%%%%%%%%%%%%%%
    %%%%Thermal velocity % sqrt(2*T/m)
    vthi0 =  sqrt(2*tim0/(phys.mp*ipa.Am)); % ion ti(eV) for ion+1
    vthiM0 = sqrt(2*tiM0/(phys.mp*ipa.AM));
    vthe0 =  sqrt(2*tee0/phys.me);
    % figure(2)
    %     plot([vthi0, vthiM0, vthe0/60])
    %     legend ('vthi','vthiM','vthe/60')
    
    lnA_ii = 30 - log(ipa.ZM^3*sqrt(ni0)./ti.^(3/2));
    G_iM = (nM0*ipa.Zm^2*ipa.ZM^2.*(phys.e)^4.*lnA_ii)./(4*pi*phys.eps0^2*mi.^2); % 4.68
    nu_i = G_iM./vthi0.^3; % 4.63
    
    %%%%%%%%%%%%%%%%%%%%% cal pabs
    % Using g_ geometry to calculate pabs!!!
    Nz = ipa.Nz;  % point nb = 2Nz
    if isfield(model.equi,'Raxe')% to take into account Shafranov shift
        xnz  = (-1:1/(Nz-1):1);
        Raxe = interp1(model.equi.xRaxe,model.equi.Raxe,abs(xnz),'linear','extrap')';
    end
    if isfield(model.equi,'Raxe')
        r_abs = Raxe+(-1:1/(Nz-1):1)'*a;
    else
        r_abs = R0+(-1:1/(Nz-1):1)'*a;
    end
    
    %%%% Determine R, the resonance position
    btor = B0*R0./r_abs;
    % bpol = B0./R0 .*(1./[flipud(q); q(2:end)]);  % JF ARTAUD
    if isfield(model.equi,'Raxe')
        grho = a ./ gradient(r_abs,1/(Nz-1)); 
        bpol = B0*a./r_abs.*interp1([flipud(model.equi.Raxe);model.equi.Raxe(2:end)]+a*(-1:(1/(nrho-1)):1)',((-1:(1/(nrho-1)):1)'./[flipud(q); q(2:end)]),r_abs,'linear') .* grho; % TRANG
    else
        bpol = btor*a./r_abs.*interp1(R0+a*(-1:(1/(nrho-1)):1)',((-1:(1/(nrho-1)):1)'./[flipud(q); q(2:end)]),r_abs,'linear'); % TRANG
    end
    %figure(21);plot(xnz,bpol);drawnow
    Bt   = sqrt( btor .^2  + bpol.^2);
    wci = ipa.Zm * phys.e./ mi .* Bt;
    w = 2*pi*ipa.f;
    % figure; plot(r_abs,[wci, ones(2*Nz-1,1)*w])
    
    % Poloidal upshift effect
    if isfield(model.equi,'Raxe')
      vthi_Nz = interp1(model.geom.rhogauss,vthi0,(r_abs(Nz:end)-Raxe(Nz:end))/a,'linear');
    else
      vthi_Nz = interp1(model.geom.rhogauss,vthi0,(r_abs(Nz:end)-R0)/a,'linear');
    end
    
    if ipa.upshift ==1
        res = ipa.p*wci + ipa.n./r_abs.*[flipud(vthi_Nz(2:end)); vthi_Nz(1:end)]; % (C.4)
    else
        res = ipa.p*wci;
    end

    % Detect bad icrh frequency w
    if w > max(res) || w < min(res)
        fprintf('Cant find the Resonance location. Need to change the ICHCD frequency. \n')
        % Fill core_sources with zeros
        core_sources = idsfill_core_sources_icrh(zeros(length(model.geom.rhogauss),1),zeros(length(model.geom.rhogauss),1),zeros(length(model.geom.rhogauss),1),model.geom.rhogauss,geop.Vp,model);
        return
    end
    
    tmp= abs(res- w); [~,idxabs] = min(tmp); % while w == res
    R = r_abs(idxabs); % Resonance position
    R_right = (idxabs>Nz); 
    idx = max(1,abs(idxabs-(Nz-1))); 
    delR = ipa.n/w*vthi_Nz(idx); % C.5)
    delZ = 2*a*eps1*ipa.fz; % C.6
 
    %%%% R in rho_METIS coordinate
    if isfield(model.equi,'Raxe')
        tmp = abs(abs(R-model.equi.Raxe)/a - model.geom.rhogauss);[~,Ridx_rho] = min(tmp);
    else
        tmp = abs(abs(R-R0)/a - model.geom.rhogauss);[~,Ridx_rho] = min(tmp);
    end
    % % Average!!!!
    % Vint = 2*pi*R*delR*delZ; % C.7
    % pabs_mean = (P_RF - P_loss)/Vint; % C.8 %%% pabs at rho = R
   
    %%%%%%%%%%% Profile of pabs with 2 gaussians
    z= (-1:1/(Nz-1):1)'*a*eps1;
    if isfield(model.equi,'Raxe')
        r = Raxe+(-1:1/(Nz-1):1)'*a;
    else
        r = R0+(-1:1/(Nz-1):1)'*a;
    end
    %%% profile for each axis
    Pz = exp(-(z).^2/(delZ^2/4)); % 4sigma = delZ central resonance (1/2 Gaussian)
    Pr = exp(-(r-R).^2/(delR^2)); % 2sigma = delR
    
    Pmean = (P_RF - P_loss)/(pi^2*R*delR*delZ);% = =A in Remi note
    P2D = Pmean*Pz*Pr';%*sqrt(2) if (delR^2/2)
    P2D_int = 2 .* pi .* P2D .* (r * ones(size(z'))) .* mean(diff(r)) .* mean(diff(z));

    %%%% Test
    % P_RF = 2*pi*trapz(r,r.*trapz(z,P2D,1)');
    %[rr,zz] = meshgrid(r,z);
    %mesh(rr,zz,P2D)
    
    %%%%%% Profile pabs 1D [Wm^-3]
    % from note: P1D ~ P2D(rho) or int(P2D,dr,dz) == int(P1D,rho)
    % In pratice, we base on the eq.15: P1D = dP/dV
    % integrate from the center R=0 <=> rho=0
    Prho = zeros(Nz,1); %%% P(rho)
    %Prho_alt = zeros(Nz,1); %%% P(rho)
    for i = 1:Nz-1
        intnum = Nz-i:Nz+i;
        %Prho(i+1) = 2*pi*trapz(r(intnum),r(intnum).*trapz(z(intnum),P2D(intnum,intnum),1)');
        Prho(i+1) = sum(sum( P2D_int(intnum,intnum)));
    end
    
    %%% P_RF == Prho(end)
    rhoNz = linspace(0,1,Nz)';
    Vrho  = 2*pi^2*R*a^2*eps1*rhoNz.^2;
    dVrho = 4*pi^2*R*a^2*eps1*rhoNz;
    P1D = diff(Prho)./diff(Vrho);
    P1D = interp1(rhoNz(1:end-1),P1D,rhoNz,'linear',0);
    % normalisattion (correction imprecision numerique)
    P1D = P1D .*  (P_RF - P_loss) ./ max(eps,trapz(rhoNz,dVrho .* P1D));
    pabs = interp1(rhoNz,P1D,model.geom.rhogauss,'linear','extrap');
    pabs(pabs < eps) = eps;
    
    % figure(1); plot(model.geom.rhogauss,[ipa.picrh pabs],rhoNz,P1D,'Linewidth',2)
    % legend('pabs_{METIS}','pabs_{interp}','pabs_{gauss}')
    % ylabel('Absorbed power density')
    % xlabel('\rho')
    
    % % P_RF - 4*pi^2*R*a^2*eps1*trapz(model.geom.rhogauss,model.geom.rhogauss.*pabs)
    % % figure; plot(model.geom.rhogauss,[pabs, ipa.picrh])
    
    %% Cal from 4.65-4.71 with the species are electron and i_M (except i_n which is heated)
    if isfield(model.equi,'Raxe')
      rho = model.equi.Raxe-(-1)^R_right*model.geom.rhogauss*a; % normalize a
    else
      rho = R0-(-1)^R_right*model.geom.rhogauss*a; % normalize a
    end
    Bt0 = interp1(r_abs,Bt,rho,'linear','extrap');
    %%%%%%%%%%%
    kpar = ipa.n/(R);
    %     kpar = w/vthe(idx); % from Stix
    npar = kpar * phys.c/w;
        %%%%% NEW cal from Stix62 %%%%%%%%%%%%%%%%%%%%%%%
    % %     Pulsation (eq III.3)
    wpm = sqrt(nm0.*ipa.Zm^2*phys.e^2/(phys.eps0*ipa.Am*phys.mp));
    wpM = sqrt(nM0.*ipa.ZM*phys.e^2/(phys.eps0*ipa.AM*phys.mp));
    %     wpi2 = (wpm.^2+wpM.^2);
    wpe = sqrt(ne0.*phys.e^2/(phys.eps0*phys.me));
    
    % Cyclotronique pulsation (eq III.4)
    wci = Bt0*ipa.Zm*phys.e/(ipa.Am*phys.mp);
    wcM = Bt0*ipa.ZM*phys.e/(ipa.AM*phys.mp);
    wce = -Bt0*phys.e/(phys.me);
    
    % Dielectric tensor elements
    e_perp = 1- (wpm.^2./(w.^2-wci.^2)+ wpM.^2./(w.^2-wcM.^2)+ wpe.^2./(w.^2-wce.^2));
    e_xy = wpm.^2.*wci./w./(w.^2-wci.^2) + wpM.^2.*wcM./w./(w.^2-wcM.^2)+ wpe.^2.*wce./w./(w.^2-wce.^2);
    
    RS = e_perp + e_xy;
    LS = e_perp - e_xy;
    SS = 1/2*(RS+LS); % S in Stix
    
    k_bot2 = -w^2/phys.c^2*(npar^2-RS).*(npar^2-LS)./(npar^2-SS); % (C.10)
    k_bot = sqrt(abs(k_bot2)); %%%
    
    
    % computation of transmission factor taking into account conversion mode losses
    % as computed in PION [L.G. Eriksson and T. Hellsten, Physica Scripta vol. 52, 70-79, 1995]
    % using Budden formula [ref 24]
    % position of the cut off
    cof = npar^2-LS;
    [cof_min,ind_co_min] = min(cof);
    [cof_max,ind_co_max] = max(cof);
    
    if sign(cof_min*cof_max) < 0
      R_cut_off = (rho(ind_co_max) .* cof_max - rho(ind_co_min) .* cof_min) ./ (cof_max - cof_min);
      %figure(21);plot(rho,(npar^2-RS),'b',rho,(npar^2-LS),'r',R,0,'*k',R_cut_off,0,'ok',rho,1000 .* pabs ./ max(pabs),'g');drawnow
      % distance between cut off and resonnance
      if ipa.tau_cut_off>0
	tau_cut_off = ipa.tau_cut_off;
      else	
        tau_cut_off  = R_cut_off - R;
      end
      
      % transmission coefficient of energy
      Transmission = exp(-pi .* abs(k_bot(Ridx_rho) * tau_cut_off));
      % there to part before and after the cut off
      ftrans = Transmission .* (r_abs < R_cut_off) + (r_abs >= R_cut_off);
      
      %figure(21);plot(r_abs,ftrans,'b',r_abs,Pr,'r');drawnow
      P_loss_conv = (P_RF - P_loss) .* trapz(r_abs,(1 - ftrans ) .* Pr) ./ max(eps,trapz(r_abs,Pr));
      % Ploss is removed externally in METIS
      %pabs = interp1(rhoNz,P1D_conv,model.geom.rhogauss,'linear','extrap');
      %pabs(pabs < eps) = eps;

    else
      P_loss_conv = 0;   
    end
     
    %%%%% u grid for integral
    %%%%% From here, all quantities in u dimension
    % u0 = (1e-3:1/ipa.Nu:1)*phys.c/vthi0(idx); % Uniform from 0 to inf
    % u0 = logspace(0,log10(phys.c/vthi0(idx)),ipa.Nu); % from 0 to inf
    u0 = logspace(-3,log10(phys.c/vthi0(Ridx_rho)),ipa.Nu); % from 0 to inf    
    %%%%%% repmat for all variables
    u = repmat(u0,nrho,1); % same u for each value of v_thermal
    vthi = repmat(vthi0,1,ipa.Nu);
    vthiM = repmat(vthiM0,1,ipa.Nu);
    vthe = repmat(vthe0,1,ipa.Nu);
    tiM = repmat(tiM0,1,ipa.Nu);
    tim = repmat(tim0,1,ipa.Nu);
    tee = repmat(tee0,1,ipa.Nu);
    nM = repmat(nM0,1,ipa.Nu);
    nm = repmat(nm0,1,ipa.Nu);
    ne = repmat(ne0,1,ipa.Nu);
    pabs_u = repmat(pabs,1,ipa.Nu);
    Bt =  repmat(Bt0,1,ipa.Nu);
    k_perp = repmat(k_bot,1,ipa.Nu);
    wci = repmat(wci,1,ipa.Nu);
    %%%%%%%% cal Dw
    up = u;
    v_bot = u.*vthi; %% v absolute
    
    %     if ipa.p == 1
    bess0 = (besselj(ipa.p-1,k_perp.*v_bot./wci));
    bess2 = (besselj(ipa.p+1,k_perp.*v_bot./wci));
    EE = abs((npar^2-LS)./(npar^2-RS));% (C.11)
    %EE = EE*0;  % Add EE => Dw is not "good"
    
    int_Dp = up.^3.*abs(bess0 + repmat(EE,1,ipa.Nu).*bess2).^2.*exp(-up.^2); % C.14
    res_u = repmat((trapz(u0,int_Dp'))',1,ipa.Nu);
    Dp = pabs_u/8./(nm.*tim)./res_u;
    
    %EE = EE*0;  % Add EE => Dw is not "good"    
    %int_Dp = up.^3.*abs(bess0 + repmat(EE,1,ipa.Nu).*bess2).^2.*exp(-up.^2); % C.14
    %res_u = repmat((trapz(u0,int_Dp'))',1,ipa.Nu);
    %Dp_EE0 = pabs_u/8./(nm.*tim)./res_u;    
    %figure(21);plot(rho,Dp,rho,Dp_EE0);drawnow
    
    %     Dp = pabs_u/4./(nm.*tim); %~ Dp Quite good! (C.15)
    Dw0 = Dp.*abs(bess0 + repmat(EE,1,ipa.Nu).*bess2).^2;
%     elseif ipa.p == 2
%         % error('not yet implemented');        
%     end
    
    %figure(41);plot(u0,Dw0);drawnow
    %mesh(u0,rho-a,Dw)
    %find(u0==1)
    % xlabel('index u')
    % ylabel('index rho')
    % title('Dw')
    % plot(u0,Dw(idx,:))
    
    %% Cal alpha, beta, gamma from Demont's note
    %%% cal u_beta (4.72)
    u_iM = vthi./vthiM.*u;  % for Majority ions
    u_e = vthi./vthe.*u;    % for electron
    
    % %% cal vi/beta_vi (4.71)
    viM_vi = 1; % for beta species ==M
    vie_vi = ne./(nM*ipa.ZM^2); %for electron
    
    Fu_M = - viM_vi.*tim./tiM.*Psi_func(u_iM); % 4.58
    Fu_e = - vie_vi.*tim./tee.*Psi_func(u_e);
    % Fu = Fu_M + Fu_e;
    
    Duu_M = 1./2./u.*(-Fu_M.*tiM./tim);  % 4.57
    Duu_e = 1./2./u.*(-Fu_e.*tee./tim);
    Duu = Duu_M + Duu_e;
    
    gam_M = 1./(2*u).*viM_vi.*Theta(u_iM); % 4.59
    gam_e = 1./(2*u).*vie_vi.*Theta(u_e);
    gam_c = gam_M + gam_e;
    
    %%%%% Cal alpha, beta, gamma (C.31)
    du_Duu_M  = -1./u.*Duu_M + 1/2./u.*viM_vi.*(vthi./vthiM).*(-2./u_iM.*Psi_func(u_iM)+ 4*exp(-u_iM.^2)/pi^0.5);
    du_Duu_e  = -1./u.*Duu_e + 1/2./u.*vie_vi.*(vthi./vthe) .*(-2./u_e .*Psi_func(u_e) + 4*exp(-u_e.^2)/pi^0.5);
    du_Duu = du_Duu_M + du_Duu_e;
    
    % du_test =  diff(Duu')./diff(u');
    % du_test1 = [du_test',du_test(end,:)'];
    
    alp_M = Fu_M + 1./u.^2.* (2*u.*Duu_M + u.^2.*du_Duu_M); % du(u^2*Duu); % D.31
    % ~Fu_M
    
    alp_e = Fu_e +1./u.^2.* (2*u.*Duu_e + u.^2.*du_Duu_e );
    alp = alp_M + alp_e;
    %%%% ~~ Fu_e
    
    beta = 2*Duu; beta_M = 2*Duu_M; beta_e = 2*Duu_e;
    gam = 4*gam_c; gam_M = 4*gam_M; gam_e = 4*gam_e; % problem 4
    % gam = gam_c; gam_M = gam_M; gam_e = gam_e;
        
    %%%%%%%%% Cal Fbot C.33 the distribution func
    % suppose u == u_bot = up (perpendicular)
    
    %% 2016 Sept 26: Try iteration loop for coefficient \Xi = Plin /pabs; or Dw_up == \Xi* Dw
    %%%%  Computing Fbot
    count = 0; % to set max iterative loop number!!!
    pabs_loop = pabs*1;
    differ = ipa.eps_Xi*max(pabs)+1;
    differ_mem = differ;
    sign_differ_flag = 1;
    count_change = NaN;
    count_max = (20 + 80 .* (nargin < 10));
    %%
    while   (differ >= ipa.eps_Xi*max(pabs) || differ <= -ipa.eps_Xi*max(pabs)) && (count < count_max)
        if sign_differ_flag  && (count  < (count_max / 2))
             Xi = pabs ./ pabs_loop .* Xi;    
        elseif (count  < (count_max / 2))
             Xi = sqrt(pabs./pabs_loop) .* Xi;
        else
             Xi = (0.7 + 0.3 .* sqrt(pabs./pabs_loop)) .* Xi;
        end
        Dw_up = repmat(Xi,1,ipa.Nu).*Dw0;
        fbot = ( -4.*alp.*up + gam + 2*beta + 2*up*2.*du_Duu)./(2*up.*beta + 4*up.*Dw_up);
        %Fbot_norm_old = exp(-int_u(u,fbot)); % C.33 without CF
        Fbot_norm = exp(-cumtrapz(u(1,:),fbot,2)); % C.33 without CF
        %figure(110);plot(u0,Fbot_norm - Fbot_norm_old); title('Fbotnorm')
        %set(gca,'ylim',[1e-10,10])
        %drawnow
        %    figure(110);semilogy(u0,Fbot_norm); title('Fbotnorm')
        %    set(gca,'ylim',[1e-10,10])
       
        CF = nm0./(2*pi.*vthi0.^3.*(trapz(u0,up.*Fbot_norm,2))); % CF from density
        
        %     du_dF = diff(Fbot')./diff(up');
        du_dF = diff(Fbot_norm')./diff(up');
        duFbot = [du_dF',du_dF(end,:)'];
        pabs_int = up.^2.*Dw_up.*duFbot;
        pabs_loop = -4*pi*nu_i.*vthi0.^3.*tim0.*CF.*(trapz(u0,pabs_int,2));
        % pabs_loop == plin
        
        %figure(500);plot([pabs pabs_loop]); legend('plin','ploop');drawnow
        %     differ = abs((mean(pabs-pabs_loop)));
        % differ = max(pabs_loop)- max(pabs);
        differ = pabs_loop - pabs;
        indmax_differ = find(abs(differ) == max(abs(differ)),1);
        differ = differ(indmax_differ);
        if count > 5
            if sign(differ_mem) ~= sign(differ)
                sign_differ_flag = 0;
                count_change = count;
            end
        end
        differ_mem = differ;
        count = count +1;
    end
    if count >= count_max
        % non convergence
        flag_non_conv = 1;
	%figure(500);plot([pabs pabs_loop]); legend('plin','ploop');drawnow
	%disp([count,count_change]);
        %keyboard
    end
    
    %%%%%%% Cal power profiles
    % P_toInt = up.*Fbot_norm.*(alp.*up -gam/4 + beta/2); % OLD
    %P_toInt = up.*Fbot_norm.*(alp.*up -gam/4 - beta/2 - up.*du_Duu) -1/2*up.^2.*beta.*duFbot; % No CF
    %Pcoll   = 4*pi*nu_i.*vthi(:,Ridx_rho).^3.*tim(:,idx).*CF.*(trapz(u0,abs(P_toInt')))'; % C.39
    %figure(22);plot([Pcoll pabs]);drawnow
    %%% recomputing CF for Pcoll == pabs
    % CF = pabs./Pcoll_norm;
    % Fbot = repmat(CF,1,ipa.Nu).*Fbot_norm; % C.33 CF == B0
    
    %%% Cal pcoll_e and pcoll_i
    PM_toInt = up.*Fbot_norm.*(alp_M.*up-gam_M/4 - beta_M/2-up.*du_Duu_M)-1/2*up.^2.*beta_M.*duFbot;
    Pe_toInt = up.*Fbot_norm.*(alp_e.*up-gam_e/4 - beta_e/2-up.*du_Duu_e)-1/2*up.^2.*beta_e.*duFbot;
    
    PM = 4*pi*nu_i.*vthi(:,Ridx_rho).^3.*tim(:,Ridx_rho) .*CF.*(trapz(u0,abs(PM_toInt')))';
    Pe = 4*pi*nu_i.*vthi(:,Ridx_rho).^3.*tim(:,Ridx_rho).*CF.*(trapz(u0,abs(Pe_toInt')))';
    
    %     Pcoll = PM + Pe;
    % jic = (Pm+PM-Pe)*ipa.tau_jic;
    
    %    figure(2);plot((rho-R0)/a,[Pcoll Pe PM],(rho-R0)/a,pabs,'--','LineWidth',2)
    %    legend('Pcoll','P_e','P_M','pabs_{qlin}')
    %    title(['P_R_F =',num2str(P_RF/1e6),' [MW]'])
    
    %%%% Useful params
    %%%% only 2/3 of energy is perpendicular
    tmp =  (2/3) .* 1/2.*mi.*(up.*vthi).^2.*up.*Fbot_norm;
    Wperp = 2*pi.*vthi(:,Ridx_rho).^3.*CF.*(trapz(u0,abs(tmp')))'; % D.34
    
    v_gam = sqrt(2/mi*14.810*te*phys.e.*(2*ipa.Am^0.5./ne0.*ipa.ZM^2.*nM0).^2/3);
    Tpar = mi/4./(up.*vthi).*repmat(v_gam.^3,1,ipa.Nu); % D.36
    tmp = up.*Tpar.*Fbot_norm .* ((up.*vthi) > repmat(v_gam .* 4 ^(-1/3),1,ipa.Nu));
    Wpar = 2*pi.*vthi(:,Ridx_rho).^3.*CF.*(trapz(u0,abs(tmp')))';  % D.37
    Tpar_out = trapz(u0,abs(tmp'))' ./ trapz(u0,abs(up.*Fbot_norm)')' ./ phys.e;
    
    %% compute distribution function for thermal (Pabs =0). Thermal contribution must be removed in power sources and energy contents
    fbot_th = ( -4.*alp.*up + gam + 2*beta + 2*up*2.*du_Duu)./(2*up.*beta + 4*up.*Dw_up .* 0);
    Fbot_norm_th = exp(-cumtrapz(u(1,:),fbot_th,2)); % C.33 without CF
    CF_th = nm0./(2*pi.*vthi0.^3.*(trapz(u0,up.*Fbot_norm_th,2))); % CF from density
    du_dF_th = diff(Fbot_norm_th')./diff(up');
    duFbot_th = [du_dF_th',du_dF_th(end,:)'];
    %P_toInt_th = up.*Fbot_norm_th.*(alp.*up -gam/4 - beta/2 - up.*du_Duu) -1/2*up.^2.*beta.*duFbot_th; % No CF
    %Pcoll_th   = 4*pi*nu_i.*vthi(:,Ridx_rho).^3.*tim(:,idx).*CF_th.*(trapz(u0,abs(P_toInt')))'; % C.39
   
    %%% Cal pcoll_e and pcoll_i
    PM_toInt_th = up.*Fbot_norm_th.*(alp_M.*up-gam_M/4 - beta_M/2-up.*du_Duu_M)-1/2*up.^2.*beta_M.*duFbot_th;
    Pe_toInt_th = up.*Fbot_norm_th.*(alp_e.*up-gam_e/4 - beta_e/2-up.*du_Duu_e)-1/2*up.^2.*beta_e.*duFbot_th;
    
    PM_th = 4*pi*nu_i.*vthi(:,Ridx_rho).^3.*tim(:,Ridx_rho) .*CF_th.*(trapz(u0,abs(PM_toInt_th')))';
    Pe_th = 4*pi*nu_i.*vthi(:,Ridx_rho).^3.*tim(:,Ridx_rho).*CF_th.*(trapz(u0,abs(Pe_toInt_th')))';
    
    %%%% Useful params
    %%%% only 2/3 of energy is perpendicular
    tmp =  (2/3) .* 1/2.*mi.*(up.*vthi).^2.*up.*Fbot_norm_th;
    Wperp_th = 2*pi.*vthi(:,Ridx_rho).^3.*CF_th.*(trapz(u0,abs(tmp')))'; % D.34
    
    %tmp = up.*(tim0 * ones(1,size(Tpar,2))).*Fbot_norm_th;
    tmp = up.*Tpar.*Fbot_norm_th .* ((up.*vthi) > repmat(v_gam .* 4 ^(-1/3),1,ipa.Nu));
    Wpar_th = 2*pi.*vthi(:,Ridx_rho).^3.*CF_th.*(trapz(u0,abs(tmp')))';  % D.37
					     
    % substraction of thermal contribution
    % prevent small negative value due to numerical limited precision
    PM    = max(0,PM - PM_th);
    Pe    = max(0,Pe - Pe_th);
    Wperp = max(0,Wperp - Wperp_th);
    Wpar  = max(0,Wpar - Wpar_th);

    PM = PM + max(0,pabs - Pe - PM); 
    pel_st = trapz(model.equi.xRaxe,geop.Vp .* Pe);
    pion   = trapz(model.equi.xRaxe,geop.Vp .* PM);
    pabs_int  = trapz(model.equi.xRaxe,geop.Vp .* pabs);

    pel   = pel_st;
    pth   = max(0,pion + pel);

    pconv = min(pion,P_loss_conv);
    P_loss_conv = P_loss_conv - pconv;
    pel  = pel  + pconv;
    pion = pion - pconv;

    % normalisation of output
    pe_int   = trapz(model.equi.xRaxe,geop.Vp .* Pe);
    pM_int   = trapz(model.equi.xRaxe,geop.Vp .* PM);
    Pe = Pe .* ((pel  ./ max(1,pe_int)));
    PM = PM .* ((pion ./ max(1,pM_int)));

    jic = Pe*ipa.tau_jic;
    core_sources = idsfill_core_sources_icrh(jic,Pe,PM,model.geom.rhogauss,geop.Vp,model);


function phys = icrh_phys
% Constants
phys.c     = 2.99792458e8;        % vitesse de la lumiere dans le vide (m/s)
phys.me    =   9.10938188e-31;    % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp    =   1.6726485e-27;     % masse au repos du proton (kg)
phys.e     = 1.602176462e-19;
phys.eps0  = 8.854187817e-12;
phys.mu    = 4*pi*1e-7;
return

function out = Psi_func(x) % \Psi (4.69)
erfprime = 2/sqrt(pi).*exp(-x.^2);
out = (erf(x) - x.*erfprime)./x.^2;
return

function out = Theta(x) % (4.70)
erfprime = 2/sqrt(pi).*exp(-x.^2);
out = 1./x.^2.*((x.^2 -1/2).*erf(x) + x/2.*erfprime);
return

function out = int_u(u,fu)
out = zeros(size(u));
x = (u(1,:));
for i = 2: length(u)
    out(:,i) = trapz(x(1:i),(fu(:,1:i)),2);
end
return

function core_sources = idsfill_core_sources_icrh(jic,Pe,PM,rho,volume,model)

%% ------------------------------------
%% FILL THE CORE_SOURCES ONLLY FOR ICRH
%% ------------------------------------

%% Generation
core_sources=ids_gen('core_sources');

%% IDS properties
core_sources.ids_properties.comment = 'IDS generated from RAPTOR output';
core_sources.ids_properties.homogeneous_time = 1;

%% Vacuum toroidal field
core_sources.vacuum_toroidal_field.b0 = model.equi.B0;
core_sources.vacuum_toroidal_field.r0 = model.equi.R0;

% Fill identifier first (this list can be expanded in the future)
nsource = 1;

core_sources.source{1}.identifier.name = 'ic';
core_sources.source{1}.identifier.index = 5;

for kk = 1:nsource
  core_sources.source{kk}.profiles_1d{1}.grid.rho_tor_norm         = rho;
  core_sources.source{kk}.profiles_1d{1}.grid.volume               = volume;
end

core_sources.source{1}.profiles_1d{1}.electrons.energy         = Pe;
core_sources.source{1}.profiles_1d{1}.total_ion_energy         = PM;
core_sources.source{1}.profiles_1d{1}.j_parallel               = jic;

return
