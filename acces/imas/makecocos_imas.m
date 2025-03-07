% appel de la fonction qui regle le bon cocos
function [data_zerod,profil0d,sigma_B0_out,sigma_bvac_r,factor_two_pi] = makecocos_imas(data_zerod,profil0d,z0dstruct,COCOS_out,factor_two_pi)

% par defaut
sigma_B0_eff = 1;
% input
sigma_Ip_in = sign(mean(data_zerod.ip));
sigma_B0_in = sign(mean(profil0d.fdia(:)));

%  Btor		Ip		Psi		Phi		safety factor
%  positive 	positive 	decreasing 	increasing 	negative
%  positive 	negative 	increasing 	increasing 	positive
%  negative 	positive 	decreasing 	decreasing 	positive
%  negative 	negative 	increasing 	decreasing 	negative
% consitency tests
% by choice in METIS ip > 0 and sign of Btor is given by option.signe
s_psi = sign(sum(sign(profil0d.psi(:,1) - profil0d.psi(:,end))));
s_phi = sign(sum(sign(profil0d.phi(:,1) - profil0d.phi(:,end))));
s_q   = sign(sum(sign(profil0d.qjli(:))));
switch z0dstruct.z0dinput.option.signe
case -1
      if sigma_Ip_in > 0
	if s_psi < 0
	  	error('Unconsistent Psi variation : Psi must be decreasing from magnetic axis to LCFS');
	end
	if s_phi < 0
	  	error('Unconsistent Phi variation : Phi must be decreasing from magnetic axis to LCFS');
	end
	if s_q < 0
	  	error('Unconsistent safety factor sign :  safety factor must be positive');
	end
      else
	if s_psi > 0
	  	error('Unconsistent Psi variation : Psi must be decreasing from magnetic axis to LCFS');
	end
	if s_phi > 0
	  	error('Unconsistent Phi variation : Phi must be decreasing from magnetic axis to LCFS');
	end
	if s_q < 0
	  	error('Unconsistent safety factor sign :  safety factor must be positive');
	end
      end
otherwise
      if sigma_Ip_in > 0
	if s_psi < 0
	  	error('Unconsistent Psi variation : Psi must be decreasing from magnetic axis to LCFS');
	end
	if s_phi > 0
	  	error('Unconsistent Phi variation : Phi must be increasing from magnetic axis to LCFS');
	end
	if s_q > 0
	  	error('Unconsistent safety factor sign :  safety factor must be negative');
	end
      else
 	if s_psi > 0
	  	error('Unconsistent Psi variation : Psi must be decreasing from magnetic axis to LCFS');
	end
	if s_phi < 0
	  	error('Unconsistent Phi variation : Phi must be increasing from magnetic axis to LCFS');
	end
	if s_q > 0
	  	error('Unconsistent safety factor sign :  safety factor must be negative');
	end     
      end
end
% CoCos de METIS
sigma_Bp_in = sigma_Ip_in .* sign(mean(sign(profil0d.psi(:,end) - profil0d.psi(:,1))));
s_dpdspi    = sign(sum(sign((profil0d.ptot(:,end) - profil0d.ptot(:,1)) ./ (profil0d.psi(:,end) - profil0d.psi(:,1)))));
if s_dpdspi ~= (-sigma_Bp_in * sigma_Ip_in)
	error('consisency error :sign(dP/dPsi) is not - sigma_Bp_in * sigma_Ip_in');
end
sigma_rhothetaphi_in = sigma_Ip_in .* sigma_B0_in .* s_q;
% we assume coordinate (R,phi,Z) right-handed
if  sigma_Bp_in > 0
    if sigma_rhothetaphi_in > 0
        COCOS_in = 1;
    else
        COCOS_in = 5;
    end
else
    if sigma_rhothetaphi_in > 0
        COCOS_in = 7;
    else
        COCOS_in = 3;
    end
end
COCOS_in = COCOS_in + 10;

% COCOS METIS
try 
    % from Olivier cocos Git
    % cocos_struct.exp_Bp:             0 or 1 depending if psi is already divided by 2pi or not, respectively
    % cocos_struct.sigma_Bp:          +1 or -1, depending if psi is increasing or decreasing with Ip and B0 positive
    % cocos_struct.sigma_RphiZ:       +1 or -1 depending if (R, phi, Z) is right-handed or (R, Z, phi), resp.
    % cocos_struct.sigma_rhothetaphi: +1 or -1 depending if (rho, theta, phi) is right-handed or (rho, phi, theta), resp.
    % cocos_struct.sign_q_pos:        +1 or -1 depending if q is positive or negative with Ip and B0 positive
    % cocos_struct.sign_pprime_pos:   +1 or -1 depending if dp/dpsi is positive or negative with Ip and B0 positive
    %                                 Note that the expected sign(dpsi) for Ip and B0 positive is -sign_pprime_pos
    cocos_struct            = cocos(COCOS_in);
    Kexp_Bp_in              = cocos_struct.exp_Bp;
    Ksigma_Bp_in            = cocos_struct.sigma_Bp;
    Ksigma_RphiZ_in         = cocos_struct.sigma_RphiZ;
    Ksigma_rhothetaphi_in   = cocos_struct.sigma_rhothetaphi;
    Ksign_q_pos_in          = cocos_struct.sign_q_pos;
    Ksign_pprime_pos_in     = cocos_struct.sign_pprime_pos;
catch
    % previous version
    [Kexp_Bp_in,Ksigma_Bp_in,Ksigma_RphiZ_in,Ksigma_rhothetaphi_in,Ksign_q_pos_in,Ksign_pprime_pos_in] = cocos(COCOS_in);
end
% COCOS cible
try 
    % from Olivier cocos Git
    % cocos_struct.exp_Bp:             0 or 1 depending if psi is already divided by 2pi or not, respectively
    % cocos_struct.sigma_Bp:          +1 or -1, depending if psi is increasing or decreasing with Ip and B0 positive
    % cocos_struct.sigma_RphiZ:       +1 or -1 depending if (R, phi, Z) is right-handed or (R, Z, phi), resp.
    % cocos_struct.sigma_rhothetaphi: +1 or -1 depending if (rho, theta, phi) is right-handed or (rho, phi, theta), resp.
    % cocos_struct.sign_q_pos:        +1 or -1 depending if q is positive or negative with Ip and B0 positive
    % cocos_struct.sign_pprime_pos:   +1 or -1 depending if dp/dpsi is positive or negative with Ip and B0 positive
    %                                 Note that the expected sign(dpsi) for Ip and B0 positive is -sign_pprime_pos
    cocos_struct             = cocos(COCOS_out);
    Kexp_Bp_out              = cocos_struct.exp_Bp;
    Ksigma_Bp_out            = cocos_struct.sigma_Bp;
    Ksigma_RphiZ_out         = cocos_struct.sigma_RphiZ;
    Ksigma_rhothetaphi_out   = cocos_struct.sigma_rhothetaphi;
    Ksign_q_pos_out          = cocos_struct.sign_q_pos;
    Ksign_pprime_pos_out     = cocos_struct.sign_pprime_pos;
catch
    [Kexp_Bp_out,Ksigma_Bp_out,Ksigma_RphiZ_out,Ksigma_rhothetaphi_out,Ksign_q_pos_out,Ksign_pprime_pos_out] = cocos(COCOS_out);
end

% verifications
if any(sign(profil0d.qjli(:).* Ksigma_rhothetaphi_in .* sigma_Ip_in .* sigma_B0_in)<= 0)
    fprintf('WARNING: sign(q) is not consistent with COCOS_in= %d\n',COCOS_in)
    fprintf('qedge = %g\n',mean(profil0d.qjli(:,end)))
    fprintf('sig_rhothetaphi*sign(Ip)*sign(B0) = %d * %d *%d = %d\n', ...
            Ksigma_rhothetaphi_in, sigma_Ip_in,sigma_B0_in, ...
            Ksigma_rhothetaphi_in*sigma_Ip_in*sigma_B0_in);
end 
if any(sign(profil0d.fdia(:) .* sigma_B0_in)<= 0)
    fprintf('WARNING: Signs of F and B0 are not consistent\n');
end
if any(sign((profil0d.psi(:,end) - profil0d.psi(:,1)).*Ksigma_Bp_in.*sigma_Ip_in) <= 0)
    if  any(sign(profil0d.psi(:,end) - profil0d.psi(:,1)) <= 0)
      fprintf('WARNING: psi should be decreasing with : sign(Ip)= %d and %d  for COCOS= %d\n',sigma_Ip_in,sigma_Ip_in,COCOS_in);
    else
      fprintf('WARNING: psi should be increasing with : sign(Ip)= %d and %d  for COCOS=%d\n',sigma_Ip_in,sigma_Ip_in,COCOS_in);
    end
end

% for cpo scenario
sigma_bvac_r = sigma_B0_in;

if COCOS_in == COCOS_out 
    sigma_B0_out = sigma_B0_in * sigma_B0_eff;
	return
end

%  
%  Define effective variables: sigma_Ip_eff, sigma_B0_eff, sigma_Bp_eff, exp_Bp_eff as in Appendix C
%  sign(Ip) in output:
%     
KIPsign_out = sigma_Ip_in;
KB0sign_out = sigma_B0_in;
sigma_RphiZ_eff  = Ksigma_RphiZ_out * Ksigma_RphiZ_in;
%sigma_IP_eff = sigma_RphiZ_eff
sigma_Ip_eff = sigma_Ip_in * KIPsign_out;
sigma_Ip_out = sigma_Ip_in * sigma_Ip_eff;
%sigma_B0_eff = Ksigma_RphiZ_in * Ksigma_RphiZ_out
sigma_B0_eff = sigma_B0_in * KB0sign_out;
sigma_B0_out = sigma_B0_in * sigma_B0_eff;
sigma_Bp_eff = Ksigma_Bp_out * Ksigma_Bp_in;
exp_Bp_eff = Kexp_Bp_out - Kexp_Bp_in;
sigma_rhothetaphi_eff  = Ksigma_rhothetaphi_out * Ksigma_rhothetaphi_in;
fact_psi = sigma_Ip_eff * sigma_Bp_eff * (2.*pi) ^ exp_Bp_eff;
fact_q = sigma_Ip_eff * sigma_B0_eff * sigma_rhothetaphi_eff;
factor_two_pi = factor_two_pi .* fact_psi;
% transformation
% changement normalisation psi
profil0d.psi            = profil0d.psi    .* fact_psi;
profil0d.dpsidt         = profil0d.dpsidt .* fact_psi;
profil0d.qjli           = profil0d.qjli   .* fact_q;  
profil0d.fdia           = profil0d.fdia .* sigma_B0_eff;
profil0d.phi            = profil0d.phi .* sigma_B0_eff;
profil0d.dphidx         = profil0d.dphidx .* sigma_B0_eff;
profil0d.jeff           = profil0d.jeff .* sigma_Ip_eff;
profil0d.jboot          = profil0d.jboot .* sigma_Ip_eff;
profil0d.jnbicd         = profil0d.jnbicd .* sigma_Ip_eff;
profil0d.jlh            = profil0d.jlh   .* sigma_Ip_eff;
profil0d.jeccd          = profil0d.jeccd .* sigma_Ip_eff;
profil0d.jfwcd          = profil0d.jfwcd .* sigma_Ip_eff;
profil0d.jfus           = profil0d.jfus .* sigma_Ip_eff;
profil0d.jrun           = profil0d.jrun  .* sigma_Ip_eff; 
profil0d.jni            = profil0d.jni .* sigma_Ip_eff; 
profil0d.jfusshape      = profil0d.jfusshape .* sigma_Ip_eff; 
profil0d.epar           = profil0d.epar .* sigma_Ip_eff;
profil0d.omega          = profil0d.omega *  sigma_RphiZ_eff;
profil0d.utheta         = profil0d.utheta * sigma_rhothetaphi_eff  *  sigma_RphiZ_eff ./ fact_psi;
profil0d.vtheta         = profil0d.vtheta * sigma_rhothetaphi_eff  *  sigma_RphiZ_eff;
profil0d.vtor           = profil0d.vtor    *  sigma_RphiZ_eff;
profil0d.rot_nbi        = profil0d.rot_nbi *  sigma_RphiZ_eff;
profil0d.rot_n0         = profil0d.rot_n0  *  sigma_RphiZ_eff;
profil0d.rot_lh         = profil0d.rot_lh  *  sigma_RphiZ_eff;

data_zerod.q0           = data_zerod.q0   .* fact_q;  
data_zerod.q95          = data_zerod.q95  .* fact_q; 
data_zerod.qa           = data_zerod.qa   .* fact_q; 
data_zerod.qmin         = data_zerod.qmin .* fact_q; 
data_zerod.phiplasma    = data_zerod.phiplasma .* sigma_B0_eff;
data_zerod.ifus         = data_zerod.ifus  .* sigma_Ip_eff;
data_zerod.jxfus        = data_zerod.jxfus .* sigma_Ip_eff;
data_zerod.j0fus        = data_zerod.j0fus .* sigma_Ip_eff;
data_zerod.iohm         = data_zerod.iohm  .* sigma_Ip_eff;
data_zerod.vloop        = data_zerod.vloop .* sigma_Ip_eff;
data_zerod.ip           = data_zerod.ip    .* sigma_Ip_eff;
data_zerod.ifwcd        = data_zerod.ifwcd .* sigma_Ip_eff;
data_zerod.ieccd        = data_zerod.ieccd .* sigma_Ip_eff;
data_zerod.inbicd       = data_zerod.inbicd.* sigma_Ip_eff;
data_zerod.iboot        = data_zerod.iboot .* sigma_Ip_eff;
data_zerod.wrad         = data_zerod.wrad  .* sigma_RphiZ_eff;
data_zerod.irun         = data_zerod.irun  .* sigma_Ip_eff;
data_zerod.sn0fr        = data_zerod.sn0fr .* sigma_RphiZ_eff;
data_zerod.ilh          = data_zerod.ilh   .* sigma_Ip_eff;
data_zerod.icd          = data_zerod.icd   .* sigma_Ip_eff;
data_zerod.qeff         = data_zerod.qeff  .* fact_q; 
data_zerod.vmes         = data_zerod.vmes  .* sigma_Ip_eff;
data_zerod.ipar         = data_zerod.ipar  .* sigma_Ip_eff;
data_zerod.ini          = data_zerod.ini   .* sigma_Ip_eff;
data_zerod.edgeflux     = data_zerod.edgeflux .* sigma_Ip_eff * sigma_Bp_eff;

% for Bphi . grad(phi)
sigma_bvac_r = sigma_B0_eff .* sigma_RphiZ_eff .* sigma_bvac_r;

