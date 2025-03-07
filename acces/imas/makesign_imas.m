% appel de la fonction qui regle le bon cocos
function [data_zerod,profil0d] = makesign_imas(data_zerod,profil0d,z0dstruct,sigma_Ip_in,sigma_B0_in)

% transformation
%  Btor		Ip		Psi		Phi		safety factor
%  positive 	positive 	decreasing 	increasing 	negative
%  positive 	negative 	increasing 	increasing 	positive
%  negative 	positive 	decreasing 	decreasing 	positive
%  negative 	negative 	increasing 	decreasing 	negative
%
s_psi = sign(sum(sign(profil0d.psi(:,1) - profil0d.psi(:,end))));
fact_psi =  sigma_Ip_in .* s_psi;
s_phi = sign(sum(sign(profil0d.phi(:,1) - profil0d.phi(:,end))));
s_q   = sign(sum(sign(profil0d.qjli(:))));
fact_q = s_q .* sigma_B0_in .* sigma_Ip_in;
sigma_rhothetaphi_in = sigma_Ip_in .* sigma_B0_in .* s_q;
sigma_RphiZ_in = 1;
sigma_Bp_in = sigma_Ip_in .* s_psi;

% changement normalisation psi
profil0d.psi            = profil0d.psi    .* fact_psi;
profil0d.dpsidt         = profil0d.dpsidt .* fact_psi;
profil0d.qjli           = profil0d.qjli   .* fact_q;  
profil0d.fdia           = profil0d.fdia .* sigma_B0_in;
profil0d.phi            = profil0d.phi .* sigma_B0_in;
profil0d.dphidx         = profil0d.dphidx .* sigma_B0_in;
profil0d.jeff           = profil0d.jeff .* sigma_Ip_in;
profil0d.jli            = profil0d.jli .* sigma_Ip_in;
profil0d.jboot          = profil0d.jboot .* sigma_Ip_in;
profil0d.jnbicd         = profil0d.jnbicd .* sigma_Ip_in;
profil0d.jlh            = profil0d.jlh   .* sigma_Ip_in;
profil0d.jeccd          = profil0d.jeccd .* sigma_Ip_in;
profil0d.jfwcd          = profil0d.jfwcd .* sigma_Ip_in;
profil0d.jfus           = profil0d.jfus .* sigma_Ip_in;
profil0d.jrun           = profil0d.jrun  .* sigma_Ip_in; 
profil0d.jni            = profil0d.jni .* sigma_Ip_in; 
profil0d.jfusshape      = profil0d.jfusshape .* sigma_Ip_in; 
profil0d.epar           = profil0d.epar .* sigma_Ip_in;
profil0d.omega          = profil0d.omega *  sigma_RphiZ_in;
profil0d.utheta         = profil0d.utheta * sigma_rhothetaphi_in  *  sigma_RphiZ_in ./ fact_psi;
profil0d.vtheta         = profil0d.vtheta * sigma_rhothetaphi_in  *  sigma_RphiZ_in;
profil0d.vtor           = profil0d.vtor    *  sigma_RphiZ_in;
profil0d.rot_nbi        = profil0d.rot_nbi *  sigma_RphiZ_in;
profil0d.rot_n0         = profil0d.rot_n0  *  sigma_RphiZ_in;
profil0d.rot_lh         = profil0d.rot_lh  *  sigma_RphiZ_in;

data_zerod.q0           = data_zerod.q0   .* fact_q;  
data_zerod.q95          = data_zerod.q95  .* fact_q; 
data_zerod.qa           = data_zerod.qa   .* fact_q; 
data_zerod.qmin         = data_zerod.qmin .* fact_q; 
data_zerod.phiplasma    = data_zerod.phiplasma .* sigma_B0_in;
data_zerod.ifus         = data_zerod.ifus  .* sigma_Ip_in;
data_zerod.jxfus        = data_zerod.jxfus .* sigma_Ip_in;
data_zerod.j0fus        = data_zerod.j0fus .* sigma_Ip_in;
data_zerod.iohm         = data_zerod.iohm  .* sigma_Ip_in;
data_zerod.vloop        = data_zerod.vloop .* sigma_Ip_in;
data_zerod.ip           = data_zerod.ip    .* sigma_Ip_in;
data_zerod.ifwcd        = data_zerod.ifwcd .* sigma_Ip_in;
data_zerod.ieccd        = data_zerod.ieccd .* sigma_Ip_in;
data_zerod.inbicd       = data_zerod.inbicd.* sigma_Ip_in;
data_zerod.iboot        = data_zerod.iboot .* sigma_Ip_in;
data_zerod.wrad         = data_zerod.wrad  .* sigma_RphiZ_in;
data_zerod.irun         = data_zerod.irun  .* sigma_Ip_in;
data_zerod.sn0fr        = data_zerod.sn0fr .* sigma_RphiZ_in;
data_zerod.ilh          = data_zerod.ilh   .* sigma_Ip_in;
data_zerod.icd          = data_zerod.icd   .* sigma_Ip_in;
data_zerod.qeff         = data_zerod.qeff  .* fact_q; 
data_zerod.vmes         = data_zerod.vmes  .* sigma_Ip_in;
data_zerod.ipar         = data_zerod.ipar  .* sigma_Ip_in;
data_zerod.ini          = data_zerod.ini   .* sigma_Ip_in;
data_zerod.edgeflux     = data_zerod.edgeflux .* sigma_Ip_in * sigma_Bp_in;

