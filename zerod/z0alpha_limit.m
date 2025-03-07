function  [kie,kie_itb,kii,kii_itb] = z0alpha_limit(A,profil)

phys = cphys;

vt   = ones(size(profil.psi,1),1);
ve   = ones(1,size(profil.psi,2));
shear = pdederive(profil.xli,profil.qjli,1,2,2,1) ./ profil.qjli .* (vt * profil.xli) ;
shear(:,1) = 0;
shear_limited = max(-0.6,shear);
% alpha critique from Freidberg
%alpha_crit = (0.5 + 1.39 .* shear_limited .^ 2) ./ (1 + 0.83 .* shear_limited);
alpha_crit = salpha_lookup(shear_limited);
% for the transport alpha_crit is always limited (to prevent convergence
% problème)
%alpha_crit = min(2,alpha_crit);
% compute pressure gradient (from dP/dr to dP/dx -> drop one r)
dptotdx_crit = -alpha_crit .* profil.bpol .^ 2 ./ (2 .* (4*pi*1e-7) .*  profil.epsi);
%
denom = - cumtrapz(profil.xli, (profil.source_el + profil.source_ion) .* profil.vpr,2);
press_dens = phys.e .* (profil.tep .* pdederive(profil.xli,profil.nep,0,2,2,1) + ...
                        profil.tip .* pdederive(profil.xli,profil.nip,0,2,2,1));
%  some limitation to prevent non convergence                  
beta = min(10,max(0.1,profil.xii ./ max(1e-3,profil.xie)));
xie  = denom ./ profil.vpr ./ profil.grho2 ./ (1 + beta) ./ ...
       min(-1,dptotdx_crit - press_dens);
xie(:,1) = xie(:,2);

% we look only for the shape
xie = xie ./ (mean(xie,2) * ve);
% to be compatible with other METIS mode
frs  = max(0.1,min(10,profil.xieshape_itb ./ max(eps,profil.xieshape)));
% others data
xii  = xie;
xie_itb  = frs .* xie;
xii_itb  = xie_itb;

% valeur centrale
xie(:,1)     = xie(:,2);
xii(:,1)     = xii(:,2);
xie_itb(:,1) = xie_itb(:,2);
xii_itb(:,1) = xii_itb(:,2);

% sorties en K
ne_shape   = profil.nep ./ max(1,max(profil.nep,[],2) * ve);
ni_shape   = profil.nip ./ max(1,max(profil.nip,[],2) * ve);
kie        = profil.xieshape .* xie     .* ne_shape;
kie_itb    = profil.xieshape .* xie_itb .* ne_shape; 
kii        = profil.xieshape .* xii     .* ni_shape;
kii_itb    = profil.xieshape .* xii_itb .* ni_shape;

%figure(21);plot(profil.xli,xie,'b',profil.xli,kie,'r',profil.xli,profil.xieshape,'k');drawnow
%keyboard

