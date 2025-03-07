% modele de coefficients de transport stiff of METIS
function [kie,kie_itb,kii,kii_itb] = z0stiff(meff,profil)

% calcul de rho_star (de l'ion majoritaire)
bphi   = profil.fdia .* profil.ri;
ve     = ones(size(profil.xli));
rhomax = profil.rmx(:,end) * ve;


rho_e   = 1.07e-4 .*  sqrt(max(13.6,profil.tep) ./ 1e3) ./ sqrt((profil.fdia .* profil.ri) .^ 2 + profil.bpol .^2);
rho_b   = rho_e .* profil.qjli ./ sqrt(profil.epsi);
rho_pot = (rho_e .^ 2 .* profil.qjli .^ 2 .* profil.Raxe) .^ (1/3);
rho_b   = rho_b .* (rho_b < profil.rmx) + rho_pot .* (rho_b >= profil.rmx);
rho_b(:,1) = rho_pot(:,1);
% debyes
rho_d  = 2.35e5  .* sqrt(max(13.6,profil.tep) ./ max(1e13,profil.nep)./ 1e3);
%figure(51);clf;plot(cons.temps,rho_b,'r',cons.temps,rho_d,'b');drawnow
rho_bd = max(rho_b,rho_d);
% cas du modele stiff2 avec particule passantes et piegees
rho_bd = profil.ftrap .* rho_bd + (1 - profil.ftrap) .* rho_e;
%
frs  = max(0.1,min(10,profil.xieshape_itb ./ max(eps,profil.xieshape)));


% les coefficients de transport
xie  = max(rho_bd(:)) ./ rho_bd;
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


