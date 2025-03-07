% modele de coefficients de transport Bohm/Gyro Bohm (JET/JETTO)
function [kie,kie_itb,kii,kii_itb] = z0bgbs(meff,profil)

% calcul de rho_star (de l'ion majoritaire)
bphi   = profil.fdia .* profil.ri;
ve     = ones(size(profil.xli));
rhomax = profil.rmx(:,end) * ve;
if meff > 3
    rho_star = 1.02e-4 .* sqrt((meff * ve) .* profil.tep) ./ (2 .* rhomax .* bphi);
else
    rho_star = 1.02e-4 .* sqrt((meff * ve) .* profil.tep) ./ (1 .* rhomax .* bphi);
end
ghe = 2;	% try to be in line with JETTO 'old' Bohm-gyroBohm
ghi = 4;	% try to be in line with JETTO 'old' Bohm-gyroBohm
%
frs  = max(0.1,min(10,profil.xieshape_itb ./ max(eps,profil.xieshape)));
% Xie bohm
ee   = 0.1602176462e-18;
pel    = ee .* profil.tep .* profil.nep;
ped1   = pdederive(profil.xli,pel,0,2,2,1);
terme1 = abs(ped1) ./ pel;
xib  = profil.tep ./ bphi .* terme1 .* profil.qjli .^ 2 .* ghi;
xeb  = profil.tep ./ bphi .* terme1 .* profil.qjli .^ 2 .* ghe;
xib_itb  =  frs .* xib;
xeb_itb  =  frs .* xeb;

% Xie gyro bohm
ted1    = abs(pdederive(profil.xli,profil.tep,0,2,2,1));
xigb    = rho_star  ./ bphi .* ted1 .* (ted1 > 0);

% les coefficients de transport
xie  = 7e-2 .* xigb + 8e-5 .* xeb;
xii  = 0.25 .* 7e-2 .* xigb + 2 .* 8e-5 .* xib;
xie_itb  = 7e-2 .* xigb + 8e-5 .* xeb_itb ;
xii_itb  = 0.25 .* 7e-2 .* xigb + 2 .* 8e-5 .* xib_itb ;

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

