function grid = mapprofiles1d_grid_imas(profil0d,indice,grid)

if nargin < 3
    %warning('Sub structure grid is not initialize');
    grid = [];
end

if nargin == 1
    grid.rho_tor_norm = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
    grid.rho_tor      = profil0d.rmx;
    grid.rho_phi_norm = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
    grid.rho_phi      = profil0d.rmx;
    grid.psi          = profil0d.psi;
    grid.rho_pol_norm = sqrt((profil0d.psi - profil0d.psi(:,1) * ones(1,size(profil0d.psi,2))) ./  ...
        ((profil0d.psi(:,end) - profil0d.psi(:,1))*ones(1,size(profil0d.psi,2))));
    % for COCOS test
    grid.psi_magnetic_axis = grid.psi(1);
    grid.psi_boundary = grid.psi(end);
    grid.volume       = cumtrapz(profil0d.xli,profil0d.vpr,2);
    grid.surface      = profil0d.vpr ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
    grid.area         = cumtrapz(profil0d.xli,profil0d.spr,2);
else
    grid.rho_tor_norm = profil0d.rmx(indice,:) ./ profil0d.rmx(indice,end);
    grid.rho_tor      = profil0d.rmx(indice,:);
    grid.rho_phi_norm = profil0d.rmx(indice,:) ./ profil0d.rmx(indice,end);
    grid.rho_phi      = profil0d.rmx(indice,:);
    grid.psi          = profil0d.psi(indice,:);
    grid.rho_pol_norm = sqrt((profil0d.psi(indice,:) - profil0d.psi(indice,1)) ./  ...
        (profil0d.psi(indice,end) - profil0d.psi(indice,1)));
    % for COCOS test
    grid.psi_magnetic_axis = grid.psi(1);
    grid.psi_boundary = grid.psi(end);
    grid.volume       = cumtrapz(profil0d.xli,profil0d.vpr(indice,:),2);
    grid.surface      = profil0d.vpr(indice,:) ./ profil0d.rmx(indice,end);
    grid.area         = cumtrapz(profil0d.xli,profil0d.spr(indice,:),2);
end