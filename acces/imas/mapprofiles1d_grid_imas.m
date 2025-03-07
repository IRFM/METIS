function grid = mapprofiles1d_grid_imas(profil0d,indice,grid)

if nargin < 3
   warning('Sub structure grid is not initialize');
   grid = [];
end

if nargin == 1
	grid.rho_tor_norm = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
	grid.rho_tor      = profil0d.rmx;
	grid.psi          = profil0d.psi;
	grid.volume       = cumtrapz(profil0d.xli,profil0d.vpr,2);
	grid.area         = cumtrapz(profil0d.xli,profil0d.spr,2);
else
	grid.rho_tor_norm = profil0d.rmx(indice,:) ./ profil0d.rmx(indice,end);
	grid.rho_tor      = profil0d.rmx(indice,:);
	grid.psi          = profil0d.psi(indice,:);
	grid.volume       = cumtrapz(profil0d.xli,profil0d.vpr(indice,:),2);
	grid.area         = cumtrapz(profil0d.xli,profil0d.spr(indice,:),2);
end