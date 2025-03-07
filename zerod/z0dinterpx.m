% interpolation of the radial  coordinate of profiles data of METIS
% nbpx = number of radial points
% mode = 'psi'= linear in normalised poloidal flux, 'sqrt(psi)' =  linear in square root of normalised poloidal flux, 'standard' = linear on Lao radial coordinate (xli)
function profil0d_out = z0dinterpx(profil0d_in,nbpx,mode)

% list of fields
noms = fieldnames(profil0d_in);
noms(strmatch('temps',noms)) = [];
noms(strmatch('xli',noms)) = [];
noms(strmatch('Rsepa',noms)) = [];
noms(strmatch('Zsepa',noms)) = [];

% reservation memoire
profil0d_out.temps = profil0d_in.temps;
for k=1:length(noms)
  profil0d_out.(noms{k}) = NaN .* ones(length(profil0d_in.temps),nbpx);
end
switch mode
case {'psi','sqrt(psi)'}
    profil0d_out.xli = NaN .* ones(length(profil0d_in.temps),nbpx);
end

% boucle sur les temps
for k=1:length(profil0d_in.temps)
      % choix des coordonnees
      switch mode
      case 'psi'
          x_in  = profil0d_in.psi(k,:);
          x_in  = x_in - x_in(end);
          x_in  = 1 - (x_in - x_in(end)) ./ (x_in(1) - x_in(end));
          x_out = linspace(0,1,nbpx);
          profil0d_out.xli(k,:) =  interp1(x_in,profil0d_in.xli,x_out,'pchip','extrap');

      case 'sqrt(psi)'
          x_in  = profil0d_in.psi(k,:);
          x_in  = x_in - x_in(end);
          x_in  = sqrt(1 - (x_in - x_in(end)) ./ (x_in(1) - x_in(end)));
          x_out = linspace(0,1,nbpx);
          profil0d_out.xli(k,:) =  interp1(x_in,profil0d_in.xli,x_out,'pchip','extrap');

      otherwise
          x_in  = profil0d_in.xli;
          x_out = linspace(0,1,nbpx);
          if k == 1
	      profil0d_out.xli = x_out;
          end
     end

      % boucle sur les champs
      for l = 1:length(noms)
          y_in = profil0d_in.(noms{l})(k,:);
          if all(isfinite(y_in))
	     profil0d_out.(noms{l})(k,:) =  interp1(x_in,y_in,x_out,'pchip','extrap');
          elseif any(isfinite(y_in))
	     profil0d_out.(noms{l})(k,:) =  interp1(x_in,y_in,x_out,'nearest','extrap');
	  end
      end
       
end


if isfield(profil0d_in,'Rsepa') && isfield(profil0d_in,'Zsepa')
  profil0d_out.Rsepa = profil0d_in.Rsepa;
  profil0d_out.Zsepa = profil0d_in.Zsepa;
end
