% sample of user defined scaling for metis
% the name of the script must be : metis_user_scaling
% available variables are :
%  nbar     = line averaged electron density (10^20 m^-3)
%  ip       = plasma current (MA)
%  Bt       = toroidal magnetic field (T)
%  ploss    = loss power (MW)
%  pin      = input power (MW)
%  R        = major radius
%  K        = plasa elongation
%  a        = minor radius
%  meff     = effective number of ion mass 
%  zeff     = effective line averaged plasma charge 
%  Vp       = plasma volume
%  Sp       = poloidal plasma surface
%  q95      = safety factor a 95 % of flux
%  sext     = external plasma surface
%  pohm     = ohmic inpur power 
%  nem      = volume averaged elctron density (m^-3)
%  eta      = volume averaged Ti/Te
%  ae       = volume averaged sum(ion_densities) / electron density
%  rap      = transition from OH (1) to L mode (0) 
taul  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
	R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s       

% ITERH-98P(y,2)        
tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
	R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s    

% OH debut de choc
tauoh =  0.5 .* Vp ./ (2 .* pi .^ 2) .* nbar .* sqrt(qcyl);
% transition d'apres J.G. Cordey rapport JET-P(85)28
rap   = min(1,max(0,max(pohm,ip ./ 1e3) ./ max(pin,pohm)));
tauthl = rap .* tauoh + (1-rap) .* taul;

