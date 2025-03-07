% calcul de la limite balooning du piedestal
function ppedmax = z0maxpped(temps,profli)

% si pas de caclcul possible
ppedmax  = Inf .* ones(size(temps));

if isfield(profli,'qjli') 
	mu0      =   4*pi*1e-7;
	% calcul de alphac (Onjun PoP 2002) + Wesson
	ppedmax  = 0.8 .* profli.rmx(:,end - 1) .*  profli.fdia(:,end - 1) .^ 2 .* profli.ri(:,end-1) .^ 3 ./ ...
	           (4 .* mu0 .* profli.qjli(:,end - 1) .^ 3) .*  ...
		   (1 + profli.kx(:,end - 1) .^ 2 .* (1 + 5 .* profli.dx(:,end - 1) .^ 2)) .*  ...
		   (profli.qjli(:,end) - profli.qjli(:,end - 1));
	ppedmax  = max(1,ppedmax) + profli.ptot(:,end);

	% ajout de limite en courant 
	% ref : S. Yu. Medvedev , PPCF 48 (2006) p 927-
	% calcul de j//_ped / <J//>
	jeff_moy = trapz(profli.xli, profli.spr .* profli.jeff,2) ./  trapz(profli.xli, profli.spr,2);
	dpped    = profli.ptot(:,end-1) -  profli.ptot(:,end);
	ppedmaxj = dpped .* 1.0 ./ max(eps,profli.jboot(:,end - 1) ./ jeff_moy) +  profli.ptot(:,end);
	ppedmax  = min(ppedmaxj,ppedmax);
end