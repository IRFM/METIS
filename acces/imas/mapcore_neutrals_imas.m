% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coreneutrals = mapcore_neutrals_itm(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coreneutrals,sigma_B0_eff)

% isotopic composition for option.gaz == 5
if z0dstruct.z0dinput.option.gaz == 5
    nHe3onD = real(z0dstruct.z0dinput.cons.iso);
    nTonD   = imag(z0dstruct.z0dinput.cons.iso);
    warning('nHe3onD & nTonD not yet used !');
else
    nHe3onD = zeros(size(z0dstruct.z0dinput.cons.iso));
    nTonD   = real(z0dstruct.z0dinput.cons.iso);
end
z0dstruct.z0dinput.cons.iso = real(z0dstruct.z0dinput.cons.iso);


% coreneutrals.datainfo est le meme que celui de scenario
if ~isfield(coreneutrals,'datainfo')
	coreneutrals.datainfo = datainfo_empty_imas;
end
coreneutrals.time			= profil0d.temps;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coreneutrals.rho_tor_norm 		= xli;
coreneutrals.rho_tor       		= profil0d.rmx;

switch z0dstruct.z0dinput.option.gaz
case 1
	coreneutrals.composition.atomlist.amn = 1;
	coreneutrals.composition.atomlist.zn  = 1;
	coreneutrals.composition.neutallist.ncomp = 1; 
	coreneutrals.composition.neutallist.tatm = 1;
	coreneutrals.composition.neutallist.multatm = 1; 
	coreneutrals.composition.typelist.ntype  = 2; 
	coreneutrals.composition.typelist.type = [0 1]; 
case 2
	coreneutrals.composition.atomlist.amn = 2;
	coreneutrals.composition.atomlist.zn  = 1;
	coreneutrals.composition.neutallist.ncomp = 1; 
	coreneutrals.composition.neutallist.tatm = 1;
	coreneutrals.composition.neutallist.multatm = 1; 
	coreneutrals.composition.typelist.ntype  = 2; 
	coreneutrals.composition.typelist.type = [0 1]; 
case 3
	coreneutrals.composition.atomlist.amn = [2 3];
	coreneutrals.composition.atomlist.zn  = [1 1];
	coreneutrals.composition.neutallist.ncomp = [1 1]; 
	coreneutrals.composition.neutallist.tatm = [1 2];
	coreneutrals.composition.neutallist.multatm = [1 1];
	coreneutrals.composition.typelist.ntype  = [2 2]; 
	coreneutrals.composition.typelist.type = [0 1;0 1]; 
case 4
	coreneutrals.composition.atomlist.amn = 2;
	coreneutrals.composition.atomlist.zn  = 4;
	coreneutrals.composition.neutallist.ncomp = 1; 
	coreneutrals.composition.neutallist.tatm = 1;
	coreneutrals.composition.neutallist.multatm = 1; 
	coreneutrals.composition.typelist.ntype  = 2; 
	coreneutrals.composition.typelist.type = [0 1]; 
otherwise 
	error('plasma compostion not yet implemented in METIITM');
end
compo_void = coreneutrals.composition.atomlist;
compo_void.zion = compo_void.zn;
coreneutrals.compositions = copy_composition_to_compositionstype(compo_void);

%
n0a = interp1_imas(data_zerod.temps,data_zerod.n0a,profil0d.temps,'pchip','extrap');
%
switch z0dstruct.z0dinput.option.gaz
case 3
	iso = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.cons.iso,profil0d.temps,'pchip','extrap');
	iso = iso * ones(size(profil0d.xli));
	coreneutrals.profiles.n0.value            = NaN .* ones(size(profil0d.n0,1),size(profil0d.n0,2),2,2);
	coreneutrals.profiles.n0.value(:,:,1,1)   = reshape(profil0d.n0m ./ (1+iso),size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.value(:,:,1,2)   = reshape(profil0d.n0  ./ (1+iso),size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.value(:,:,2,1)   = reshape(profil0d.n0m ./ (1+iso) .* iso,size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.value(:,:,2,2)   = reshape(profil0d.n0  ./ (1+iso) .* iso,size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.flux             = NaN .* coreneutrals.profiles.n0.value;
	coreneutrals.profiles.n0.flux(:,:,1,1)    = reshape(cumtrapz(profil0d.xli,profil0d.s0m .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho) ./ (1+iso), ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.flux(:,:,1,2)    = reshape(cumtrapz(profil0d.xli,profil0d.s0 .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho) ./ (1+iso), ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.flux(:,:,2,1)    = reshape(cumtrapz(profil0d.xli,profil0d.s0m .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho) ./ (1+iso) .* iso, ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.flux(:,:,2,2)    = reshape(cumtrapz(profil0d.xli,profil0d.s0 .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho) ./ (1+iso) .* iso, ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.boundary.value   = zeros(size(profil0d.n0,1),3,2,2);
	coreneutrals.profiles.n0.boundary.value(:,1,1,1)   = n0a ./ (1+iso(:,1));
	coreneutrals.profiles.n0.boundary.value(:,1,1,2)   =0;
	coreneutrals.profiles.n0.boundary.value(:,1,2,1)   = n0a ./ (1+iso(:,1)) .* iso(:,1);
	coreneutrals.profiles.n0.boundary.value(:,1,2,2)   =0;
	coreneutrals.profiles.n0.boundary.type    = NaN .* ones(size(profil0d.n0,1),2,2);
	coreneutrals.profiles.n0.boundary.type(:)    = 4;  
	coreneutrals.profiles.n0.boundary.rho_tor = NaN .* ones(size(profil0d.n0,1),2,2);
	coreneutrals.profiles.n0.boundary.rho_tor(:,1,1) = coreneutrals.rho_tor(:,end);
	coreneutrals.profiles.n0.boundary.rho_tor(:,1,2) = coreneutrals.rho_tor(:,end);
	coreneutrals.profiles.n0.boundary.rho_tor(:,2,1) = coreneutrals.rho_tor(:,end);
	coreneutrals.profiles.n0.boundary.rho_tor(:,2,2) = coreneutrals.rho_tor(:,end);
otherwise
	coreneutrals.profiles.n0.value            = NaN .* ones(size(profil0d.n0,1),size(profil0d.n0,2),1,2);
	coreneutrals.profiles.n0.value(:,:,1,1)   = reshape(profil0d.n0m,size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.value(:,:,1,2)   = reshape(profil0d.n0,size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.flux             = NaN .* coreneutrals.profiles.n0.value;
	coreneutrals.profiles.n0.flux(:,:,1,1)    = reshape(cumtrapz(profil0d.xli,profil0d.s0m .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho), ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.flux(:,:,1,2)    = reshape(cumtrapz(profil0d.xli,profil0d.s0 .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho), ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.boundary.value   = zeros(size(profil0d.n0,1),3,1,2);
	coreneutrals.profiles.n0.boundary.value(:,1,1,1)   =n0a;
	coreneutrals.profiles.n0.boundary.value(:,1,1,2)   =0;
	coreneutrals.profiles.n0.boundary.type    = NaN .* ones(size(profil0d.n0,1),1,2);
	coreneutrals.profiles.n0.boundary.type(:)    = 4;  
	coreneutrals.profiles.n0.boundary.rho_tor = NaN .* ones(size(profil0d.n0,1),1,2);
	coreneutrals.profiles.n0.boundary.rho_tor(:,1,1) = coreneutrals.rho_tor(:,end);
	coreneutrals.profiles.n0.boundary.rho_tor(:,1,2) = coreneutrals.rho_tor(:,end);
end

