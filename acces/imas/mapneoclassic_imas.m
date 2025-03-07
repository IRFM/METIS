% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function neoclassic = mapneoclassic_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,neoclassic,vtor,vpol,sigma_B0_eff)

% neoclassic.datainfo est le meme que celui de scenario
if ~isfield(neoclassic,'datainfo')
	neoclassic.datainfo = datainfo_empty_imas;
end
neoclassic.time			        = profil0d.temps;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
neoclassic.rho_tor_norm 		= profil0d.xli;
neoclassic.rho_tor       		= profil0d.rmx;
rb0 =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
r0  = mean(z0dstruct.z0dinput.geo.R);
%b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ r0;
% les courants ont déja le bon signe
b0  = rb0 ./ r0;
%
neoclassic.jboot = (profil0d.jrun + profil0d.jboot) .* ((b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
neoclassic.sigma = 1./max(1e-307,profil0d.eta);
neoclassic.vpol  = vpol;
neoclassic.er    = profil0d.er;
neoclassic.ne_neo.vconv_eff = profil0d.ware;


