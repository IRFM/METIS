% mapping du CPO coreprof (seule les variables d'interet sont remplies, un model vide doit etre fourni)
function coresource = mapcore_sources_generic_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coresource,sigma_B0_eff)


coresource.ids_properties.comment = 'METIS sources';
coresource.time		    = profil0d.temps;
rb0 =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
coresource.vacuum_toroidal_field.r0          = mean(z0dstruct.z0dinput.geo.R);
coresource.vacuum_toroidal_field.b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ coresource.vacuum_toroidal_field.r0;

