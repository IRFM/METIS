k = 150;
nbth = 65;
profil0d = post.profil0d;
z0dstruct = post

   % donnees de metis 
    prof.x    = profil0d.xli;
    prof.kx   = profil0d.kx(k,:);     
    prof.dx   = profil0d.dx(k,:);      
    prof.rmx  = profil0d.rmx(k,:);     
    prof.Raxe = profil0d.Raxe(k,:);
    prof.psi  = profil0d.psi(k,:);
    prof.fdia = profil0d.fdia(k,:);
    prof.jmoy = profil0d.jli(k,:);
    prof.ptot = profil0d.ptot(k,:);
    %
    geo_input.a     = interp1(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.a ,profil0d.temps(k),'pchip','extrap');
    geo_input.R     = interp1(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps(k),'pchip','extrap');
    geo_input.K     = interp1(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.K ,profil0d.temps(k),'pchip','extrap');
    geo_input.d     = interp1(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.d ,profil0d.temps(k),'pchip','extrap');
    geo_input.b0    = interp1(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.b0 ,profil0d.temps(k),'pchip','extrap'); 
    geo_input.z0    = interp1(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.z0 ,profil0d.temps(k),'pchip','extrap');
    geo_input.sp    = interp1(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.sp ,profil0d.temps(k),'pchip','extrap');
    geo_input.vp    = interp1(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.vp ,profil0d.temps(k),'pchip','extrap');
    geo_input.sext  = interp1(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.sext ,profil0d.temps(k),'pchip','extrap');
    if isfield(profil0d,'Rsepa') &&isfield(profil0d,'Zsepa')
    	geo_input.Rsepa = profil0d.Rsepa(k,:);
    	geo_input.Zsepa = profil0d.Zsepa(k,:);
    else
    	geo_input.Rsepa = [];
    	geo_input.Zsepa = []; 
	z0dstruct.z0dinput.option.morphing = 0;
    end
    %
    phys.mu0 = 4 .* pi .* 1e-7;
    [profil,deuxd,moment,scalaire,facteur] = metis2equi1t(prof,geo_input,phys,nbth,1, ...
                                             z0dstruct.z0dinput.option.morphing,0);
