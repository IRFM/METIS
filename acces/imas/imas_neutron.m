function [t,ntot,ndd,ndd_th,ndd_nbi_th,ndd_nbi_nbi,ndt,ndt_th,ndt_nbi_th,ndt_nbi_nbi,ndt_nbi_icrh,ntt,ntt_th,ntt_nbi_th,ntt_nbi_nbi] = imas_neutron(post)

% compatibility with zerodevolution data structure
if ~isfield(post,'profil0d') 
    post.profil0d = post.profil;
end

% initialisation
zz           = zeros(size(post.z0dinput.cons.temps));
%  ntot         = zz;
%  ndd          = zz;
%  ndd_th       = zz;
%  ndd_nbi_th   = zz;
%  ndd_nbi_nbi  = zz;
ndt          = zz;
ndt_th       = zz;
ndt_nbi_th   = zz;
ndt_nbi_nbi  = zz;
ndt_nbi_icrh = zz;
ntt          = zz;
ntt_th       = zz;
ntt_nbi_th   = zz;
ntt_nbi_nbi  = zz;

% DD neutrons
t            = post.z0dinput.cons.temps;
ndd          = post.zerod.ndd;
ndd_th       = post.zerod.ndd_th;
ndd_nbi_th   = post.zerod.ndd_nbi_th;
ndd_nbi_nbi  = post.zerod.ndd_nbi_nbi;

% switch between nuclear and non nuclear tokamak
if post.z0dinput.option.gaz == 3
    if length(post.profil0d.temps) == length(post.zerod.temps)
        zs   = post.zerod;
        geo  = post.z0dinput.geo;
        cons = post.z0dinput.cons;
    else
        noms = fieldnames(post.zerod);
        temps = post.zerod.temps;
        for k=1:length(noms)
            nomc = noms{k};
            var  = post.zerod.(nomc);
            if length(var) == length(temps)
                zs.(nomc) = interp1(temps,var,post.profil0d.temps,'nearest','extrap');
            else
                zs.(nomc) = var;
            end
        end
        zs.temps = post.profil0d.temps;
        noms = fieldnames(post.z0dinput.cons);
        temps = post.z0dinput.cons.temps;
        for k=1:length(noms)
            nomc = noms{k};
            var  = post.z0dinput.cons.(nomc);
            if length(var) == length(temps)
                cons.(nomc) = interp1(temps,var,post.profil0d.temps,'nearest','extrap');
            else
                cons.(nomc) = var;
            end
        end
        cons.temps = post.profil0d.temps;
        noms = fieldnames(post.z0dinput.geo);
        for k=1:length(noms)
            nomc = noms{k};
            var  = post.z0dinput.geo.(nomc);
            if length(var) == length(temps)
                geo.(nomc) = interp1(temps,var,post.profil0d.temps,'nearest','extrap');
            else
                geo.(nomc) = var;
            end
        end
    end
    [~,salpha,~,~,~,~,~,~,~,~,~,~,~,splustd,splusdt,splusff,splusicrh] = ...
        zfus0tae(zs.nDm,zs.nTm,zs.tem,zs.nem,zs.zeff,zs.tite,geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,zs.vp,zs.sp, ...
        zs.pnbi_th,zs.taus_nbi,zs.ecrit_nbi,post.z0dinput.option.einj,post.z0dinput.option.einj2 ,cons.ftnbi, ...
        zs.pion_icrh,zs.taus_icrh,zs.ecrit_nbi,zs.einj_icrh,post.profil0d.temps,cons.pnbi,zs.d0,zs.qa,zs.qmin, ...
        zs.te0,zs.nebord,zs.tebord,zs.pped,zs.nim,zs.wth, ...
        post.z0dinput.option.tae,post.z0dinput.option.nb_nbi,post.z0dinput.option.fspot,post.z0dinput.option.e_shielding,post.profil0d, ...
        post.z0dinput.option.fpolarized, post.z0dinput.option.forced_H_NBI);
    
    % resmplaing if needed
    if length(post.profil0d.temps) ~= length(post.zerod.temps)
        salpha    = interp1(post.profil0d.temps,salpha,post.zerod.temps,'nearest','extrap');
        splustd   = interp1(post.profil0d.temps,splustd,post.zerod.temps,'nearest','extrap');
        splusdt   = interp1(post.profil0d.temps,splusdt,post.zerod.temps,'nearest','extrap');
        splusff   = interp1(post.profil0d.temps,splusff,post.zerod.temps,'nearest','extrap');
        splusicrh = interp1(post.profil0d.temps,splusicrh,post.zerod.temps,'nearest','extrap');
    end
    ndt          = salpha;
    ndt_th       = salpha - (splustd + splusdt + splusff + splusicrh);
    ndt_nbi_th   = splustd + splusdt;
    ndt_nbi_nbi  = splusff;
    ndt_nbi_icrh = splusicrh;
    
    [ntt,ntt_th,ntt_nbi_th,ntt_nbi_nbi,pttfus,proton_tt,picrh_nbi_tt,einj_tt] = ...
        z0neutron_tt(post.z0dinput.option,cons,zs,post.profil0d);
    
    % resmplaing if needed
    if length(post.profil0d.temps) ~= length(post.zerod.temps)
        ntt            = interp1(post.profil0d.temps,ntt,post.zerod.temps,'nearest','extrap');
        ntt_th         = interp1(post.profil0d.temps,ntt_th,post.zerod.temps,'nearest','extrap');
        ntt_nbi_th     = interp1(post.profil0d.temps,ntt_nbi_th,post.zerod.temps,'nearest','extrap');
        ntt_nbi_nbi    = interp1(post.profil0d.temps,ntt_nbi_nbi,post.zerod.temps,'nearest','extrap');
        pttfus         = interp1(post.profil0d.temps,pttfus,post.zerod.temps,'nearest','extrap');
        proton_tt      = interp1(post.profil0d.temps,proton_tt,post.zerod.temps,'nearest','extrap');
        picrh_nbi_tt   = interp1(post.profil0d.temps,picrh_nbi_tt,post.zerod.temps,'nearest','extrap');
        einj_tt        = interp1(post.profil0d.temps,einj_tt,post.zerod.temps,'nearest','extrap');
   end
    
    
    
    
end
%
ntot = ndt + ndd + ntt;

