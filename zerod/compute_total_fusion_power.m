% compute total fusion power from METIS data.
% do not include energy released in breeding blankets:
% see for example: https://iopscience.iop.org/article/10.1088/1742-6596/1493/1/012003/pdf
function [pfusion_total,pfusion_DT,pfusion_DD,pfusion_TT,neutron_total_tt, ...
          frac_dt_thermal,frac_dt_beam_plasma,frac_dt_beam_beam,frac_dt_icrh, ...
          frac_dd_thermal,frac_dd_beam_plasma,frac_dd_beam_beam, ...
          frac_tt_thermal,frac_tt_beam_plasma,frac_tt_beam_beam]  =  ...
          compute_total_fusion_power(post)

% various used
if ~isfield(post,'profil0d')
    post.profil0d = post.profil;
end
if ~isfield(post,'zerod')
    post.zerod = post.zs;
end
      
% retreive data for particle energies
[dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion(post.zerod.tite .* post.zerod.tem);
%
pfusion_DT = max(0,(dt.n + dt.he4) ./ dt.he4  * (post.zerod.pfus   - post.zerod.pddfus - post.zerod.pttfus + post.zerod.pfus_loss));
%
pfusion_DD = max(0,post.zerod.pddfus  + post.zerod.ndd .* dd_n.n .* 1.602176462e-19); % (proton and He3 in pddfus, should add neutron )
%
% TT neutron
%pttfus    = phys.e .* 1.6e6 .* neutron_total ./ 2;
%ntt         = 2 .* post.zerod.pttfus ./ 1.6e6 ./ 1.602176462e-19;
%pfusion_TT  = post.zerod.pttfus  + ntt .* (tt.n/2) .* 1.602176462e-19;  
% tt.n is the total released energy by the reaction, including 2 neutrons and He4
% as it is a reaction with 3 products, the spectrum in energy of each
% spectrum is complex:
% https://www.sciencedirect.com/science/article/pii/0029558265900404
%
pfusion_TT  = max(0,post.zerod.pttfus ./ 1.6e6 .* tt.n); 
%
pfusion_total = pfusion_DT + pfusion_DD + pfusion_TT;

% now compute contribution from thermal, beam-plasma, beam-beam
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
    %
    pfusion_total = interp1(temps,pfusion_total,post.profil0d.temps,'nearest','extrap');
    pfusion_DT    = interp1(temps,pfusion_DT,post.profil0d.temps,'nearest','extrap');
    pfusion_DD    = interp1(temps,pfusion_DD,post.profil0d.temps,'nearest','extrap');
    pfusion_TT    = interp1(temps,pfusion_TT,post.profil0d.temps,'nearest','extrap');

end
%
[~,salpha,~,~,~,~,~,~,~,~,~,~,~,splustd,splusdt,splusff,splusicrh] = ...
    zfus0tae(zs.nDm,zs.nTm,zs.tem,zs.nem,zs.zeff,zs.tite,geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,zs.vp,zs.sp, ...
    zs.pnbi_th,zs.taus_nbi,zs.ecrit_nbi,post.z0dinput.option.einj,post.z0dinput.option.einj2 ,cons.ftnbi, ...
    zs.pion_icrh,zs.taus_icrh,zs.ecrit_nbi,zs.einj_icrh,post.profil0d.temps,cons.pnbi,zs.d0,zs.qa,zs.qmin, ...
    zs.te0,zs.nebord,zs.tebord,zs.pped,zs.nim,zs.wth, ...
    post.z0dinput.option.tae,post.z0dinput.option.nb_nbi,post.z0dinput.option.fspot,post.z0dinput.option.e_shielding, ...
    post.profil0d,post.z0dinput.option.fpolarized,post.z0dinput.option.forced_H_NBI);
%
[neutron_total_tt,neutron_th_tt,neutron_nbi_th_tt,neutron_nbi_nbi_tt,pttfus,proton_tt,picrh_nbi_tt,einj_tt] = ...
    z0neutron_tt(post.z0dinput.option,cons,zs,post.profil0d);

% total source DT
dt_total        = salpha .* (zs.pfus + zs.pfus_loss) ./ zs.pfus;
frac_dt_thermal = max(0,(salpha .* (zs.pfus + zs.pfus_loss) ./ zs.pfus -  ...
                  (splustd + splusdt + splusff + splusicrh)) ./ max(eps,dt_total));
frac_dt_beam_plasma = max(0,(splustd + splusdt) ./ dt_total);
frac_dt_beam_beam   = max(0,splusff ./ dt_total);
frac_dt_icrh        = max(0,splusicrh ./ dt_total);

% DD
frac_dd_thermal     = max(0,zs.ndd_th ./ max(eps,zs.ndd));
frac_dd_beam_plasma = max(0,zs.ndd_nbi_th./ max(eps,zs.ndd));
frac_dd_beam_beam   = max(0,zs.ndd_nbi_nbi./ max(eps,zs.ndd));

% TT
frac_tt_thermal     = max(0,neutron_th_tt ./ max(eps,neutron_total_tt));
frac_tt_beam_plasma = max(0,neutron_nbi_th_tt ./ max(eps,neutron_total_tt));
frac_tt_beam_beam   = max(0,neutron_nbi_nbi_tt ./ max(eps,neutron_total_tt));


if nargout == 0
    h = findobj(0,'type','figure','tag','fusion_total');
    if isempty(h)
        h=figure('tag','fusion_total');
    else
        figure(h);
    end
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
        'defaultlinelinewidth',1,'color',[1 1 1])
    ax(1) = subplot(4,1,1);
    semilogy(zs.temps,pfusion_total/1e6,zs.temps,pfusion_DT/1e6,'-.',zs.temps,pfusion_DD/1e6,'--',zs.temps,pfusion_TT/1e6,'-.');
    legend('Total','DT','DD','TT');
    title('Fusion power (including neutrons contribution)');
    %xlabel('time (s)');
    ylabel('MW');
    z0loglin(gca);
    %
    ax(2) = subplot(4,1,2);
    semilogy(zs.temps,pfusion_DT/1e6, ...
         zs.temps,frac_dt_thermal .* pfusion_DT/1e6, ...
         zs.temps,frac_dt_beam_plasma .* pfusion_DT/1e6, ...
         zs.temps,frac_dt_beam_beam .* pfusion_DT/1e6, ...ntt
         zs.temps,frac_dt_icrh .* pfusion_DT/1e6);
    legend('Total','Thermal','Beam-plasma','Beam-beam','ICRH induced');
    title('Fusion power (DT reactions)');
    %xlabel('time (s)');
    ylabel('MW');
    z0loglin(gca);
    %
    ax(3) = subplot(4,1,3);
    semilogy(zs.temps,pfusion_DD/1e6, ...
         zs.temps,frac_dd_thermal .* pfusion_DD/1e6, ...
         zs.temps,frac_dd_beam_plasma .* pfusion_DD/1e6, ...
         zs.temps,frac_dd_beam_beam .* pfusion_DD/1e6);
    legend('Total','Thermal','Beam-plasma','Beam-beam');
    title('Fusion power (DD reactions)');
    %xlabel('time (s)');
    ylabel('MW');
    z0loglin(gca);
    %
    ax(4) = subplot(4,1,4);
    semilogy(zs.temps,pfusion_TT/1e6, ...
         zs.temps,frac_tt_thermal .* pfusion_TT/1e6, ...
         zs.temps,frac_tt_beam_plasma .* pfusion_TT/1e6, ...
         zs.temps,frac_tt_beam_beam .* pfusion_TT/1e6);
    legend('Total','Thermal','Beam-plasma','Beam-beam');
    title('Fusion power (TT reactions)');
    xlabel('time (s)');
    ylabel('MW');
    z0loglin(gca);
    % 
    linkaxes(ax,'x');
    %
    edition2
    warndlg(sprintf('%s\n%s\n%s\n%s', ...
           'This founction computes the total fusion power from METIS data.', ...
           'Energy released in breeding blankets is not included', ...
           'For details about energy released in breeding blankets, see for example:', ... 
           'https://iopscience.iop.org/article/10.1088/1742-6596/1493/1/012003/pdf'));
end


    
