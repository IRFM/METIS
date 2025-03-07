rep =load('testimas');
index = 3;
el = 0;
ion = 0;
ps = 0;
disp(' ')
for k=1:length(rep.exported_data.core_sources.source)
    sl = rep.exported_data.core_sources.source{k};
    gq = sl.global_quantities{index};
    fprintf('%s:\t%6.3g electron;\t %6.3g ion \tand %6.3g total (sum = %6.3g)\n', ...
        sl.identifier.name,gq.electrons.power,gq.total_ion_power, ...
        gq.power, gq.electrons.power+gq.total_ion_power);
    switch sl.identifier.name
        case 'total'
            disp('not total');
        otherwise
            el = el + gq.electrons.power;
            ion = ion + gq.total_ion_power;
            ps = ps + gq.power;
    end
end
fprintf('sum:\t%6.3g electron;\t %6.3g ion \tand %6.3g total (sum = %6.3g)\n', ...
    el,ion,ps,el+ion);

% for Ohmic shot is OK on electron and ion
metis_tot = post.zerod.pohm(index) - post.zerod.prad(index) - post.zerod.pbrem(index) - post.zerod.pcyclo(index) - post.zerod.pioniz(index) + ...
            post.zerod.plh(index) + post.zerod.picrh(index) + real(post.zerod.pnbi(index)) + imag(post.zerod.pnbi(index)) + ...
            post.zerod.pfus(index) + post.zerod.pecrh(index);
        
metis_ion = post.zerod.pion_fus(index) +   post.zerod.pion_icrh(index) - post.zerod.pioniz_i(index) + real(post.zerod.pion_nbi(index)) + imag(post.zerod.pion_nbi(index));
fprintf('METIS  el: %6.3g\n',metis_tot - metis_ion - post.zerod.pei(index));
fprintf('METIS  ion: %6.3g\n',metis_ion + post.zerod.pei(index));
fprintf('METIS  tot: %6.3g\n',metis_tot);


% recompute metis source from profiles and compare to 0d
%fptot    = max(1./(Vp*ve),fplh + fpecrh + fpfweh + fpnbi + fpfus + fpicrh + fpoh - fprad - fpcyclo - fpbrem - fpioniz);
%profil.source_ion = profil.pnbi_ion + profil.picrh_ion + profil.pfus_ion - profil.pioniz_i;
%profil.source_el  = fptot - profil.source_ion;

x = post.profil0d.xli;
vpr = post.profil0d.vpr(index,:);
fprintf('Pohm: %g & %g\n',trapz(x,post.profil0d.pohm(index,:) .* vpr,2),post.zerod.pohm(index));
fprintf('Prad: %g & %g\n',trapz(x,post.profil0d.prad(index,:) .* vpr,2),post.zerod.prad(index));
fprintf('Pbrem: %g & %g\n',trapz(x,post.profil0d.pbrem(index,:) .* vpr,2),post.zerod.pbrem(index));
fprintf('Pcyclo: %g & %g\n',trapz(x,post.profil0d.pcyclo(index,:) .* vpr,2),post.zerod.pcyclo(index));
fprintf('Pioniz: %g & %g\n',trapz(x,post.profil0d.pioniz(index,:) .* vpr,2),post.zerod.pioniz(index));
fprintf('Pei: %g & %g\n',post.profil0d.qei(index,end) ,post.zerod.pei(index));
fprintf('PLH: %g & %g\n',trapz(x,post.profil0d.plh(index,:) .* vpr,2),post.zerod.plh(index));
fprintf('PICRH: %g & %g\n',trapz(x,post.profil0d.picrh(index,:) .* vpr,2)+trapz(x,post.profil0d.pfweh(index,:) .* vpr,2),post.zerod.picrh(index));
fprintf('PNBI: %g & %g\n',trapz(x,post.profil0d.pnbi(index,:) .* vpr,2),real(post.zerod.pnbi(index)) + imag(post.zerod.pnbi(index)));
fprintf('PFUS: %g & %g (dd = %g)\n',trapz(x,post.profil0d.pfus(index,:) .* vpr,2),post.zerod.pfus(index),post.zerod.pddfus(index));
fprintf('PECRH: %g & %g\n',trapz(x,post.profil0d.pecrh(index,:) .* vpr,2),post.zerod.pecrh(index));

ptotfrom1d = trapz(x,post.profil0d.pohm(index,:) .* vpr,2) - trapz(x,post.profil0d.prad(index,:) .* vpr,2) - ...
             trapz(x,post.profil0d.pbrem(index,:) .* vpr,2) - trapz(x,post.profil0d.pcyclo(index,:) .* vpr,2) - ...
             trapz(x,post.profil0d.pioniz(index,:) .* vpr,2) + trapz(x,post.profil0d.picrh(index,:) .* vpr,2) + ...
             trapz(x,post.profil0d.pfweh(index,:) .* vpr,2) + trapz(x,post.profil0d.pnbi(index,:) .* vpr,2) + ...,        
             trapz(x,post.profil0d.pfus(index,:) .* vpr,2) + trapz(x,post.profil0d.pecrh(index,:) .* vpr,2) + ...
             trapz(x,post.profil0d.plh(index,:) .* vpr,2);
%
ptotrom1d_ion =  trapz(x,post.profil0d.pnbi_ion(index,:) .* vpr,2) +  trapz(x,post.profil0d.picrh_ion(index,:) .* vpr,2)  + ...
                 trapz(x,post.profil0d.pfus_ion(index,:) .* vpr,2) - trapz(x,post.profil0d.pioniz_i(index,:) .* vpr,2);
             
ptotrom1d_el =   ptotfrom1d -  ptotrom1d_ion;
fprintf('from profil0d:\t%6.3g electron;\t %6.3g ion \tand %6.3g total \n', ...
    ptotrom1d_el - post.profil0d.qei(index,end),ptotrom1d_ion + post.profil0d.qei(index,end), ptotfrom1d);

 fprintf('from source_* :\t%6.3g electron;\t %6.3g ion \tand %6.3g total \n', ...
    trapz(x,post.profil0d.source_el(index,:) .* vpr,2)- post.profil0d.qei(index,end), ...
    trapz(x,post.profil0d.source_ion(index,:) .* vpr,2) + post.profil0d.qei(index,end), ...
    trapz(x,post.profil0d.source_ion(index,:) .* vpr,2) + trapz(x,post.profil0d.source_el(index,:) .* vpr,2) );


             



