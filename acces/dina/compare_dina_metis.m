% script comparing DINA and METIS P' and FF'
clear
shot = 105092
uri_dina = sprintf('imas:hdf5?path=/home/ITER/artaudj1/PDS/pds/scenario_configs/%d_metis_predictive_nice_inverse_from_dina/tmp/data/%d_dina_in',shot,shot);
uri_metis = sprintf('imas:hdf5?path=/home/ITER/artaudj1/PDS/pds/scenario_configs/%d_metis_predictive_nice_inverse_from_dina/tmp/data/%d_metis_in',shot,shot);
idx = imas_open(uri_dina);
equi_dina = ids_get(idx,'equilibrium');
imas_close(idx);
idx = imas_open(uri_metis);
equi_metis = ids_get(idx,'equilibrium');
imas_close(idx);

tdina = equi_dina.time(:)';
for k=1:length(tdina)
    dpdpsi_dina(k,:) = equi_dina.time_slice{k}.profiles_1d.dpressure_dpsi;
    df2dpsi_dina(k,:) = 2 .* equi_dina.time_slice{k}.profiles_1d.f_df_dpsi;
    psi_dina(k,:) = equi_dina.time_slice{k}.profiles_1d.psi;
    jtor_dina(k,:) = equi_dina.time_slice{k}.profiles_1d.j_phi;
    psin_dina(k,:) = (psi_dina(k,:) - psi_dina(k,1)) ./ (psi_dina(k,end) - psi_dina(k,1));
end

tmetis = equi_metis.time(:)';
for k=1:length(tmetis)
    dpdpsi_metis(k,:) = equi_metis.time_slice{k}.profiles_1d.dpressure_dpsi;
    df2dpsi_metis(k,:) = 2 .* equi_metis.time_slice{k}.profiles_1d.f_df_dpsi;
    psi_metis(k,:) = equi_metis.time_slice{k}.profiles_1d.psi;
    jtor_metis(k,:) = equi_metis.time_slice{k}.profiles_1d.j_phi;
    psin_metis(k,:) = (psi_metis(k,:) - psi_metis(k,1)) ./ (psi_metis(k,end) - psi_metis(k,1));
end


figure;
subplot(2,2,1)
hold on
zplotprof(gca,tdina,psin_dina,dpdpsi_dina,'color','r');
zplotprof(gca,tmetis,psin_metis,dpdpsi_metis,'color','b');
xlabel('PSI_N')
ylabel('dPdpsi')
title(sprintf('shot = %d',shot))
legend({'DINA','METIS'})
subplot(2,2,2)
hold on
zplotprof(gca,tdina,psin_dina,df2dpsi_dina,'color','r');
zplotprof(gca,tmetis,psin_metis,df2dpsi_metis,'color','b');
xlabel('PSI_N')
ylabel('dF2dpsi')
title(sprintf('shot = %d',shot))
legend({'DINA','METIS'})
subplot(2,2,3)
hold on
zplotprof(gca,tdina,psin_dina,psi_dina,'color','r');
zplotprof(gca,tmetis,psin_metis,psi_metis,'color','b');
xlabel('PSI_N')
ylabel('Psi')
title(sprintf('shot = %d',shot))
legend({'DINA','METIS'})
subplot(2,2,4)
hold on
zplotprof(gca,tdina,psin_dina,jtor_dina,'color','r');
zplotprof(gca,tmetis,psin_metis,jtor_metis,'color','b');
xlabel('PSI_N')
ylabel('J_{phi}')
title(sprintf('shot = %d',shot))
legend({'DINA','METIS'})





