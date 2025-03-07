function flux_data_output_ = poloidal_flux_data(post,couplage)

% z0plotflux.m
NOPLUTFLUX_HELIOS = 1;
z0plotflux;

% create data structure
liste_of_fields_ = who;
for kkk_ = 1:length(liste_of_fields_)
    flux_data_output_.(liste_of_fields_{kkk_}) = eval(liste_of_fields_{kkk_});
end
