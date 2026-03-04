function option = set_missing_al_parameters(option)

% add codeparam to options
info = metis4imas(1);
option_default = info.valeur;
noms = fieldnames(option_default);
for k=1:length(noms)
    if ~isfield(option,noms{k})
        option.(noms{k}) = option_default.(noms{k});
    end
end
% secure value
option = secure_option_imas(option);
