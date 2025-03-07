root = fileparts(which('metis'));
rep_ = dir(fullfile(root,'certification','metis','*.mat'));

for k_ = 1:length(rep_)

	name_  = fullfile(root,'certification','metis',rep_(k_).name);
	metis_load(name_)
        drawnow
	if length(post.zerod.temps) == length(post.profil0d.temps)
		% calcul complet
		metis_run;
	else
		% calcul rapide
		metis_fast;
	end
        drawnow
	metis_save(name_);
        drawnow
        clear name_
        munlock
        clear functions
end




