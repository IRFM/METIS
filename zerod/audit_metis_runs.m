% test d'existence
if exist('post','var') && exist('jeux1','var')
    if isfield(post,'z0dinput') && isfield(jeux1,'post') && isfield(jeux1.post,'z0dinput')
	% comparaison des structures
	if length(post.profil0d.temps) ==  length(jeux1.post.profil0d.temps)
	    fprintf('difference in zerod data :\n');
	    zcompstruct(post.zerod,jeux1.post.zerod,max(1e-3,post.z0dinput.option.tol0d));
	    fprintf('difference in profiles :\n');
	    zcompstruct(post.profil0d,jeux1.post.profil0d,max(1e-3,post.z0dinput.option.tol0d));
	end
	zplotstruct(post.zerod,jeux1.post.zerod,'zerod data',sprintf('METIS : %s@%d/COMPARAISON zerod ',post.z0dinput.option.machine,post.z0dinput.option.shot));
	fullscreen = get(0,'ScreenSize');
        set(gcf,'Position',fullscreen);
        drawnow
	if length(post.profil0d.temps) ==  length(jeux1.post.profil0d.temps)
	      zplotstruct(post.profil0d,jeux1.post.profil0d,'profiles',sprintf('METIS : %s@%d/COMPARAISON profiles ',post.z0dinput.option.machine,post.z0dinput.option.shot));
	      fullscreen = get(0,'ScreenSize');
	      set(gcf,'Position',fullscreen);
	      drawnow
	end
    end
end