function post = copy_ual_parameters_2_output(post,z0dinput)

info = metis4imas;
noms = fieldnames(info.valeur);

for k=1:length(noms)
  nomc = noms{k};
  switch info.section.(nomc)
  case {'UAL','Occurrence UAL'}
      if isfield(z0dinput.option,nomc)
	  post.z0dinput.option.(nomc) = z0dinput.option.(nomc);
      end
  end
end
