% recompute ndd to prevent effect of transient during convergence on Ti
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
[out.ndd,out.ndd_th,out.ndd_nbi_th,out.ndd_nbi_nbi,out.pddfus,void_proton_dd,out.pnbi_icrh,out.einj_nbi_icrh] = ...
z0neutron_dd(post.z0dinput.option,cons,zs,post.profil0d);

noms = fieldnames(out);
for k=1:length(noms)
  post.zerod.(noms{k}) = interp1(zs.temps,out.(noms{k}),post.zerod.temps,'nearest','extrap');
end
