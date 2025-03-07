% script servant a prolonger un run metis
tmax = input('tmax ? (s)');
dt = input(' time step ? (s)');
tadd = ((max(z0dinput.cons.temps) + dt):dt:tmax)';
z0dinput_mem = z0dinput;
va   = ones(size(tadd));

noms = fieldnames(z0dinput.cons);
for l = 1:length(noms)
	switch noms{l}
	case 'temps'	
		z0dinput.cons.temps = cat(1,z0dinput.cons.temps,tadd);
	otherwise
		z0dinput.cons.(noms{l}) = cat(1,z0dinput.cons.(noms{l}),z0dinput.cons.(noms{l})(end) .* va);
	end
end
noms = fieldnames(z0dinput.geo);
for l = 1:length(noms)
	switch noms{l}
	case 'temps'	
		z0dinput.geo.temps = cat(1,z0dinput.geo.temps,tadd);
	otherwise
		z0dinput.geo.(noms{l}) = cat(1,z0dinput.geo.(noms{l}),z0dinput.geo.(noms{l})(end) .* va);
	end
end

 noms = fieldnames(z0dinput.exp0d);
for l = 1:length(noms)
	switch noms{l}
	case 'temps'	
		z0dinput.exp0d.temps = cat(1,z0dinput.exp0d.temps,tadd);
	otherwise
		z0dinput.exp0d.(noms{l}) = cat(1,z0dinput.exp0d.(noms{l}),va * z0dinput.exp0d.(noms{l})(end,:));
	end
end
 




