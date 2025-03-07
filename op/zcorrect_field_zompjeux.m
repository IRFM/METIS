% script de modification des donnees avant appel de zcompjeux
try
	jeux1.param.memoire.neo.data.neodir = param.memoire.neo.data.neodir;
end
try
	jeux1.param.from.sample.signal.defaut.inf = param.from.sample.signal.defaut.inf;
end
try
	jeux1.param.from.sample.groupe.defaut.inf = param.from.sample.groupe.defaut.inf;
end
try
	jeux1.param.edit.currentfile = param.edit.currentfile;
end
try
	jeux1.param.gene.date_exec = param.gene.date_exec;
end
try
	jeux1.param.gene.file = param.gene.file;
end
try
	jeux1.param.gene.origine = param.gene.origine;
end
try
	jeux1.param.gene.rapsauve = param.gene.rapsauve;
end
try
	jeux1.data.gene.datation = data.gene.datation;
end
try
	jeux1.param.gene.computer = param.gene.computer;
end
try
	jeux1.param.gene.memrapsauve = param.gene.memrapsauve;
end
try
	jeux1.data.gene.neutralite = data.gene.neutralite;
end
% pour les derivees de ae qui sont souvent quasi nulles
try
	jeux1.data.prof.aed1 = jeux1.data.prof.aed1 + data.prof.ae ./ mean(diff(param.gene.x)) ./ 1e-6;
end
try
	jeux1.data.prof.aed2 = jeux1.data.prof.aed2 + data.prof.ae ./ mean(diff(param.gene.x)) .^ 2 ./ 1e-6;
end
try
	data.prof.aed1 = data.prof.aed1 + data.prof.ae ./ mean(diff(param.gene.x)) ./ 1e-6;
end
try
	data.prof.aed2 = data.prof.aed2 + data.prof.ae ./ mean(diff(param.gene.x)) .^ 2 ./ 1e-6;
end
try
	jeux1.data.equi.mhd.mercier = data.equi.mhd.mercier;
end
try
	jeux1.data.equi.mhd.ideal = data.equi.mhd.ideal;
end

