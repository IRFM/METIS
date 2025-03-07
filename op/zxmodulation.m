function zxmodulation(data,param,post);

% les donnees
cronos.paddfce  = data.gene.paddfce;
cronos.temps    = data.gene.temps;
cronos.tece     = post.ece.tece;
cronos.Rece     = post.ece.Rece;
cronos.rsuraece = post.ece.rsuraece;

% creation du repertoire
racine  = sprintf('/usr/drfc/cgc/matlab5/zineb/data/zxmodulation_%s',param.from.machine);
if ~isdir(racine)
	unix(sprintf('mkdir %s',racine));
end
% nom du fichier
[pp,ff] = fileparts(param.gene.file);
ff      = strrep(ff,'zineb','modulation');
save(fullfile(racine,ff),'cronos');
