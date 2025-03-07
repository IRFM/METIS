function zxrema(data,param,temps,riponoff);

if nargin < 4 
   riponoff = 0;
elseif isempty(riponoff)
   riponoff =0;
end

% extraction des donnees a 1 temps
if length(data.gene.temps) > 1
   ind = min(find(temps <= data.gene.temps));
   if isempty(ind)
      disp('Aucun temps ne correspond');
      return
   end
   datak = zget1t(data,ind);
else
   datak =data;
end

% les donnees
BR     = double(squeeze(datak.equi.BR));
BZ     = double(squeeze(datak.equi.BZ));
BPHI   = double(squeeze(datak.equi.BPHI));
R      = double(squeeze(datak.equi.R));
Z      = double(squeeze(datak.equi.Z));
RHORZ  = double(squeeze(datak.equi.rhoRZ));
if 1>2
	NE     = pchip(param.gene.x.* RHORZ(end),datak.prof.ne,RHORZ);
	TE     = pchip(param.gene.x.* RHORZ(end),datak.prof.te,RHORZ);
	ZEFF   = pchip(param.gene.x.* RHORZ(end),datak.prof.zeff,RHORZ);
	RHO    = pchip(param.gene.x.* RHORZ(end),param.gene.x,RHORZ);
	RSA    = pchip(param.gene.x .* RHORZ(end),datak.equi.a ./datak.equi.a(end),RHORZ);
	VPR    = pchip(param.gene.x .* RHORZ(end),datak.equi.vpr,RHORZ);
else
	NE     = datak.prof.ne;
	TE     = datak.prof.te;
	ZEFF   = datak.prof.zeff;
	RHO    = param.gene.x;
	RSA    = datak.equi.a ./datak.equi.a(end);
	VPR    = datak.equi.vpr;
end
% creation du repertoire
if riponoff ==1
   racine  = sprintf('%s/zxrema_ripple_%s_%d@%g',fileparts(param.gene.file),param.from.machine,fix(param.from.shot.num),temps);
else
   racine  = sprintf('%s/zxrema_%s_%d@%g',fileparts(param.gene.file),param.from.machine,fix(param.from.shot.num),temps);
end
if ~isdir(racine)
	unix(sprintf('mkdir %s',racine));
end
% nom du fichier
save(fullfile(racine,'data'),'NE','TE','ZEFF','RHO','RSA','BR','BZ','BPHI','R','Z','RHORZ','VPR','-V4');
save(fullfile(racine,'NE'),'NE','-V4');
save(fullfile(racine,'TE'),'TE','-V4');
save(fullfile(racine,'ZEFF'),'ZEFF','-V4');
save(fullfile(racine,'RHO'),'RHO','-V4');
save(fullfile(racine,'RSA'),'RSA','-V4');
save(fullfile(racine,'BR'),'BR','-V4');
save(fullfile(racine,'BZ'),'BZ','-V4');
save(fullfile(racine,'BPHI'),'BPHI','-V4');
save(fullfile(racine,'R'),'R','-V4');
save(fullfile(racine,'Z'),'Z','-V4');
save(fullfile(racine,'RHORZ'),'RHORZ','-V4');
save(fullfile(racine,'VPR'),'VPR','-V4');


fprintf('creation de %s\n',racine);