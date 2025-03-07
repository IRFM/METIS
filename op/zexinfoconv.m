function out = zexinfoconv(file)

if nargin < 1
   error('il faut donner le nom du fichier')
end

[lfile,rmfile,cr]=zgzip(file,'uncompress');
load(lfile)
if ~isempty(rmfile)
	[voids,voidt]=unix(['rm -f ',rmfile,' >& /dev/null']);
end
% test du contenu
if ~isstruct(param) | ~isstruct(data)
	error('Format de fichier incorrect !')
else
   % decompactage
   data = zreduit(param,data,'uncompact');
end

% recupere les donnees utiles
out.cn       = param.gene.cn;
out.amorti   = param.gene.amorti;
out.x        = param.gene.x;
out.t        = data.gene.temps;
out.deltat   = data.gene.dt;

cpu =real(data.gene.cputime);
ind = find(isfinite(cpu));
dcpudt = zdxdt(cpu(ind),data.gene.temps(ind));
cpu = cumtrapz(data.gene.temps(ind),dcpudt .* (dcpudt >=0));
out.cpu      = cpu;

dd =real(data.gene.datation);
ind = find(isfinite(dd));
ddddt = zdxdt(dd(ind),data.gene.temps(ind));
dd = cumtrapz(data.gene.temps(ind),ddddt .* (ddddt >=0));
out.date     = dd;

out.memory   = data.gene.memory;
out.conv     = data.gene.conv;
out.nbsplit  = data.gene.nbsplit;

out.psi      = data.prof.psi;
out.ne       = data.prof.ne;
out.ae       = data.prof.ae;
out.pe       = data.prof.pe;
out.pion     = data.prof.pion;
out.qdiff    = data.prof.q;
out.qeq      = data.equi.q;
out.jmoy     = data.prof.jmoy;
out.jmoyeq   = data.equi.jmoy;


[dt,tmin,tmax] = zgetdt(file);
out.dt       = dt;
out.tmin     = tmin;
out.tmax     = tmax;


out.kmin     = param.gene.kmin;
out.kmax     = param.gene.k;
out.ipexp    = data.exp.ip;
out.ip       = data.gene.ip;
out.ipeq     = data.equi.ip;

out.vexp     = data.exp.vloop;
out.vres     = data.gene.vres;
out.vloop    = data.gene.vloop;
out.vsurf    = data.gene.vsurf;

out.liexp    = data.exp.li;
out.li       = data.gene.li;
out.licr     = data.gene.licr;
out.lieq     = data.equi.li;
out.liip     = data.gene.liip;

out.qaexp    = data.exp.qa;
out.qa       = data.prof.q(:,end);
out.qaeq     = data.equi.q(:,end);


