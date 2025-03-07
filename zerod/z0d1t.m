% script de preparation d'un temps pour un simulation demo
info = z0dacces;
opc = info.valeur;
[cr,data,param] = z0dacces('''',opc,post);
indc = [];
tsscenario;
while isempty(indc)
	tc   = input('temps ? (s) ');
	indc = min(find(data.gene.temps >= tc));
end
param.gene.kmin    = indc;
param.gene.k       = indc;
param.gene.reprise = 0;
param.gene.rapsauve = '';
param.gene.file = '';
param.gene.origine = '';
param.gene.nbsauve  = Inf;
param.gene.corrae   = 0;
data.mode.fci      = ones(size(data.mode.fci));
data.mode.hyb      = ones(size(data.mode.fci));
data.mode.idn      = ones(size(data.mode.fci));
data.mode.fus      = ones(size(data.mode.fci));
data.mode.fce      = ones(size(data.mode.fci));
data.mode.impur    = ones(size(data.mode.fci));
data.mode.bord     = ones(size(data.mode.fci));
data.mode.brem     = ones(size(data.mode.fci));
data.mode.n0       = zeros(size(data.mode.fci));
data.mode.rip      = zeros(size(data.mode.fci));
data.mode.rad      = ones(size(data.mode.fci));
data.mode.cyclo    = ones(size(data.mode.fci));
data.mode.ohm      = ones(size(data.mode.fci));
data.mode.zeff     = ones(size(data.mode.fci));
data.mode.ae       = ones(size(data.mode.fci));

[cr,data,param]=zpremiertemps(data,param);
datak    = zget1t(data,indc);
figure
plot(double(squeeze(datak.equi.R))',double(squeeze(datak.equi.Z))');


figure
subplot(2,2,1)
plot(param.gene.x,squeeze(datak.impur.impur),'or',param.gene.x,squeeze(dataka.impur.impur),'b')

subplot(2,2,2)
plot(param.gene.x,datak.prof.ni,'ro',param.gene.x,dataka.prof.ni,'b')

subplot(2,2,3)
plot(param.gene.x,datak.prof.ae,'ro',param.gene.x,dataka.prof.ae,'b', ...
     param.gene.x,datak.impur.ae,'mo',param.gene.x,dataka.impur.ae,'c')

subplot(2,2,4)
plot(param.gene.x,datak.prof.ne,'ro',param.gene.x,dataka.prof.ne,'b')



