ind = 4
x   = param.gene.x;
ve  = ones(size(x));
vt  = ones(length(ind),1);
qmin = min(data.equi.q(ind,:));
Bt   = sqrt(data.equi.b2(ind,end));
a    = data.equi.a(ind,end);
R    = data.equi.rmoy(ind,end);
pfus = data.gene.paddfus(ind);
indg   = find(param.compo.z == 2 & param.compo.a == 4);
s       = squeeze(data.source.n.fus(ind,:,indg));
zeff = data.gene.zeffm(ind);
spr    = data.equi.spr(ind,:) .* (data.equi.rhomax(ind) * ve);
esupfus = (3/2) .* data.equi.rhomax(ind) .* trapz(x,data.equi.vpr(ind,:) .* data.source.fus.psupra(ind,:),2);

valpha  = sqrt(2 .*  3.56e6 .* 1.602176462e-19 ./  1.66053873e-27 ./ 4);
wca     = 2 .* 1.602176462e-19 .* Bt ./ 1.66053873e-27 ./ 4;
dp0R    = (2 .* (max(1,qmin)*ve) .* valpha ./ (wca * ve) ./ (R*ve)) .^ (2/3);
% calcul de l'effet du courant de retour
GZ      = (1.55 + 0.85./(zeff*ve)) .* sqrt((a./R) * x) - (0.2 + 1.55./(zeff*ve)) .* ((a./R) * x);
jfus    = 0.5 .* 2 .* 1.602176462e-19 .* valpha .* dp0R .* s .* (1 - 2 .* (1-GZ) ./ (zeff*ve));
mask    = (jfus == (max(jfus,[],2)*ve));
xfus    = sum((vt*x) .* mask,2) ./sum(mask,2); 
jxfus   = sum(jfus .* mask,2) ./sum(mask,2);  
j0fus   = jfus(:,1);  
%dfus    = sqrt(trapz(x,(vt*x) .^ 2 .* jfus,2) ./trapz(x,jfus,2) - xfus.^ 2);
ifus0d   = trapz(x,spr .* jfus,2);
ifus0d  = ifus0d .* esupfus ./ max(1,pfus);
ifus    = trapz(x,spr .* data.source.fus.j(ind,:),2);
logCoulomb = 15;
tslow = 6.27e8*4*(data.prof.te(ind,1).^1.5)/(data.prof.ne(ind,1)*1e-6*4*logCoulomb); % estimation simple du temps de slowing down
fspot   = ifus ./ ifus0d
