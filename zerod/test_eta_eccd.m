% test le calcul de l'efficacite de ECCD dans METIS avec les donnees CRONOS	
aece = z0dinput.option.angle_ece./180*pi;
mutm  = sqrt(data.equi.a .* (1 + cos(aece))  ./ (data.equi.raxe + data.equi.a .* cos(aece)));
etaecem = 1e20 ./ (1 + 100 ./ (data.prof.te./1000)) .* (1 - (1 + (5+data.prof.zeff) ./ 3 ./ (1+data.prof.zeff)) .* ...
		(sqrt(2) .* mutm) .^ ((5+data.prof.zeff)./ (1+data.prof.zeff))) .* 6 ./ (1 + 4 .*(1- data.equi.ftrap) + data.prof.zeff);
jeccd = 2 .* pi .*  data.source.fce.el .* etaecem ./ data.prof.ne;
ieccd = data.equi.rhomax .* trapz(param.gene.x,jeccd .* data.equi.spr,2);

figure
plot(data.gene.temps,data.gene.ifce,data.gene.temps,ieccd);
figure
zplotprof(gca,data.gene.temps,param.gene.x,data.source.fce.j)
zplotprof(gca,data.gene.temps,param.gene.x,jeccd,'color','r')

