% calcul 0d
pelc = data.equi.rhomax .* trapz(param.gene.x,data.equi.vpr .* data.source.fci.el,2);
pionc = data.equi.rhomax .* trapz(param.gene.x,data.equi.vpr .* data.source.fci.ion,2);
cmin  = mean(trapz(param.gene.x,data.equi.vpr .* squeeze(data.impur.impur(:,:,2)),2) ./ trapz(param.gene.x,data.equi.vpr .* squeeze(data.impur.impur(:,:,1)),2));
zs.picrh = pelc + pionc;
cons = z0dinput.cons;
cons.picrh = pelc +pionc;
option = z0dinput.option;
option.transitoire =0;
option.rip =0;
option.cmin = cmin;
cmin = input(sprintf('cmin [%g]? ',option.cmin));
if ~isempty(cmin)
   option.cmin = cmin;
end

[pth,pel,pion,esupra,einj,egamma,teff0,taus] = z0icrh(zs,z0dinput.geo,cons,option);

esup_icrh = (3/2) .* data.equi.rhomax .*trapz(param.gene.x,data.source.fci.psupra .* vpr,2);

figure(19);
subplot(3,1,1);
plot(zs.temps,pel,'r',zs.temps,pion,'b',data.gene.temps,pelc,'m',data.gene.temps,pionc,'c',zs.temps,zs.picrh,'k',zs.temps,pth,'g');
subplot(3,1,2);
plot(zs.temps,esupra,'r',data.gene.temps,esup_icrh,'m',zs.temps,einj,'b:',zs.temps,egamma,'g:');
subplot(3,1,3);
plot(zs.temps,taus,'r',data.gene.temps,esup_icrh ./ pelc .* 2,'m')
