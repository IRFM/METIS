% script pour le test de z0icrh
picrh_ion = data.equi.rhomax .* trapz(param.gene.x,data.equi.vpr .* data.source.fci.ion,2);
picrh_el  = data.equi.rhomax .* trapz(param.gene.x,data.equi.vpr .* data.source.fci.el,2);
esupra_icrh = data.equi.rhomax .* trapz(param.gene.x,data.equi.vpr .* ((3/2) .* data.source.fci.psupra+data.source.fci.paniso),2);
indh        = find(param.compo.a == 1 & param.compo.z == 1)
nmino       =  data.equi.rhomax .* trapz(param.gene.x,data.equi.vpr .* squeeze(data.impur.impur(:,:,indh)),2);
nd          =  data.equi.rhomax .* trapz(param.gene.x,data.equi.vpr .* squeeze(data.impur.impur(:,:,1)),2);

%[zs.picrh_th,zs.pel_icrh,zs.pion_icrh,zs.esup_icrh,zs.einj_icrh,zs.ecrit_icrh,teff0,zs.taus_icrh,zs.frloss_icrh] = z0icrh(zs,geo,cons,option);
zs   = post.zerod;
cons = post.z0dinput.cons;
geo  = post.z0dinput.geo;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;

figure;
subplot(2,2,1)
plot(t,picrh_ion,'r',t,zs.pion_icrh,'b',t,picrh_el,'m',t,zs.pel_icrh,'c')
subplot(2,2,2)
plot(t,picrh_el+picrh_ion,'r',t,zs.pel_icrh+zs.pion_icrh,'b',t,cons.picrh,'k')
subplot(2,2,3)
plot(t,esupra_icrh,'r',t,zs.esup_icrh,'b')
subplot(2,2,4)
plot(t,nmino,'r',t,zs.vp .* zs.nDm .* post.z0dinput.option.cmin,'b',t,nd,'m',t,zs.vp .* zs.nDm ,'c')





