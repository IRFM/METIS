% script de test de zsauter0d vs nclass dans cronos
tsscenario
[temps,void] = ginput(1);
indt = min(find(temps <= data.gene.temps));
datak = zget1t(data,indt);
x    = param.gene.x;
tep  = datak.prof.te;
tip  = datak.prof.ti;
nep  = datak.prof.ne;
nip  = datak.prof.ni;
Raxe = datak.equi.raxe;
ftrap = datak.equi.ftrap;
fdia = datak.equi.F;
dpsidx = datak.prof.psid1;
epsi   = datak.equi.a ./ datak.equi.raxe;
qp     = datak.equi.q;
zeffp  = datak.prof.zeff;
zion   = 1;
modeh  = 1;
jboot          = zsauter0d(x,tep,tip,nep,nip,qp,zeffp,zion,Raxe,ftrap,epsi,dpsidx,fdia,modeh) ./ datak.geo.b0;
jboot_koh  = zsauter0d_Koh(x,tep,tip,nep,nip,qp,zeffp,zion,Raxe,ftrap,epsi,dpsidx,fdia,modeh) ./ datak.geo.b0;
eta    = zeta0(tep,nep,qp,zeffp,datak.equi.a,Raxe,ftrap,epsi);

[nustare,nustari,L31,L32,L34,alpha,Isau,eta2,etasptz,IsauZi]=sauter(data,param);
[cr,neo,memoire] = zneoclass(param.fonction.neo,param.cons.neo,param.gene,param.phys,param.compo,datak,datak.neo,0,param.memoire.impur.etatcharge ,param.memoire.neo);

figure
subplot(2,1,1)
semilogy(x,neo.eta,'r',x,eta,'ob',x,eta2(indt,:),'m+',x,etasptz(indt,:),'cx');
legend('NClass','METIS','Sauter CRONOS','Spitzer');
subplot(2,1,2)
plot(x,neo.jboot,'r',x,jboot,'ob',x,Isau(indt,:),'m+',x,IsauZi(indt,:),'cx',x,jboot_koh,'g*');
legend('NClass','METIS','Sauter CRONOS Zion','Sauter CRONOS Zeff','Sauter modified Koh');

