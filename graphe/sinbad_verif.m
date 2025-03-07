% plot de comparaison  sinbad /transp
strchoc = int2str(param.from.shot.num);
load(strcat('/usr/drfc/cgc/matlab5/tcron/JET/',strchoc,'/transp',strchoc));
load(strcat('/usr/drfc/cgc/matlab5/tcron/JET/',strchoc,'/efit',strchoc));
transrho;
t=data.gene.temps;
x=param.gene.x;
figure
plotprof(gca,t,x,data.source.idn.j,'color','r');
plotprof(gca,ttransp,xefit,jnbtrx,'color','c')
