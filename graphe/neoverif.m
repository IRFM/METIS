% script de verification des donnees neoclassique
num =param.from.shot.num;
if fix(num) == num
  num = num + 0.1;
end

if ~exist('boot','var');
    [boot,s] = cgcgettrait(num,'tbootstrap');
end
if ~exist('bile','var');
     [bile,s] = cgcgettrait(num,'tprof');
end

test =data.equi.rmoy;

figure('color',[1 1 1]);
warning off
t=data.gene.temps;
x=data.equi.a ./ (data.equi.a(:,end)*ones(1,size(data.equi.a,2)));

subplot(3,3,1)
plotprof(gca,boot.times,boot.rhofit,boot.j1,'linestyle','+');
plotprof(gca,boot.times,boot.rhofit,boot.j2,'linestyle','x');
plotprof(gca,t,x,data.neo.jboot,'linestyle','-','color','r');
xlabel('r/a')
ylabel('Jboot (A)')
hl=legend('Ti = a Te','Ti = a Ne','Cronos');
hs=subplot(3,3,2);
set(hl,'position',get(hs,'position'));
delete(hs);

subplot(3,3,3)
plotprof(gca,boot.times,boot.rhofit,1./boot.resp1,'linestyle','+');
plotprof(gca,boot.times,boot.rhofit,1./boot.resp2,'linestyle','x');
plotprof(gca,t,x,1./data.neo.eta,'linestyle','-','color','r');
xlabel('r/a')
ylabel('eta (Ohm^-1 * m^-1)')
%legend('Ti = a Te','Ti = a Ne','Cronos')

subplot(3,3,4)
plotprof(gca,boot.times,boot.rhofit,-boot.qe1,'linestyle','+');
plotprof(gca,boot.times,boot.rhofit,-boot.qe2,'linestyle','x');
plotprof(gca,t,x,data.neo.flux.qe.*test,'linestyle','-','color','r');
xlabel('r/a')
%legend('Ti = a Te','Ti = a Ne','Cronos')
ylabel('Qe W/m^2')

subplot(3,3,5)
plotprof(gca,boot.times,boot.rhofit,-boot.qi1,'linestyle','+');
plotprof(gca,boot.times,boot.rhofit,-boot.qi2,'linestyle','x');
plotprof(gca,t,x,data.neo.flux.qion.*test,'linestyle','-','color','r');
xlabel('r/a')
ylabel('Qi W/m^2')
%legend('Ti = a Te','Ti = a Ne','Cronos')

subplot(3,3,6)
plotprof(gca,boot.times,boot.rhofit,-boot.fe1,'linestyle','+');
plotprof(gca,boot.times,boot.rhofit,-boot.fe2,'linestyle','x');
plotprof(gca,t,x,data.neo.flux.ne.*test,'linestyle','-','color','r');
xlabel('r/a')
ylabel('Ge m^-2')
%legend('Ti = a Te','Ti = a Ne','Cronos')

subplot(3,3,7)
plotprof(gca,boot.times,boot.rhofit,-boot.fi1,'linestyle','+');
plotprof(gca,boot.times,boot.rhofit,-boot.fi2,'linestyle','x');
plotprof(gca,t,x,data.neo.flux.nion.*test,'linestyle','-','color','r');
xlabel('r/a')
ylabel('Gi m^-2')
%legend('Ti = a Te','Ti = a Ne','Cronos')

subplot(3,3,8)
plotprof(gca,bile.times,bile.rhofit,bile.pel,'linestyle','o','color','r');
plotprof(gca,bile.times,bile.rhofit,bile.pion,'linestyle','+','color','b');
plotprof(gca,t,x,data.prof.pe,'linestyle','-','color','r');
plotprof(gca,t,x,data.prof.pion,'linestyle','-.','color','b');
xlabel('r/a')
ylabel('Pa')
hl=legend('Pe prof','Pion prof','Pe cronos','Pion cronos');
hs=subplot(3,3,9);
set(hl,'position',get(hs,'position'));
delete(hs);
