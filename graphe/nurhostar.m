function nurhostar(data,param,jeux1,temps,option)
%constantes
c=param.phys.c;qe=param.phys.e;eps0=param.phys.epsi0;depi=2*pi;
me=param.phys.me;mp=param.phys.mp;muo=param.phys.mu0;
comp = 0;
if isfield(jeux1,'data')
  x1      = jeux1.param.gene.x;
  choc1   = jeux1.param.from.shot.num;
  times1  = jeux1.data.gene.temps;
  nl1     = 2*jeux1.data.gene.nbar.*jeux1.data.geo.a;
  ploss1  = jeux1.data.gene.ploss;
  betap1  = jeux1.data.gene.betap;
  wdia1   = jeux1.data.gene.wdia;
  tefit1  = jeux1.data.prof.te;
  Ai1     = jeux1.param.compo.a(1);
  zi1     = jeux1.param.compo.z(1);
  mi1     = Ai1*mp;
  R1      = jeux1.data.geo.r0*ones(size(x1));
  rmin1   = jeux1.data.geo.a*ones(size(x1));
  Bo1     = jeux1.data.geo.b0*ones(size(x1));
  z1      = jeux1.data.prof.zeff;
  ray1    = jeux1.data.equi.rhomax*x1;
  tex1    = abs(jeux1.data.prof.te);
  nex1    = jeux1.data.prof.ne;
  tix1    = jeux1.data.prof.ti;
  nix1    = jeux1.data.prof.ni;
  Betax1  = 2*muo*(jeux1.data.prof.pion+jeux1.data.prof.pe)./(Bo1.*Bo1);
  epsilon1= ray1./R1;
  lnc1    = 31.474 + log(nex1.^(-.5).*tex1);
  ftrap1  = jeux1.data.equi.ftrap;

%determination de nue* et rhostar

  vte1    = ((2*qe*tex1)/me).^(.5);
  tei1    = 3*eps0*eps0*depi^1.5*sqrt(me)/qe^4*((qe*tex1).^(3/2))./lnc1./nex1;
  wpe1    = sqrt(qe^2*nex1/eps0/me);
  wce1    = qe*Bo1./me;
  wci1    = qe*Bo1./mi1;

  rhoe1   = vte1./wce1;
  rhostare1=rhoe1./rmin1;
  q1      = jeux1.data.prof.q;
  nustar1 = sqrt(2)*R1.*q1./vte1./tei1./epsilon1.^(1.5);
  nustar1(1)=nustar1(2);

  rhoi1   = (2*qe*tix1*mi1).^0.5./(qe*zi1.*Bo1);
  rhostari1= rhoi1./rmin1;

  tps1    = temps;
  itps1   = iround(times1,tps1);
  kx1     = 5:length(x1)-5;
  comp    = 1;
end

choc   = param.from.shot.num;
times  = data.gene.temps;
nl     = 2*data.gene.nbar.*data.geo.a;
ploss  = data.gene.ploss;
betap  = data.gene.betap;
wdia   = data.gene.wdia;
tefit  = data.prof.te;

if 1>2
figure
subplot(511)
plot(times,nl)
ylabel('nl')
title(['# ', int2str(choc)])

subplot(512)
plot(times,ploss)
ylabel('ptot')

subplot(513)
plot(times,betap)
ylabel('betap')

subplot(514)
plot(times,wdia./ploss)
ylabel('tauE')

subplot(515)
plot(times,tefit(:,1))
ylabel('Te0')
end

qen=qe*sqrt(1/(muo*sqrt(me)));
Ai=param.compo.a(1);
zi=param.compo.z(1);
mi=Ai*mp;

x=param.gene.x;
R=data.geo.r0*ones(size(x));
rmin=data.geo.a*ones(size(x));
Bo=data.geo.b0*ones(size(x));
z=data.prof.zeff;
ray=data.equi.rhomax*x;
tex=abs(data.prof.te);nex=data.prof.ne;
tix=data.prof.ti;
nix=data.prof.ni;
Betax=2*muo*(data.prof.pion+data.prof.pe)./(Bo.*Bo);
epsilon=ray./R;

lnc=31.474 + log(nex.^(-.5).*tex);
ftrap=data.equi.ftrap;

%determination de nue* et rhostar

vte=((2*qe*tex)/me).^(.5);
tei=3*eps0*eps0*depi^1.5*sqrt(me)/qe^4*((qe*tex).^(3/2))./lnc./nex;
wpe=sqrt(qe^2*nex/eps0/me);
wce=qe*Bo./me; wci = qe*Bo./mi;

rhoe=vte./wce;
rhostare=rhoe./rmin;
q = data.prof.q;
nustar=sqrt(2)*R.*q./vte./tei./epsilon.^(1.5);
nustar(1)=nustar(2);

rhoi=(2*qe*tix*mi).^0.5./(qe*zi.*Bo);
rhostari=rhoi./rmin;

tps = temps;
itps=iround(times,tps);
kx=5:length(x)-5;


subplot(221)
plot(x(kx),q(itps,kx),'r')
if comp == 1
  hold on
  if ~isnan(itps1)
    plot(x1(kx1),q1(itps1,kx1),'r--')
  end
  hold off
end
ylabel('q')
xlabel('r/a')
if comp == 0
  title(['# ', int2str(choc),  '  @ t = ', num2str(tps), ' B(T) = ', num2str(Bo(itps,1),2)])
else
  title(['- : # ', int2str(choc),' -- : #',  int2str(choc1), '  @ t = ', num2str(tps)])
  text1=['-  B(T) = ', num2str(Bo(itps,1),2)];
  text2=['-- B(T) = ', num2str(Bo1(itps,1),2)];
  text('units','normalized','position',[0.6 0.9],'string',text1)
  text('units','normalized','position',[0.6 0.8],'string',text2)
end
grid

subplot(222)

plot(x(kx),rhostari(itps,kx),'r')
if comp == 1
  if ~isnan(itps1)
    hold on
    plot(x1(kx1),rhostari1(itps1,kx1),'r--')
    hold off
  end
end
ylabel('rhostari')
xlabel('r/a')
grid

subplot(223)

plot(x(kx),nustar(itps,kx),'r',x(kx),rhostare(itps,kx)*1e4,'b')
if comp == 1
  if ~isnan(itps1)
    hold on
    plot(x1(kx1),nustar(itps1,kx1),'r--',x1(kx1),rhostare(itps1,kx1)*1e4,'b--')
    hold off
  end
end
xlabel('r/a')
ylabel('nustar & rhostre (1e-4)')
grid

subplot(224)

plot(x(kx),Betax(itps,kx),'r')
if comp == 1
  if ~isnan(itps1)
    hold on
    plot(x1(kx1),Betax1(itps1,kx1),'r--')
    hold off
  end
end
xlabel('r/a')
ylabel('betat')
grid
