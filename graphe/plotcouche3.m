function rcouche=plotcouche3(data,param,temps,nant_char)
%  [rcouche]=plotcouche2(data,param,temps,nant)
%  position des couches fci H, D et He3
%  temps = temps d'analyse
%  nant  = num�o de l'antenne
%  rcouche(harmonique,minoritaire)
%
equi=data.equi;
nant  = str2num(nant_char(6));
freq  = param.cons.fci.frequence;

ind    = iround(data.gene.temps,temps);
frmin  = squeeze(data.impur.impur(ind,1,2)) / squeeze(data.impur.impur(ind,1,1)) * 100; 
frboot = data.gene.iboot(ind) / data.gene.ip(ind) * 100; 
ne0    = data.prof.ne(ind,1)/1e19;
a      = data.geo.a(ind);
ip     = data.gene.ip(ind)/1e6;
R0     = data.geo.r0(ind);
R1     = squeeze(double(data.equi.R(ind,:,:)));
Z      = squeeze(double(data.equi.Z(ind,:,:)));
bR     = squeeze(double(data.equi.BR(ind,:,:)));
bZ     = squeeze(double(data.equi.BZ(ind,:,:)));
bphi   = squeeze(double(data.equi.BPHI(ind,:,:)));
green  = data.gene.nbar(ind) / (data.gene.ip(ind)/pi/data.geo.a(ind)^2/1e6*10) / 1e19 * 100;
bp     = data.gene.betap(ind);
bN     = data.gene.beta(ind) /(a*data.geo.b0(ind)/(data.gene.ip(ind)/1e6));

%
% effet du ripple
%
if strcmp(param.from.machine,'TS')
  Ritor         = data.geo.r0;
  itor          = Ritor .* data.geo.b0 ./ (param.phys.mu0 .* 18 .* 2028  ./ 2 ./ pi);
  itor          = itor(ind);
  [btr,brr,bzr] = bcronoscall(R1,Z,10/180*pi,itor);
  btot          = sqrt((bphi+btr).^2 + (bR+brr).^2 + (bZ+bzr).^2);
  [btr,brr,bzr] = bcronoscall(R1,Z,0,itor);
  btotplus      = sqrt((bphi+btr).^2 + (bR+brr).^2 + (bZ+bzr).^2);
else
  btot          = sqrt(bphi.^2 + bR.^2 + bZ.^2);
  btotplus      = btot;
end
bmag  = data.geo.b0(ind);
axe   = R1(1,1);
B0    = btot(1,1);
x     = param.gene.x;
Te    = data.prof.te(ind,:);
ne    = data.prof.ne(ind,:);
Ti    = data.prof.ti(ind,:);
ni    = data.prof.ni(ind,:);
Pcons = sum(abs(data.cons.fci(ind,:)))/1e6;
Pel   = data.source.fci.el(ind,:);
Pion  = data.source.fci.ion(ind,:);
Psup  = data.source.fci.psupra(ind,:);
convm = real(data.source.fci.err(ind));
renorm = abs(1-imag(data.source.fci.err(ind)));
totion =  zintvol(Pion,x,data.equi.vpr(ind,:),data.equi.rhomax(ind))/1e6;
totel  =  zintvol(Pel,x,data.equi.vpr(ind,:),data.equi.rhomax(ind))/1e6;
wth    = data.gene.wth(ind,:)/1e6;
totsup =  zintvol(Psup,x,data.equi.vpr(ind,:),data.equi.rhomax(ind))/1e6;
Pabs   = data.gene.paddfci/1e6;
nmoy   = data.gene.nemoy(ind)/1e19;
nbar   = data.gene.nbar(ind)/1e19;
nl     = nbar*2*a;
zeff   = data.gene.zeffm(ind)*2;
ptot   = data.gene.paddtot(ind)/1e6;
ploss  = data.gene.ploss(ind);
if totsup ~= 0
  wdia   = data.gene.wdia(ind)*(1-totsup/(totion+totel))/1e6;
else
  wdia   = data.gene.wdia(ind)/1e6;
end
dwth   = gradient(data.gene.wth,data.gene.temps);
elong  = data.geo.e1(ind);
M      = param.compo.a(1);
%Wrlw   = zrlw(ip,R0,a,nmoy,bmag,zeff,ptot-dwth(ind)/1e6,elong)*1e6;
%Hrlw   = data.gene.we(ind) / Wrlw;
if strcmp(param.from.machine,'TS')
  fit(1)      = scallaw('ts',ip,R0,a,elong,nl,bmag,1,ptot-dwth(ind)/1e6)*ploss/1e6/wdia;
else
  M           = param.compo.a(1);
  fit(1)      = IPB98y2(ip,bmag,nbar,ptot-dwth(ind)/1e6,R0,elong,a,M);
  fit(2)      = H97(ip,bmag,nbar,ptot-dwth(ind)/1e6,R0,elong,a,M);
end
%
%
%
mp                   = 1.6726e-27;
el                   = 1.6022e-19;
me                   = 9.1095e-31;
%
% champ B
%
btot(1,1) = btot(1,2);
btotplus(1,1) = btotplus(1,2);
B      = btot;
fp     = el*B/2/pi/mp/1e6;
fd     = el*B/2/pi/mp/2/1e6;
fhe    = 2*el*B/2/pi/mp/3/1e6;
ft     = el*B/2/pi/mp/3/1e6;
if strcmp(param.from.machine,'TS')
  fpplus     = el*btotplus/2/pi/mp/1e6;
  fdplus     = el*btotplus/2/pi/mp/2/1e6;
  fheplus    = 2*el*btotplus/2/pi/mp/3/1e6;
  ftplus    = el*btotplus/2/pi/mp/3/1e6;
end
theta  = linspace(0,2*pi,200);

subplot(2,2,1)
plot(R1(1:10:end,:)',Z(1:10:end,:)','r')


axis('equal');
hold on
xlabel('R (m)');
ylabel('Z (m)');


f    = freq(nant);
nh   = 1;
nd   = 1;
nhe  = 1;
nt   = 1;
tex1 = [' '];
tex2 = [' '];
tex3 = [' '];
tex4 = [' '];

rcouche(4,4) = 0;
for k=1:4

  cc     = contour(R1,Z,fp,[f/k f/k]);
  if ~isempty(cc)
     ind         = find(cc(1,:) == f/k);
	  [maxcc,jnd] = max(cc(2,ind));
	  deb         = ind(jnd)+1;
     Rc  = cc(1,deb:cc(2,ind(jnd))+deb-1);
	  Zc  = cc(2,deb:cc(2,ind(jnd))+deb-1);
     if k == 1
       tex1=[' 1H = ',num2str(mean(Rc),3),' m'];
     end
     if k == 2
       tex2=[' 2H = ',num2str(mean(Rc),3),' m'];
     end
	  knd           = find(Zc>0);
	  knd           = knd(1);
     rcouche(k,1) = Rc(knd);
     if strcmp(param.from.machine,'TS')
	    cc     = contour(R1,Z,fpplus,[f/k f/k]);
	  end
     text('units','data','position',[rcouche(k,1) a-0.05],'string',[int2str(k),' H'])
     nh = nh+1;
   end

  cc     = contour(R1,Z,fd,[f/k f/k]);
  if ~isempty(cc)
     ind         = find(cc(1,:) == f/k);
	  [maxcc,jnd] = max(cc(2,ind));
	  deb         = ind(jnd)+1;
     Rc  = cc(1,deb:cc(2,ind(jnd))+deb-1);
	  Zc  = cc(2,deb:cc(2,ind(jnd))+deb-1);
     if k == 3
       tex3=[' 3D = ',num2str(mean(Rc),3),' m'];
     end
	  knd           = find(Zc>0);
	  knd           = knd(1);
     rcouche(k,2)  = Rc(knd);
     if strcmp(param.from.machine,'TS')
	    cc     = contour(R1,Z,fdplus,[f/k f/k]);
	  end
     nd = nd+1;
   end
  clear cc

  cc     = contour(R1,Z,fhe,[f/k f/k]);
  if ~isempty(cc)
     ind         = find(cc(1,:) == f/k);
     [maxcc,jnd] = max(cc(2,ind));
     if maxcc > 2
       deb         = ind(jnd)+1;
       Rc  = cc(1,deb:cc(2,ind(jnd))+deb-1);
       Zc  = cc(2,deb:cc(2,ind(jnd))+deb-1);
       if k == 1
         tex4=[' 1He3 = ',num2str(mean(Rc),3),' m'];
       end
       knd            = find(Zc>0);
       knd            = knd(1);
       rcouche(k,3)   = Rc(knd);
       if strcmp(param.from.machine,'TS')
	  cc     = contour(R1,Z,fheplus,[f/k f/k]);
       end
       if k < 3 & strcmp(param.from.machine,'ITER')==0
         text('units','data','position',[rcouche(k,3) a-.25*k],'string',[int2str(k),' He3'])
       end
       nhe = nhe+1;
     end
  end
  cc     = contour(R1,Z,ft,[f/k f/k]);
  if ~isempty(cc)
     ind         = find(cc(1,:) == f/k);
	  [maxcc,jnd] = max(cc(2,ind));
	  if maxcc > 2
	    deb         = ind(jnd)+1;
       Rc  = cc(1,deb:cc(2,ind(jnd))+deb-1);
	    Zc  = cc(2,deb:cc(2,ind(jnd))+deb-1);
       if k == 2
         tex4=[' 2T = ',num2str(mean(Rc),3),' m'];
       end
	    knd            = find(Zc>0);
	    knd            = knd(1);
       rcouche(k,4)   = Rc(knd);
       if strcmp(param.from.machine,'TS')
	      cc     = contour(R1,Z,fheplus,[f/k f/k]);
	    end
       if k < 3 & strcmp(param.from.machine,'ITER')
         text('units','data','position',[rcouche(k,3) a-.25*k],'string',[int2str(k),' T'])
       end
       nt = nt+1;
	  end

   end
      end
text('units','normalized','position',[0.1 0.95],'string',tex1);
text('units','normalized','position',[0.4 0.95],'string',tex2);
text('units','normalized','position',[0.7 0.95],'string',tex3);
text('units','normalized','position',[1 0.95],'string',tex4);

title(['ICRH layer, ant. ',int2str(nant),', temps =',num2str(temps,4),' s'])
plot(axe,0,'r*')

axis([min(R1(:)) max(R1(:))+0.1 -a a])
hold off
%
%
%
subplot(2,2,2)
Pfci = abs(data.cons.fci)/1e6;
if size(Pfci,2) == 3
plot(data.gene.temps,Pabs,'ro',data.gene.temps,Pfci(:,1),...
     data.gene.temps,Pfci(:,2),...
     data.gene.temps,Pfci(:,3))
     legend('Pabs',['ant 1 : ' int2str(freq(1)) ' MHz'],['ant 2 : ' int2str(freq(2)) ' MHz'],['ant 3 : ' int2str(freq(3)) ' MHz'])
end
if size(Pfci,2) == 4
plot(data.gene.temps,Pabs,'ro',data.gene.temps,Pfci(:,1),...
     data.gene.temps,Pfci(:,2),...
     data.gene.temps,Pfci(:,3),...
     data.gene.temps,Pfci(:,3))
     leg1=['ant 1 : ' int2str(freq(1)) ' MHz'];
     leg2=['ant 1 : ' int2str(freq(2)) ' MHz'];
     leg3=['ant 1 : ' int2str(freq(3)) ' MHz'];
     leg4=['ant 1 : ' int2str(freq(4)) ' MHz'];

 legend('Pabs',leg1,leg2,leg3,leg4)
end
hold on
Pmax = max(max(Pfci(:),max(Pabs)));
v=axis;
axis([v(1) v(2) 0 Pmax*2])
plot([temps temps],[0 v(4)],'--')
hold off
ylabel('MW')
xlabel('time')
if strcmp(param.from.machine,'TS')
%  title(['P. / ant., Hrlw =',num2str(Hrlw,2),' Hscalts =',num2str(fit(1),2)])
  title(['P. / ant.,  Hscalts =',num2str(fit(1),2)])
else
  title(['P. / ant.'])
end
grid

subplot(2,2,3)
indc=find(rcouche > 0);
nat = ['H' 'H' 'H' 'H';'D' 'D' 'D' 'D';'X' 'X' 'X' 'X';'T' 'T' 'T' 'T']';
for k=1:length(indc)
  if rcouche(indc(k))  > axe
    Rcor = R1(:,1);
  else
    Rcor = R1(:,34);
  end
  pos(k) = interp1(Rcor,x,rcouche(indc(k)),'nearest');
  ch(k)  = rem(indc(k),4);
  minvb(k) = nat(indc(k));
  if strcmp(minvb(k), 'X')
    mino(k,:) = ['He3'];
  else
    mino(k,:) = [minvb(k) '  '];
  end

end
plot(x,Te,x,Ti,[pos;pos],[0*pos;ones(size(pos))*max(Te)])
legend('Te','Ti',2)
var1H=0
if rcouche(1) ~=0 & rcouche(6) ~=0

  var1H=1;

end
ylabel('eV')
for k=1:length(indc)
  if indc(k) ~= 6 & var1H == 1
    text('units','data','string',[int2str(ch(k)) mino(k,:),'= ',num2str(pos(k),3)],'position',[pos(k)*1.01 max(Te)*0.9])
  end
end

%ylabel('eV')
%for k=1:length(indc)
%text('units','data','string',[int2str(ch(k)) mino(k,:),'= ',num2str(pos(k),3)],'position',[pos(k)*1.01 max(Te)*0.9])
%end
if strcmp(param.cons.fci.rip,'Yes')
  title(['B0=',num2str(bmag,2),' T, R0=',num2str(R0,2),' m, a=',num2str(a,2),' m, + ripple effect'])
else
  title(['B0=',num2str(bmag,2),' T, R0=',num2str(R0,2),' m, a=',num2str(a,2),' m, without ripple'])
end
subplot(2,2,4)

plot(x,Pel/1e6,x,Pion/1e6,x,(Pel+Pion)/1e6,x,Psup/1e6,'--',[pos;pos],[0*pos;ones(size(pos))*max(Pel/1e6)])
for k=1:length(indc)
text('units','data','string',[int2str(ch(k)) mino(k,:)],'position',[pos(k)*1.01 max(Pel/1e6)*0.9])
end
ylabel('MW/m^3')
legend('el','ion','tot','sup (MJ)')
text('units','normalized','position',[0.6 0.7],'string',['Pel       = ',num2str(totel,2),' MW'])
text('units','normalized','position',[0.6 0.6],'string',['Pion      = ',num2str(totion,2),' MW'])
fracion = totion / (totel + totion) * 100;
fracsup = totsup / wth * 100;
fracrip = 100 - (totel+totion) / Pcons * 100;
text('units','normalized','position',[0.6 0.5],'string',['ion fraction     = ',num2str(fracion,2),' %'])
text('units','normalized','position',[0.6 0.4],'string',['supra fraction     = ',num2str(fracsup,2),' %'])
if strcmp(param.cons.fci.rip,'Yes')
  text('units','normalized','position',[0.6 0.3],'string',['ripple fraction     = ',num2str(fracrip,2),' %'])
end
text('units','normalized','position',[0.6 0.2],'string',['mode convertion = ',num2str(convm*100,2),' %'])
text('units','normalized','position',[0.6 0.1],'string',['deviation pion/cons     = ',num2str(renorm*100,2),' %'])
%text('units','normalized','position',[0.6 0.01],'string',['bN       = ',num2str(bN,2)])

if rcouche(1,1) < 2.9 & rcouche(1,1) > 1.7
  title(['nh/nd=',num2str(frmin,3),' %, frboot=',num2str(frboot,2),' %, ne(0)=',num2str(ne0),' 10^1^9 m^-^3','bp=',num2str(bp,2),' gw=',num2str(green,3),' %'])
end
if rcouche(1,3) < 3.0 & rcouche(1,3) > 1.9
  title(['nhe3/nd=',num2str(frmin,3),' %, frboot=',num2str(frboot,2),' %, ne(0)=',num2str(ne0),' 10^1^9 m^-^3','bp=',num2str(bp,2),' gw=',num2str(green,3),' %'])
end

% integralle surfacique (d'une grandeur independante de theta)
%  s = integrale de surface
%  e = valeur a integree
%  x = coordonnees normalisee
%  sp = datak.equi.sp
%  rhomax = datak.equi.rhomax
function s=zintsurf(e,x,sp,rhomax)

    s = rhomax .* trapz(x,sp .* e,2);
% integralle volumique
%  s = integrale de volume
%  e = valeur a integree
%  x = coordonnees normalisee
%  vpr = datak.equi.vpr
%  rhomax = datak.equi.rhomax
function s=zintvol(e,x,vpr,rhomax)

  s = rhomax.*trapz(x,vpr .* e,2);
