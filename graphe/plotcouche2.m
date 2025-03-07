function rcouche=plotcouche2(data,param,temps,nant)
%  [rcouche]=plotcouche2(data,param,temps,nant)
%  position des couches fci H, D et He3
%  temps = temps d'analyse
%  nant  = num�ro de l'antenne
%  rcouche(harmonique,minoritaire)
%


freq  = param.cons.fci.frequence;
if nargin < 4
  if length(freq) == 3
    tex1 = ['antenne 1, f=',num2str(freq(1),3),' MHz [defaut]'];
    tex2 = ['antenne 2, f=',num2str(freq(2),3),' MHz'];
    tex3 = ['antenne 3, f=',num2str(freq(3),3),' MHz'];
    nant = menuTS('num�ro de l''antenne FCI',tex1,tex2,tex3);
  else
    tex1 = ['antenne 1, f=',num2str(freq(1),3),' MHz [defaut]'];
    tex2 = ['antenne 2, f=',num2str(freq(2),3),' MHz'];
    tex3 = ['antenne 3, f=',num2str(freq(3),3),' MHz'];
    tex4 = ['antenne 4, f=',num2str(freq(4),3),' MHz'];
    nant = menuTS('num�ro de l''antenne FCI',tex1,tex2,tex3,tex4);
  end
  if isempty(nant)
    nant = 1;
  end
end
if nargin < 3
  temps = 10;
  tex0  = [num2str(param.gene.tdeb,3),' s < t < ',num2str(param.gene.tfin,3),' s'];
  temps = inputd(['temps d''analyse [ ',tex0,' ]'],temps);
  if temps < param.gene.tdeb
    temps = param.gene.tdeb;
  end
  if temps > param.gene.tfin
    temps = param.gene.tfin;
  end
end
ind    = iround(data.gene.temps,temps);
frmin  = squeeze(data.impur.impur(ind,1,2)) / squeeze(data.impur.impur(ind,1,1)) * 100; 
frboot = data.gene.iboot(ind) / data.gene.ip(ind) * 100; 
ne0    = data.prof.ne(ind,1)/1e19;
a      = data.geo.a(ind);
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
  [btr,brr,bzr] = bmex(R1,Z,10/180*pi,itor);
  btot          = sqrt((bphi+btr').^2 + (bR+brr').^2 + (bZ+bzr').^2);
  [btr,brr,bzr] = bmex(R1,Z,0,itor);
  btotplus      = sqrt((bphi+btr').^2 + (bR+brr').^2 + (bZ+bzr').^2);
else
  btot          = sqrt(bphi.^2 + bR.^2 + bZ.^2);
end
bmag  = data.geo.b0(ind);
axe   = R1(1,1);
B0    = btot(1,1);
x     = param.gene.x;
Te    = data.prof.te(ind,:);
Ti    = data.prof.ti(ind,:);
Pel   = data.source.fci.el(ind,:);
Pion  = data.source.fci.ion(ind,:);
Psup  = data.source.fci.psupra(ind,:);
totion =  zintvol(Pion,x,data.equi.vpr(ind,:),data.equi.rhomax(ind))/1e6;
totel  =  zintvol(Pel,x,data.equi.vpr(ind,:),data.equi.rhomax(ind))/1e6;
totsup =  zintvol(Psup,x,data.equi.vpr(ind,:),data.equi.rhomax(ind))/1e6;
Pabs   = data.gene.paddfci/1e6;

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
if strcmp(param.from.machine,'TS')
  fpplus     = el*btotplus/2/pi/mp/1e6;
  fdplus     = el*btotplus/2/pi/mp/2/1e6;
  fheplus    = 2*el*btotplus/2/pi/mp/3/1e6;
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
tex1 = [' '];
tex2 = [' '];
tex3 = [' '];
tex4 = [' '];

rcouche(4,3) = 0;
for k=1:4
 
  cc     = contour(R1,Z,fp,[-f/k f/k]);
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
     rcouche(nh,1) = Rc(knd);
     if strcmp(param.from.machine,'TS')
	    cc     = contour(R1,Z,fpplus,[-f/k f/k]);
	  end
     text('units','data','position',[rcouche(nh,1) a-0.05],'string',[int2str(k),' H'])
     nh = nh+1;
   end

  cc     = contour(R1,Z,fd,[-f/k f/k]);
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
     rcouche(nd,2) = Rc(knd);
     if strcmp(param.from.machine,'TS')
	    cc     = contour(R1,Z,fdplus,[-f/k f/k]);
	  end
     nd = nd+1;
   end
  clear cc
  cc     = contour(R1,Z,fhe,[-f/k f/k]);
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
       rcouche(nhe,3) = Rc(knd);
       if strcmp(param.from.machine,'TS')
	      cc     = contour(R1,Z,fheplus,[-f/k f/k]);
	    end
       if k < 3
         text('units','data','position',[rcouche(nhe,3) a-.25*nhe],'string',[int2str(k),' He3'])
       end
       nhe = nhe+1;
	  end
   end



    end

text('units','normalized','position',[0.1 0.95],'string',tex1);
text('units','normalized','position',[0.4 0.95],'string',tex2);
text('units','normalized','position',[0.7 0.95],'string',tex3);
text('units','normalized','position',[1 0.95],'string',tex4);

title(['couches FCI, ant. ',int2str(nant),', temps =',num2str(temps,4),' s'])
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
	 % legend('Pabs',['ant 1 : ' int2str(freq(1)) ' MHz'],['ant 2 : ' int2str(freq(2)) ' MHz'],['ant 3 : ' int2str(freq(3)) ' MHz'],['ant 4 : 'int2str(freq(4)) ' MHz'])
end	
hold on
if max(Pfci(:) > 0)
  Pmax = max(max(Pfci(:),max(Pabs)));
  v=axis;
  axis([v(1) v(2) 0 Pmax*2])
  plot([temps temps],[0 v(4)],'--')
end
hold off
ylabel('MW')
xlabel('temps')
title('Puissance / antenne')
grid
  
subplot(2,2,3)
indc=find(rcouche > 0);
nat = ['H' 'H' 'H' 'H';'D' 'D' 'D' 'D';'X' 'X' 'X' 'X']';
for k=1:length(indc)
  if rcouche(indc(k)  > axe)
    Rcor = R1(:,1);
  else
    Rcor = R1(:,34);
  end
  pos(k) = interp1(Rcor,x,rcouche(indc(k)),'nearest');
  ch(k)  = rem(indc(k),4);
  minf(k) = nat(indc(k));
  if strcmp(minf(k), 'X')
    mino(k,:) = ['He3'];
  else
    mino(k,:) = [minf(k) '  '];
  end
  
end
plot(x,Te,x,Ti,[pos;pos],[0*pos;ones(size(pos))*max(Te)])
legend('Te','Ti',2)
ylabel('eV')
for k=1:length(indc)
text('units','data','string',[int2str(ch(k)) mino(k,:),'= ',num2str(pos(k),3)],'position',[pos(k)*1.01 max(Te)*0.9])
end
title(['B0=',num2str(bmag,3),' T, R0=',num2str(R0,4),' m, a=',num2str(a,4),' m'])

subplot(2,2,4)

plot(x,Pel/1e6,x,Pion/1e6,x,(Pel+Pion)/1e6,[pos;pos],[0*pos;ones(size(pos))*max(Pel/1e6)])
for k=1:length(indc)
text('units','data','string',[int2str(ch(k)) mino(k,:)],'position',[pos(k)*1.01 max(Pel/1e6)*0.9])
end
ylabel('MW/m^3')
legend('el','ion','tot')
text('units','normalized','position',[0.6 0.6],'string',['Pel       = ',num2str(totel,2),' MW'])
text('units','normalized','position',[0.6 0.5],'string',['Pion      = ',num2str(totion,2),' MW'])
fracion = totion / (totel + totion) * 100;
fracsup = totsup / (totel + totion) * 100;
text('units','normalized','position',[0.6 0.4],'string',['frion     = ',num2str(fracion,2),' %'])
text('units','normalized','position',[0.6 0.3],'string',['frsup     = ',num2str(fracsup,2),' %'])
text('units','normalized','position',[0.6 0.2],'string',['greenwald = ',num2str(green,2),' %'])
text('units','normalized','position',[0.6 0.1],'string',['betap     = ',num2str(bp,2)])
%text('units','normalized','position',[0.6 0.01],'string',['bN       = ',num2str(bN,2)])

if rcouche(1,1) < 2.9 & rcouche(1,1) > 1.7
  title(['nh/nd=',num2str(frmin,3),' %, frboot=',num2str(frboot,2),' %, ne(0)=',num2str(ne0),' 10^1^9 m^-^3'])
end
if rcouche(1,3) < 2.9 & rcouche(1,3) > 1.7
  title(['nhe3/nd=',num2str(frmin,3),' %, frboot=',num2str(frboot,2),' %, ne(0)=',num2str(ne0),' 10^1^9 m^-^3'])
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
