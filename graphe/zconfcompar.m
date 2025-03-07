if ~exist('jeux1')
  jeux1 = [];
end
Ipf    = data.gene.ip/1e6;
bt     = data.geo.b0;
nbarf  = data.gene.nbar/1e19;
if ~isempty(jeux1)
  nIpf   = jeux1.data.gene.ip/1e6;
  nbt    = jeux1.data.geo.b0;
  nnbarf = jeux1.data.gene.nbar/1e19;
end
M=1;
%
%  plth = puissance totale absorbee (on enleve la puissance perdue dans le ripple) 
%
% indice 0 : ploss (prise en compte de dW/dt)
% indice 1 : paddtot (sans ripple)
%
plth0    = data.gene.ploss/1e6;
plth1   = (sum(abs(data.cons.fci'))' + ...
             sum(abs(data.cons.hyb'))' + ...
	     data.gene.paddohm)/1e6;
%plth1    = (data.gene.paddfci+data.gene.paddhyb+data.gene.paddohm)/1e6;
tauth0   = data.gene.wth./plth0/1e6;
tauth1   = data.gene.wth./plth1/1e6;
R0       = data.geo.r0;
af       = data.geo.a;
if ~isempty(jeux1)
  nplth0   = jeux1.data.gene.ploss/1e6;
  nplth1   = (jeux1.data.gene.paddfci+jeux1.data.gene.paddhyb+jeux1.data.gene.paddohm)/1e6;
  nplth1   = (sum(abs(jeux1.data.cons.fci'))' + ...
             sum(abs(jeux1.data.cons.hyb'))' + ...
	     jeux1.data.gene.paddohm)/1e6;
  ntauth0  = jeux1.data.gene.wth./nplth0/1e6;
  ntauth1  = jeux1.data.gene.wth./nplth1/1e6;
  nR0      = jeux1.data.geo.r0;
  naf      = jeux1.data.geo.a;
end
%
% pinj0 = puissance injectee
%
pinj0    = (sum(abs(data.cons.fci'))'+sum(abs(data.cons.hyb'))'+data.gene.paddohm)/1e6;
tauinj0  = data.gene.wth./pinj0/1e6;
temps    = data.gene.temps;
tauL     = zeros(size(temps));
tauL1    = tauL;
if ~isempty(jeux1)
  npinj0    = (sum(abs(jeux1.data.cons.fci'))'+sum(abs(jeux1.data.cons.hyb'))'+jeux1.data.gene.paddohm)/1e6;
  ntauinj0  = jeux1.data.gene.wth./npinj0/1e6;
  ntemps    = jeux1.data.gene.temps;
  temps1    = ntemps;
  ntauL     = zeros(size(ntemps));
  ntauL1    = ntauL;
end

for k=param.gene.kmin:param.gene.k

  tauL(k)  = H97(Ipf(k),bt(k),nbarf(k),plth0(k),R0(k),1,af(k),M);
  tauL1(k) = H97(Ipf(k),bt(k),nbarf(k),plth1(k),R0(k),1,af(k),M);

end
if ~isempty(jeux1)
  for k=jeux1.param.gene.kmin:jeux1.param.gene.k

    ntauL(k)  = H97(nIpf(k),nbt(k),nnbarf(k),nplth0(k),nR0(k),1,naf(k),M);
    ntauL1(k) = H97(nIpf(k),nbt(k),nnbarf(k),nplth1(k),nR0(k),1,naf(k),M);
  end
end
oui=1;
if oui
  if ~exist('tauf','var')
    [tauf,ttauf]=tsbase(fix(param.from.shot.num),'staudia');
    if ~isempty(jeux1)
      [ntauf,nttauf]=tsbase(fix(jeux1.param.from.shot.num),'staudia');
    end
  end
  inderr = find(diff(ttauf)==0);
  if ~isempty(inderr) 
     ttauf(inderr) = [];
     tauf(inderr) = [];
   end
  if ~isempty(jeux1)
    inderr = find(diff(nttauf)==0);
    if ~isempty(inderr) 
      nttauf(inderr) = [];
      ntauf(inderr) = [];
    end
  end
  if ~exist('wdia','var')
    [wdia,tdia]=tsbase(fix(param.from.shot.num),'swdia');
    if ~isempty(jeux1)
      [wdia1,tdia1]=tsbase(fix(jeux1.param.from.shot.num),'swdia');
    end
  end
  inderr = find(diff(tdia)==0);
  if ~isempty(inderr) 
     tdia(inderr) = [];
     wdia(inderr) = [];
  end
  nwdia = interp1(tdia,wdia*1e6,data.gene.temps);
  if ~isempty(jeux1)
    inderr = find(diff(ntdia)==0);
    if ~isempty(inderr) 
      ntdia(inderr) = [];
      nwdia(inderr) = [];
    end
    nwdia1 = interp1(tdia1,wdia1*1e6,jeux1.data.gene.temps);
  end
%
% tauexp : taudia mesure a partir de Wdia et ploss
%
  tauexp   = interp1(ttauf,tauf,data.gene.temps);
  if ~isempty(jeux1)
    ntauexp  = interp1(nttauf,ntauf,jeux1.data.gene.temps);
  end
%
% tauexp : taudia corrig� des pertes ripple
%
%  tauexp  = tauexp .* plth1 ./ plth;
%  ntauexp = ntauexp .* nplth1 ./ nplth;
%
end
h = findobj(0,'type','figure','tag','conf1');
if isempty(h)
       h=figure('tag','conf1');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
	
%
% protection pour le graphique
%
ver           = (plth1 < 15 & plth1 > 1) | tauL < 1;
ver(ver==0)   = -1;
if ~isempty(jeux1)
  ver1          = (nplth1 < 15 & nplth1 > 1) | ntauL < 1;
  ver1(nver==0) = -1;
end
%
% puissance de suprafci
%
for k=1:length(temps)
  psupfci(k,:)  =trapz(param.gene.x,data.source.fci.psupra(k,:).*data.equi.vpr(k,:)).*data.equi.rhomax(k);
end
wdiac                      = data.gene.wdia;
indc                       = find(wdiac > mean(wdiac)*3);
wdiac(indc)                = NaN;
sup                        = smooth(abs(data.gene.wdia-data.gene.wth)./data.gene.wth,5);
sup(data.gene.paddfci<1e5) = 0;
sup(sup>0.5)               = 0;
rap0                       = tauexp./tauL;
rap1                       = tauexp./tauL1;

if ~isempty(jeux1)
  for k=1:length(ntemps)
    npsupfci(k,:)  =trapz(param.gene.x,jeux1.data.source.fci.psupra(k,:).*jeux1.data.equi.vpr(k,:)).*jeux1.data.equi.rhomax(k);
  end  
  nwdiac                            = jeux1.data.gene.wdia;
  nindc                             = find(nwdiac > mean(nwdiac)*3);
  nwdiac(nindc)                      = NaN;
  sup1                              = smooth(abs(nwdiac-jeux1.data.gene.wth)./jeux1.data.gene.wth,5);
  sup1(jeux1.data.gene.paddfci<1e5) = 0;
  sup1(sup1>0.5)                    = 0;
  nrap0                             = ntauexp./ntauL;
  nrap1                             = ntauexp./ntauL1;
end


ver(rap0>2)      = -1;
ver(rap0<0.1)    = -1;
if ~isempty(jeux1)
  ver1(nrap0>2)    = -1;
  ver1(nrap0<0.1)  = -1;
end

subplot(2,2,1)

if ~isempty(jeux1)
  plot(temps,(tauexp.*(1-sup)./tauL1).*ver,'b',...
      ntemps,(ntauexp.*(1-sup1)./ntauL1).*ver1,'r--')
  axis([param.gene.tdeb param.gene.t 0 2])
  legend('Cronos 1 ','Cronos 2 ')
  grid
  title(['scaling L mode (sup effect)'])
  ylabel('H97')
  xlabel('time (s)')
else
  plot(temps,(tauexp.*(1-sup)./tauL1).*ver,'b')
  axis([param.gene.tdeb param.gene.t 0 2])
  legend(['shot ',int2str(param.from.shot.num)])
  grid
  title(['scaling L mode (sup effect)'])
  ylabel('H97')
  xlabel('time (s)')
end

subplot(2,2,2)

if ~isempty(jeux1)
  plot(temps,(tauexp./tauL1).*ver,'b',...
      temps1,(ntauexp./ntauL1).*ver1,'r--')
  axis([param.gene.tdeb param.gene.t 0 2])
  legend('Cronos 1','Cronos 2')
  grid
  title(['scaling L mode (no sup)'])
  ylabel('H97')
  xlabel('time (s)')
else
  plot(temps,(tauexp./tauL).*ver,'b')
  axis([param.gene.tdeb param.gene.t 0 2])
  legend(['shot ',int2str(param.from.shot.num)],['shot ',int2str(jeux1.param.from.shot.num)])
  grid
  title(['scaling L mode (no sup.)'])
  ylabel('H97')
  xlabel('time (s)')
end


subplot(2,2,3)

plot(temps,sup*100,'b',temps1,sup1*100,'r.-')
ylabel('frsup%')
axis([param.gene.tdeb param.gene.t 0 40])


subplot(2,2,4)
frrip=(sum(abs(data.cons.fci'))'-data.gene.paddfci)./sum(abs(data.cons.fci'))'*100;
frrip1=(sum(abs(jeux1.data.cons.fci'))'-jeux1.data.gene.paddfci)./sum(abs(jeux1.data.cons.fci'))'*100;
frrip(frrip>50) = 0;
frrip(frrip<0) = 0;
frrip1(frrip1>50) = 0;
frrip1(frrip1<0) = 0;
plot(temps,frrip,'b',temps1,frrip1,'r.-')
axis([param.gene.tdeb param.gene.t 0 40])

ylabel('frrip %')

