function sor = zbatchconf(data,param)

h = findobj(0,'type','figure','tag','conf');
if isempty(h)
       h=figure('tag','conf');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

tps = data.gene.temps;
if strcmp(param.from.machine,'TS')
  [tauv,ttauv]=tsbase(fix(param.from.shot.num),'staudia');
  tauv = interp1(ttauv,tauv,data.gene.temps);
else
  tauv = data.gene.taue;
end
if size(data.cons.fci,2) > 1
  pconsfci = sum(abs(data.cons.fci'))';
else
  pconsfci = abs(data.cons.fci);
end
indh = find(pconsfci>0.5e6 & ~isnan(data.gene.paddfci) & tps > param.gene.tdeb & tps < param.gene.t);
hachcons = mean(abs(diff(pconsfci(indh))));
hachfci = mean(abs(diff(data.gene.paddfci(indh))));
if hachfci < hachcons*10
padd = data.gene.paddfci+data.gene.paddhyb+data.gene.paddfus+data.gene.paddfce+data.gene.paddohm;
else
padd = pconsfci+data.gene.paddhyb+data.gene.paddfus+data.gene.paddfce+data.gene.paddohm;
end
ploss_rec=padd-gradient(data.gene.wdia,tps(2)-tps(1));
%betap 
betap=data.gene.betap; 
Ipla=data.gene.ip/1e6;
a = data.equi.a;
Btor = data.geo.b0;
%betan 
betan=(4*betap.*Ipla)./(a(:,101).*Btor);
Iboot = data.gene.iboot;
Iidn = data.gene.iidn;
Ilh = data.gene.ihyb;
nbar = data.gene.nbar/1e19;
ngreen = 10 * Ipla ./pi ./a(:,end) ./a(:,end);
fr = nbar ./ ngreen;
sor.ngreen = ngreen;


rlwcr = rlw(Ipla,data.equi.rmoy(:,end),a(:,end),data.gene.nemoy/1e19,Btor,data.gene.zeffm,ploss_rec/1e6,data.equi.e(:,end));
Hrlw = data.gene.we./rlwcr/1e6;

%Hfactors 
M=2

%scaling ITER L-mode - tau L thermique - NF39(1999)2206
tauL=0.023*(data.geo.r0).^1.89*M^0.2.*(data.cons.ip/1e6).^0.96.*(data.geo.b0).^0.03.*(data.equi.e(:,101)).^0.64.*...
(data.geo.a).^-0.06.*(data.gene.nbar./2e19./data.geo.a).^0.4.*(ploss_rec/1e6).^-0.73;
tauL(find(isnan(tauL)==1))=0;
tauL(find(iscomplex(tauL)==1))=0;

%scaling ITER H-mode - IPB98(y,2) - NF39(1999)2208
tauH=0.0562*(data.geo.r0).^1.97*M^0.19.*(data.cons.ip/1e6).^0.93.*(data.geo.b0).^0.15.*(data.equi.e(:,101)).^0.78.*...
(data.geo.a./data.geo.r0).^0.58.*(data.gene.nbar./2e19./data.geo.a).^0.41.*(ploss_rec/1e6).^-0.69;
tauH(find(isnan(tauH)==1))=0;
tauH(find(iscomplex(tauH)==1))=0;

%scaling ITER H-mode - ITER89-P - NF37(1997)1303
tauH2=0.038*(data.geo.r0).^1.5*M^0.5.*(data.cons.ip/1e6).^0.85.*(data.geo.b0).^0.2.*(data.equi.e(:,101)).^0.5.*...
(data.geo.a./data.geo.r0).^0.3.*(data.gene.nbar./2e19./data.geo.a).^0.1.*(ploss_rec/1e6).^-0.5;
tauH2(find(isnan(tauH2)==1))=0;
tauH2(find(iscomplex(tauH2)==1))=0;


tstart = param.gene.tdeb;
tend   = param.gene.t;
x = param.gene.x;
subplot(2,2,1)
if strcmp(param.from.machine,'TS')
  plot(tps,tauv./tauL,'r',tps,Hrlw,'b');
  legend('H-L','Hrlw')
  axis([tstart tend 0.5 2]) 
else
  plot(tps,tauv./tauL,'r',tps,tauv./tauH,'g',tps,tauv./tauH2,'b'); 
  legend('H-L','H98th','H89-P')
  axis([tstart tend 0.5 6]) 
end


title('confinement, M = 2')
xlabel('t (s)')
sor.Hhl = tauv./tauL;
sor.H98th = tauv./tauH;
sor.H89_P=tauv./tauH2;
sor.tauL = tauL;
sor.tauH=tauH;
sor.tauH2=tauH2;

subplot(2,2,2)
plot(tps,betap,tps,betan,tps,fr)
title('beta')
legend('b_p','b_N','fgreen')
axis([tstart tend 0 ceil(max(max(betan),max(betap)))])
xlabel('t (s)')
grid
sor.betaN=betan;
sor.betap = betap;
sor.fgreen = fr;

subplot(2,2,3)
plot(tps,Ipla,tps,Iboot/1e6,tps,Iidn/1e6,tps,Ilh/1e6)
legend('Itot','Iboot','Inbi','Ilh')
ylabel('MA')
xlabel('t (s)')
sor.fboot = Iboot./Ipla;
sor.fini = (Iboot + Iidn + Ilh) ./ Ipla;
axis([tstart tend 0 ceil(max(Ipla))])
title('current')
subplot(4,2,6)
[qm,iqm]=min(data.prof.q,[],2);
plot(tps,qm)
title('q_m_i_n')
grid
ind = find(tps > tstart & tps < tend);
axis([tstart tend 0 ceil(max(qm(ind)))])
subplot(4,2,8)
plot(tps,x(iqm))
title('r q_m_i_n')

axis([tstart tend 0 1])
ylabel('x')
xlabel('t (s)')


