function sor = zbatchconf(data,param)

% Fortement corrige (FI, Novembre 2007)
% Est maintenant adapte a verifier les sorties de zkiautonew (memes definitions)


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
  ind=find(diff(ttauv)==0);
  if ~isempty(ind)
     ttauv(ind)=[];
     tauv(ind)=[];
  end   
  ind=find(ttauv>0& tauv>0);
  tauv = tauv(ind);
  ttauv = ttauv(ind);
  tauv = interp1(ttauv,tauv,data.gene.temps);
else
  tauv = data.gene.tauth;   %taue avant. Completement faux ???  Recalcule plus bas !!!
end
if strcmp(param.from.machine,'DIIID')
  racine = strcat(getenv('HOME'),'/zineb/data/diiid');
  racine = [racine,'/',int2str(param.from.shot.num)];
  eval(['load ',racine,'/diiidtemp']);
end
if size(data.cons.fci,2) > 1
  pconsfci = sum(abs(data.cons.fci'))';
else
  pconsfci = abs(data.cons.fci);
end
if isnan(pconsfci)
  pconsfci = 0 * size(pconsfci);
  hachfci  = 0;
  hachcons = 1;
else
  indh = find(pconsfci>0.5e6 & ~isnan(data.gene.paddfci) & tps > param.gene.tdeb & tps < param.gene.t); 
  hachcons = mean(abs(diff(pconsfci(indh))));
  hachfci = mean(abs(diff(data.gene.paddfci(indh))));
end
%if hachfci < hachcons*10
%padd = data.gene.paddfci+data.gene.paddhyb+data.gene.paddfus+data.gene.paddfce+data.gene.paddohm+data.gene.paddidn;


%%%%%%%%%%%%%% CALCUL DE PADD (PLOSS) : COMME DANS CHIAUTO   %%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(param.fonction.coefa,'zkiautonew')
   inte = data.source.totale.el + data.source.totale.ion + (1 - param.cons.coefa.frad) .* data.source.prad;
else
   inte = data.source.totale.el + data.source.totale.ion + data.source.prad;
end
ind    = find(~isfinite(inte));
   if ~isempty(ind)
      inte(ind) = 0;
   end
padd  = data.equi.rhomax .* trapz(param.gene.x,data.equi.vpr .* inte,2);   % same definition of Ploss as in zkiautonew
disp('Ok')
%else
%padd = pconsfci+data.gene.paddhyb+data.gene.paddfus+data.gene.paddfce+data.gene.paddohm+data.gene.paddidn;
%disp('Pas Ok')
%end
%ploss_rec=padd-gradient(data.gene.wdia,tps(2)-tps(1));
ploss_rec = padd - data.gene.dwthdt;
tauv = data.gene.wth ./ ploss_rec;       % Ceci est le tauTH de la simu Cronos, tenant compte du dW/dt

%betap 
betap=data.gene.betap; 
Ipla=data.gene.ip/1e6;
a = data.equi.a;
Btor = data.geo.b0;
%betan 
betan=data.gene.beta./Ipla.*a(:,101).*Btor.*100;
Iboot = data.gene.iboot;
Iidn = data.gene.iidn;
Ifce = data.gene.ifce;
Ilh = data.gene.ihyb;
nbar = data.gene.nbar/1e19;
ngreen = 10 * Ipla ./pi ./a(:,end) ./a(:,end);
fr = nbar ./ ngreen;
sor.ngreen = ngreen;
q95 = data.prof.q(:,95);

rlwcr = rlw(Ipla,data.equi.rmoy(:,end),a(:,end),data.gene.nemoy/1e19,Btor,data.gene.zeffm,ploss_rec/1e6,data.equi.e(:,end));
Hrlw = data.gene.we./rlwcr/1e6;

%Hfactors 

% effective mass
x = param.gene.x;
nions   = data.impur.impur;  % size(nions) = [nbg,nbrho] !
M    = trapz(x, squeeze(sum(nions(:,:,1) .* param.compo.a(1) + nions(:,:,2) .* param.compo.a(2) + nions(:,:,3) .* param.compo.a(3)  + ...
                    nions(:,:,4) .* param.compo.a(4) + nions(:,:,5) .*param.compo.a(5),3) ./  sum(nions(:,:,1) .* param.compo.z(1) +  ...
		    nions(:,:,2) .* param.compo.z(2) + nions(:,:,3) .* param.compo.z(3)  + ...
                    nions(:,:,4) .* param.compo.z(4) + nions(:,:,5) .*param.compo.z(5),3)),2);



%scaling ITER L-mode - tau L thermique - NF39(1999)2206
tauL=0.023*(data.equi.raxe(:,end)).^1.89.*M.^0.2.*(data.cons.ip/1e6).^0.96.*(data.geo.b0).^0.03.*(data.equi.e(:,101)).^0.64.*...
(data.equi.a(:,end)./data.equi.raxe(:,end)).^-0.06.*(trapz(x,data.prof.ne,2) ./ 1e19).^0.4.*(padd/1e6).^-0.73;
tauL(find(isnan(tauL)==1))=0;
tauL(find(iscomplex(tauL)==1))=0;

%scaling ITER H-mode - IPB98(y,2) - NF39(1999)2208
tauH=0.0562*(data.equi.raxe(:,end)).^1.97.*M.^0.19.*(data.cons.ip/1e6).^0.93.*(data.geo.b0).^0.15.*(data.equi.e(:,101)).^0.78.*...
(data.equi.a(:,end)./data.equi.raxe(:,end)).^0.58.*(trapz(x,data.prof.ne,2) ./ 1e19).^0.41.*(padd/1e6).^-0.69;
tauH(find(isnan(tauH)==1))=0;
tauH(find(iscomplex(tauH)==1))=0;
G = betan.*(tauv./tauH)./q95./q95;

%scaling ITER H-mode - ITER89-P - NF37(1997)1303
tauH2=0.038*(data.equi.raxe(:,end)).^1.5.*M.^0.5.*(data.cons.ip/1e6).^0.85.*(data.geo.b0).^0.2.*(data.equi.e(:,101)).^0.5.*...
(data.equi.a(:,end)./data.equi.raxe(:,end)).^0.3.*(trapz(x,data.prof.ne,2) ./ 1e19).^0.1.*(padd/1e6).^-0.5;
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


title('confinement plot')
xlabel('t (s)')

sor.Hhl = tauv./tauL;
sor.H98th = tauv./tauH;
sor.H89_P=tauv./tauH2;

set(gca,'xlim',[tstart tend])

sor.tauL = tauL;
sor.tauH=tauH;
sor.tauH2=tauH2;

subplot(2,2,2)



if strcmp(param.from.machine,'DIIID')
  plot(tps,betap,tps,betan,tps,fr,tps,G,diiidtemp.tbp,diiidtemp.bp,'r--')
  title(['beta (-- ',diiidtemp.nomefit,')'])
  legend('b_p','b_N','fgreen','G')
else
  plot(tps,betap,tps,betan,tps,fr,tps,G)
  title('beta')
  legend('b_p','b_N','fgreen','G')
end

axis([tstart tend 0 ceil(max(max(betan),max(betap)))])
xlabel('t (s)')
grid
sor.betaN=betan;
sor.betap = betap;
sor.fgreen = fr;

subplot(2,2,3)
if Iidn > 0
plot(tps,Ipla,tps,Iboot/1e6,tps,Iidn/1e6,tps,Ilh/1e6,tps,Ifce/1e6)
legend('Itot','Iboot','Inbi','Ilh','Ieccd')
else
plot(tps,Ipla,tps,Iboot/1e6,tps,Ilh/1e6,tps,Ifce/1e6)
legend('Itot','Iboot','Ilh','Ieccd')
end
ylabel('MA')
xlabel('t (s)')
sor.fboot = Iboot./Ipla;
sor.fini = (Iboot + Iidn + Ilh) ./ Ipla;
axis([tstart tend 0 ceil(max(Ipla))])
title('current')
subplot(4,2,6)
[qm,iqm]=min(data.prof.q(:,10:end),[],2);
if strcmp(param.from.machine,'DIIID')
  plot(tps,qm,diiidtemp.tqmin,diiidtemp.qmin,'r--')
  title(['q_m_i_n (-- ',diiidtemp.nomefit,')'])
else
  plot(tps,qm)
  title('q_m_i_n')
end

grid
ind = find(tps > tstart & tps < tend);
axis([tstart tend floor(min(qm(ind))) ceil(max(qm(ind)))])
subplot(4,2,8)
plot(tps,x(iqm))
title('r q_m_i_n')

axis([tstart tend 0 1])
ylabel('x')
xlabel('t (s)')


