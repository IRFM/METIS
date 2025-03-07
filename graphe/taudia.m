function tau = taudia(data,param,jeux1)



[taud,ttaud]=tsbase(param.from.shot.num,'staudia');

ip     = data.cons.ip/1e6;
rgeo   = data.geo.r0;
amin   = data.geo.a;
kappa  = data.equi.e(:,1);
nel    = data.gene.nemoy/1e19;
bt     = data.geo.b0;
meff   = 2;
plth   = data.gene.paddtot/1e6;

tps    = data.gene.temps;
ploss_rec  = data.gene.ploss+data.gene.dwdiadt-gradient(data.gene.wdia,tps(2)-tps(1));
gene.temps = data.gene.temps;
gene.x = param.gene.x;
[frac,Eloss,Ilost,Emean,frlost,iso]=...
ripple_fci(data.cons.fci,data.geo,data.prof,data.impur,param.phys,gene,param.compo,param.cons.fci,data.equi);
tau.Emean = Emean;
tau.Eloss = Eloss;
tau.frlost = frlost;
tau.Ilost = Ilost;
tau.iso = iso;
%scaling ITER L-mode - tau L thermique - NF39(1999)2206
tauL=scallaw('iter0',ip,rgeo,amin,kappa,nel,bt,meff,plth);
tauL(find(isnan(tauL)==1))=0;
tauL(find(iscomplex(tauL)==1))=0;
tau.tauL = tauL;
if nargin == 3
  [taudn,ttaudn]=tsbase(jeux1.param.from.shot.num,'staudia');
  gene.temps = jeux1.data.gene.temps;
  gene.x = param.gene.x;
  [fracn,Elossn,Ilostn,Emeann,frlostn,ison]=...
  ripple_fci(jeux1.data.cons.fci,jeux1.data.geo,jeux1.data.prof,jeux1.data.impur,param.phys,gene,jeux1.param.compo,jeux1.param.cons.fci,jeux1.data.equi);

  tpsn       = jeux1.data.gene.temps;
  ploss_recn = jeux1.data.gene.ploss+jeux1.data.gene.dwdiadt-gradient(jeux1.data.gene.wdia,tpsn(2)-tpsn(1));
  ipn        = jeux1.data.cons.ip/1e6;
  rgeon      = jeux1.data.geo.r0;
  aminn      = jeux1.data.geo.a;
  kappan     = jeux1.data.equi.e(:,1);
  neln       = jeux1.data.gene.nemoy/1e19;
  btn        = jeux1.data.geo.b0;
  meffn      = 2;
  plthn      = jeux1.data.gene.paddtot/1e6;

%scaling ITER L-mode - tau L thermique - NF39(1999)2206
  tauLn=scallaw('iter0',ipn,rgeon,aminn,kappan,neln,btn,meff,plthn);

  tauLn(find(isnan(tauLn)==1))=0;
  tauLn(find(iscomplex(tauLn)==1))=0;

  tau.Emeann = Emeann;
  tau.Elossn = Elossn;
  tau.frlostn = frlostn;
  tau.Ilostn = Ilostn;
  tau.ison = ison;

  tau.tauLn = tauLn
end

h = findobj(0,'type','figure','tag','taudia');
if isempty(h)
       h=figure('tag','taudia');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

if nargin == 3
    subplot(2,1,1)
    plot(ttaud,smooth(taud,3),'r',tps,tau.tauL,'b--',...
         ttaudn,smooth(taudn,3),'m',tpsn,tau.tauLn,'g--' )
    legend(['Mesure #',int2str(param.from.shot.num)],...
            'scaling',['Ref.#',int2str(jeux1.param.from.shot.num)],'scal Ref.')
    ylabel('s')
    xlabel('time (s)')
    axis([min(tps) max(tps) 0 0.4])
    subplot(2,1,2)
    tauexp = smooth(interp1(ttaud,taud,tps,'nearest'),3);
    tauexpn = smooth(interp1(ttaud,taud,tpsn,'nearest'),3);
    ecart  = (tauexp-tau.tauL)./tau.tauL*100;
    ecart1 = (tauexpn-tau.tauLn)./tau.tauLn*100;
    indfci = data.gene.paddfci/1e6>0.4;
    indfcin = jeux1.data.gene.paddfci/1e6>0.4;
    if sum(indfci) > 2 & sum(indfcin) > 2
      plot(ecart(indfci),tau.tauL(indfci),'ro',...
           ecart(~indfci),tau.tauL(~indfci),'m+',...
           ecart1(indfcin),tau.tauLn(indfcin),'b*',...
           ecart1(~indfcin),tau.tauLn(~indfcin),'gs')
      legend(['#',int2str(param.from.shot.num),' FCI'],...
             'Ohmic phase',...
              ['#',int2str(jeux1.param.from.shot.num),' Ref. (FCI)'],'Ref. (Ohmic)')
    else
      plot(ecart,tau.tauL,'ro',ecart1,tau.tauLn,'b*')
    end
     ylabel('tau (s)')
    xlabel('ecart %')
    title('ecart taudia p/r aux mesures')
end
