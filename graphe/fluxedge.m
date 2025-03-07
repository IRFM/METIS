% figure pour verification de la variation temporelle du flux magnetique au
% bord
h = findobj(0,'type','figure','tag','flux');
if isempty(h)
       h=figure('tag','flux');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
if param.gene.t == param.gene.tdeb
    param.gene.t = param.gene.tfin;
    param.gene.k = param.gene.kmax;
end
numchoc=param.from.shot.num;
if strcmp(param.from.machine,'TS')
  [flux,tflux]=tsbase(fix(numchoc),'GFLUXNOY%4');
end
if strcmp(param.from.machine,'JET')
  if all(isnan(data.cons.flux))  % can happen if the simulation has been created via JAMS --> need to reconstruct the flux as during the creation of a "normal" Cronos file
    clear choix     % prepare dummy options for using zjetprofile
    choix = zjetacces;
    choix.efit = 0;
    choix.exfile = 0;
    choix.fitTe = 3;
    choix.creux = 1;
    choix.modezeff = 0;
    choix.normnl = 1;
    choix.helena = 2; % -> derniere surface de flux a 0.995
    choix.transp = 1;
    choix.centre = 2; % -> utilisation de zmag 
    racine = ['/home/gsem/basiuk/zineb/data/JET/',int2str(param.from.shot.num)];
    clear optionex
    optionex.te = 0;
    optionex.ne = 0;
    optionex.ti = 0;
    optionex.rf = 0;
    optionex.lh = 0;
    optionex.zeff = 0;
    optionex.nbi = 0;
    rti = 1;  % dummy
    [temps,jet] = zprofilejet(param.from.shot.num,param.phys,data.gene.temps',param.gene.nbrho,2,choix,optionex,rti,racine);
     flux = 2.*pi.*param.gene.signe.ip*(jet.psip(end)+(jet.fluxbord-jet.fluxbord(1))');
  else
     flux = data.cons.flux*2*pi;
  end
  tflux = data.gene.temps;
end
if strcmp(param.from.machine,'DIIID')
  flux = data.cons.flux*2*pi;
  tflux = data.gene.temps;
end
if strcmp(param.from.machine,'HL2A')
  flux = data.cons.flux*2*pi;
  tflux = data.gene.temps;
end

if ~exist('jeux1')
  %chargedata
  jeux1.data.prof=[];  
end
if isempty(jeux1)
  jeux1.data.prof=[];
end
if isfield(jeux1.data.prof,'psi')
  subplot(2,1,1)
  fluxc = data.prof.psi(:,end)*2*pi;
  tfluxc = data.gene.temps;
  indf = min(find(tflux > tfluxc(1)));
  corf = flux(indf)-fluxc(1);
  plot(tflux,flux-corf,'k.',data.gene.temps,fluxc,'r');
  hold on
  fluxc1  = jeux1.data.prof.psi(:,end)*2*pi;
  tfluxc1 = jeux1.data.gene.temps;
  indf1   = min(find(tflux > tfluxc1(1)));
  corf1 = flux(indf1)-fluxc1(1);

  plot(jeux1.data.gene.temps,fluxc1,'b--');
  xlabel('time (s)')
  ylabel('Flux (Wb)')
  fmax = ceil(max([max(data.prof.psi(:,end)*2*pi),...
             max(jeux1.data.prof.psi(:,end)*2*pi),...
             max(flux-corf)]));
  fmin = floor(min([min(data.prof.psi(:,end)*2*3.14159),...
             min(jeux1.data.prof.psi(:,end)*2*3.14159),...
             min(flux-corf)]));
  axis([param.gene.tdeb param.gene.t fmin fmax])
  legend('Exp.','Cronos','Ref.')
  title(['poloidal edge flux, #',num2str(param.from.shot.num)])
  grid
  subplot(2,2,3)
  ind1  = param.gene.kmin:param.gene.k;
  ind2  = jeux1.param.gene.kmin:jeux1.param.gene.k;
  val1 = interp1(tflux,flux-corf,data.gene.temps);
  valc1 = data.prof.psi(ind1,end)*2*pi;
  val2 = interp1(tflux,flux-corf,jeux1.data.gene.temps);
  valc2 = jeux1.data.prof.psi(ind2,end)*2*pi;
  fluxf = val1(end);
  if sign(val1(1))*sign(val1(end)) < 0
    fluxe = val1(end);
    val1  = val1-fluxe;
    val2  = val2-fluxe;
    valc1 = valc1 - fluxe;
    valc2 = valc2 - fluxe;
  end
  errc1 = abs((val1(ind1)-valc1)./val1(ind1))*100;
  errc2 = abs((val2(ind2)-valc2)./val2(ind2))*100;

  plot(data.gene.temps(ind1),errc1,'r',jeux1.data.gene.temps(ind2),errc2,'b--')
  fmin = floor(min([min(errc1),min(errc2)]));
  fmax = ceil(max([max(errc1),max(errc2)]));
  axis([param.gene.tdeb param.gene.t 0 min(100,fmax)])
  legend('Cronos','Ref.')
  title('deviation')
  xlabel('time (s)')
  ylabel('%')
  grid;
  subplot(2,2,4)
  val1 = smooth(gradient(flux-corf)./gradient(tflux),5);
  plot(tflux,val1,'k+',data.gene.temps,gradient(fluxc)./gradient(tfluxc),'r',...
                       jeux1.data.gene.temps,gradient(fluxc1)./gradient(tfluxc1),'b--');
  title('dpsi/dt')
  indg=find(tflux>param.gene.tdeb & tflux<param.gene.t );
  val1 = gradient(flux(indg)-corf)./gradient(tflux(indg));
  axis([param.gene.tdeb param.gene.t min(val1) max(val1)]);

else
  subplot(2,1,1)
  fluxc = data.prof.psi(:,end)*2*pi;
  tfluxc = data.gene.temps;
  indf = min(find(tflux > tfluxc(1)));
  corf = flux(indf)-fluxc(1);
  plot(tflux,flux-corf,'k+',data.gene.temps,fluxc,'r');

  xlabel('time (s)')
  ylabel('Flux (Wb)')
  fmax = ceil(max([max(fluxc),...
             max(flux-corf)]));
  fmin = floor(min([max(fluxc),...
             min(flux-corf)]));
  axis([param.gene.tdeb param.gene.t fmin fmax])
  legend('Exp.','Cronos')
  grid
  subplot(2,2,3)
  val1 = interp1(tflux,flux-corf,data.gene.temps);
  ind1  = param.gene.kmin:param.gene.k;
  errc1 = (val1(ind1)-fluxc(ind1));
  plot(data.gene.temps(ind1),errc1,'r')
  fmin = floor(min(errc1));
  fmax = ceil(max(errc1));
  axis([param.gene.tdeb param.gene.t fmin fmax])
  title('deviation')
  xlabel('time (s)')
  ylabel('Wb')
  grid;
  subplot(2,2,4)
  plot(tflux,smooth(gradient(flux-corf)./gradient(tflux),5),'k+',data.gene.temps,gradient(fluxc)./gradient(tfluxc),'r');
  title('dpsi/dt')
  indg=find(tflux>param.gene.tdeb & tflux<param.gene.t );
  val1 = gradient(flux(indg)-corf)./gradient(tflux(indg));
  axis([param.gene.tdeb param.gene.t min(val1) max(val1)]);

end
